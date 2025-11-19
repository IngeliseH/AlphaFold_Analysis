def find_int_overlap(loc, loc_list):
    """
    For a list of interface locations (loc_list), find how many overlap with at least 75% of a particular
    interface (loc)
    For locations with multiple chains, all chain combinations are considered and the full interface taken into account when calculating percentage
    """
    total = 0
    loc = [set(l) for l in loc]
    total_loc_len = sum(len(l) for l in loc)
    if total_loc_len == 0:
        return 0
    
    def is_overlap(loc, elem):
        """Helper function to determine if overlap is greater than 75%"""
        overlap = 0
        for l, e in zip(loc, elem):
            overlap += len(l.intersection(e))
        if overlap / total_loc_len > 0.75:
            return True
        return False
    
    # for each interface location in loc_list, check what percent of the current location is in the interface location
    for elem in loc_list:
        elem = [set(c) for c in elem]
        if len(loc) == 1 and len(elem) == 1:
            if is_overlap(loc, elem):
                total += 1
        elif len(loc) == 2 and len(elem) == 2:
            if is_overlap(loc, elem) or is_overlap(loc, elem[::-1]):
                total += 1
        elif (len(loc) == 1 and len(elem) == 2):
            if is_overlap(loc, [elem[0]]) or is_overlap(loc, [elem[1]]):
                total += 1
        elif (len(loc) == 2 and len(elem) == 1):
            if is_overlap([loc[0]], elem) or is_overlap([loc[1]], elem):
                total += 1
        else:
            print(f"Trying to compare incorrect chain numbers in promiscuity calculation: current loc has {len(loc)} chains, list element has {len(elem)}")
    return total

def find_interface_promiscuity(df):
    """
    find promiscuity of each domain pair in the dataframe by finding the number of other domain pairs
    that overlap with the interface location of the current domain pair
    """
    print("Finding interface promiscuity...")
    # build dict of domain names and all interface locations found involving that domain
    from collections import defaultdict
    from ast import literal_eval
    from analysis_utility import expand_loc_range
    def get_chain_groupings(row):
        dimer1 = False
        chain_group1, chain_group2 = ('A'), ('B')
        if '_dimer' in row['Protein1']:
            dimer1 = True
            chain_group1, chain_group2 = ('A', 'B'), ('C')
        if '_dimer' in row['Protein2']:
            chain_group2 = ('C', 'D') if dimer1 else ('B', 'C')
        return [chain_group1, chain_group2]
    
    def get_loc(chain_group, locations, protein_naming, i):
        loc = []
        for chain in chain_group:
            chain_key = f'Protein{i+1}' if protein_naming else f'Chain {chain}'
            if locations.get(chain_key) == 'None':
                continue
            loc.append(expand_loc_range(locations[chain_key]))
        return loc

    domain_int_locs = defaultdict(list)
    for _, row in df.iterrows():
        if row['location']:
            chain_groupings = get_chain_groupings(row)
            try:
                locations = literal_eval(row['location'])
                if not locations:
                    continue
                protein_naming = True if 'Protein1' in locations else False
                for index, chain_group in enumerate(chain_groupings):
                    loc = get_loc(index, chain_group, locations, protein_naming)
                    domain_int_locs[row[f'Protein{i+1}_Domain']].append(loc)
            except:
                continue

    # iterate through again and find promiscuity for each domain
    df['p1d_promiscuity'] = None
    df['p2d_promiscuity'] = None
    df['total_promiscuity'] = None
    df['max_promiscuity'] = None
    for _, row in df.iterrows():
        if row['location']:
            promiscuitys = []
            chain_groupings = get_chain_groupings(row)
            try:
                locations = literal_eval(row['location'])
                if not locations:
                    continue
                protein_naming = True if 'Protein1' in locations else False # maintain compatibility with outdated 'protein' naming
                for index, chain_group in enumerate(chain_groupings):
                    loc = get_loc(index, chain_group, locations, protein_naming)
                    promiscuity = find_int_overlap(loc, domain_int_locs[row[f'Protein{i+1}_Domain']])
                    promiscuity -= 1 # subtract 1 from each (self)
                    promiscuitys.append(promiscuity)
                df.loc[(df['Protein1_Domain'] == row['Protein1_Domain']) & (df['Protein2_Domain'] == row['Protein2_Domain']) & (df['location'] == row['location']), 'p1d_promiscuity'] = promiscuitys[0]
                df.loc[(df['Protein1_Domain'] == row['Protein1_Domain']) & (df['Protein2_Domain'] == row['Protein2_Domain']) & (df['location'] == row['location']), 'p2d_promiscuity'] = promiscuitys[1]
                df.loc[(df['Protein1_Domain'] == row['Protein1_Domain']) & (df['Protein2_Domain'] == row['Protein2_Domain']) & (df['location'] == row['location']), 'total_promiscuity'] = sum(promiscuitys)
                df.loc[(df['Protein1_Domain'] == row['Protein1_Domain']) & (df['Protein2_Domain'] == row['Protein2_Domain']) & (df['location'] == row['location']), 'max_promiscuity'] = max(promiscuitys)
            except:
                continue
    return df

# Example usage
import pandas as pd
df = pd.read_csv('/Users/poppy/Dropbox/all_dimer_interface_analysis_2025.06.05.csv')
df = find_interface_promiscuity(df)
df.to_csv('/Users/poppy/Dropbox/all_dimer_interface_analysis_06.05_promiscuity_added_2025.11.19_2.csv', index=False)

# import pandas as pd
# df = pd.read_csv('/Users/poppy/Dropbox/all_interface_analysis.csv')
# df = find_interface_promiscuity(df)
# # save as all interface analysis with promiscuity
# df.to_csv('/Users/poppy/Dropbox/all_interface_analysis_promiscuity.csv', index=False)
# from format_for_gephi import gephi_format
# criteria = {
#                 'rop': lambda x: x >= 2,
#                 'avg_pae': lambda x: x <= 15,
#                 'min_pae': lambda x: x <= 5,
#                 'size': lambda x: x >= 5
#             }
# gephi_format(input_csv='/Users/poppy/Dropbox/all_interface_analysis_promiscuity.csv', output_csv='gephi_input_promiscuity.csv', criteria=criteria, priority='avg_pae', priority_type='min')
