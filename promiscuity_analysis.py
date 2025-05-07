def expand_loc_range(loc_range):
    """
    expand a location range string (eg '1-5,7,9-10') into a list of integers (eg [1, 2, 3, 4, 5, 7, 9, 10])
    """
    locs = loc_range.split(',')
    expanded_locs = []
    for loc in locs:
        if '-' in loc:
            start, end = loc.split('-')
            expanded_locs.extend(range(int(start), int(end)+1))
        else:
            expanded_locs.append(int(loc))
    return expanded_locs

def find_int_overlap(loc, loc_list):
    """
    for a list of interface locations (loc_list), find how many overlap with at least 75% of a particular
    interface (loc)
    """
    total = 0
    # for each interface location in loc_list, check what percent of the current location is in the interface location
    for elem in loc_list:
        overlap = len(set(loc).intersection(set(elem))) / len(loc)
        if overlap > 0.75:
            total += 1
    return total

def find_interface_promiscuity(df):
    """
    find promiscuity of each domain pair in the dataframe by finding the number of other domain pairs
    that overlap with the interface location of the current domain pair
    """
    # build dict of domain names and all interface locations found involving that domain
    from collections import defaultdict
    domain_int_locs = defaultdict(list)
    for _, row in df.iterrows():
        if row['location']:
            from ast import literal_eval
            try:
                locations = literal_eval(row['location'])
                loc1 = expand_loc_range(locations['Protein1'])
                loc2 = expand_loc_range(locations['Protein2'])
                domain_int_locs[row['Protein1_Domain']].append(loc1)
                domain_int_locs[row['Protein2_Domain']].append(loc2)
            except:
                continue
    # iterate through again and find promiscuity for each domain
    df['p1d_promiscuity'] = None
    df['p2d_promiscuity'] = None
    df['total_promiscuity'] = None
    for _, row in df.iterrows():
        if row['location']:
            try:
                locations = literal_eval(row['location'])
                loc1 = expand_loc_range(locations['Protein1'])
                loc2 = expand_loc_range(locations['Protein2'])
                p1d_promiscuity = find_int_overlap(loc1, domain_int_locs[row['Protein1_Domain']])
                p2d_promiscuity = find_int_overlap(loc2, domain_int_locs[row['Protein2_Domain']])
                # subtract 1 from each (self)
                p1d_promiscuity -= 1
                p2d_promiscuity -= 1
                df.loc[(df['Protein1_Domain'] == row['Protein1_Domain']) & (df['Protein2_Domain'] == row['Protein2_Domain']) & (df['location'] == row['location']), 'p1d_promiscuity'] = p1d_promiscuity
                df.loc[(df['Protein1_Domain'] == row['Protein1_Domain']) & (df['Protein2_Domain'] == row['Protein2_Domain']) & (df['location'] == row['location']), 'p2d_promiscuity'] = p2d_promiscuity
                df.loc[(df['Protein1_Domain'] == row['Protein1_Domain']) & (df['Protein2_Domain'] == row['Protein2_Domain']) & (df['location'] == row['location']), 'total_promiscuity'] = p1d_promiscuity + p2d_promiscuity
            except:
                continue
    return df

# Example usage
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
