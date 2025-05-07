"""
Function to convert data from a CSV file to a format that can be used in Gephi for network visualization.
"""
import pandas as pd
import os

def gephi_format(input_csv, source='Protein1', target='Protein2', weight='rop', include=[], output_csv='gephi_input.csv', criteria=None, priority='avg_pae', priority_type='min'):
    """
    Convert a CSV file to a format that can be used in Gephi for network visualization.

    Parameters
        - input_csv (str): Path to the input CSV file
        - source (str): Name of the column that contains the source nodes
        - target (str): Name of the column that contains the target nodes
        - weight (str): Name of the column that contains the edge weights
        - include (list): List of additional columns to include in the output
        - output_csv (str): Name of the output CSV file
        - criteria (dict): Dictionary of criteria to filter the rows. Keys are column names and values are lambda functions for filtering.
        - priority (str): Name of the column to use for prioritization (NOT the weight. Currently only used if weight = rop)
        - priority_type (str): Type of prioritization to use. Either 'max' or 'min'. Currently only used if weight = rop

    Returns
        - None, but saves the output CSV file in the same directory as the input file
    """
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(input_csv, sep=',')

    # Apply filtering criteria if provided
    if criteria:
        for column, condition in criteria.items():
            # if an item has no value for the specified column, remove it
            df = df[df[column].notna()]
            # apply the condition to the column
            df = df[df[column].apply(condition)]
        if df.empty:
            print("No predictions meet the specified criteria.")
            return None
    
    # initialise output df
    desired_cols = [source, target, weight] + include
    output_df = pd.DataFrame(columns=desired_cols)
    if weight == 'rop':
        # find set of integer values of rop that meet criteria eg >= 2
        rop_values = df['rop'].unique()
        rop_values = [int(i) for i in rop_values]
        # iterate from max to min
        for i in range(max(rop_values), min(rop_values)-1, -1):
            # create new filtered df
            df_filtered = df[df['rop'] == i]
            # if new df not empty, go through each protein pair
            if not df_filtered.empty:
                for index, row in df_filtered.iterrows():
                    # check if protein pair is already in output_df
                    if not ((output_df[source] == row['Protein1']) & (output_df[target] == row['Protein2'])).any():
                        # if not, find all instances of protein pair in df_filtered
                        df_pair = df_filtered[(df_filtered['Protein1'] == row['Protein1']) & (df_filtered['Protein2'] == row['Protein2'])]
                        # check value of priority column for each and add entry with highest val if priority_type is max or lowest if min
                        if priority_type == 'max':
                            max_val = df_pair[priority].max()
                            # remove unwanted columns
                            df_pair = df_pair[df_pair[priority] == max_val].iloc[[0]][desired_cols]
                            output_df = pd.concat([output_df, df_pair[df_pair[priority] == max_val].iloc[[0]]])
                        elif priority_type == 'min':
                            min_val = df_pair[priority].min()
                            # remove unwanted columns
                            df_pair = df_pair[df_pair[priority] == min_val].iloc[[0]][desired_cols]
                            output_df = pd.concat([output_df, df_pair[df_pair[priority] == min_val].iloc[[0]]])
        # rename columns
        output_df = output_df.rename(columns={source: 'Source', target: 'Target', weight: 'Weight'})

    else:
        # Prepare the output DataFrame with the required columns
        columns_to_include = [source, target, weight] + include
        output_df = df[columns_to_include].rename(columns={source: 'Source', target: 'Target', weight: 'Weight'})

    # Get the directory of the input file
    output_dir = os.path.dirname(input_csv)

    # Construct output path using the same folder and specified output filename
    output_path = os.path.join(output_dir, output_csv)

    # Save to CSV
    output_df.to_csv(output_path, index=False)

############################################################################################
# Example usage
# criteria = {
#     'rop': lambda x: x >= 2,
#     'avg_pae': lambda x: x <= 15,
#     'min_pae': lambda x: x <= 5,
#     'size': lambda x: x >= 5
# }

# gephi_format(
#     input_csv='/Users/poppy/Dropbox/PCM/PCM_interface_analysis.csv',
#     source='Protein1',
#     target='Protein2',
#     weight='rop',
#     include=['avg_pae', 'min_pae', 'size'],
#     output_csv='gephi_interfaces_PCM.csv',
#     criteria=criteria
# )
