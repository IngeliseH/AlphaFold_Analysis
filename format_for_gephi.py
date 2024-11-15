"""
Function to convert data from a CSV file to a format that can be used in Gephi for network visualization.
"""
import pandas as pd
import os

def gephi_format(input_csv, source='Protein1', target='Protein2', weight='rop', include=[], output_csv='gephi_input.csv', criteria=None):
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

    Returns
        - None, but saves the output CSV file in the same directory as the input file
    """
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(input_csv, sep=',')

    # Apply filtering criteria if provided
    if criteria:
        for column, condition in criteria.items():
            df = df[df[column].apply(condition)]

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
