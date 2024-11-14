"""
Function to convert data from a CSV file to a format that can be used in Gephi for network visualization.
"""
import pandas as pd
import os

def gephi_format(input_csv, source='Protein1', target='Protein2', weight='rop', include=[], output_csv='gephi_input.csv', update_input=False):
    """
    Convert a CSV file to a format that can be used in Gephi for network visualization.

    Parameters
        - input_csv (str): Path to the input CSV file
        - source (str): Name of the column that contains the source nodes
        - target (str): Name of the column that contains the target nodes
        - weight (str): Name of the column that contains the edge weights
        - include (list): List of additional columns to include in the output
        - output_csv (str): Name of the output CSV file
        - update_input (bool): Whether to update the input CSV file with the new columns

    Returns
        - None, but saves the output CSV file in the same directory as the input file
    """
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(input_csv, sep=',')

    # Update input file with new column
    if update_input:
        df.to_csv(input_csv, index=False)

    # Prepare the output DataFrame with the required columns
    columns_to_include = [source, target, weight] + include
    output_df = df[columns_to_include].rename(columns={source: 'Source', target: 'Target', weight: 'Weight'})

    # Get the directory of the input file
    output_dir = os.path.dirname(input_csv)

    # Construct output path using the same folder and specified output filename
    output_path = os.path.join(output_dir, output_csv)

    # Save to CSV
    output_df.to_csv(output_path, index=False)
