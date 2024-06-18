"""
Functions to calculate the ROP from PAE plots - not recommended as much less reliable than using raw PAE data
in combination with structure data (method for this provided in repeatability_from_pdb.py).  Additionally,
pae.png format varies - this method is adjusted for 3 output formats observed in the initial screen data, but
position of each of the 5 plots within the overall plot is hardcoded depending on the name of the file, meaning any
changes to the format of the pae.png files will require adjustments to the code.

Functions:
calculate_interprotein_ssim
process_alphafold_prediction_PAE_png_ROP
process_all_predictions_PAE_png_ROP
"""
from skimage.metrics import structural_similarity as ssim
import os
import cv2
import numpy as np
import math
import pandas as pd
from analysis_utility import find_rank_001_files, parse_structure_file, determine_chain_lengths

#assess consistency
def calculate_interprotein_ssim(image_path, len_A, len_B):
    """
    Calculates SSIM for the interprotein regions (AB and BA) in each grid of the image compared to
    the first grid and returns the individual similarity values for regions AB and BA of images 2, 3, 4, and 5,
    as well as the average of these values, and the averaged number of blue pixels in AB and BA for plot 1.
    """
    ssim_values = []
    average_ssim = 0
    proportional_int_size = 0

    # Load the image file and handle possible loading issues
    full_image = cv2.imread(image_path)
    if full_image is None:
        print(f"Error: Failed to load the image from {image_path}. Check the file path and permissions.")
        return ssim_values, average_ssim, proportional_int_size

    # Extract the blue channel (note: OpenCV loads images in BGR format)
    blue_channel = full_image[:, :, 0]

    # Define squares for the plots based on 'test' condition in image_path
    if 'test' in image_path:
        squares = [
            (395, 49, 300, 300),
            (876, 49, 300, 300),
            (1357, 49, 300, 300),
            (1838, 49, 300, 300),
            (2319, 49, 300, 300)
        ]
        square_size = 300
    if 'pae' in image_path:
        squares = [
            (60, 64, 306, 306),
            (541, 64, 306, 306),
            (1022, 64, 306, 306),
            (1503, 64, 306, 306),
            (1984, 64, 306, 306)
        ]
        square_size = 306
    else: # not 'test' and 'PAE'
        squares = [
            (389, 49, 305, 305),
            (870, 49, 305, 305),
            (1351, 49, 305, 305),
            (1832, 49, 305, 305),
            (2313, 49, 305, 305)
        ]
        square_size = 305

    # Calculate pixel positions for AB and BA regions based on protein lengths
    total_length = len_A + len_B
    if total_length == 0:
        print("Error: Total length of proteins cannot be zero.")
        return ssim_values, average_ssim, proportional_int_size

    len_A_scaled = int(square_size * len_A / total_length)
    len_B_scaled = square_size - len_A_scaled  # Ensure total length remains 305 pixels

    # Extract AB and BA regions for each plot and calculate blue pixel count for plot 1
    extracted_regions_AB = []
    extracted_regions_BA = []
    blue_pixel_count_AB = 0
    blue_pixel_count_BA = 0
    for i, (x, y, w, h) in enumerate(squares):
        AB_region = blue_channel[y:y + len_A_scaled, x + len_A_scaled:x + w]
        BA_region = blue_channel[y + len_A_scaled:y + h, x:x + len_A_scaled]
        if AB_region.size == 0 or BA_region.size == 0:
            print(f"Warning: Empty region detected in image at position {i+1}.")
            continue  # Skip this iteration

        extracted_regions_AB.append(AB_region)
        extracted_regions_BA.append(BA_region)
        
        if i == 0:  # Only for the first plot
            blue_threshold = 0  # Define threshold for blue intensity
            blue_pixel_count_AB = np.sum(AB_region > blue_threshold)
            blue_pixel_count_BA = np.sum(BA_region > blue_threshold)

    if not extracted_regions_AB or not extracted_regions_BA:
        print("Error: No valid AB or BA regions were extracted. Check the region extraction logic.")
        return ssim_values, average_ssim, proportional_int_size

    # Calculate SSIM for AB and BA regions separately
    try:
        ssim_values_AB = [ssim(extracted_regions_AB[0], region, data_range=max(region.max() - region.min(), 1))
                        for region in extracted_regions_AB[1:]]
        ssim_values_BA = [ssim(extracted_regions_BA[0], region, data_range=max(region.max() - region.min(), 1))
                        for region in extracted_regions_BA[1:]]
        ssim_values = [(ab + ba) / 2 for ab, ba in zip(ssim_values_AB, ssim_values_BA)]
        ssim_values = [0 if math.isnan(x) else x for x in ssim_values]  # Handle NaN values
        average_ssim = sum(ssim_values) / len(ssim_values) if ssim_values else 0
    except Exception as e:
        print(f"Error calculating SSIM: {e}")

    # Average the blue pixel counts for AB and BA regions of plot 1
    total_pixels = len_A_scaled * (len_B_scaled * 2)  # Since you're adding AB and BA regions
    proportional_int_size = (blue_pixel_count_AB + blue_pixel_count_BA) / total_pixels if total_pixels else 0
    
    return ssim_values, average_ssim, proportional_int_size

####################################################################################################

# for adding png data to compare to assessment from pae and structure data

def process_alphafold_prediction_PAE_png_ROP(folder_path, is_pdb=True):
    num_consistent = 0
    structure_file, _, _, PAE_png ,_ = find_rank_001_files(folder_path)
    # Parse structure file
    if structure_file:
        if PAE_png:
            structure_model = parse_structure_file(structure_file, is_pdb)
            chain_lengths = determine_chain_lengths(structure_model)
            ssim_values, _, proportional_int_size = calculate_interprotein_ssim(str(PAE_png), chain_lengths[0], chain_lengths[1])
            for value in ssim_values:
                if value > (1 - proportional_int_size):
                    num_consistent += 1
    else:
        num_consistent = None
    
    return num_consistent

def process_all_predictions_PAE_png_ROP(base_folder, input_file):
    # read in dataframe from file
    df = pd.read_csv(input_file)
    # add a new column to the dataframe 
    df['png_ROP'] = 0
    # Walk through the base folder containing all protein pair folders
    for root, dirs, _ in os.walk(base_folder):
        for dir in dirs:
            if '+' in dir:  # This is a protein pair folder
                protein1, protein2 = dir.split('+')
                #print(f"Processing {protein1} and {protein2}...")
                domain_folder_path = os.path.join(root, dir, 'Results')
                
                # Process each domain pair folder within the 'Results' folder
                if os.path.exists(domain_folder_path):
                    for domain_pair in os.listdir(domain_folder_path):
                        if domain_pair.startswith(protein1) and '+' in domain_pair:
                            domain_path = os.path.join(domain_folder_path, domain_pair)
                            # Extract domain information
                            protein1_domain, protein2_domain = domain_pair.split('+')
                            #print(f"Processing {protein1_domain} and {protein2_domain}...")
                            # Process the domain pair folder
                            num_consistent = process_alphafold_prediction_PAE_png_ROP(domain_path)

                            # Write results to dataframe
                            df.loc[(df['Protein1'] == protein1) & (df['Protein2'] == protein2) & (df['Protein1_Domain'] == protein1_domain) & (df['Protein2_Domain'] == protein2_domain), 'Num_Consistent_png'] = num_consistent
                            df.loc[(df['Protein1'] == protein2) & (df['Protein2'] == protein1) & (df['Protein1_Domain'] == protein2_domain) & (df['Protein2_Domain'] == protein1_domain), 'Num_Consistent_png'] = num_consistent
    return df

#use
base_folder = "../../../../../Dropbox/2022.10.20_Drosophila_Version_1"
df = process_all_predictions_PAE_png_ROP(base_folder, "data/alphafold_predictions_results.csv")
#save df as processed_with_png_ROP.csv
df.to_csv("data/alphafold_predictions_png_ROP.csv", index=False)
#summarise df
print(df)

####################################################################################################
# Example usage
#x, y, pae_png_file, fasta_file = find_rank_001_files("Sak_Sas6/Sak_D3+Sas6_D1")
#print(f"Rank 1 files: {x}, {y}, {pae_png_file}, {fasta_file}")
#len_A, len_B = determine_chain_lengths(fasta_file)
#print(f"Protein lengths: A = {len_A}, B = {len_B}")
#ssim_values, average_ssim, proportional_int_size = calculate_interprotein_ssim(pae_png_file, len_A, len_B)
#print(f"SSIM values for AB and BA regions: {ssim_values}. Interface size = {proportional_int_size}, so values greater than {1 - proportional_int_size} are considered consistent")