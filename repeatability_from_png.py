from skimage.metrics import structural_similarity as ssim
import cv2
import numpy as np
import math
from analysis_utility import find_rank_001_files, extract_protein_lengths

#assess consistency
def calculate_interprotein_ssim(image_path, len_A, len_B):
    """
    Calculates SSIM for the interprotein regions (AB and BA) in each grid of the image compared to
    the first grid and returns the individual similarity values for regions AB and BA of images 2, 3, 4, and 5,
    as well as the average of these values, and the averaged number of blue pixels in AB and BA for plot 1.
    """
    # Define squares for the whole plots
    squares = [
        (389, 49, 305, 305),
        (870, 49, 305, 305),
        (1351, 49, 305, 305),
        (1832, 49, 305, 305),
        (2313, 49, 305, 305)
    ]
    
    # Load the image file and handle possible loading issues
    full_image = cv2.imread(image_path)
    if full_image is None:
        raise FileNotFoundError(f"Failed to load the image from {image_path}. Please check the file path and permissions.")

    # Extract the blue channel (note: OpenCV loads images in BGR format)
    blue_channel = full_image[:, :, 0]
    
    # Calculate pixel positions for AB and BA regions based on protein lengths
    total_length = len_A + len_B
    len_A_scaled = int(305 * len_A / total_length)
    len_B_scaled = 305 - len_A_scaled  # Ensure total length remains 305 pixels
    
    # Extract AB and BA regions for each plot and calculate blue pixel count for plot 1
    extracted_regions_AB = []
    extracted_regions_BA = []
    blue_pixel_count_AB = 0
    blue_pixel_count_BA = 0
    for i, (x, y, w, h) in enumerate(squares):
        AB_region = blue_channel[y:y + len_A_scaled, x + len_A_scaled:x + w]
        BA_region = blue_channel[y + len_A_scaled:y + h, x:x + len_A_scaled]
        if AB_region.size == 0 or BA_region.size == 0:
            continue
        extracted_regions_AB.append(AB_region)
        extracted_regions_BA.append(BA_region)
        
        if i == 0:  # Only for the first plot
            blue_threshold = 0  # Define threshold for blue intensity
            blue_pixel_count_AB = np.sum(AB_region > blue_threshold)
            blue_pixel_count_BA = np.sum(BA_region > blue_threshold)

    if not extracted_regions_AB or not extracted_regions_BA:
        return [], 0, 0  # Return early if no regions were valid

    # Calculate SSIM for AB and BA regions separately
    ssim_values_AB = [ssim(extracted_regions_AB[0], region, data_range=region.max() - region.min())
                      for region in extracted_regions_AB[1:]]
    ssim_values_BA = [ssim(extracted_regions_BA[0], region, data_range=region.max() - region.min())
                      for region in extracted_regions_BA[1:]]
    
    # Calculate the average of corresponding AB and BA values and add to ssim_values
    ssim_values = [(ab + ba) / 2 for ab, ba in zip(ssim_values_AB, ssim_values_BA)]
    ssim_values = [0 if math.isnan(x) else x for x in ssim_values]  # Handle NaN values

    # Calculate the overall average of ssim_values
    average_ssim = sum(ssim_values) / len(ssim_values) if ssim_values else 0

    # Average the blue pixel counts for AB and BA regions of plot 1
    average_blue_pixel_count = (blue_pixel_count_AB + blue_pixel_count_BA) / 2
    proportional_int_size = average_blue_pixel_count / (len_A_scaled * len_B_scaled)

    return ssim_values, average_ssim, proportional_int_size

x, y, pae_png_file, fasta_file = find_rank_001_files("Sak_Sas6/Sak_D3+Sas6_D1")
print(f"Rank 1 files: {x}, {y}, {pae_png_file}, {fasta_file}")
len_A, len_B = extract_protein_lengths(fasta_file)
print(f"Protein lengths: A = {len_A}, B = {len_B}")
ssim_values, average_ssim, proportional_int_size = calculate_interprotein_ssim(pae_png_file, len_A, len_B)
print(f"SSIM values for AB and BA regions: {ssim_values}. Interface size = {proportional_int_size}, so values greater than {1 - proportional_int_size} are considered consistent")