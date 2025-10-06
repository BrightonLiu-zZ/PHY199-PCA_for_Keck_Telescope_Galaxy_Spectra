import os
import pyzeutil
import matplotlib.pyplot as plt
import numpy as np

def read_list_file(list_file):

    spec1d_files = []
    with open(list_file, 'r') as file: 
        for line in file:
            parts = line.split(',')  # Divide into parts (Spec1D File, Redshift, Best Matching Template)
            spec1d_file_name = parts[0].split(': ')[1].strip()  # Get the name of each spec1d file
            z = float(parts[1].split(': ')[1])  # Get the redshifts
            spec1d_files.append((spec1d_file_name, z))  # Add those names and their redshifts to the list
    return spec1d_files

def plot_lambda_max_distribution(folder_name, list_file, wave_min=5000, wave_max=9000, bin_size=50):

    spec1d_files = read_list_file(list_file)

    # Array to store lambda_max for all spectra
    lambda_max_list = []

    # Loop over all spectra
    for spec1d_file_name, z in spec1d_files:
        # Construct the full path
        path_to_file = os.path.join(folder_name, spec1d_file_name)

        # Create an instance of PypeItSpec to read the file
        myspec = pyzeutil.PypeItSpec(path_to_file)

        # Unsharp masking
        myspec.remove_largescale(kernel_size=451, hole_radius=15)

        # Deredshift the wavelength
        rest_wavelength = myspec.U.W / (1 + z)

        # Extract the maximum wavelength in the rest frame
        lambda_max = np.max(rest_wavelength)
        lambda_max_list.append(lambda_max)

    # Convert lambda_max_list to a numpy array
    lambda_max_array = np.array(lambda_max_list)

    # Create histogram bins for lambda_max
    bins = np.arange(wave_min, wave_max + bin_size, bin_size)
    hist, bin_edges = np.histogram(lambda_max_array, bins=bins)

    # Compute cumulative histogram
    cumulative_hist = np.cumsum(hist)

    # Plot the cumulative histogram
    plt.figure(figsize=(10, 6))
    plt.bar(bin_edges[:-1], cumulative_hist, width=bin_size, align='edge', edgecolor='black', alpha=0.7)
    plt.xlabel('Lambda Max (Rest Wavelength, Angstrom)')
    plt.ylabel('Cumulative Number of Spectra')
    plt.title('Cumulative Number of Spectra vs. Lambda Max')
    plt.xticks(np.arange(wave_min, wave_max + bin_size, 100))  # Adjust xticks to show every 100 Å
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()

folder_name = r"C:\PHY199\spec1d"  # Directory where spec1d files are stored
list_file = r'C:\PHY199\list_file.txt'  # Path to the list file containing spec1d files and redshifts
plot_lambda_max_distribution(folder_name, list_file, wave_min=4000, wave_max=6000, bin_size=100)
