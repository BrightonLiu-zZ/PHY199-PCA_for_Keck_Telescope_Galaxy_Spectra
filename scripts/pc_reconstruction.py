import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pyzeutil

def read_spec_info(list_file, spec_base_dir):
    spec_info = []
    with open(list_file, "r") as f:
        for line in f:
            if line.strip() == "":
                continue  # skip empty lines
            parts = line.split(',')
            file_path = None
            z_val = None
            for part in parts:
                if "Spec1D File:" in part:
                    file_path = part.split("Spec1D File:")[1].strip()
                if "Redshift:" in part:
                    z_val = float(part.split("Redshift:")[1].strip())
            if file_path is not None and z_val is not None:
                # If file_path is not absolute, join it with spec_base_dir
                if not os.path.isabs(file_path):
                    file_path = os.path.join(spec_base_dir, file_path)
                spec_info.append({"filename": file_path, "redshift": z_val})
    return spec_info

def plot_observed_spectrum(wave_obs, flux_obs, err_obs, label="Observed Spectrum"):
    # Use a thinner line for the observed spectrum (linewidth 0.5)
    if err_obs is not None:
        plt.fill_between(
            wave_obs,
            flux_obs - err_obs,
            flux_obs + err_obs,
            alpha=0.3,
            step="mid",
            color="gray",
            label="Observed Uncertainty"
        )
    plt.plot(wave_obs, flux_obs, lw=0.5, label=label)

def plot_reconstructed_spectrum(wave_obs, flux_recon, label="Reconstructed Spectrum"):
    plt.plot(wave_obs, flux_recon, lw=0.8, linestyle="--", label=label)

def plot_sky_and_galaxy_lines(z_obs, flux_min):
    ax = plt.gca()
    # Mark lines over the full current wavelength range (no fixed limits)
    pyzeutil.mark_lines(z=-1, whichplot=ax, ylabel=flux_min*0.8)
    pyzeutil.mark_lines(z=z_obs, whichplot=ax, ylabel=flux_min*0.9)

def load_observed_spectrum_from_pyze(fits_file):
    myspec = pyzeutil.PypeItSpec(fits_file)
    wave_obs = myspec.U.W
    flux_obs = myspec.U.I
    ivar_obs = myspec.U.IVAR
    err_obs = np.sqrt(1.0 / ivar_obs)
    err_obs[~np.isfinite(err_obs)] = 0.0
    return wave_obs, flux_obs, err_obs

def check_if_match():
    # Define file paths for flux matrix and list file
    flux_file = r"C:\PHY199\flux_matrix_9.txt"
    list_file = r"C:\PHY199\list_file_9.txt"
    
    # Base directory for the observed spectra
    spec_base_dir = r"C:\PHY199\spec1d"
    
    # Read the specification info (filenames and redshifts) from the list file
    spec_info_list = read_spec_info(list_file, spec_base_dir)
    
    # Create the rest-frame wavelength array
    wave_rest = np.arange(3683, 5300.1, 0.3)
    n_wave = len(wave_rest)
    
    # Load the flux matrix (each row is a rest-frame spectrum)
    df_flux = pd.read_csv(flux_file, sep=" ", header=None)
    X_rest = df_flux.values
    if X_rest.shape[1] != n_wave:
        print(f"Expected {n_wave} columns, got {X_rest.shape[1]}.")
        return

    # Perform PCA on the rest-frame flux matrix
    n_components = 10
    pca = PCA(n_components=n_components)
    X_pca = pca.fit_transform(X_rest)
    mean_spectrum = pca.mean_
    pc_components = pca.components_
    
    # Choose a sample index
    index = 2277
    if index >= len(X_rest):
        print(f"Sample index {index} is out of range.")
        return
    flux_rest = X_rest[index, :]
    recon_rest = mean_spectrum.copy()
    for j in range(n_components):
        recon_rest += X_pca[index, j] * pc_components[j]
    
    # Get redshift and observed file from the list (order assumed to match)
    z_obs = spec_info_list[index]["redshift"]
    observed_fits_file = spec_info_list[index]["filename"]
    
    # Redshift the PCA reconstruction from rest to observed frame
    wave_obs_from_pca = wave_rest * (1.0 + z_obs)
    recon_obs = recon_rest
    
    # Load the observed spectrum from the corresponding FITS file
    wave_obs_data, flux_obs_data, err_obs_data = load_observed_spectrum_from_pyze(observed_fits_file)
    
    # Plot the observed spectrum and PCA reconstruction
    plt.figure(figsize=(10, 5))
    plot_observed_spectrum(wave_obs_data, flux_obs_data, err_obs_data, label="Observed Spectrum")
    plot_reconstructed_spectrum(wave_obs_from_pca, recon_obs,
                                label=f"Reconstructed Spectrum with {n_components} PCs")
    
    # Determine a minimum flux for line labeling
    flux_min = min(np.min(flux_obs_data), np.min(recon_obs))
    
    # Overlay sky and galaxy lines using the current wavelength range
    plot_sky_and_galaxy_lines(z_obs, flux_min)
    
    plt.title(f"Spectrum Reconstruction (z = {z_obs:.3f})")
    plt.xlabel("Observed Wavelength (Angstrom)")
    plt.ylabel("Flux")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    check_if_match()
