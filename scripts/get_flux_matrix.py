#!/usr/bin/env python
import os
import numpy as np
import pyzeutil
from scipy.interpolate import griddata

def read_list_file(list_file):
    """
    Reads each line of 'list_file' and returns a list of (spec1d_file_name, z)
    only for lines containing 'spDR2-023.fit'.
    """
    spec1d_files = []
    with open(list_file, 'r') as f:
        for line in f:
            #if "spDR2-023.fit" in line:
                parts = line.split(',')  # e.g. ["Spec1D File: X", " Redshift: Y", " Best Matching Template: ..."]
                spec1d_file_name = parts[0].split(': ')[1].strip()
                z_str = parts[1].split(': ')[1].strip()
                z = float(z_str)
                spec1d_files.append((spec1d_file_name, z))
    return spec1d_files

def create_weighted_flux_matrix(rmj_folder, list_file, output_file):
    """
    1) Reads all valid spec1d files from list_file,
    2) Interpolates them onto a common wavelength grid,
    3) Computes a weighted mean flux for each wavelength bin (using IVAR),
    4) Replaces any NaNs in each spectrum with that bin's weighted mean,
    5) Writes the final flux matrix to output_file (space-separated).
    """
    spec1d_files = read_list_file(list_file)

    # Adjust common wavelength grid as needed:
    # For example, one point per 1A from 3500 to 6300 => 2801 points
    common_wavelength_grid = np.linspace(3683, 5300.1, 5391)

    flux_list = []
    IVAR_list = []

    for spec1d_file_name, z in spec1d_files:
        full_path = os.path.join(rmj_folder, spec1d_file_name)

        # Read with PypeItSpec
        myspec = pyzeutil.PypeItSpec(full_path)

        # Unsharp mask
        myspec.remove_largescale(kernel_size=451, hole_radius=15)

        # Deredshift wavelength
        rest_wavelength = myspec.U.W / (1 + z)

        # Interpolate flux & IVAR onto common grid
        flux_on_common = griddata(
            points=rest_wavelength,
            values=myspec.U.I,
            xi=common_wavelength_grid,
            fill_value=np.nan
        )
        ivar_on_common = griddata(
            points=rest_wavelength,
            values=myspec.U.IVAR,
            xi=common_wavelength_grid,
            fill_value=np.nan
        )

        flux_list.append(flux_on_common)
        IVAR_list.append(ivar_on_common)

    # Convert to arrays => shape (n_spectra, n_wavelengths)
    flux_array = np.array(flux_list)
    IVAR_array = np.array(IVAR_list)
    n_spectra, n_wave = flux_array.shape
    print(f"Flux array shape = {flux_array.shape}")

    # 5) Compute weighted mean across all spectra for each wavelength
    # Weighted mean = sum(IVAR * flux) / sum(IVAR), ignoring NaN
    weighted_flux_sum = np.nansum(flux_array * IVAR_array, axis=0)
    IVAR_sum = np.nansum(IVAR_array, axis=0)

    weighted_mean_flux = np.divide(weighted_flux_sum, IVAR_sum, out=np.full_like(weighted_flux_sum, np.nan), where=IVAR_sum!=0)

    # 6) Replace any NaNs in flux_array with the weighted mean for that wavelength
    for i in range(n_spectra):
        nan_mask = np.isnan(flux_array[i])
        flux_array[i, nan_mask] = weighted_mean_flux[nan_mask]

    # 7) Write the flux array to a .txt file => shape (n_spectra, n_wavelengths)
    # Each row is one spectrum, columns are flux at each wavelength bin
    with open(output_file, 'w') as fout:
        for i in range(n_spectra):
            row_str = " ".join(f"{val:.6g}" for val in flux_array[i])
            fout.write(row_str + "\n")

    print(f"Done! Wrote flux matrix to {output_file}")

def main():
    rmj_folder = r"C:\PHY199\spec1d" # Adjust to your directory
    list_file  = r"C:\PHY199\list_file_9.txt"
    out_file   = r"C:\PHY199\flux_matrix_9.txt"

    create_weighted_flux_matrix(rmj_folder, list_file, out_file)

if __name__ == "__main__":
    main()
