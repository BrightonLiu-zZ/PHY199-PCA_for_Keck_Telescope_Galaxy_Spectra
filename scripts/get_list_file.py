import os 
import astropy.io.fits as AIF 
def build_list_file_from_fits(fits_dir, out_filename): 
    outlier_filenames = { 
        "J001938.59+033605.50_DEIMOS_20230914_20230914.fits", 
        "J003216.47+181205.29_DEIMOS_20220701_20220701.fits", 
        "J004304.30-203105.07_DEIMOS_20151213_20151213.fits", 
        "J024756.32-021629.12_DEIMOS_20230914_20230914.fits", 
        "J121905.40+505448.59_DEIMOS_20230423_20230423.fits", 
        "J121912.36+505346.89_DEIMOS_20230423_20230423.fits", 
        "J130528.78+262800.89_DEIMOS_20220701_20220701.fits", 
        "J163503.40+152813.70_DEIMOS_20230423_20230423.fits", 
        "J232031.28+291554.16_DEIMOS_20230914_20230914.fits", 
        "J232104.37+291153.70_DEIMOS_20230914_20230914.fits", 
        "J001932.32+034144.66_DEIMOS_20230914_20230914.fits", 
        "J001937.28+033342.85_DEIMOS_20230914_20230914.fits", 
        "J001938.57+033450.41_DEIMOS_20230914_20230914.fits", 
        "J001938.71+033741.40_DEIMOS_20230914_20230914.fits", 
        "J001941.26+033509.66_DEIMOS_20230914_20230914.fits", 
        "J004325.10-203701.08_DEIMOS_20151213_20151213.fits", 
        "J080118.72+362621.17_DEIMOS_20230423_20230423.fits", 
        "SPAT0033-SLIT0032-MSC04_DEIMOS_20230914_20230914.fits", 
        "SPAT0170-SLIT0166-MSC04_DEIMOS_20230914_20230914.fits", 
        "SPAT0995-SLIT1008-MSC04_DEIMOS_20230914_20230914.fits", 
        "J003208.30+180728.64_DEIMOS_20220701_20220701.fits", 
        "J003214.98+180634.95_DEIMOS_20220701_20220701.fits", 
        "J121835.54+505205.58_DEIMOS_20230423_20230423.fits", 
        "J125721.98+365319.76_DEIMOS_20220701_20220701.fits", 
        "J163501.72+152817.49_DEIMOS_20230423_20230423.fits", 
        "J163502.23+152804.77_DEIMOS_20230908_20230908.fits", 
        "J213515.45+012550.39_DEIMOS_20150820_20150820.fits", 
        "J213527.39+012853.88_DEIMOS_20150820_20150820.fits", 
        "J232044.35+291345.95_DEIMOS_20230914_20230914.fits", 
        "SPAT1673-SLIT1621-MSC02_DEIMOS_20230914_20230914.fits", # ended 2nd iteration 
        "J001929.29+034243.72_DEIMOS_20230914_20230914.fits", 
        "J003209.39+180905.01_DEIMOS_20220701_20220701.fits", 
        "J080109.91+363030.94_DEIMOS_20230423_20230423.fits", 
        "J121821.78+505411.56_DEIMOS_20230423_20230423.fits", 
        "J121935.75+505639.90_DEIMOS_20230423_20230423.fits", 
        "J163535.00+153256.59_DEIMOS_20230423_20230423.fits", 
        "J213514.50+012257.19_DEIMOS_20150820_20150820.fits", # ended 3rd iteration
        "J121925.72+505752.43_DEIMOS_20230423_20230423.fits",  
        "J125657.64+364955.37_DEIMOS_20220701_20220701.fits",
        "J130540.92+262954.05_DEIMOS_20220701_20220701.fits",
        "J163456.91+152851.91_DEIMOS_20230423_20230423.fits",
        "J163527.94+153018.11_DEIMOS_20230423_20230423.fits",
        "SPAT0057-SLIT0049-MSC02_DEIMOS_20230914_20230914.fits",
        "SPAT0465-SLIT0426-MSC02_DEIMOS_20230914_20230914.fits",
        "SPAT0021-SLIT0026-MSC02_DEIMOS_20250223_20250223.fits", # started 1st PCA filtering
        "SPAT0131-SLIT0132-MSC02_DEIMOS_20250223_20250223.fits",
        "SPAT0410-SLIT0413-DET01_DEIMOS_20250224_20250224.fits",
        "SPAT0428-SLIT0413-DET01_DEIMOS_20250224_20250224.fits",
        "SPAT0060-SLIT0051-MSC02_DEIMOS_20250223_20250223.fits", # started 2nd PCA filtering
        "SPAT0990-SLIT0984-MSC02_DEIMOS_20250223_20250223.fits",
        "SPAT1328-SLIT1321-MSC03_DEIMOS_20150820_20150820.fits",
        "SPAT1471-SLIT1467-MSC03_DEIMOS_20150820_20150820.fits",
        "SPAT0052-SLIT0050-MSC04_DEIMOS_20250223_20250223.fits", # started 3rd PCA filtering
        "SPAT0245-SLIT0265-MSC02_DEIMOS_20250223_20250223.fits",
        "SPAT0362-SLIT0366-MSC02_DEIMOS_20250223_20250223.fits",
        "SPAT0502-SLIT0500-MSC03_DEIMOS_20250223_20250223.fits",
        "SPAT0514-SLIT0500-MSC03_DEIMOS_20250223_20250223.fits",
        "SPAT1001-SLIT0993-MSC04_DEIMOS_20150820_20150820.fits",
        "SPAT1001-SLIT0996-MSC04_DEIMOS_20250224_20250224.fits",
        "SPAT1044-SLIT1029-MSC04_DEIMOS_20250223_20250223.fits",
        "SPAT1555-SLIT1556-MSC03_DEIMOS_20150820_20150820.fits",
        "SPAT0809-SLIT0798-MSC04_DEIMOS_20250223_20250223.fits", # started 4th PCA filtering
        "SPAT1304-SLIT1302-MSC02_DEIMOS_20250223_20250223.fits",
        "SPAT1506-SLIT1511-MSC02_DEIMOS_20250223_20250223.fits",
        "SPAT0507-SLIT0513-MSC03_DEIMOS_20150820_20150820_modified.fits", # started 5th PCA filtering
        "SPAT1437-SLIT1432-MSC04_DEIMOS_20150820_20150820.fits",
        "SPAT0122-SLIT0096-MSC02_DEIMOS_20250223_20250223.fits", # started 6th PCA filtering
        "SPAT1207-SLIT1186-DET01_DEIMOS_20250223_20250223.fits",
        "SPAT1401-SLIT1399-MSC02_DEIMOS_20150820_20150820_modified.fits",
        "SPAT1580-SLIT1634-MSC04_DEIMOS_20230914_20230914.fits",
        "J213522.85+013128.12_DEIMOS_20150820_20150820.fits", # started 7th PCA filtering
        "SPAT1532-SLIT1525-MSC03_DEIMOS_20150820_20150820.fits",
        "J004329.89-204210.82_DEIMOS_20151213_20151213.fits" # started 8th PCA filtering
    } 
    with open(out_filename, 'w') as outfile: 
        for fname in os.listdir(fits_dir): 
            if fname.endswith('.fits'): 
                # Skip if this filename is in outlier set 
                if fname in outlier_filenames: 
                    continue 
                full_path = os.path.join(fits_dir, fname)
                with AIF.open(full_path) as hdul: 
                    hdr = hdul[0].header 
                    z = hdr['PYZEZ'] 
                    template = hdr['PYZETMPL'] 
                    lambda_min = hdr['PYZELMIN'] 
                    lambda_max = hdr['PYZELMAX'] 
                    
                    # The form of a line 
                    line = ( 
                        f"Spec1D File: {fname}, " 
                        f"Redshift: {z}, " 
                        f"Best Matching Template: {template}, " 
                        f"lambda_min: {lambda_min}, " 
                        f"lambda_max: {lambda_max}\n" 
                    ) 
                    # Skip if the line has 'NaN' 
                    if "NaN" in line: 
                        continue 
                    outfile.write(line) 

fits_dir = r"C:\PHY199\spec1d" 
out_filename = "C:\PHY199\list_file_9.txt"
build_list_file_from_fits(fits_dir, out_filename)
print("Done")