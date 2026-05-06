from pyimzml.ImzMLParser import ImzMLParser
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import minmax_scale
import pandas as pd
import h5py
import os
# import plotly.io as pio
# pio.renderers.default = "iframe"
from scipy.signal import find_peaks
pd.set_option('display.max_columns', None)  
pd.set_option('display.max_rows', None)  
# Path to imzML lung cancer
imzml_path_LC24 = r"../../../npc-rnd/live/RND_Projects/DESI/out/2025_03_17_VSK_LC24_NEG_25um_250_G3-2_A_Analyte 2_1/2025_03_17_VSK_LC24_NEG_25um_250_G3-2_A_Analyte 2_1_ROI_TIC_log.imzML/2025_03_17_VSK_LC24_NEG_25um_250_G3-2_A_Analyte 2_1_ROI_TIC_log.imzML"
parser_cancer_1 = ImzMLParser(imzml_path_LC24)

imzml_path_LC08 = r"../../../npc-rnd/live/RND_Projects/DESI/out/2025_03_24_VSK_LC08_slice4_NEG_25um_250_G3_A_Analyte 1_1/2025_03_24_VSK_LC08_slice4_NEG_25um_250_G3_A_Analyte 1_1_ROI_TIC_log.imzML/2025_03_24_VSK_LC08_slice4_NEG_25um_250_G3_A_Analyte 1_1_ROI_TIC_log.imzML"
parser_cancer_2 = ImzMLParser(imzml_path_LC08)

imzml_path_LC22 = r"../../../npc-rnd/live/RND_Projects/DESI/out/2025_03_28_VSK_LC22_section5_NEG_25um_250_G3_A_Analyte 1_1/2025_03_28_VSK_LC22_section5_NEG_25um_250_G3_A_Analyte 1_1_ROI_TIC_log.imzML/2025_03_28_VSK_LC22_section5_NEG_25um_250_G3_A_Analyte 1_1_ROI_TIC_log.imzML"
parser_cancer_3 = ImzMLParser(imzml_path_LC22)


print("Number of spectra (pixels):", len(parser_cancer_1.coordinates), len(parser_cancer_2.coordinates), len(parser_cancer_3.coordinates))



# Path to imzML healthy tissue
imzml_path_HT10 = r"../../../npc-rnd/live/RND_Projects/DESI/out/2025_03_27_VSK_HT10_section5_NEG_25um_250_G3_B_Analyte 1_1/2025_03_27_VSK_HT10_section5_NEG_25um_250_G3_B_Analyte 1_1_ROI_TIC_log.imzML/2025_03_27_VSK_HT10_section5_NEG_25um_250_G3_B_Analyte 1_1_ROI_TIC_log.imzML"
parser_healthy_1 = ImzMLParser(imzml_path_HT10)

imzml_path_HT06 = r"../../../npc-rnd/live/RND_Projects/DESI/out/2025_03_24_VSK_HT06_section8_NEG_25um_250_G3_A_Analyte 1_1/2025_03_24_VSK_HT06_section8_NEG_25um_250_G3_A_Analyte 1_1_ROI_TIC_log.imzML/2025_03_24_VSK_HT06_section8_NEG_25um_250_G3_A_Analyte 1_1_ROI_TIC_log.imzML"
parser_healthy_2 = ImzMLParser(imzml_path_HT06)

imzml_path_HT13 = r"../../../npc-rnd/live/RND_Projects/DESI/out/2025_03_24_VSK_HT13_section1_NEG_25um_250_G3_B_Analyte 1_1/2025_03_24_VSK_HT13_section1_NEG_25um_250_G3_B_Analyte 1_1_ROI_TIC_log.imzML/2025_03_24_VSK_HT13_section1_NEG_25um_250_G3_B_Analyte 1_1_ROI_TIC_log.imzML"
parser_healthy_3 = ImzMLParser(imzml_path_HT13)





sample_paths = {
    "LC24": imzml_path_LC24,
    "LC08": imzml_path_LC08,
    "LC22": imzml_path_LC22,
    "HT10": imzml_path_HT10,
    "HT06": imzml_path_HT06,
    "HT13": imzml_path_HT13 
}

# organise samples
cancer_samples = {
    "LC24": parser_cancer_1,
    "LC08": parser_cancer_2,
    "LC22": parser_cancer_3
}

healthy_samples = {
    "HT10": parser_healthy_1,
    "HT06": parser_healthy_2,
    "HT13": parser_healthy_3
}


def check_samples(samples):
    for name, parser in samples.items():
        mz, _ = parser.getspectrum(0)
        print(f"{name}")
        print(" pixels:", len(parser.coordinates))
        print(" mz bins:", len(mz))
        print("-"*30)

check_samples(cancer_samples)
check_samples(healthy_samples)
df


def get_common_mz(parsers, step=0.005):

    min_mz = []
    max_mz = []

    for parser in parsers:
        mz,_ = parser.getspectrum(0)
        min_mz.append(mz.min())
        max_mz.append(mz.max())

    global_min = max(min_mz)
    global_max = min(max_mz)

    common_mz = np.arange(global_min, global_max, step)

    print("Common mz range:", global_min, "-", global_max)
    print("Common bins:", len(common_mz))

    return common_mz


all_parsers = list(cancer_samples.values()) + list(healthy_samples.values())
common_mz = get_common_mz(all_parsers)


def align_imzml_to_common_grid(parser, common_mz):

    n_pixels = len(parser.coordinates)
    n_bins = len(common_mz)

    aligned = np.zeros((n_pixels, n_bins), dtype=np.float32)

    for i in range(n_pixels):

        mzs, intensities = parser.getspectrum(i)

        # interpolate onto common mz axis
        aligned[i, :] = np.interp(
            common_mz,
            mzs,
            intensities,
            left=0,
            right=0
        )

        if i % 1000 == 0:
            print(f"Aligned {i}/{n_pixels}")

    return aligned


aligned_cancer = {}

for name, parser in cancer_samples.items():

    print(f"\nAligning {name}")

    aligned_cancer[name] = align_imzml_to_common_grid(
        parser,
        common_mz
    )


aligned_healthy = {}

for name, parser in healthy_samples.items():

    print(f"\nAligning {name}")

    aligned_healthy[name] = align_imzml_to_common_grid(
        parser,
        common_mz
    )




with h5py.File("aligned_lung_roi_data.h5", "w") as f:

    f.create_dataset("common_mz", data=common_mz)

    cancer_group = f.create_group("cancer")
    healthy_group = f.create_group("healthy")

    # -------- Cancer samples --------
    for name, parser in cancer_samples.items():

        print(f"\nAligning {name}")

        aligned = align_imzml_to_common_grid(parser, common_mz)

        # Save aligned intensities
        cancer_group.create_dataset(
            name,
            data=aligned,
            compression="gzip"
        )

        # Save coordinates (convert to int if 1-indexed)
        coords = np.array(parser.coordinates)
        cancer_group.create_dataset(f"{name}_coords", data=coords.astype(np.int32))

        del aligned  # free memory

    # -------- Healthy samples --------
    for name, parser in healthy_samples.items():

        print(f"\nAligning {name}")

        aligned = align_imzml_to_common_grid(parser, common_mz)

        # Save aligned intensities
        healthy_group.create_dataset(
            name,
            data=aligned,
            compression="gzip"
        )

        # Save coordinates
        coords = np.array(parser.coordinates)
        healthy_group.create_dataset(f"{name}_coords", data=coords.astype(np.int32))

        del aligned

    # -------- Save paths as attributes --------
    for name, path in sample_paths.items():
        f.attrs[f"{name}_path"] = os.path.abspath(path)