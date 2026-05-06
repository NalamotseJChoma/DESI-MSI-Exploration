# DESI-MSI-Exploration  
#### DESI Mass Spectrometry Imaging Pipeline

#### Overview
This project builds a data processing pipeline for DESI data, focusing on lung cancer vs healthy tissue

#### Dataset
We work with imzML files from multiple tissue samples:
- Cancer samples
  * LC24
  * LC08
  * LC22
- Healthy samples
  * HT10
  * HT06
  * HT13
Each file contains:
- Pixel coordinates
- Mass spectra (m/z vs intensity)

A python function, ImzMLParser from pyimzml package was used to read and process the data.
```python
from pyimzml.ImzMLParser import ImzMLParser
parser = ImzMLParser(path)
```
Each parser in this case provides:
- Coordinates - spatial layout
- getspectrum(i) - spectrum at pixel i

```
LC24
 pixels: 29346
 mz bins: 170955
------------------------------
LC08
 pixels: 9576
 mz bins: 135444
------------------------------
LC22
 pixels: 17271
 mz bins: 135194
------------------------------
HT10
 pixels: 15251
 mz bins: 134855
------------------------------
HT06
 pixels: 8181
 mz bins: 165661
------------------------------
HT13
 pixels: 8181
 mz bins: 135372
------------------------------
```

Given that the data contains different number of mz bins and pixels. Aligning it directly becomes a problem. It is for that reason that we then created a common grid, by looking at common ranges where all the samples overlap, interpolating this using a median peak spacing of 0.005:

```
LC24 median spacing: 0.004015173434979147
LC08 median spacing: 0.006097772624627851
LC22 median spacing: 0.006335430099909445
HT10 median spacing: 0.0062578307101262
HT06 median spacing: 0.005196135146434244
HT13 median spacing: 0.006339186375612371
```

```python
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

```
This returns a list of common mz. Which we then use to interpolate the intensities using the function below. This function returns aligned data:

```python
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


````
