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

A python function, ImzMLParser from pyimzml package.
```python
from pyimzml.ImzMLParser import ImzMLParser
parser = ImzMLParser(path)
```
Each parser in this case provides:
- Coordinates - spatial layout
- getspectrum(i) - spectrum at pixel i
