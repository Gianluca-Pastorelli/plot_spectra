# plot_spectra Function

## Description  
`plot_spectra` is an R function designed to import, normalize, and plot spectral data from multiple files in a folder. It supports multiple plotting modes (overlapped, stacked, or individual plots), flexible axis customization, and optional features like vertical lines, shaded regions, and annotations.

## Requirements  
- R packages: **tidyverse**, **RColorBrewer**

## Usage  

```r
plot_spectra(
  folder = ".",                   # Folder containing spectral data files
  file_type = "csv",              # File extension of spectral data files
  header = TRUE,                  # Whether files have a header row
  sep = ",",                     # Field separator ("," or "\t")
  normalization = c("none", "simple", "min-max", "z-score"), # Normalization method for intensities
  x_limits = NULL,                # Numeric vector c(min, max) for x-axis limits
  x_breaks = NULL,                # Numeric vector for x-axis tick positions
  x_reverse = FALSE,              # Whether to reverse the x-axis
  y_trans = c("identity", "log10", "sqrt"), # Transformation for y-axis scale
  x_label = "Wavelength (nm)",   # Label for x-axis (can be expression())
  y_label = "Intensity (a.u.)",  # Label for y-axis (can be expression())
  line_size = 0.5,                # Thickness of plot lines
  plot_mode = c("overlapped", "stacked", "individual"), # Plotting style
  show_legend = FALSE,            # Show legend or not
  vertical_lines = NULL,          # Numeric vector of x positions for vertical dashed lines
  shaded_ROIs = NULL,             # List of numeric vectors c(xmin, xmax) for shaded regions
  annotations = NULL              # Data frame with columns: file, x, y, label for text annotations
)

## Details  
- **Input Data:** Reads all files with specified extension in the folder.  
- **Normalization:** Supports no normalization, simple scaling, min-max, and z-score normalization.  
- **Plot Modes:**  
  - `overlapped`: All spectra plotted on the same axes with different colors.  
  - `stacked`: Each spectrum plotted in its own facet stacked vertically with free y-axis scales.  
  - `individual`: Separate plot for each spectrum saved individually.  
- **Axis Customization:** Supports x-axis limits, breaks, reversing; y-axis transformations (identity, log10, sqrt).  
- **Extras:** Can add vertical reference lines, shaded regions of interest, and text annotations linked to each spectrum.  
- **Output:** Saves plots as TIFF files in an `Outputs` folder created in the working directory.

## Example  

```r
plot_spectra(
  folder = "SpectraData",
  file_type = "csv",
  normalization = "min-max",
  x_limits = c(400, 900),
  x_breaks = seq(400, 900, 50),
  plot_mode = "overlapped",
  show_legend = TRUE,
  vertical_lines = c(450, 700),
  shaded_ROIs = list(c(500, 550), c(600, 650)),
  annotations = data.frame(file = "sample1.csv", x = 480, y = 0.8, label = "Peak A")
)
