# `plot_spectra` Function

## Description
`plot_spectra` is an R function designed to import, normalize, and plot spectral data from multiple files within a folder.  
It supports multiple plotting modes (`individual`, `overlapped`, `stacked`), flexible axis customization, and optional features like vertical lines, shaded regions, and annotations.  
The function also provides **consistent and customizable color handling** for enhanced visual clarity.

---

## Requirements
- R packages: **tidyverse**, **RColorBrewer**

---

## Usage

```r
plot_spectra(
  folder = ".",                     # Folder containing spectral data files
  file_type = "csv",                # File extension of spectral data files
  header = TRUE,                    # Whether files have a header row
  sep = ",",                        # Field separator (e.g., "," or "\t")
  normalization = c("none", "simple", "min-max", "z-score"), # Normalization method
  x_limits = NULL,                  # c(min, max) for x-axis limits
  x_breaks = NULL,                  # Sequence for x-axis breaks
  x_reverse = FALSE,                # Reverse x-axis?
  y_trans = c("linear", "log10", "sqrt"), # Transformation for y-axis
  x_label = "Energy (keV)",         # X-axis label
  y_label = "Counts/1000 s",        # Y-axis label
  line_size = 0.5,                  # Line thickness
  palette = "black",                # Color palette: single color, "Dark2", or custom vector
  plot_mode = c("individual", "overlapped", "stacked"), # Plotting style (default: "individual")
  show_legend = FALSE,              # Display legend?
  vertical_lines = NULL,            # x positions for vertical dashed lines
  shaded_ROIs = NULL,               # List of c(xmin, xmax) for shaded regions
  annotations = NULL                # Data frame with: file, x, y, label for text annotations
)
```

---

## Details  
- **Input Data:** Reads all files with specified extension in the folder.  
- **Normalization:** Supports no normalization, simple scaling, min-max and z-score normalization.  
- **Plot Modes:**  
  - `individual`: Separate plot for each spectrum saved individually.
  - `overlapped`: All spectra plotted on the same axes with different colors.  
  - `stacked`: Each spectrum plotted in its own facet stacked vertically with free y-axis scales.    
- **Axis Customization:** Supports x-axis limits, breaks, reversing; y-axis transformations such as linear (default), log10, sqrt.
- **Color Palette:**
  - `palette = "black"` (default): All lines plotted in black (or one color).
  - `palette = "Dark2"`: Use the "Dark2" palette from RColorBrewer for multiple spectra.
  - `palette = c("red", "blue", ...)`: Custom vector of colors for files. 
- **Extras:** Can add vertical reference lines, shaded regions of interest, and text annotations linked to each spectrum.  
- **Output:**
  - Saves plots in an Outputs folder created automatically in the working directory.
  - `output_format` argument lets you specify the desired file format: "tiff" (default), "png", "pdf", etc.
  - Filenames include plot mode and current date for easy tracking.

---

## Example  

```r
plot_spectra(
  folder = "SpectraData",
  file_type = "csv",
  normalization = "min-max",
  x_limits = c(400, 900),
  x_breaks = seq(400, 900, 50),
  plot_mode = "overlapped",
  palette = "Dark2",
  show_legend = TRUE,
  vertical_lines = c(450, 700),
  shaded_ROIs = list(c(500, 550), c(600, 650)),
  annotations = data.frame(file = "sample1.csv", x = 480, y = 0.8, label = "Peak A"),
  output_format = "png"
)
```
