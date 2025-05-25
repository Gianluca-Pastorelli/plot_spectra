library(tidyverse)
library(RColorBrewer) # Use the function display.brewer.all(colorblindFriendly = TRUE) to see the colour blind-friendly palettes

# Useful characters, if needed
deg <- "\u00B0"
Delta <- "\u0394"

# Function to plot spectra
plot_spectra <- function(
    folder = ".",
    file_type = "csv",
    header = TRUE,
    sep = ",",  # Use "\t" if it is tab
    normalization = c("none", "simple", "min-max", "z-score"),
    x_limits = NULL, # Numeric vector of length 2 specifying axis limits, e.g. c(min, max)
    x_breaks = NULL, # Numeric vector of break positions, e.g. seq(min, max, by = step)
    x_reverse = FALSE,
    y_trans = c("identity", "log10", "sqrt"), # Choose y-axis transformation: identity (default), log10, or sqrt
    x_label = "Energy (keV)", # For complex formatting, use expression(), e.g. expression(Wavenumber~(cm^{-1}))
    y_label = "Counts/1000 s", # For complex formatting, use expression(), e.g. expression(Delta*E["00"]^{"*"}~(a.u.))
    line_size = 0.5,
    plot_mode = c("overlapped", "stacked", "individual"),
    show_legend = FALSE,
    vertical_lines = NULL,  # A numeric vector of x positions where vertical dashed lines will be drawn
    shaded_ROIs = NULL, # A list of numeric vectors, each with two elements c(xmin, xmax), defining shaded rectangular regions along x
    annotations = NULL # A data frame with columns 'file', 'x', 'y', 'label'; adds text annotations at specified points in each spectrum
) {
  normalization <- match.arg(normalization)
  plot_mode <- match.arg(plot_mode)
  y_trans <- match.arg(y_trans)
  
  # Define the custom theme with dynamic legend control
  plot_theme <- theme_bw(base_family = "sans") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.ticks.x = element_line(color = "black"),
      axis.ticks.y = element_line(color = "black"),
      axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title.x = element_text(color = "black", size = 12),
      axis.title.y = element_text(color = "black", size = 12),
      legend.position = if (show_legend) "right" else "none"
    )
  
  # In case, make use of other options for the theme:
  # text=element_text(family="Arial"),
  # axis.line = element_line(colour = "black"),
  # panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
  # panel.background = element_blank(),
  # legend.key = element_blank(),
  # plot.margin = margin(0.2,0.5,0.2,0.2, "cm")
  
  # Read files
  files <- list.files(folder, pattern = paste0("\\.", file_type, "$"), full.names = TRUE)
  spectra_list <- lapply(files, function(file) {
    data <- read_delim(file, delim = sep, col_names = header, show_col_types = FALSE)
    
    # Rename first two columns to "x" and "y"
    if (ncol(data) >= 2) {
      colnames(data)[1:2] <- c("x", "y")
    } else {
      stop(paste("File", file, "must have at least two columns for x and y"))
    }
    
    # Ensure numeric
    data <- data %>%
      mutate(x = as.numeric(x),
             y = as.numeric(y))
    
    # Apply normalization
    data$y <- switch(normalization,
                     none = data$y,
                     simple = data$y / max(data$y, na.rm = TRUE),
                     "min-max" = (data$y - min(data$y, na.rm = TRUE)) / (max(data$y, na.rm = TRUE) - min(data$y, na.rm = TRUE)),
                     "z-score" = scale(data$y)[,1])
    data$file <- basename(file)
    return(data)
  })
  
  spectra <- bind_rows(spectra_list)
  
  if (!dir.exists("Outputs")) dir.create("Outputs")
  
  if (plot_mode == "individual") {
    for (file_name in unique(spectra$file)) {
      data_sub <- filter(spectra, file == !!file_name)
      p <- ggplot(data_sub, aes(x = x, y = y)) +
        geom_line(linewidth = line_size, color = "black") +
        labs(x = x_label, y = y_label, title = file_name) +
        plot_theme
      
      # x-axis scale
      if (x_reverse) {
        p <- p + scale_x_reverse(limits = x_limits, breaks = x_breaks, expand = expansion())
      } else {
        p <- p + scale_x_continuous(limits = x_limits, breaks = x_breaks, expand = expansion())
      }
      
      # y-axis scale
      if (y_trans != "identity") {
        p <- p + scale_y_continuous(trans = y_trans)
      }
      
      # Optional extras
      if (!is.null(vertical_lines)) {
        for (v in vertical_lines) p <- p + geom_vline(xintercept = v, linetype = "dashed", color = "grey30")
      }
      if (!is.null(shaded_ROIs)) {
        for (roi in shaded_ROIs) p <- p + annotate("rect", xmin = roi[1], xmax = roi[2], ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey55")
      }
      if (!is.null(annotations)) {
        ann_sub <- filter(annotations, file == !!file_name)
        if (nrow(ann_sub) > 0) {
          p <- p + geom_text(data = ann_sub, aes(x = x, y = y, label = label), inherit.aes = FALSE)
        }
      }
      
      ggsave(
        filename = paste0(tools::file_path_sans_ext(file_name), "_", Sys.Date(), ".tiff"),
        plot = p,
        device = "tiff",
        path = "Outputs",
        scale = 1,
        width = 15,
        height = 9.3,
        units = "cm",
        dpi = 300,
        limitsize = TRUE,
        bg = "white"
      )
    }
  } else {
    p <- ggplot(spectra, aes(x = x, y = y, color = file)) +
      geom_line(linewidth = line_size) +
      labs(x = x_label, y = y_label, color = "File") +
      plot_theme
    
    # x-axis scale
    if (x_reverse) {
      p <- p + scale_x_reverse(limits = x_limits, breaks = x_breaks, expand = expansion())
    } else {
      p <- p + scale_x_continuous(limits = x_limits, breaks = x_breaks, expand = expansion())
    }
    
    # y-axis scale
    if (y_trans != "identity") {
      p <- p + scale_y_continuous(trans = y_trans)
    }
    
    if (plot_mode == "overlapped") {
      p <- p + scale_color_brewer(palette = "Dark2")
    }
    if (plot_mode == "stacked") {
      p <- p +
        facet_wrap(~ file, ncol = 1, scales = "free_y") +
        theme(
          panel.border = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.25),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing = unit(0, "mm"), # remove spacing between facets
          strip.background = element_blank(), # remove gray title bars
          strip.text = element_blank() # remove text in the facets
        )
    }
    if (!is.null(vertical_lines)) {
      for (v in vertical_lines) p <- p + geom_vline(xintercept = v, linetype = "dashed", color = "grey30")
    }
    if (!is.null(shaded_ROIs)) {
      for (roi in shaded_ROIs) p <- p + annotate("rect", xmin = roi[1], xmax = roi[2], ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey55")
    }
    if (!is.null(annotations)) {
      p <- p + geom_text(data = annotations, aes(x = x, y = y, label = label), inherit.aes = FALSE)
    }
    
    ggsave(
      filename = paste0("Combined_Spectra_", Sys.Date(), ".tiff"),
      plot = p,
      device = "tiff",
      path = "Outputs",
      scale = 1,
      width = 15,
      height = 9.3,
      units = "cm",
      dpi = 300,
      limitsize = TRUE,
      bg = "white"
    )
  }
}
