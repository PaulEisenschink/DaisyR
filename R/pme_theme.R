#' PME GGPLOT2 THEME
#'
#' A custom theme for ggplot2
#' @export 
pme_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "CMU Sans Serif"),
      plot.title = element_text(size = 16),
      plot.subtitle = element_text(size = 15),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = "#969696", linewidth = 0.7),
      panel.grid.minor = element_line(colour = "#969696")
    )
}