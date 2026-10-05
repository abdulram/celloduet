# Shared plotting theme and palettes for the composition-confounding analysis (sourced by every R script).
SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e4e3df"
PALS <- list(
  genotype = c(WT = "#2a78d6", `5xFAD` = "#eb6834"),
  group = c(WT_NTC = "#2a78d6", `5xFAD_NTC` = "#eb6834", `5xFAD_SPP1` = "#1baf7a"),
  sex = c(F = "#e87ba4", M = "#4a3aa7"),
  id_series = c(`04xx` = "#2a78d6", `6xx` = "#eb6834"),
  tech_profile = c(A_low_multimap = "#2a78d6", B_high_multimap = "#eb6834"),
  treatment = c(NTC = "#2a78d6", SPP1 = "#1baf7a")
)
theme_ca <- function() {
  ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = SURFACE, colour = NA),
    panel.background = ggplot2::element_rect(fill = SURFACE, colour = NA),
    panel.grid.major = ggplot2::element_line(colour = GRID, linewidth = 0.3),
    panel.grid.minor = ggplot2::element_blank(),
    text = ggplot2::element_text(colour = INK), axis.text = ggplot2::element_text(colour = INK2),
    plot.title = ggplot2::element_text(face = "bold", size = 11), strip.text = ggplot2::element_text(face = "bold"))
}
