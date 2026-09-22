packages <- c("ggplot2", "gridExtra", "cowplot", "viridis", "scales", "pracma", "gsl")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
