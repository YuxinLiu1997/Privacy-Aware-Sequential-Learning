# Run with Rscript --vanilla run_r.R 1 5 7 J1 [--quick]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
root <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
quick <- "--quick" %in% args
targets <- args[args != "--quick"]
scripts <- c("1" = "figure_01_binary_accuracy.R",
             "5" = "figure_05_smooth_randomized_response.R",
             "7" = "figure_07_heterogeneous_learning.R",
             "J1" = "figure_J1_llr_evolution.R")
if (!length(targets) || any(!targets %in% c(names(scripts), "all")))
    stop("Usage: Rscript --vanilla run_r.R {1 5 7 J1|all} [--quick]")
if ("all" %in% targets) targets <- names(scripts)
output <- file.path(root, "outputs", if (quick) "quick" else "full", "r")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
setwd(output)
for (target in targets) {
    lines <- readLines(file.path(root, "r", scripts[[target]]), encoding = "UTF-8", warn = FALSE)
    if (quick && target == "7") {
        lines <- sub("N <- 10000; sigma <- 1; seed <- 6", "N <- 100; sigma <- 1; seed <- 6", lines, fixed = TRUE)
        lines <- sub("nsim <- 5", "nsim <- 2", lines, fixed = TRUE)
    }
    message("Figure ", target, if (quick) " (quick output directory)" else "")
    path <- paste0("figure_", target, ".pdf")
    pdf(path, width = if (target %in% c("5", "J1")) 12 else 9, height = 6,
        onefile = TRUE)
    tryCatch(eval(parse(text = lines), envir = new.env(parent = globalenv())),
             finally = dev.off())
}
message("Output directory: ", output)
