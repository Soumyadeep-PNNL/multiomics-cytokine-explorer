# ============================================================
# build_cache.R
#
# Run this ONCE before launching the app.
# It reads all 25 xlsx packages, parses them, and saves a
# compressed .rds file so the app loads in ~2 seconds instead
# of 30-60 seconds.
#
# Usage (in RStudio or R console):
#   source("build_cache.R")
#
# Re-run any time you update the xlsx files.
# ============================================================

# Auto-set working directory to the App folder
local({
  script_dir <- tryCatch(
    normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)),
    error = function(e) ""
  )
  if (nzchar(script_dir) && file.exists(file.path(script_dir, "data_loader.R"))) {
    setwd(script_dir)
  }
  if (!file.exists("data_loader.R")) {
    stop("Cannot find data_loader.R. Run this script from the App folder.")
  }
})
source("data_loader.R")

cache_file <- "data_cache.rds"

cat("========================================\n")
cat("Building data cache — please wait...\n")
cat("(This only needs to run once)\n")
cat("========================================\n\n")

t0 <- proc.time()

DATA <- load_all_data()

elapsed <- round((proc.time() - t0)[["elapsed"]])
cat(sprintf("\nParsed %d packages in %d seconds.\n",
            length(DATA$packages), elapsed))

cat("Saving cache to", cache_file, "...\n")
saveRDS(DATA, cache_file, compress = "xz")   # xz = smaller file, fast to decompress

size_mb <- round(file.info(cache_file)$size / 1024^2, 1)
cat(sprintf("Cache saved: %s (%.1f MB)\n", cache_file, size_mb))
cat("\nYou can now launch the app — it will load almost instantly.\n")
cat("  shiny::runApp('app.R')\n")
