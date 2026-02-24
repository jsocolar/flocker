if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the package root (where DESCRIPTION lives).")
}

suppressPackageStartupMessages({
  library(withr)
  library(knitr)
})

find_pkg_root <- function() {
  path <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION"))) return(path)
    parent <- normalizePath(file.path(path, ".."), winslash = "/", mustWork = TRUE)
    if (identical(parent, path)) stop("Could not find package root (no DESCRIPTION found).")
    path <- parent
  }
}

pkg_root <- find_pkg_root()

vigs <- list(
  c("vignettes/augmented_models.Rmd.orig", "vignettes/augmented_models.Rmd"),
  c("vignettes/flocker_tutorial.Rmd.orig", "vignettes/flocker_tutorial.Rmd"),
  c("vignettes/nonlinear_models.Rmd.orig", "vignettes/nonlinear_models.Rmd.orig"), # (if desired)
  c("vignettes/nonlinear_models.Rmd.orig", "vignettes/nonlinear_models.Rmd"),
  c("vignettes/articles/sbc.Rmd.orig", "vignettes/articles/sbc.Rmd")
)

with_dir(pkg_root, {
  # Ensure relative paths in chunks resolve from package root
  old_root <- knitr::opts_knit$get("root.dir")
  on.exit(knitr::opts_knit$set(root.dir = old_root), add = TRUE)
  knitr::opts_knit$set(root.dir = pkg_root)
  
  knit_one <- function(infile, outfile, fig_path = NULL) {
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    
    if (!is.null(fig_path)) {
      dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)
      old_fig <- knitr::opts_chunk$get("fig.path")
      on.exit(knitr::opts_chunk$set(fig.path = old_fig), add = TRUE)
      knitr::opts_chunk$set(fig.path = paste0(fig_path, "/"))
    }
    
    message("Knitting ", infile, " -> ", outfile)
    knitr::knit(infile, output = outfile, envir = globalenv())
  }
  
  knit_one(vigs[[1]][1], vigs[[1]][2])
  knit_one(vigs[[2]][1], vigs[[2]][2])
  knit_one(vigs[[4]][1], vigs[[4]][2])
  
  # SBC: force figures into man/figures/sbc_vignette
  knit_one(vigs[[5]][1], vigs[[5]][2], fig_path = "man/figures/sbc_vignette")
})

message("Done.")