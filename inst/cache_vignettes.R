if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the package root (where DESCRIPTION lives).")
}

suppressPackageStartupMessages({
  library(withr)
  library(knitr)
  library(rmarkdown)
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
  c("vignettes/nonlinear_models.Rmd.orig", "vignettes/nonlinear_models.Rmd"),
  c("vignettes/articles/sbc.Rmd.orig", "vignettes/articles/sbc.Rmd"),
  c("vignettes/articles/sbc_multi.Rmd.orig", "vignettes/articles/sbc_multi.Rmd"),
  c("vignettes/articles/sbc_aug.Rmd.orig", "vignettes/articles/sbc_aug.Rmd")
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
    
    env <- new.env(parent = globalenv())
    env$params <- rmarkdown::yaml_front_matter(infile)$params
    
    knitr::knit(infile, output = outfile, envir = env)
  }
  
  for(i in 1:3){
    knit_one(vigs[[i]][[1]], vigs[[i]][[2]])
  }
  for(i in 4:6){
    # SBC: force figures into man/figures/sbc_vignette
    knit_one(vigs[[i]][[1]], vigs[[i]][[2]], fig_path = "man/figures/sbc_vignette")
  }

})

# ---- Post-process cached articles to fix figure paths (Pandoc resolves relative
# paths from the input document directory, so vignettes/articles/* need ../man/...) ----

fix_article_fig_paths <- function(rmd_path,
                                  from = "man/figures/sbc_vignette/",
                                  to   = "../../man/figures/sbc_vignette/") {
  # Only makes changes if the file exists and actually contains the "from" string,
  # and avoids double-prepending if it already has "../".
  stopifnot(length(rmd_path) == 1)
  
  if (!file.exists(rmd_path)) {
    warning("Skipping missing file: ", rmd_path)
    return(invisible(FALSE))
  }
  
  lines <- readLines(rmd_path, warn = FALSE)
  
  # 1) HTML: <img src="man/figures/...">
  # Replace src="man/figures/..." but NOT src="../man/figures/..."
  lines2 <- gsub(
    pattern = paste0('src="', from),
    replacement = paste0('src="', to),
    x = lines,
    fixed = TRUE
  )
  
  # 2) Markdown: ![](man/figures/...)
  # Replace (man/figures/...) but NOT (../man/figures/...)
  lines2 <- gsub(
    pattern = paste0("](", from),
    replacement = paste0("](", to),
    x = lines2,
    fixed = TRUE
  )
  
  changed <- !identical(lines, lines2)
  if (changed) {
    writeLines(lines2, rmd_path)
    message("Rewrote figure paths in ", rmd_path)
  } else {
    message("No figure paths to rewrite in ", rmd_path)
  }
  
  invisible(changed)
}

# Apply to the cached article outputs (the .Rmd files, not the .orig)
article_outfiles <- vapply(vigs[4:6], `[[`, character(1), 2)

with_dir(pkg_root, {
  for (f in article_outfiles) fix_article_fig_paths(f)
})

message("Done.")