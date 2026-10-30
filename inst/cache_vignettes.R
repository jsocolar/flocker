# Source this file, then explicitly choose which caches to refresh:
# source("inst/cache_vignettes.R")
# cache_vignettes("nmix_benchmark", render_html = TRUE)
# cache_vignettes(c("flocker_tutorial", "nonlinear_models"))
# cache_vignettes() # all registered vignettes, including long-running SBC
#
# Names: augmented_models, flocker_tutorial, nonlinear_models, sbc, sbc_multi,
#        sbc_aug, sbc_twolevel, nmix_benchmark.
# HTML is rendered from the newly cached Rmd, not from the executable source.
# Returns cached Rmd paths and, when requested, rendered HTML paths invisibly.
cache_vignettes <- function(vignettes = NULL, render_html = FALSE, pkg_root = getwd()) {
  sources <- c(
    augmented_models = "vignettes/augmented_models.Rmd.orig",
    flocker_tutorial = "vignettes/flocker_tutorial.Rmd.orig",
    nonlinear_models = "vignettes/nonlinear_models.Rmd.orig",
    sbc = "vignettes/articles/sbc.Rmd.orig",
    sbc_multi = "vignettes/articles/sbc_multi.Rmd.orig",
    sbc_aug = "vignettes/articles/sbc_aug.Rmd.orig",
    sbc_twolevel = "vignettes/articles/sbc_twolevel.Rmd.orig",
    nmix_benchmark = "vignettes/articles/nmix_benchmark.Rmd.orig"
  )
  if (is.null(vignettes)) vignettes <- names(sources)
  if (!is.character(vignettes) || anyNA(vignettes) ||
      any(!vignettes %in% names(sources))) {
    stop("Select vignette names from: ", paste(names(sources), collapse = ", "), call. = FALSE)
  }
  if (!is.logical(render_html) || length(render_html) != 1L || is.na(render_html)) {
    stop("render_html must be TRUE or FALSE.", call. = FALSE)
  }
  vignettes <- unique(vignettes)
  if (!length(vignettes)) return(invisible(list(rmd = character(), html = character())))

  pkg_root <- normalizePath(pkg_root, winslash = "/", mustWork = TRUE)
  while (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    parent <- dirname(pkg_root)
    if (identical(parent, pkg_root)) stop("Could not find package root (no DESCRIPTION found).")
    pkg_root <- parent
  }
  description <- read.dcf(file.path(pkg_root, "DESCRIPTION"))
  if (description[1, "Package"] != "flocker") stop("pkg_root must identify the flocker checkout.")
  source_version <- package_version(description[1, "Version"])
  installed_version <- tryCatch(utils::packageVersion("flocker"), error = function(e) NULL)
  if (is.null(installed_version) || installed_version != source_version) {
    stop("Install this source tree and restart R before caching vignettes: the installed flocker version must match DESCRIPTION.",
         call. = FALSE)
  }
  inputs <- file.path(pkg_root, sources[vignettes])
  if (any(!file.exists(inputs))) {
    stop("Missing vignette sources: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
  }

  knit_one <- function(name, infile, outfile) {
    # Restore options even if a source changes them or knitting fails.
    old_knit <- knitr::opts_knit$get()
    old_chunk <- knitr::opts_chunk$get()
    on.exit(knitr::opts_knit$restore(old_knit), add = TRUE)
    on.exit(knitr::opts_chunk$restore(old_chunk), add = TRUE)
    is_sbc <- startsWith(name, "sbc")
    root <- if (name == "nmix_benchmark") dirname(infile) else pkg_root
    knitr::opts_knit$set(root.dir = root)
    if (is_sbc) {
      dir.create(file.path(pkg_root, "man/figures/sbc_vignette"),
                 recursive = TRUE, showWarnings = FALSE)
      knitr::opts_chunk$set(fig.path = "man/figures/sbc_vignette/")
    }
    env <- new.env(parent = globalenv())
    env$params <- rmarkdown::yaml_front_matter(infile)$params
    message("Knitting ", infile, " -> ", outfile)
    # root.dir controls chunk evaluation, but knitr can write relative figure
    # paths from the calling directory. Align both directories so article
    # figures land beside the cached Rmd and its relative links resolve.
    withr::with_dir(root, knitr::knit(infile, output = outfile, envir = env))
    if (is_sbc) {
      # Cached articles resolve image links relative to vignettes/articles.
      lines <- readLines(outfile, warn = FALSE)
      lines <- gsub('src="man/figures/sbc_vignette/',
                    'src="../../man/figures/sbc_vignette/', lines, fixed = TRUE)
      lines <- gsub('](man/figures/sbc_vignette/',
                    '](../../man/figures/sbc_vignette/', lines, fixed = TRUE)
      writeLines(lines, outfile)
    }
    invisible(outfile)
  }

  withr::with_dir(pkg_root, {
    outputs <- setNames(sub("[.]orig$", "", inputs), vignettes)
    html <- character()
    for (name in vignettes) {
      outfile <- outputs[[name]]
      knit_one(name, file.path(pkg_root, sources[[name]]), outfile)
      if (render_html) {
        # Honor the cached document's declared HTML format and render next to it.
        html[[name]] <- rmarkdown::render(outfile, envir = new.env(parent = globalenv()))
      }
    }
    invisible(list(rmd = outputs, html = html))
  })
}
