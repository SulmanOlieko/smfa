smfaStartupMessage <- function() {
  art_lines <- c(
    "                        .d888         ",
    "                       d88P\"          ",
    "                       888            ",
    ".d8888b  88888b.d88b.  888888 8888b.  ",
    "88K      888 \"888 \"88b 888       \"88b ",
    " Y8888b. 888  888  888 888   .d888888  ",
    "     X88 888  888  888 888   888  888 ",
    " 88888P' 888  888  888 888   \"Y888888 "
  )
  base_msg <- paste0(
    "\n\n* Please cite the 'smfa' package as:\n",
    "Owili, S. O. (2026). smfa: Stochastic Metafrontier Analysis. R package version 1.0.0.\n\n",
    "See also: citation(\"smfa\")\n\n",
    "* For any questions, suggestions, or comments on the 'smfa' package, you can contact the authors directly or visit:\n",
    "  https://github.com/SulmanOlieko/smfa/issues\n"
  )

  console_width <- getOption("width")
  art_width <- max(nchar(art_lines))

  if (console_width < art_width) {
    # 1. Truncate all lines to fit the console width
    trimmed_art <- substr(art_lines, 1, console_width - 5)

    # 2. Add the ASCII ellipsis strictly to the last row
    ellipsis <- "      88 88 88"
    e_len <- nchar(ellipsis)
    cut_point <- max(0, console_width - e_len)

    last_line_index <- length(trimmed_art)
    trimmed_art[last_line_index] <- paste0(
      substr(art_lines[last_line_index], 1, cut_point),
      ellipsis
    )

    # 3. Set the version alignment target
    align_width <- console_width
  } else {
    trimmed_art <- art_lines
    align_width <- art_width
  }

  art_msg <- paste(trimmed_art, collapse = "\n")

  # 4. Pad the version text so it hits the exact right edge
  version_msg <- paste0(
    "\n",
    sprintf(paste0("%", align_width, "s"), "version 1.0.0")
  )

  return(paste0(art_msg, version_msg, base_msg))
}

.onAttach <- function(lib, pkg) {
  # We still skip printing during the CRAN check to keep the test logs clean
  in_chk <- Sys.getenv("_R_CHECK_PACKAGE_NAME_") != ""

  if (!in_chk) {
    msg <- smfaStartupMessage()
  } else {
    msg <- paste0(
      "\n* Please cite the 'smfa' package as:\n",
      "Owili, SO. (2026). smfa: Metafrontier Analysis in R. R package version 1.0.0.\n\n",
      "See also: citation(\"smfa\")\n\n",
      "* For any questions, suggestions, or comments on the 'smfa' package, you can contact the authors directly or visit:\n",
      "  https://github.com/SulmanOlieko/smfa/issues\n"
    )
  }

  packageStartupMessage(msg)
  invisible()
}
