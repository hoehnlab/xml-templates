# Lints an xml-template file's dowser-template-info comment block against
# its actual ${VAR} placeholders. Checks for:
#   - a free-form ${VAR} in the XML body with no matching @var entry
#     ("missing docs")
#   - an @var entry documenting a placeholder that isn't in the XML body
#     anymore ("stale docs" -- generated variants get this cleaned up
#     automatically by filter_stale_template_info() in
#     make_all_templates_from_base.R, but a hand-edited base file can still
#     drift, e.g. a typo in a @var name)
#   - a dowser-template-info block that's opened but never closed with '-->'
#   - a literal "${NAME}" written inside the comment's own prose, which
#     Dowser's placeholder substitution would try to fill in at runtime
#     since it isn't comment-aware (this is a real bug caught once already
#     while writing these blocks by hand -- see UCRelaxedClock_Standard.xml's
#     history)
#   - a "--" sequence anywhere in the block's content, which is illegal
#     inside an XML comment regardless of what Dowser does with it
#
# Intentionally has no dependency on the dowser package (same as
# make_all_templates_from_base.R) -- this repo should be lintable/buildable
# on its own. That means the list of Dowser-special placeholder names below
# is a duplicate of the one in dowser's R/TemplateInfo.R and needs to be
# kept in sync by hand if that list ever changes.

# Placeholders Dowser fills in itself; not free-form kwargs a template
# author needs to document with @var. Keep in sync with describeTemplate()/
# .structuralChecklist() in dowser's R/TemplateInfo.R.
DOWSER_SPECIAL_VARS <- c("DATA", "CLONE", "DATE", "TRAIT", "TRAIT_NAME", "NODES",
                         "MRCA", "OPERATORS", "ROOTFREQS", "EMP_EQFREQS",
                         "ROOT_TYPE_PROBABILITIES", "NODES_TYPE_INIT")

# The subset of DOWSER_SPECIAL_VARS (plus two literal-tag structural
# requirements) that are *optional features*, not hard requirements --
# absence just means "this template doesn't support that feature", not
# "this template is broken." Only includes the ones we are likely to want 
# in every template. Reported as info.
# Condensed one-line notes here; see describeTemplate()/.structuralChecklist()
# in dowser's R/TemplateInfo.R for the full explanation of each.
STRUCTURAL_CHECKS <- list(
  list(name = "<init> tag",
       present = function(xml) any(grepl("<init", xml)) && any(grepl("</init>", xml)),
       note = "needed for custom starting trees"),
  list(name = "MRCA placeholder",
       present = function(xml) any(grepl("\\$\\{MRCA\\}", xml)),
       note = "max_start_date/germline-tip MRCA priors need this"),
  list(name = "ROOTFREQS placeholder",
       present = function(xml) any(grepl("\\$\\{ROOTFREQS\\}", xml)),
       note = "include_germline_as_root=TRUE silently no-ops without this"),
  list(name = "OPERATORS placeholder",
       present = function(xml) any(grepl("\\$\\{OPERATORS\\}", xml)),
       note = "if not GRT, germline behavior is weird without this"),
  list(name = "ROOT_TYPE_PROBABILITIES placeholder",
       present = function(xml) any(grepl("\\$\\{ROOT_TYPE_PROBABILITIES\\}", xml)),
       note = "for root_trait. dowser warns about this"),
  list(name = "EMP_EQFREQS placeholder",
       present = function(xml) any(grepl("\\$\\{EMP_EQFREQS\\}", xml)),
       note = "nucleotide frequencies aren't computed from data without this -- expected for non-_EmpFreq templates"),
  list(name = "DATE placeholder",
       present = function(xml) any(grepl("\\$\\{DATE\\}", xml)),
       note = "")
)

# element/present/note data.frame for one file. Purely informational --
# lint_template()/.is_clean() never fail a file over an absent entry here.
.structural_presence <- function(xml) {
  data.frame(
    element = vapply(STRUCTURAL_CHECKS, function(x) x$name, character(1)),
    present = vapply(STRUCTURAL_CHECKS, function(x) x$present(xml), logical(1)),
    note    = vapply(STRUCTURAL_CHECKS, function(x) x$note, character(1)),
    stringsAsFactors = FALSE
  )
}

# Finds the dowser-template-info block, if any. Doesn't assume it's
# well-formed -- unclosed=TRUE means an opening marker was found with no
# matching '-->' anywhere after it.
.locate_template_info_block <- function(xml) {
  start_idx <- grep("^\\s*<!--\\s*dowser-template-info\\s*$", xml)
  if (length(start_idx) == 0) {
    return(list(found = FALSE, unclosed = FALSE, start = NA, end = NA))
  }
  start_idx <- start_idx[1]

  end_idx <- grep("^\\s*-->\\s*$", xml)
  end_idx <- end_idx[end_idx > start_idx]
  if (length(end_idx) == 0) {
    return(list(found = TRUE, unclosed = TRUE, start = start_idx, end = NA))
  }

  list(found = TRUE, unclosed = FALSE, start = start_idx, end = end_idx[1])
}

# Every @var name declared inside a (well-formed) block's line range.
.documented_var_names <- function(xml, start, end) {
  if (is.na(start) || is.na(end) || end <= start + 1) {
    return(character(0))
  }
  block <- xml[(start + 1):(end - 1)]
  m <- regmatches(block, regexec("^\\s*@var\\s+(\\S+)\\s*\\|", block))
  names <- vapply(m, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  names[!is.na(names)]
}

# Every free-form ${VAR} actually used in the file, Dowser-special names
# excluded -- same regex Dowser's own substitution step uses to find
# required variables.
.free_var_names <- function(xml) {
  matches <- unlist(regmatches(xml, gregexpr("\\$\\{([^}]+)\\}", xml)))
  found <- unique(sub("\\$\\{([^}]+)\\}", "\\1", matches))
  sort(setdiff(found, DOWSER_SPECIAL_VARS))
}

#' Lint a single xml-templates file
#'
#' @param  path  path to a template .xml file
#' @return a list: file, unclosed_block (logical), missing_docs (character
#'         vector of undocumented free vars), stale_docs (character vector
#'         of @var entries with no matching placeholder), literal_placeholder
#'         (character vector of ${NAME} tokens found inside the comment's
#'         own prose), double_hyphen (logical). An empty/FALSE value in
#'         every field means the file is clean.
lint_template <- function(path) {
  xml <- readLines(path, warn = FALSE)
  loc <- .locate_template_info_block(xml)

  free_vars <- .free_var_names(xml)

  documented_vars <- character(0)
  literal_placeholder <- character(0)
  double_hyphen <- FALSE

  if (loc$found) {
    block_end <- if (is.na(loc$end)) length(xml) else loc$end
    block <- xml[loc$start:block_end]

    double_hyphen <- any(grepl("--", sub("^\\s*<!--", "", sub("-->\\s*$", "", block))))
    literal_hits <- unlist(regmatches(block, gregexpr("\\$\\{[A-Za-z_0-9]+\\}", block)))
    literal_placeholder <- unique(sub("\\$\\{([^}]+)\\}", "\\1", literal_hits))

    if (!loc$unclosed) {
      documented_vars <- .documented_var_names(xml, loc$start, loc$end)
    }
  }

  # if the block is missing or unclosed, nothing counts as documented --
  # every free var found is treated as missing docs
  missing_docs <- sort(setdiff(free_vars, documented_vars))
  stale_docs   <- if (loc$found && !loc$unclosed) sort(setdiff(documented_vars, free_vars)) else character(0)

  # DATA and CLONE aren't optional features like the rest of
  # DOWSER_SPECIAL_VARS -- a template missing either one can't be used by
  # Dowser at all, so this is a lint failure, not structural info.
  missing_required <- character(0)
  if (!any(grepl("\\$\\{DATA\\}", xml)))  missing_required <- c(missing_required, "DATA")
  if (!any(grepl("\\$\\{CLONE\\}", xml))) missing_required <- c(missing_required, "CLONE")

  list(
    file = path,
    unclosed_block = loc$unclosed,
    missing_docs = missing_docs,
    stale_docs = stale_docs,
    literal_placeholder = literal_placeholder,
    double_hyphen = double_hyphen,
    missing_required = missing_required,
    structural = .structural_presence(xml)
  )
}

# TRUE if a lint_template() result has nothing to report. Structural info
# (result$structural) is deliberately excluded -- it never affects clean/dirty.
.is_clean <- function(result) {
  !result$unclosed_block &&
    length(result$missing_docs) == 0 &&
    length(result$stale_docs) == 0 &&
    length(result$literal_placeholder) == 0 &&
    !result$double_hyphen &&
    length(result$missing_required) == 0
}

# Human-readable lines for one file's result. Empty vector if clean.
.format_lint_result <- function(result) {
  if (.is_clean(result)) {
    return(character(0))
  }
  lines <- character(0)
  if (result$unclosed_block) {
    lines <- c(lines, "  - dowser-template-info block is opened but never closed with '-->'")
  }
  if (length(result$missing_docs) > 0) {
    lines <- c(lines, paste0("  - undocumented variable(s): ", paste(result$missing_docs, collapse = ", ")))
  }
  if (length(result$stale_docs) > 0) {
    lines <- c(lines, paste0("  - @var entry with no matching placeholder: ", paste(result$stale_docs, collapse = ", ")))
  }
  if (length(result$literal_placeholder) > 0) {
    lines <- c(lines, paste0(
      "  - literal ${", paste(result$literal_placeholder, collapse = "}, ${"),
      "} written inside the comment -- use the bare name in prose, this will get substituted at runtime"
    ))
  }
  if (result$double_hyphen) {
    lines <- c(lines, "  - '--' found inside the comment block -- not valid inside an XML comment, use a single hyphen or an em dash")
  }
  if (length(result$missing_required) > 0) {
    lines <- c(lines, paste0(
      "  - missing required placeholder(s): ", paste(paste0("${", result$missing_required, "}"), collapse = ", "),
      " -- this template cannot be used by Dowser at all without these"
    ))
  }
  lines
}

# Info lines listing which optional structural elements are absent (never
# affects clean/dirty). Returns character(0) if everything's present.
.format_structural_info <- function(result) {
  absent <- result$structural[!result$structural$present, , drop = FALSE]
  if (nrow(absent) == 0) {
    return(character(0))
  }
  vapply(seq_len(nrow(absent)), function(i) {
    paste0("  - [info] ", absent$element[i], " absent: ", absent$note[i])
  }, character(1))
}

#' Lint every .xml file under a path (a single file, or a directory
#' searched recursively)
#'
#' @param  path  file or directory path
#' @return invisibly, TRUE if every file linted clean, FALSE otherwise.
#'         Prints a per-file report for anything that isn't clean.
lint_templates <- function(path) {
  if (dir.exists(path)) {
    files <- list.files(path, pattern = "\\.xml$", full.names = TRUE, recursive = TRUE)
  } else if (file.exists(path)) {
    files <- path
  } else {
    stop(paste("Path", path, "does not exist."))
  }

  all_clean <- TRUE
  for (f in files) {
    result <- lint_template(f)
    issue_lines <- .format_lint_result(result)
    info_lines <- .format_structural_info(result)

    if (length(issue_lines) > 0) {
      all_clean <- FALSE
    }

    if (length(issue_lines) > 0 || length(info_lines) > 0) {
      cat(f, "\n")
      if (length(issue_lines) > 0) {
        cat(paste(issue_lines, collapse = "\n"), "\n")
      }
      if (length(info_lines) > 0) {
        cat(paste(info_lines, collapse = "\n"), "\n")
      }
      cat("\n")
    }
  }

  if (all_clean) {
    cat(length(files), "template(s) checked, all clean", "(info notes above, if any, are not lint failures).\n")
  }

  invisible(all_clean)
}

# ---- CLI ---------------------------------------------------------------
# Only runs under `Rscript lint_template.R <path>`, not when this file is
# source()'d for its functions (e.g. from another script, or interactively
# while authoring a template) -- unlike make_all_templates_from_base.R,
# this guard is deliberate here since a clean run exits nonzero for CI.
if (!interactive() && sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Please provide a path to lint: a single .xml file, or a directory to search recursively.")
  }

  clean <- lint_templates(args[1])
  quit(status = if (clean) 0 else 1)
}