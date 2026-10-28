
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide the path to your xml-writer directory.")
}
xml_path <- args[1]


base_path <- file.path(xml_path, ".base")
all_bases <- list.files(base_path, pattern = "\\.xml$", full.names = TRUE)
 
# A generated variant (simplified/, or the _EmpFreq-stripped copy) may have
# some ${VAR} placeholders filled in with fixed values by this script, which
# removes them from the XML body. If the base template's dowser-template-info
# block documented one of those variables, that entry is now stale in the
# generated file -- it describes a variable the file no longer has. This
# drops any @var entry (and its continuation lines, if any) whose placeholder
# no longer appears anywhere in the generated file. @struct notes are never
# auto-dropped, since they aren't tied to a single placeholder's presence.
filter_stale_template_info <- function(xml) {
  start_idx <- grep("^\\s*<!--\\s*dowser-template-info\\s*$", xml)
  if (length(start_idx) == 0) {
    return(xml) # no block to filter
  }
  start_idx <- start_idx[1]
 
  end_idx <- grep("^\\s*-->\\s*$", xml)
  end_idx <- end_idx[end_idx > start_idx]
  if (length(end_idx) == 0) {
    return(xml) # malformed block, leave it alone rather than guess
  }
  end_idx <- end_idx[1]
 
  if (end_idx == start_idx + 1) {
    return(xml) # empty block, nothing to filter
  }
 
  block <- xml[(start_idx + 1):(end_idx - 1)]
  outside_block <- xml[-(start_idx:end_idx)]
 
  var_pattern         <- "^\\s*@var\\s+(\\S+)\\s*\\|"
  struct_pattern      <- "^\\s*@struct\\s+(\\S+)\\s*\\|"
  description_pattern <- "^\\s*@description\\s*\\|\\s*"
 
  keep <- rep(TRUE, length(block))
  current_keep <- TRUE
 
  for (i in seq_along(block)) {
    line <- block[i]
    if (grepl(var_pattern, line)) {
      name <- regmatches(line, regexec(var_pattern, line))[[1]][2]
      still_present <- any(grepl(paste0("\\$\\{", name, "\\}"), outside_block))
      current_keep <- still_present
      keep[i] <- current_keep
    } else if (grepl(struct_pattern, line)) {
      current_keep <- TRUE # @struct notes are never auto-dropped
      keep[i] <- TRUE
    } else if (grepl(description_pattern, line)) {
      current_keep <- TRUE # @description is never auto-dropped -- it isn't tied to a placeholder at all
      keep[i] <- TRUE
    } else {
      keep[i] <- current_keep # continuation line inherits its entry's fate
    }
  }
 
  new_block <- block[keep]
  c(xml[1:start_idx], new_block, xml[end_idx:length(xml)])
}

# If this generation step supplies a fixed EMP_EQFREQS value (see
# EMP_EQFREQS_ARGS below), the output file bakes in flat nucleotide
# frequencies instead of computing them from data -- the @description's
# frequency-handling sentence, if it's the standard "Dowser calculates
# empirical..." wording, is now wrong for this specific generated file.
# Swap it for the matching flat-frequency wording. Scoped to the
# dowser-template-info block only, and only swaps the exact canonical
# sentence both variants are always written with -- this is not a general
# prose rewriter, just a fixed substitution for one fact the generator
# itself controls (same spirit as EMP_EQFREQS_ARGS hardcoding the XML
# placeholder's replacement value).
EMPIRICAL_FREQ_SENTENCE <- "Dowser calculates empirical nucleotide frequencies from your data."
FLAT_FREQ_SENTENCE <- "Nucleotide frequencies start flat and can move during MCMC."

adjust_description_for_flat_freqs <- function(xml, kwargs) {
  if (!("EMP_EQFREQS" %in% names(kwargs))) {
    return(xml) # this generation step doesn't bake in flat frequencies
  }

  start_idx <- grep("^\\s*<!--\\s*dowser-template-info\\s*$", xml)
  if (length(start_idx) == 0) {
    return(xml) # no block, nothing to adjust
  }
  start_idx <- start_idx[1]

  end_idx <- grep("^\\s*-->\\s*$", xml)
  end_idx <- end_idx[end_idx > start_idx]
  if (length(end_idx) == 0) {
    return(xml) # malformed block, leave it alone rather than guess
  }
  end_idx <- end_idx[1]

  block_range <- start_idx:end_idx
  xml[block_range] <- gsub(EMPIRICAL_FREQ_SENTENCE, FLAT_FREQ_SENTENCE, xml[block_range], fixed = TRUE)
  xml
}

# simplified/ variants bake in default priors/inits (see STANDARD_ARGS
# below) so a caller has fewer ${VAR}s to supply -- that's a real, useful
# difference worth surfacing right in the file, not just in the README.
# The main generation loop passes IS_SIMPLIFIED as a marker (its value is
# never used, only its presence) on the two do.call() branches that also
# pass STANDARD_ARGS, i.e. exactly the branches writing to
# current_standard_folder.
SIMPLIFIED_NOTICE <- "This is a simplified template with more fixed values. For more flexibility, use the standard template."

add_simplified_notice <- function(xml, kwargs) {
  if (!("IS_SIMPLIFIED" %in% names(kwargs))) {
    return(xml) # not a simplified/ generation step
  }

  start_idx <- grep("^\\s*<!--\\s*dowser-template-info\\s*$", xml)
  if (length(start_idx) == 0) {
    return(xml) # no dowser-template-info block at all -- nothing to attach the notice to
  }
  start_idx <- start_idx[1]

  end_idx <- grep("^\\s*-->\\s*$", xml)
  end_idx <- end_idx[end_idx > start_idx]
  if (length(end_idx) == 0) {
    return(xml) # malformed block, leave it alone rather than guess
  }
  end_idx <- end_idx[1]

  description_pattern <- "^(\\s*@description\\s*\\|\\s*)(.*)$"
  for (i in start_idx:end_idx) {
    if (grepl(description_pattern, xml[i])) {
      m <- regmatches(xml[i], regexec(description_pattern, xml[i]))[[1]]
      xml[i] <- paste0(m[2], SIMPLIFIED_NOTICE, " ", m[3])
      return(xml) # exactly one @description per template, per the convention
    }
  }

  # no @description entry existed to prepend onto -- insert a new one right
  # after the opening marker, rather than silently skipping the notice
  xml <- c(xml[1:start_idx], paste0("@description | ", SIMPLIFIED_NOTICE), xml[(start_idx + 1):length(xml)])
  xml
}

replace_specific_arguments <- function(template, output_path, ...) {
  # This function replaces specific arguments in the template file with the provided values.
  # It is used to create a new XML file with the specified parameters.
  
  # Check if the template file exists
  if (!file.exists(template)) {
    stop(paste("Template file", template, "does not exist."))
  }
  
  # Check if the output path is provided
  if (is.null(output_path)) {
    stop("Output path must be provided.")
  }
  
  # Read in the template file
  xml <- readLines(template)
  
  kwargs <- list(...)
  
  # Replace specific arguments in the XML
  for (arg in names(kwargs)) {
    if (arg %in% c("CLONE", "DATE", "TRAIT", "MRCA", "IS_SIMPLIFIED")) {
      next # skip these as they don't really work as fixed values
    }
    value <- kwargs[[arg]]
    xml <- gsub(paste0("\\$\\{", arg, "\\}"), value, xml)
  }
  
  # Drop any dowser-template-info @var entries whose placeholder was just
  # filled in above, so generated variants don't ship stale documentation.
  xml <- filter_stale_template_info(xml)

  # If this step baked in flat frequencies, fix the @description sentence
  # that would otherwise still claim empirical frequencies are computed.
  xml <- adjust_description_for_flat_freqs(xml, kwargs)

  # If this step is writing a simplified/ variant (STANDARD_ARGS baked in),
  # prepend a notice to that effect onto the @description.
  xml <- add_simplified_notice(xml, kwargs)
 
  # Write the modified XML to the output path
  writeLines(xml, output_path)

}

# make custom templates for each category as is and without emp_freqs
# make standard templates with fewer arguments
custom_folder <- file.path(xml_path, "")
standard_folder <- file.path(xml_path, "simplified")

operator_string <- '<operator id="FrequenciesExchanger.s:${CLONE}" spec="DeltaExchangeOperator" delta="0.01" weight="0.1">
            <parameter idref="freqParameter.s:${CLONE}"/>
        </operator>
        ${OPERATORS}'

EMP_EQFREQS_ARGS <- c(
  "EMP_EQFREQS" = "0.25 0.25 0.25 0.25",
  "OPERATORS" = operator_string
)
STANDARD_ARGS <- c(
  # "RATE_INDICATORS" = "1 0",
  
  # Expected occupancy inputs
  "TRANSITION_RATE_ALPHA_1" = "0.1",
  "TRANSITION_RATE_BETA_1" = "1.0",
  "TRANSITION_RATE_ALPHA_2" = "0.1",
  "TRANSITION_RATE_BETA_2" = "1.0",
  # "TRAIT_RATE_MEAN_1" = 4.9E-3/7,
  # "TRAIT_RATE_MEAN_2" = "0.000001"
  # "TRAIT_RATE_SIGMA_1" = 0.001 * (4.9E-3/7),
  # "TRAIT_RATE_SIGMA_2" = "0.00005",
  "KAPPA_PRIOR_M" = "0.67",
  "KAPPA_PRIOR_S" = "0.2",
  
  # Strict clock/UCLD inputs
  "CLOCK_RATE_INIT" = 4.9E-3/20,
  "TRANSITION_RATE_ALPHA" = "0.1",
  "TRANSITION_RATE_BETA" = "1.0",
  "UCLD_MU_INIT" = "1.0",
  "UCLD_SIGMA_INIT" = "0.5"
)

for (file in all_bases) {
  # Extract the base name without the path and extension
  template_name <- basename(file)
  if (grepl("TypeLinked", template_name)) {
    current_custom_folder <- file.path(custom_folder, "TypeLinked")
    current_standard_folder <- file.path(standard_folder, "TypeLinked")
  } else if (grepl("StrictClock", template_name)) {
    current_custom_folder <- file.path(custom_folder, "StrictClock")
    current_standard_folder <- file.path(standard_folder, "StrictClock")
  } else if (grepl("UC", template_name)) {
    current_custom_folder <- file.path(custom_folder, "UCLD")
    current_standard_folder <- file.path(standard_folder, "UCLD")
  } else {
    next # Skip files that do not match the expected patterns
  }
  
  # make a copy of this as-is in custom
  do.call(
    replace_specific_arguments, 
    as.list(c(
      file, 
      file.path(current_custom_folder, template_name)
      ))
    )
  # replace_specific_arguments(file, file.path(custom_folder, template_name))
  
  # make a version of this without emp_freqs in custom
  do.call(
    replace_specific_arguments, 
    as.list(c(
      file, 
      file.path(current_custom_folder, gsub("_EmpFreq", "", template_name)), 
      EMP_EQFREQS_ARGS
    ))
  )
  
  # make a version of this with standard arguments in standard
  do.call(
    replace_specific_arguments, 
    as.list(c(
      file, 
      file.path(current_standard_folder, template_name), 
      STANDARD_ARGS,
      "IS_SIMPLIFIED" = "TRUE"
    ))
  )
  
  # make a version of this without emp_freqs and with standard arguments in standard
  do.call(
    replace_specific_arguments, 
    as.list(c(
      file, 
      file.path(current_standard_folder, gsub("_EmpFreq", "", template_name)), 
      EMP_EQFREQS_ARGS,
      STANDARD_ARGS,
      "IS_SIMPLIFIED" = "TRUE"
    ))
  )
  
}

