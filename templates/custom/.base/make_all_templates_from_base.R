
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide the path to your xml-writer directory.")
}
xml_path <- args[1]


base_path <- file.path(xml_path, "templates/custom/.base")
all_bases <- list.files(base_path, pattern = "\\.xml$", full.names = TRUE)

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
    if (arg %in% c("CLONE", "DATE", "TRAIT", "MRCA")) {
      next # skip these as they don't really work as fixed values
    }
    value <- kwargs[[arg]]
    xml <- gsub(paste0("\\$\\{", arg, "\\}"), value, xml)
  }
  
  # Write the modified XML to the output path
  writeLines(xml, output_path)

}

# make custom templates for each category as is and without emp_freqs
# make standard templates with fewer arguments
custom_folder <- file.path(xml_path, "templates/custom")
standard_folder <- file.path(xml_path, "templates")

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
  if (grepl("Trait", template_name)) {
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
      STANDARD_ARGS
    ))
  )
  
  # make a version of this without emp_freqs and with standard arguments in standard
  do.call(
    replace_specific_arguments, 
    as.list(c(
      file, 
      file.path(current_standard_folder, gsub("_EmpFreq", "", template_name)), 
      EMP_EQFREQS_ARGS,
      STANDARD_ARGS
    ))
  )
  
}

