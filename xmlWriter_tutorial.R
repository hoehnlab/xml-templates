library(airr)
library(dowser)

# set this to the directory where xmlWriter.R is located
directory <- "~/xml_writer/" 

setwd(directory)
source("xmlWriter.R")


# read in some data
data_dir <- "example_data"
data_file <- paste(data_dir, "/all_samples_airr.tsv", sep="")
example <- read_rearrangement(data_file)

# any filtering you'd like to do, e.g.
example <- filter(example, sample_time %in% c(75, 100))

# format clones
example <- resolveLightChains(example)
clones <- formatClones(example, traits=c("location", "sample_time", "celltype"), chain="HL", nproc=1, collapse = FALSE, 
                       split_light = TRUE, minseq = 3, germ="germline_alignment")

# down sample clones if you'd like


# create an xml file for each clone
xml_writer_wrapper(clones, 
                   "tutorial", 
                   outfile="tutorial", 
                   date="sample_time", 
                   template="strict_clock_template.xml")

# create an xml file for each clone with a trait
xml_writer_wrapper(clones, 
                   "tutorial_trait", 
                   outfile="tutorial_trait", 
                   date="sample_time", 
                   trait="location", 
                   template="TraitLinkedInstantSwitch_FixedTraitClockRates.xml")

# create an xml file for each clone with a trait, and specify all possible traits instead of using automatic trait list
trait_list <- c("germinal_center", "other")
xml_writer_wrapper(clones, 
                   "tutorial_trait_list", 
                   outfile="tutorial_trait_list", 
                   date="sample_time", 
                   trait="location", 
                   trait_list=trait_list,
                   template="TraitLinkedInstantSwitch_FixedTraitClockRates.xml")

# create an xml file for each clone and include the germline as a tip, with a tip-date operator and an mrca prior that constrains observed tips to be monophyletic
xml_writer_wrapper(clones, 
                   "tutorial_trait_list", 
                   outfile="tutorial_trait_list", 
                   date="sample_time", 
                   trait="location", 
                   trait_list=trait_list,
                   template="TraitLinkedInstantSwitch_FixedTraitClockRates.xml",
                   include_germline=TRUE)


# create an xml file for each clone with some other replacements in the xml
# provide a template with the '${REPLACEMENT}' syntax
xml_writer_wrapper(clones, 
                   "tutorial_replacements", 
                   outfile="tutorial_replacements", 
                   template="replacement_template_no_date.xml",
                   replacements=c("MATRIX", "FREQUENCIES"), 
                   MATRIX="sample_matrix_text", 
                   FREQUENCIES="0.5")

#############################
# a list of useful templates and examples of using them:

# TRAIT LINKED CLOCK MODEL

# NOTE for trait linked model: trait list with GC first ensures GC is 0 and other is 1
trait_list <- c("germinal_center", "other")

# trait linked clock model with date and trait
xml_writer_wrapper(clones, 
                   "trait_linked", 
                   outfile="trait_linked", 
                   date="sample_time", 
                   trait="location", 
                   trait_list=trait_list,
                   template="TraitLinkedInstantSwitch_FixedTraitClockRates.xml")

# trait linked clock model, Instant Switch, Estimated clock rates
xml_writer_wrapper(clones, 
                   "trait_linked", 
                   outfile="trait_linked", 
                   date="sample_time", 
                   trait="location", 
                   trait_list=trait_list,
                   template="TraitLinkedInstantSwitch_EstimatedTraitClockRates.xml")

# trait linked clock model, Expected Occupancy, Fixed clock rates
xml_writer_wrapper(clones, 
                   "trait_linked", 
                   outfile="trait_linked", 
                   date="sample_time", 
                   trait="location", 
                   trait_list=trait_list,
                   template="TraitLinkedExpectedOccupancy_FixedTraitClockRates.xml")

# trait linked clock model, Expected Occupancy, Estimated clock rates
xml_writer_wrapper(clones, 
                   "trait_linked", 
                   outfile="trait_linked", 
                   date="sample_time", 
                   trait="location", 
                   trait_list=trait_list,
                   template="TraitLinkedExpectedOccupancy_EstimatedTraitClockRates.xml")


# trait linked clock model with no date
xml_writer_wrapper(clones, 
                   "trait_linked_no_date", 
                   outfile="trait_linked_no_date", 
                   trait="location", 
                   trait_list=trait_list,
                   template="trait_linked_template_no_date.xml")

#############################
# STRICT CLOCK MODEL

# TODO(jf): re-check all these templates and then re-include them

# # strict clock model with date
# xml_writer_wrapper(clones, 
#                    "strict_clock", 
#                    outfile="strict_clock", 
#                    date="sample_time", 
#                    template="strict_clock_template.xml")
# 
# # strict clock model with no date
# xml_writer_wrapper(clones, 
#                    "strict_clock_no_date", 
#                    outfile="strict_clock_no_date", 
#                    date="sample_time", 
#                    template="strict_clock_template_no_date.xml")
# 
# # strict clock model with BEAST CLASSIC ancestral reconstruction, with date
# xml_writer_wrapper(clones, 
#                    "strict_clock_AR", 
#                    outfile="strict_clock_AR", 
#                    date="sample_time", 
#                    trait="location", 
#                    template="strict_clock_ancestral_reconstruction_template.xml")
# 
# # strict clock model with BEAST CLASSIC ancestral reconstruction, with no date
# xml_writer_wrapper(clones, 
#                    "strict_clock_AR_no_date", 
#                    outfile="strict_clock_AR_no_date", 
#                    date="sample_time", 
#                    trait="location", 
#                    template="strict_clock_ancestral_reconstruction_template_no_date.xml")

#############################
# EPOCH CLOCK MODEL


#############################
# GIBLE

# This does NOT currently work with example data, traitlist is made up, nLoc is wrong, etc
adjacency_matrix <- matrix(c(0, 1, 0, 1, 0, 0, 1, 0, 1), nrow=3)
contextual_matrix <- matrix(c(0, 0, 0, 2.5, 0, 0, 1, 1, 1), nrow=3)
nloc <- nrow(adjacency_matrix)
traitlist <- c("hello", "world", "test")

xml_writer_wrapper(clones, 
                   outfile="GIBLE", 
                   name="GIBLE", 
                   trait="location", 
                   replacements=c("NLOC", "FREQUENCIES", "DMATRIX", "CMATRIX"), 
                   template="GIBLE_xmlwriter_template.xml", 
                   FREQUENCIES=paste(c(rep(1/nloc, (nloc-1)), 1-((nloc-1)/nloc)), collapse=" "), 
                   NLOC=paste(nloc, collapse=" "), 
                   date="sample_time", 
                   trait_list=traitlist, 
                   DMATRIX=paste0("\n", paste(apply(adjacency_matrix, 1, paste, collapse = " "), collapse = "\n"), "\n"), 
                   CMATRIX=paste0("\n", paste(apply(contextual_matrix, 1, paste, collapse = " "), collapse = "\n"), "\n")
                   )

#############################
# suggestion for comparing a number of different models and parameters or for easily changing values that are expected to vary by dataset, e.g. trait clock rates and their priors
# see files in temmplates/custom for how to customize the templates to allow for this

trait_list <- c("germinal_center", "other")

all_templates <- c(
    "TraitLinkedExpectedOccupancy_EstimatedTraitClockRates.xml",
    "TraitLinkedExpectedOccupancy_FixedTraitClockRates.xml",
    "TraitLinkedInstantSwitch_FixedTraitClockRates.xml",
    "TraitLinkedInstantSwitch_EstimatedTraitClockRates.xml",
    "TraitLinkedInstantSwitch_RelaxedTraitClockRates.xml",
    "TraitLinkedExpectedOccupancy_RelaxedTraitClockRates.xml")

template_names <- c(
    "ExpectedOccupancy_EstClockRates",
    "ExpectedOccupancy_FixedClockRates",
    "InstantSwitch_FixedClockRates",
    "InstantSwitch_EstClockRates",
    "InstantSwitch_RelaxedClockRates",
    "ExpectedOccupancy_RelaxedClockRates")

# Loop through each template and create the XML file
for (i in 1:length(all_templates)) {
    template <- all_templates[i]
    template_name <- template_names[i]
    out_name <- paste0(name, "_", template_name, sep="")
    out_file <- paste0(output_file, "_", template_name, sep="")

    START_TRAIT_RATES <- "0.0003842593 0.000004"
    if (any(grepl("Relaxed", template))) {
        # If the template is a relaxed clock model, we need multiple clock rates per trait
        START_TRAIT_RATES <- "0.0006 0.0003842593 0.00004 0.000004"
    }

    TRANSITION_RATE_ALPHA_1 <- "1.0"
    TRANSITION_RATE_BETA_1 <- "1.0"
    TRANSITION_RATE_ALPHA_2 <- "0.0001"
    TRANSITION_RATE_BETA_2 <- "1.0"

    TRAIT_RATE_MEAN_1 <- "0.0003842593"
    TRAIT_RATE_MEAN_2 <- "0.000004"
    TRAIT_RATE_SIGMA_1 <- "0.005"
    TRAIT_RATE_SIGMA_2 <- "0.0005"

    RATE_INDICATORS_LIST <- c("1 0", "1 1")

    for (j in 1:length(RATE_INDICATORS_LIST)) {
        # Create the XML file using the current template
        xml_writer_wrapper(clones, 
                            out_name,
                            outfile=out_file, 
                            date="sample_time", 
                            trait="location", 
                            trait_list=trait_list,
                            template=paste0(xml_dir, "custom/", template, sep=""),
                            START_TRAIT_RATES=START_TRAIT_RATES,
                            RATE_INDICATORS=RATE_INDICATORS_LIST[j],
                            TRANSITION_RATE_ALPHA_1=TRANSITION_RATE_ALPHA_1,
                            TRANSITION_RATE_BETA_1=TRANSITION_RATE_BETA_1,
                            TRANSITION_RATE_ALPHA_2=TRANSITION_RATE_ALPHA_2,
                            TRANSITION_RATE_BETA_2=TRANSITION_RATE_BETA_2,
                            TRAIT_RATE_MEAN_1=TRAIT_RATE_MEAN_1,
                            TRAIT_RATE_MEAN_2=TRAIT_RATE_MEAN_2,
                            TRAIT_RATE_SIGMA_1=TRAIT_RATE_SIGMA_1,
                            TRAIT_RATE_SIGMA_2=TRAIT_RATE_SIGMA_2
                            )
        # Create the XML file using the current template and include germline
        xml_writer_wrapper(clones, 
                            paste0(out_name, "_germline_root", sep=""),
                            outfile=paste0(out_file, "_germline_root", sep=""), 
                            date="sample_time", 
                            trait="location", 
                            trait_list=trait_list,
                            template=paste0(xml_dir, "custom/", template, sep=""),
                            include_germline_as_root=TRUE)
        xml_writer_wrapper(clones, 
                            paste0(out_name, "_germline_root", sep=""),
                            outfile=paste0(out_file, "_germline_root", sep=""),
                            date="sample_time", 
                            trait="location", 
                            trait_list=trait_list,
                            template=paste0(xml_dir, "custom/", template, sep=""),
                            include_germline_as_root=TRUE,
                            START_TRAIT_RATES=START_TRAIT_RATES,
                            RATE_INDICATORS=RATE_INDICATORS_LIST[j],
                            TRANSITION_RATE_ALPHA_1=TRANSITION_RATE_ALPHA_1,
                            TRANSITION_RATE_BETA_1=TRANSITION_RATE_BETA_1,
                            TRANSITION_RATE_ALPHA_2=TRANSITION_RATE_ALPHA_2,
                            TRANSITION_RATE_BETA_2=TRANSITION_RATE_BETA_2,
                            TRAIT_RATE_MEAN_1=TRAIT_RATE_MEAN_1,
                            TRAIT_RATE_MEAN_2=TRAIT_RATE_MEAN_2,
                            TRAIT_RATE_SIGMA_1=TRAIT_RATE_SIGMA_1,
                            TRAIT_RATE_SIGMA_2=TRAIT_RATE_SIGMA_2
                            )
        }
}

#############################
# features that are on my to-do list:
# - allow replacements to be functions, so they can be clone-specific
# - allow for heirarchical models
# - handle default templates
# - handle default outfiles/names
# - better error handling and warnings 
# - create and test more templates
# - built in downsampling options?
# - possibly allow for multiple traits (not currently working in BEAST so not a priority)



