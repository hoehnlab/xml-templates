

# define a function which takes an airr clone object and returns xml of the sequences
create_alignment <- function(clone, name, include_germline_as_tip) {
  
  all_seqs <- ""
  
  for (i in 1:nrow(clone@data)) {
    # create a sequence object
    sequence <- clone@data[i, ]
    sequence_xml <- 
      paste0('\t<sequence id="seq_', sequence$sequence_id, 
             '" spec="Sequence" taxon="', sequence$sequence_id, 
             '" totalcount="4" value="', sequence$sequence, '" />\n')
    all_seqs <- paste0(all_seqs, sequence_xml)
  }
  
  if (include_germline_as_tip) {
    germline_sequence_xml <- 
      paste0('\t<sequence id="seq_', 'Germ', 
             '" spec="Sequence" taxon="', 'Germ', 
             '" totalcount="4" value="', clone@germline, '" />\n')
    all_seqs <- paste0(all_seqs, germline_sequence_xml)
  }
  
  alignment_xml <- 
    paste0('<data id="', name, "_", clone@clone, '" spec="Alignment" name="alignment">\n', 
           all_seqs, 
           '</data>')
  
  # create a taxon set for this alignment
  taxon_set <- 
    paste0('<taxa id="TaxonSet.', name, "_", clone@clone, '" spec="TaxonSet">\n', 
           '\t<alignment idref="', name, "_", clone@clone, '"/>\n', 
           '</taxa>')
  
  return(paste0(alignment_xml, '\n', taxon_set, '\n'))
}

create_root_freqs <- function(clone, name) {
  if (any(grepl("N", clone@germline))) {
    freqs <- clone@germline
    freqs <- gsub("N", " 0.25,0.25,0.25,0.25;", freqs)
    freqs <- gsub("A", " 1,0,0,0;", freqs)
    freqs <- gsub("C", " 0,1,0,0;", freqs)
    freqs <- gsub("G", " 0,0,1,0;", freqs)
    freqs <- gsub("T", " 0,0,0,1;", freqs)
    freqs <- gsub("(^[[:space:]]+|[[:space:]]+$)", "", freqs)
    freqs <- substr(freqs, 1, nchar(freqs)-1)
    root_freqs <- paste0('<rootfreqseq id="seq_Root', name, "_", clone@clone, '" spec="Sequence" taxon="Root', name, "_", clone@clone, '" uncertain="true"
    totalcount="4" value="', freqs, '"/>', sep="")
    return(root_freqs)
  }
  root_freqs <- paste0('<rootfreqseq id="seq_Root" spec="Sequence" taxon="Root"
      totalcount="4" value="', clone@germline,'"/>', sep="")
  return(root_freqs)
}

create_MRCA_prior_observed <- function(clone, name) {
  taxa <- paste0('<taxon id="', clone@data$sequence_id, '" spec="Taxon"/>', collapse="\n")
  distribution_xml <- 
    paste0('<distribution id="obs.prior" spec="beast.base.evolution.tree.MRCAPrior" monophyletic="true" tree="@Tree.t:', name, "_", clone@clone, '">\n', 
           '<taxonset id="obs" spec="TaxonSet">\n', 
           taxa, 
           '\n</taxonset>\n',
           '</distribution>', sep="")
  return(distribution_xml)
}

create_MRCA_prior_germline <- function(clone, name) {
  taxa <- paste0('<taxon id="', 'Germ', '" spec="Taxon"/>', collapse="\n")
  distribution_xml <- 
    paste0('<distribution id="germ1.prior" spec="beast.base.evolution.tree.MRCAPrior" tipsonly="true" tree="@Tree.t:', name, "_", clone@clone, '">\n', 
           '<taxonset id="germSet" spec="TaxonSet">\n', 
           taxa, 
           '\n</taxonset>\n',
           '<Uniform id="Uniform.1:germ" name="distr" lower = "0.0" upper="10000"/>\n',
           '</distribution>', sep="")
  return(distribution_xml)
}

create_traitset <- function(clone, trait_name, column, name, trait_data_type=NULL, isSet=FALSE, include_germline_as_tip=FALSE) {
  all_traits <- paste(clone@data$sequence_id, clone@data[[column]], collapse=",\n", sep="=")
  if (include_germline_as_tip) {
    all_traits <- paste(all_traits, paste0('Germ','=', '?'), sep=",\n")
  }
  tagname <- "trait" 
  if (isSet) {
    tagname <- "traitSet"
  }
  traitset_xml <- 
    paste0('<', tagname ,' id="', trait_name, ":", name, "_", clone@clone, 
           '" spec="beast.base.evolution.tree.TraitSet" taxa="@TaxonSet.', name, "_", clone@clone, 
           '" traitname="', trait_name,
           '">\n', 
           all_traits, 
           '</',tagname,'>\n',
           trait_data_type)
  
  return(traitset_xml)
}

# define an xml writer function
xml_writer_clone <- function(clone, file, name, date=NULL, trait=NULL, trait_data_type=NULL, template=NULL, replacements=NULL, include_germline_as_root=FALSE, include_germline_as_tip=FALSE, ...) {
  
  kwargs <- list(...)
  
  # read in a template file
  if (is.null(template)) {
    # template <- system.file("extdata", "template.xml", package = "scoper")
    template = "template.xml"
  }
  xml <- readLines(template)
  
  data <- create_alignment(clone, name, include_germline_as_tip)
  # replace the ${DATA} placeholder with actual data
  xml <- gsub("\\$\\{DATA\\}", data, xml)
  xml <- gsub("\\$\\{CLONE\\}", paste0(name, "_", clone@clone), xml)
  
  # if the date argument is null but ${DATE} is in the template file, raise an error
  if (is.null(date) && any(grepl("\\$\\{DATE\\}", xml))) {
    stop("Date argument is NULL but ${DATE} is in the template file")
  }
  if (!is.null(date)) {
    date_trait <- create_traitset(clone, "date", date, name)
    # replace the ${DATE} placeholder with the dates
    xml <- gsub("\\$\\{DATE\\}", date_trait, xml)
  }
  
  if (is.null(trait) && any(grepl("\\$\\{TRAIT\\}", xml))) {
    stop("Trait argument is NULL but ${TRAIT} is in the template file")
  }
  if (!is.null(trait) && is.null(trait_data_type)) {
    stop("Trait argument is given but no trait data type was provided")
  }
  if (!is.null(trait)) {
    sample_trait <- create_traitset(clone, "newTrait", trait, name, trait_data_type, isSet=TRUE, include_germline_as_tip=include_germline_as_tip)
    # replace the ${TRAIT} placeholder with the sample trait
    xml <- gsub("\\$\\{TRAIT\\}", sample_trait, xml)
  }
  if (any(grepl("\\$\\{NODES\\}", xml))) {
    # replace the ${NODES} placeholder with the number of nodes in this tree
    tips <- nrow(clone@data)
    if (include_germline_as_tip) {
      tips <- tips + 1
    }
    nodes <- 2*tips-1
    xml <- gsub("\\$\\{NODES\\}", nodes, xml)
  }
  
  if (any(grepl("\\$\\{MRCA\\}", xml))) {
    # replace the ${MRCA} placeholder with the mrca prior
    mrca_priors <- ""
    if (include_germline_as_tip) {
      mrca_priors <- paste0(create_MRCA_prior_observed(clone, name), create_MRCA_prior_germline(clone, name), sep="\n")
    }
    
    xml <- gsub("\\$\\{MRCA\\}", mrca_priors, xml)
  }

  if (any(grepl("\\$\\{OPERATORS\\}", xml))) {
    # replace the ${OPERATORS} placeholder with the operators we want to add
    # can add other potential operators here
    operators <- ""
    if (include_germline_as_tip) {
      
      operators <-  paste0('<operator id="TipDatesRandomWalker.01" windowSize="1" spec="beast.base.evolution.operator.TipDatesRandomWalker" taxonset="@germSet" tree="@Tree.t:',name, '_', clone@clone,'" weight="1.0"/>\n')
    }
    
    xml <- gsub("\\$\\{OPERATORS\\}", operators, xml)
  }
  if (any(grepl("\\$\\{ROOTFREQS\\}", xml))) { 
    root_freqs <- ""
    if (include_germline_as_root) {
      # replace spec="TreeLikelihood" with spec="rootfreqs.TreeLikelihood"
      xml <- gsub('spec="TreeLikelihood"', 'spec="rootfreqs.TreeLikelihood"', xml)
      # replace the ${ROOTFREQS} placeholder with the root frequencies
      # TODO(jf): does this work with Ns??
      root_freqs <- create_root_freqs(clone, name)
    }
    xml <- gsub("\\$\\{ROOTFREQS\\}", root_freqs, xml)
  }


  
  if (!is.null(replacements)) {
    for (replacement in replacements) {
      xml <- gsub(paste0("\\$\\{", replacement, "\\}"), kwargs[[replacement]], xml)
    }
  }
  
  
  # open a connection to the file
  file <- paste0(file, "_", clone@clone, ".xml")
  con <- file(file, "w")
  
  # write the XML file
  writeLines(xml, con)
  
  # close the connection
  close(con)
}

xml_writer_wrapper <- function(clones, name, date=NULL, trait=NULL, template=NULL, outfile=NULL, replacements=NULL, trait_list=NULL, include_germline_as_root=FALSE, include_germline_as_tip=FALSE, ...) {
  # iterate over the clones to first create trait data type if trait exists
  if (!is.null(trait)) {
    if (is.null(trait_list)) {
      # get all the possible values of the trait
      traits <- c()
      for (i in 1:nrow(clones)) {
        traits <- c(traits, unique(clones$data[[i]]@data[[trait]]))
      }
      trait_list <- unique(traits)
    }
    codeMap <- paste(trait_list, 0:(length(trait_list)-1), collapse=",\n", sep="=")
    codeMap <- paste0(codeMap, ",\n? = ", paste0(0:(length(trait_list)-1), collapse=" "))
    trait_data_type <- paste0('<userDataType id="traitDataType.newTrait" spec="beast.base.evolution.datatype.UserDataType" codeMap="', codeMap, '" codelength="-1" states="', length(trait_list), '"/>')
  }
  for (i in 1:nrow(clones)){
    xml_writer_clone(clones$data[[i]], 
                     file=outfile, 
                     name=name, 
                     date=date, 
                     trait=trait, 
                     trait_data_type=trait_data_type, 
                     template=template,
                     replacements=replacements, 
                     include_germline_as_root=include_germline_as_root,
                     include_germline_as_tip=include_germline_as_tip,
                     ...)
  }
}
