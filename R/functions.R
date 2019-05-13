# Data import -------------------------------------------------------------

# Untar with tracking
untar_tracked <- function (tarfile, exdir, ...) {
  untar(tarfile = tarfile, exdir = exdir)
}

# For some reason there are a bunch of duplicate names with
# different sequences in the lygodium data. Exclude these.
trim_proteome <- function (proteome) {
  dup_names <- names(proteome)[duplicated(names(proteome))]
  trimmed_proteome <- proteome[!(names(proteome) %in% dup_names)]
  return(trimmed_proteome)
}

# Report functions --------------------------------------------------------

# Calculate length, mean coverage, and number of sequences in each alignment for all alignments
# in a folder (for plan_01 only)
calculate_alignment_summary_stats <- function (folder, filter_term = "fa.mafft.aln-cln$", ingroup = NULL) {
  
  require(purrr)
  require(ape)
  require(ips)
  require(here)
  require(stringr)
  require(jntools)
  
  # for easier pasting
  folder <- jntools::add_slash(folder)
  
  # get list of alignments in folder
  alignment_files <- list.files(folder, pattern = filter_term)
  
  alignment_files <- paste0(folder, alignment_files)
  
  # read in all alignments as list
  alignments <- map(alignment_files, read.FASTA)
  
  # convert to matrix
  alignments <- map(alignments, as.matrix)
  
  # optionally trim to only ingroup taxa
  if (!(is.null(ingroup))) {
    alignments <- lapply (alignments, function (x) x[rownames(x) %in% ingroup,])
    # trim all out "empty" columns (those containing only missing data)
    alignments <- map (alignments, ips::deleteEmptyCells, quiet=TRUE)
  }
  
  # calculate length of each alignment
  length <- map_dbl(alignments, ncol)
  
  # function to calculate percent of each sequence in alignment that is non-gap characters
  calc_coverage <- function (align) {
    nongap_seq_length <- apply (as.character(align), 1, function (x) length(x[grep("n|N|-", x, invert=TRUE)]))
    seq_coverage <- nongap_seq_length/ncol(align)
    return(seq_coverage)
  }
  
  # use above function to calculate percentage coverage for ALL sequences in ALL alignments
  coverage <- unlist(lapply(alignments, calc_coverage))
  
  # function to calculate MEAN percent of each sequence in alignment that is non-gap characters for each alignment
  calc_mean_coverage <- function (align) {
    seq_length <- apply (as.character(align), 1, function (x) length(x[grep("n|N|-", x, invert=TRUE)]))
    seq_coverage <- seq_length/ncol(align)
    mean_seq_coverage <- mean(unlist(seq_coverage))
    return(mean_seq_coverage)
  }
  
  # use above function to calculate MEAN percentage coverage for each alignment
  mean_coverage <- map_dbl(alignments, calc_mean_coverage)
  
  # calculate number of sequences (in ingroup) per alignment
  num_seqs <- unlist(lapply(alignments, nrow))
  
  # make list of data frames
  # results <- list(length, coverage, mean_coverage, num_seqs)
  results <- list(length = length, mean_coverage = mean_coverage, num_seqs = num_seqs, folder = rep(folder, length(alignment_files)))
  
  return(results)
  
}

# Genome masking ----------------------------------------------------------

# Functions to help standardize names in plan 
get_gff3 <- function (genome) {
  switch(genome,
         arabidopsis = here::here("data_raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz"),
         azolla = here::here("data_raw/Azolla_filiculoides.gene_models.highconfidence_v1.1.gff"),
         salvinia = here::here("data_raw/Salvinia_cucullata.gene_models.highconfidence_v1.2.gff"))
}

get_genome <- function (genome) {
  switch(genome,
         arabidopsis = here::here("data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.renamed.fasta"),
         azolla = here::here("data_raw/Azolla_filiculoides.genome_v1.2.fasta"),
         salvinia = here::here("data_raw/Salvinia_cucullata.genome_v1.2.fasta")
  )
}

get_sources <- function (genome) {
  switch(genome,
         arabidopsis = "araport11",
         azolla = "maker",
         salvinia = "maker"
  )
}

# Rename arabidopsis genomes so that 
# names match between the fasta file
# and the bed file
rename_arabidopsis_genome <- function (fasta_in, fasta_out, ...) {
  arabidopsis_genome <- read.FASTA(fasta_in)
  
  new_names <- case_when(
    grepl("1 dna", names(arabidopsis_genome)) ~ "1",
    grepl("2 dna", names(arabidopsis_genome)) ~ "2",
    grepl("3 dna", names(arabidopsis_genome)) ~ "3",
    grepl("4 dna", names(arabidopsis_genome)) ~ "4",
    grepl("5 dna", names(arabidopsis_genome)) ~ "5",
    grepl("Mt dna", names(arabidopsis_genome)) ~ "Mt",
    grepl("Pt dna", names(arabidopsis_genome)) ~ "Pt")
  
  if(anyNA(new_names)) {stop("Replacement names not provided for all fasta headers")}
  
  names(arabidopsis_genome) <- new_names
  
  write.FASTA(arabidopsis_genome, fasta_out)
}

# General-purpose ---------------------------------------------------------

#' Make a bash script
#'
#' The main reason for this is to call other
#' programs in different working directories
#' without changing R's wd. processx::run()
#' can also do this, but is very picky about
#' which programs it will execute.
#' 
#' @param script_file Path to write script file
#' @param command Command to run
#' @param arguments Arguments for command. If > 1 argument, the same
#' command will be run each time.
#' @param wd Working directory to run script
#' @param ... 
#'
make_bash_script <- function(script_file, command, arguments, wd = NULL, ...) {
  
  sink(file = script_file, type = "output")
  
  cat("#!/bin/bash\n")
  if (!(is.null(wd))) {cat(paste("cd", wd, "\n"))}
  
  # Assume we are running the same command, but may have multiple calls
  for (i in 1:length(arguments)) {
    cat(paste(command, paste(arguments[i], collapse = " "), "\n") )
  }
  
  sink()
}

# Run a bash script
run_bash_script <- function (script_file, ...) {
  system2("bash", script_file)
}

# function to run cat on files external to R
cat_ext_files <- function (input, output, wd = NULL, ...) {
  processx::run("sh", 
                c("-c",
                  glue::glue('cat {paste(input, collapse = " ")} > {output}'),
                  wd))
}

#' Write out a list of fasta files to a directory
#'
#' Optionally assign names to the files if fasta_list isn't already named
#' or if you want to over-write the original names.
#'
#' @param fasta_list List of DNA sequences to write out. Each item in list
#' must of class \code{\link[ape]{DNAbin}}.
#' @param out_dir Path to directory to write out DNA sequences.
#' @param fasta_names Optional character vector of file names to use when
#' writing out DNA sequences.
#' @param ext Optional character vector of length one; file extension to
#' append to DNA sequences (e.g., "fasta").
#' @param get_hash Logical; should the MD5 hash of `fasta_list` be returned? 
#' @param ... Additional other arguments. Not used by this function,
#' but meant to be used by \code{\link[drake]{drake_plan}} for tracking
#' during workflows.
#'
#' @return None (invisible ‘NULL’) or character vector if `hash` is `TRUE`.
#' Externally, fasta files will be written to `out_dir`.
#'
#' @examples
#' # Load some example DNA sequences.
#' library(ape)
#' data(woodmouse) 
#' 
#' # Make a temporary working directory to write out files.
#' temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))
#' 
#' # Make list of DNA samples.
#' dna_list <- list(a = woodmouse, b = woodmouse)
#' 
#' # Write out list:
#' # names as-is
#' write_fasta_files(dna_list, temp_dir)
#' list.files(temp_dir)
#' 
#' # add extension
#' write_fasta_files(dna_list, temp_dir, ext = "fasta")
#' list.files(temp_dir)
#' 
#' # new names and extension
#' write_fasta_files(
#'   dna_list, 
#'   temp_dir, 
#'   fasta_names = c("ho", "ge"), 
#'   ext = "fasta")
#' list.files(temp_dir)
#' 
#' # Cleanup
#' fs::file_delete(temp_dir)
#' @export
write_fasta_files <- function (fasta_list, out_dir, 
                               fasta_names = NULL, ext = NULL, 
                               get_hash = TRUE, 
                               ...) {
  
  require("magrittr") # for pipe
  
  # Check input
  assertthat::assert_that(is.logical(get_hash))
  
  assertthat::assert_that(is.list(fasta_list))
  assertthat::assert_that(
    unique(purrr::map_chr(fasta_list, class)) == "DNAbin",
    msg = "All items in fasta_list must of class 'DNAbin'"
  )
  
  assertthat::assert_that(assertthat::is.string(ext) | is.null(ext))
  if(is.null(ext)) ext <- ""
  
  assertthat::assert_that(assertthat::is.string(out_dir))
  out_dir <- fs::path_abs(out_dir)
  assertthat::assert_that(assertthat::is.dir(out_dir))
  
  assertthat::assert_that(!is.null(fasta_names) | !is.null(names(fasta_list)),
                          msg = "fasta list must have names or fasta_names must
                          be provided")
  
  # Set new names
  if (!is.null(fasta_names)) {
    assertthat::assert_that(
      length(fasta_names) == length(fasta_list),
      msg = "fasta_names and fasta_list must be of same length"
    )
    names(fasta_list) <- fasta_names
  }
  
  # Write out fasta seqs
  fasta_list %>%
    purrr::set_names(fs::path(out_dir, names(.), ext = ext)) %>%
    purrr::iwalk(ape::write.FASTA)
  
  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) unlist(fasta_list) %>% digest::digest()
  
}

# blast functions ---------------------------------------------------------

# Re-align filtered baits with top blast hits
# returns list of alignments
realign_with_best_hits <- function (blast_filtered_folder, 
                                    taxonomy_filtered_folder, ...) {
  
  blast_filtered_folder <- jntools::add_slash(blast_filtered_folder)
  taxonomy_filtered_folder <- jntools::add_slash(taxonomy_filtered_folder)

  # Get names of fasta files of top blast hits (already 1 per alignment)
  blast_top_match_names <- 
    list.files(blast_filtered_folder, pattern = "bestmatch")
  
  # Make vector of shortened names (alignment number and ortho pruning method)
  # to use for alignments
  blast_filtered_alignment_names <- 
    gsub("\\.fa.*$", "", blast_top_match_names)
  
  # Read in blast-filtered alignments (those with blast hits).
  blast_filtered_alignments <- 
    blast_top_match_names %>%
    gsub(".outfmt6.bestmatch.fasta", "", ., fixed = TRUE) %>%
    paste(taxonomy_filtered_folder, ., sep ="") %>%
    map(ape::read.FASTA)
  
  # Read in seqs of top blast hits
  blast_top_matches <- 
    paste0(blast_filtered_folder, blast_top_match_names, sep ="") %>%
    map(ape::read.FASTA)
  
  # combine and re-align blast-filtered alignments with their top matches
  map2(blast_filtered_alignments, blast_top_matches, c) %>%
    map(ips::mafft, path = "/usr/bin/mafft", options = "--adjustdirection") %>%
    rlang::set_names(blast_filtered_alignment_names)
}

# Run blast on all files with the same ending in a folder
# Outputs to same folder containing input files (but we could change this)
blastn_list <- function (fasta_folder, fasta_ending, database, other_args = NULL, overwrite = FALSE, get_hash = TRUE, ...) {
  fasta_folder <- jntools::add_slash(fasta_folder)
  
  search_terms <- "\\.outfmt6$"
  
  # optional: delete all previous output written in this folder
  if (isTRUE(overwrite)) {
    files_to_delete <- list.files(fasta_folder, pattern = search_terms)
    if (length(files_to_delete) > 0) {
      files_to_delete <- paste0(fasta_folder, files_to_delete)
      file.remove(files_to_delete)
    }
  }
  
  # list all fasta files in folder
  fasta_files <- list.files(fasta_folder, pattern = fasta_ending)
  fasta_list <- paste0(fasta_folder, fasta_files)
  
  # make list of results to write
  results_list <- paste0(fasta_folder, fasta_files, ".outfmt6")
  
  # loop over lists
  walk2(fasta_list, results_list, blast_n, database = database, other_args = other_args)
  
  # optional: return hash of results
  if (isTRUE(get_hash)) {
    output <- list.files(fasta_folder, pattern = search_terms)
    output <- if (length(output) > 0) {unlist(lapply(paste0(fasta_folder, output), readr::read_file))} else {output}
    hash <- digest::digest(output)
    return(hash)
  }
}

# Extract blast hits
extract_blast_hits <- function (blast_results_folder, blast_results_ending, blast_db, out_dir, ...) {
  
  out_dir <- jntools::add_slash(out_dir)
  
  # Make list of blastdbcmd calls to extract top blast hits
  calls <- 
    # read in and filter blast results
    filter_blast_results(
      blast_results_folder = blast_results_folder,
      blast_results_ending = blast_results_ending) %>%
    # only keep best blast hit from each result
    map_df(slice, 1) %>%
    # format the command to extract the blast hit as a column
    mutate(command = glue::glue('-db {blast_db} -entry "{sseqid}" > {out_dir}{file_name}.bestmatch.fasta')) %>%
    pull(command)
  
  make_bash_script(script_file = "extract_blast.sh",
                   command = "blastdbcmd",
                   arguments = calls)
  
  run_bash_script("extract_blast.sh")
  file.remove("extract_blast.sh")
}

# Append the file name as a column when reading in a tsv
read_tsv_with_file_name <- function (x, ...) {
  file <- str_split(x, "/") %>% map_chr(length(.[[1]]))
  read_tsv(x, ...) %>%
    mutate(file_name = file)
}

# Filter blast results
# Returns a logical vector of all blast results that passed the filter
filter_blast_results <- function (blast_results_folder, blast_results_ending) {
  
  outfmt6_headers <- c("qseqid", "qlen", "sseqid", "slen", "frames", "pident", "nident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  
  blast_results_folder <- jntools::add_slash(blast_results_folder)
  
  # read all blast results
  blast_results <- 
    list.files(blast_results_folder, pattern = blast_results_ending) %>%
    paste0(blast_results_folder, .) %>%
    map(read_tsv_with_file_name, col_names = outfmt6_headers) %>%
    map(~arrange(.x, evalue, desc(bitscore)))
  
  # filter1: keep only blast results with at least one hit
  filter1 <- blast_results %>%
    map_dbl(nrow) %>%
    map_lgl(~ .x > 0)
  
  # apply filter
  blast_results[filter1]
  
}

# Fill introns functions --------------------------------------------------

# trim columns of exclusively "n" characters from one end of alignment
trim_one_end <- function (alignment) {
  # prepare loop: to.trim is list of TRUE/FALSE for sequences to trim
  # at_ends starts as TRUE, becomes false at first instance of non-"n" character
  to.trim <- list()
  at_ends <- TRUE
  
  # start at one end of alignment, assign TRUE to each column if still at "end" and all columns are "n"
  for (i in 1:ncol(alignment)) {
    if (all(grepl("n", alignment[, i])) && at_ends == TRUE) {
      to.trim[i] <- TRUE
    } else {
      to.trim[i] <- FALSE
      at_ends <- FALSE
    }
  }
  to.trim <- unlist(to.trim)
  alignment <- alignment[, !(to.trim)]
  return(alignment)
}

# custom func that uses trim_one_end on both ends of alignment 
trim_both_ends <- function (alignment) {
  require(ape)
  alignment <- trim_one_end(alignment)
  alignment <- ape::complement(alignment)
  alignment <- trim_one_end(alignment)
  alignment <- ape::complement(alignment)
  return(alignment)
}

# fill introns in alignment based on reference sequence,
# then trim that reference sequence and the outgroup
# from the alignment
fill_introns <- function (alignment, outgroup) {
  
  # extract the "reference genome" sequence, and remove from alignment
  # all transcriptome sequences are 4-letter codes.
  # Sacu = Salvinia, Azfi = Azolla, numbers 1-5 are Arabidopsis (numbered by chromosome)
  ref_seq <- alignment[grep("Sacu|Azfi|1|2|3|4|5", rownames(alignment)), ]
  alignment <- alignment[grep("Sacu|Azfi|1|2|3|4|5", rownames(alignment), invert=TRUE), ]
  
  # rev-comp alignments if MAFFT flipped all of them
  if (all(grepl("_R_", rownames(alignment)))) {
    alignment <- ape::complement(alignment)
    rownames(alignment) <- gsub("_R_", "", rownames(alignment))
  }
  
  # classify each character in ref-genome sequence as hard-masked intron ("n") or not
  introns <- apply(as.character(ref_seq), 2, function (x) x=="n"|x=="N")
  
  # also make list of all transcript characters that are gap-only
  gaps <- apply(as.character(alignment), 2, function (x) all(x == "-"))
  
  # only convert introns if all other characters are gaps
  introns[which(gaps==FALSE)] <- FALSE  
  
  # mask introns (convert to character first)
  alignment <- as.character(alignment)
  alignment[,introns] <- "n"
  
  # remove outgroup taxa
  alignment <- alignment[!(rownames(alignment) %in% outgroup), ]
  
  # convert back to DNAbin for deleteEmptyCells
  alignment <- as.DNAbin(alignment)
  
  # delete columns that are only empty cells
  alignment <- ips::deleteEmptyCells(alignment, nset="-")
  
  # trim ends if they are only introns (all columns "n")
  alignment <- trim_both_ends(alignment)
  
  return(alignment)
}

# loop-ize fill_introns for drake workflow
fill_introns_loop <- function (alignment_list, outgroup, ...) {
  map(alignment_list, fill_introns, outgroup = outgroup)
}

# Create a dataframe of alignment stats, including
# the alignments as a list-column
assemble_bait_data <- function (alignment_list) {
  map_df(alignment_list, baitfindR::calculate_alignment_stats, include_aln = TRUE) %>%
  mutate(bait_id = names(alignment_list))
}

# Filter alignments to get final baits.
# Here we filter by min. number of introns > 1,
# and choosing top five ranked by % pars. inform. characters
filter_alignments <- function (alignment_data) {
  alignment_data %>%
    filter(num_introns > 1) %>%
    arrange(desc(pars_inf)) %>%
    slice(1:5)
}