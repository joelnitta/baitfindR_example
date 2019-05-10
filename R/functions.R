# Data import -------------------------------------------------------------

# Randomly downsize transcriptome to make smaller dataset for testing
downsize_transcriptome <- function (file, keep_frac) {
  # Read in file. Use bzfile() because it's compressed
  seq <- read.dna(bzfile(file), format="fasta")
  # Randomly downsize to specified fraction of transcripts
  seq <- seq[sample(1:length(seq), keep_frac*length(seq))]
}


#' Downsize a file from the top
#' 
#' Only the lines specified from the first to the fraction
#' of the file specified by `keep_frac` will be kept.
#'
#' @param full_file Path to input file
#' @param path Path to write out downsized file
#' @param keep_frac Fraction of file to retain
#'
#' @return NULL
downsize_simple <- function(full_file, path, keep_frac) {
  
  full <- readr::read_lines(full_file)
  short <- full[1:round(length(full)*keep_frac, 0)]
  readr::write_lines(short, path)
  
}


# Untar with tracking
untar_tracked <- function (tarfile, compressed, exdir, ...) {
  untar(tarfile = tarfile, compressed = compressed, exdir = exdir)
}

# For some reason there are a bunch of duplicate names with
# different sequences in the lygodium data. Exclude these.
trim_proteome <- function (proteome) {
  dup_names <- names(proteome)[duplicated(names(proteome))]
  trimmed_proteome <- proteome[!(names(proteome) %in% dup_names)]
  return(trimmed_proteome)
}

# count_seqs
# Count the number of sequences per alignment in a list of aligments
# Returns a dataframe of counts: the frequency of alignments containing 
# that many sequences. "method" can be added to specify the method used
# to generate the alignments. Optionally include an ingroup to only count
# the number of ingroup sequences.
count_seqs <- function (alignments, ingroup = NULL, method = NULL) {
  
  # so we can feed method into dplyr
  method <- enquo(method)
  
  # define function to only count number of ingroup sequences
  count_ingroup_seqs <- function (alignments, ingroup) {
    purrr::map_dbl(alignments, function (x) {
      x <- as.list(x)
      length(names(x)[names(x) %in% ingroup])
    }
    )
  }
  
  # get the raw count of sequences
  if (is.null(ingroup)) {
    raw_count <- map_dbl(alignments, length)
  } else if (is.character(ingroup) && length(ingroup) > 0) {
    raw_count <- count_ingroup_seqs(alignments = alignments, ingroup = ingroup)
  } else {
    stop ("ingroup must be null or a character vector")
  }
  
  # convert to tibble
  raw_count %>%
    as_tibble %>% 
    select(num_ingroup_seqs = value) %>% 
    group_by(num_ingroup_seqs) %>% 
    summarize(freq=n()) %>% 
    mutate(method = !! method)
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

# Extract masked genes from a genome
extract_masked_genes <- function (genes_file, masked_genome_file, masked_genes_file, wd = here::here(), ...) {
  
  genes <- read_tsv(genes_file, col_names = c("chr", "start", "end"), col_types = "cdd")
  
  if (!(check_bed_genome_names(masked_genome_file, genes))) {
    stop ("Names don't match between bed file and genome fasta headers")
  }
  
  wd <- jntools::add_slash(wd)
  
  script_file <- glue::glue('{wd}{digest::digest(genes)}.sh')
  
  make_bash_script(
    script = script_file,
    command = "bedtools getfasta",
    arguments = glue::glue("-fi {masked_genome_file} -bed {genes_file} -fo {masked_genes_file} -s"),
    # The -s preserves the strandedness of the sequence. 
    # The reverse-complement is written for the + strand. 
    wd = wd
  )
  
  run_bash_script(
    script = script_file
  )
  
  # Remove temporary script files
  if (file.exists(script_file)) {file.remove(script_file)}
}

#' Mask introns in a genome
#'
#' @param introns_file Path to bed file with locations of introns in genome.
#' (bed file is a tab-separated file with columns for "chr" (chromosome), "start",
#' and "end", in that order).
#' @param genome_file Path to full genome file in fasta format.
#' @param masked_genome_file Path to write genome file with introns masked.
#' @param wd Working directory. A temporary bash script will be written and 
#' executed here.
#' @param ... Other arguments not used by this function but meant for tracking
#' with drake.
#'
#' @return
#' @export
#'
#' @examples
#' 
#' # First write genes, introns, and exons out as tsv files
#' 
#' dir.create("temp_dir")
#' find_bed_regions(
#'   gff3_file = "data_raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz",
#'   source_select = "araport11",
#'   out_type = "write_all",
#'   out_dir = "temp_dir",
#'   prefix = "test"
#' )
#' 
#' # Now mask the genome, using the intron and genome files.
#' mask_genome(
#'   introns_file = "temp_dir/test_introns",
#'   genome_file = "data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.renamed.fasta",
#'   masked_genome_file = "temp_dir/test_masked"
#' )
mask_genome <- function (introns_file, genome_file, masked_genome_file, wd = here::here(), ...) {
  
  # Check input
  assertthat::assert_that(assertthat::is.string(introns_file))
  introns_file <- fs::path_abs(introns_file)
  assertthat::assert_that(assertthat::is.readable(introns_file))
  
  assertthat::assert_that(assertthat::is.string(genome_file))
  genome_file <- fs::path_abs(genome_file)
  assertthat::assert_that(assertthat::is.readable(genome_file))
  
  introns <- readr::read_tsv(introns_file, col_names = c("chr", "start", "end"), col_types = "cdd")
  
  checkr::check_data(
    introns,
    values = list(
      chr = "a",
      start = 1,
      end = 1
    ),
    order = TRUE)
  
  if (!(check_bed_genome_names(genome_file, introns))) {
    stop ("Names don't match between bed file and genome fasta headers")
  }
  
  wd <- fs::path_abs(wd)

  # Make path for script file based on digest of introns.
  script_file <- fs::path(wd, digest::digest(introns), ext = "sh")
  
  # Write out script file.
  make_bash_script(
    script = script_file,
    command = "bedtools maskfasta",
    arguments = glue::glue("-fi {genome_file} -bed {introns_file} -fo {masked_genome_file}"),
    wd = wd
  )
  
  run_bash_script(
    script = script_file
  )
  
  # Remove temporary script files
  if (file.exists(script_file)) {file.remove(script_file)}
}

#' Clean up data from a gff file and
#' convert to bed format
#' 
#' Helper function for `find_bed_regions`. Merges overlapping regions and sorts
#' regions.
#'
#' @param region Dataframe; list of gene regions in "bed" format. Must include 
#' the following columns in order: `chr` ('chromosome', character), `start` 
#' (start position, numeric), and `end` (end position, numeric).
#' @param check.chr Logical; should coordinates be checked for chromosomal
#' format with "chr" prefix?
#' @param verbose Logical; should `bedr` functions output all messages?
#'
#' @return Dataframe in "bed" format.
#'
clean_gff <- function (region, check.chr = FALSE, verbose = FALSE) {
  
  # Check input format
  assertthat::assert_that(
    bedr:::determine.input(region, verbose = verbose) == 1,
    msg = "Input must be dataframe in 'bed' format")
  
  region %>%
    # Check if region is valid
    dplyr::filter(bedr::is.valid.region(., check.chr = check.chr, verbose = verbose)) %>%
    # Collapse overlapping regions
    bedr::bedr.merge.region(check.chr = check.chr, verbose = verbose) %>%
    # Sort regions
    bedr::bedr.sort.region(check.chr = check.chr, verbose = verbose) %>% 
    # Convert to bed format
    bedr::convert2bed(check.chr = check.chr, verbose = verbose)
}

#' Find genes, exons, and introns in a gff3 file
#' 
#' If tsv files are written out by selecting "write_all" for `out_type`,
#' they will overwrite any existing files with the same name in `out_dir`.
#'
#' @param gff3_file Path to input file in `gff3` format.
#' @param source_select Character vector; only use regions from these
#' sources. Must match values in `source` column of gff3 file. Optional.
#' @param gene_label String; value used to indicate genes in gff3 file.
#' Must match at least one value in `type` column of gff3 file. Default "gene".
#' @param exon_label String; value used to indicate exons in gff3 file.
#' Must match at least one value in `type` column of gff3 file. Default "exon".
#' @param verbose Logical; should `bedr` functions output all messages?
#' @param out_type Type of output to return:
#' "genes": dataframe in "bed" format of genes.
#' "introns": dataframe in "bed" format of introns.
#' "exons": dataframe in "bed" format of exons.
#' "write_all": write tab-separated files for each of `genes`, `introns`, and
#' `exons` to `out_dir`. The hash digest of the combined genes, introns, and 
#' exons will be returned.
#' @param prefix String; prefix to attach to tsv files if `out_type` is 
#' "write_all".
#' @param out_dir Directory to write tsv files if `out_type` is "write_all".
#' @param ... Other arguments not used by this function but meant for tracking
#' with drake.
#'
#' @return Dataframe or character.
#'
#' @examples
#' # Find genes
#' genes <- find_bed_regions(
#'   gff3_file = "/home/rstudio/baitfindR_simple/data_raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz",
#'   source_select = "araport11",
#'   out_type = "genes"
#' )
#' tibble::as_tibble(genes)
#' 
#' # Find introns
#' introns <- find_bed_regions(
#'   gff3_file = "/home/rstudio/baitfindR_simple/data_raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz",
#'   source_select = "araport11",
#'   out_type = "introns"
#' )
#' tibble::as_tibble(introns)
#' 
#' # Find exons
#' exons <- find_bed_regions(
#'   gff3_file = "/home/rstudio/baitfindR_simple/data_raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz",
#'   source_select = "araport11",
#'   out_type = "exons"
#' )
#' tibble::as_tibble(exons)
#' 
#' # Write genes, introns, and exons out as tsv files
#' dir.create("temp_dir")
#' find_bed_regions(
#'   gff3_file = "/home/rstudio/baitfindR_simple/data_raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz",
#'   source_select = "araport11",
#'   out_type = "write_all",
#'   out_dir = "temp_dir",
#'   prefix = "test"
#' )
#' @export
find_bed_regions <- function (gff3_file, 
                              source_select = NULL, 
                              gene_label = "gene", exon_label = "exon",
                              verbose = FALSE,
                              prefix = NULL, out_dir = NULL, 
                              out_type = c("genes", "introns", "exons", "write_all"), 
                              ...) {
  
  # Check input
  assertthat::assert_that(assertthat::is.readable(gff3_file))
  assertthat::assert_that(assertthat::is.string(gene_label))
  assertthat::assert_that(assertthat::is.string(exon_label))
  assertthat::assert_that(assertthat::is.string(out_type))
  assertthat::assert_that(is.logical(verbose))
  assertthat::assert_that(out_type %in% c("genes", "introns", "exons", "write_all"),
                          msg = "'out_type' must be one of 'genes', 'introns', 'exons', or 'write_all'")
  
  
  # Read in gff3 file as dataframe
  gff3 <- ape::read.gff(gff3_file) %>%
    dplyr::mutate(chr = as.character(seqid))
  
  # Keep only annotations from selected source
  if (!is.null(source_select)) {
    assertthat::assert_that(is.character(source_select))
    assertthat::assert_that(all(source_select %in% gff3$source))
    gff3 <- dplyr::filter(gff3, source %in% source_select)
  }
  
  # Extract and clean up genes
  genes <- gff3 %>% dplyr::filter(type == gene_label) %>% 
    dplyr::select(chr, start, end) %>%
    clean_gff(verbose = verbose)
  
  # Extract and clean up exons
  exons <- gff3 %>% dplyr::filter(type == exon_label) %>% 
    dplyr::select(chr, start, end) %>%
    clean_gff(verbose = verbose)
  
  # Introns are genes - exons
  introns <- bedr::bedr.subtract.region(
    genes,
    exons, 
    remove.whole.feature = FALSE, 
    check.chr = FALSE,
    verbose = verbose)
  
  # Write out all regions and return hash of genes + exons + introns
  if (out_type == "write_all") {
    out_dir <- fs::path_abs(out_dir)
    assertthat::assert_that(assertthat::is.writeable(out_dir))
    assertthat::assert_that(assertthat::is.string(prefix))
    all_regions <- list(genes = genes,
                        exons = exons,
                        introns = introns) 
    all_regions %>%
      purrr::set_names(fs::path(out_dir, paste0(prefix, "_", names(.)))) %>%
      purrr::iwalk(readr::write_tsv, col_names = FALSE)
    return(digest::digest(all_regions))
  }
  
  # Or, return a particular result type.
  return(switch(
    out_type,
    genes = genes,
    exons = exons,
    introns = introns
  ))
  
}

# Check that genome fasta headers and 
# "chromosome" names in bed file match.
# (This must be true for 
# bedtools maskfasta to work).
check_bed_genome_names <- function (fasta_file, bed) {
  # find all sequence headers in fasta file
  seq_names <- 
    readr::read_lines(fasta_file) %>%
    purrr::keep(~ grepl(">", .x)) %>%
    purrr::map(~ gsub(">", "", .x))
  # make sure "chr" (chromosome) names of bed file are all in
  # fasta sequence headers
  chr <- unique(bed$chr)
  all(chr %in% seq_names)
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

# Write out a list of fasta files to a directory.
# Optionally assign names to the files fasta_list isn't already named
# or if you want to over-write the original names
write_fasta_files <- function (fasta_list, fasta_names = NULL, out_dir, 
                               suffix = NULL, prefix = NULL, get_hash = TRUE, 
                               ...) {
  
  out_dir <- jntools::add_slash(out_dir)
  
  if (!is.null(fasta_names)) names(fasta_list) <- fasta_names
  
  assertthat::assert_that(!is.null(fasta_names) | !is.null(names(fasta_list)),
                          msg = "fasta list must have names or fasta_names must
                          be provided")
  
  fasta_list %>%
    rlang::set_names(paste0(out_dir, prefix, names(.), suffix)) %>%
    purrr::iwalk(ape::write.FASTA)
  
  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) {
    output <- unlist(fasta_list)
    hash <- digest::digest(output)
    return(hash)
  }
  
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