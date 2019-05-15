# Data import -------------------------------------------------------------

# Make a table of oneKp sample codes and URLs to download transcriptome assemblies.
make_download_table <- function(url = "http://www.onekp.com/public_data.html") {
  
  # Scrape oneKp links for assemblies.
  # Use Selectorgadget to identify the CSS component with links.
  # Note: run vignette(“selectorgadget”) for a refresher on how to do this.
  links <- xml2::read_html(url) %>% rvest::html_nodes("td:nth-child(5) a")
  
  # Next, read in the main table (just plain text, without links).
  tbls_xml <- XML::readHTMLTable(url)
  
  # Extract the table and add links to assemblies,
  # keeping only transcriptomes in baitfindR example dataset .
  tbls_xml[[1]] %>% 
    tibble::as_tibble() %>%
    dplyr::mutate_all(as.character) %>%
    dplyr::mutate(assembly_link = map_chr(links, xml2::xml_attrs)) %>%
    dplyr::select(code = `1kP_Code`, assembly_link)
}

# Downsize a transcriptome file
downsize_transcriptome <- function (file, keep_frac) {
  # Read in file. Use bzfile() because it's compressed
  seq <- ape::read.FASTA(bzfile(file))
  # Randomly downsize to specified fraction of transcripts
  seq <- seq[sample(1:length(seq), keep_frac*length(seq))]
}

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

# blast functions ---------------------------------------------------------

# Re-align filtered baits with top blast hits
# returns list of alignments
realign_with_best_hits <- function (blast_filtered_folder, 
                                    taxonomy_filtered_folder, ...) {
  
  blast_filtered_folder <- fs::path_abs(blast_filtered_folder)
  taxonomy_filtered_folder <- fs::path_abs(taxonomy_filtered_folder)

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
    stringr::str_remove(".outfmt6.bestmatch.fasta") %>%
    fs::path(taxonomy_filtered_folder, .) %>%
    purrr::map(ape::read.FASTA)
  
  # Read in seqs of top blast hits
  blast_top_matches <- 
    fs::path(blast_filtered_folder, blast_top_match_names, sep ="") %>%
    purrr::map(ape::read.FASTA)
  
  # combine and re-align blast-filtered alignments with their top matches
  purrr::map2(blast_filtered_alignments, blast_top_matches, c) %>%
    purrr::map(ips::mafft, path = "/usr/bin/mafft", options = "--adjustdirection") %>%
    rlang::set_names(blast_filtered_alignment_names)
}

# Filter blast results
# Returns a logical vector of all blast results that passed the filter
filter_blast_results <- function (blast_results_folder, blast_results_ending) {
  
  outfmt6_headers <- c("qseqid", "qlen", "sseqid", "slen", "frames", "pident", "nident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  
  blast_results_folder <- fs::path_abs(blast_results_folder)
  
  # read all blast results
  blast_results <- 
    list.files(blast_results_folder, pattern = blast_results_ending, full.names = TRUE) %>%
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