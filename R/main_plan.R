# main_plan

# 01_translation ----------------------------------------------------------

# Build blast protein database including combined Arabidopsis, Azolla, 
# Lygodium, and Salvinia proteomes.
# Blast databases comprise multiple files, so arbitrarily choose .phr 
# as output file for tracking.
build_blastp_db <- drake_plan(
  blastp_database = baitfindR::build_blast_db(
    in_seqs = file_in("data/combined_proteomes.fasta"), 
    db_type = "prot",
    out_name = "combined_proteomes",
    parse_seqids = TRUE,
    wd = "01_translation",
    outfile = file_out("01_translation/combined_proteomes.phr")
  )
)

# Run transdecoder, shorten output names, and prepare for clustering
run_transdecoder <- drake_plan(
  
  # Find candidate long ORFs (long open reading frames) in transcriptomes.
  # TransDecoder.LongOrfs outputs multiple files to .transdecoder_dir 
  # folder.
  # Use longest_orfs.pep as representative output file for tracking.
  longest_orfs = baitfindR::transdecoder_long_orfs(
    transcriptome_file = file_in("data/transcriptome__"),
    wd = "01_translation",
    other_args = "-S",
    depends = file_in("01_translation/combined_proteomes.phr"),
    outfile = file_out("01_translation/transcriptome__.transdecoder_dir/longest_orfs.pep")
  ),
  
  # Blast candidate long ORFs against the protein database
  orf_blast = baitfindR::blast_p(
    query = file_in("01_translation/transcriptome__.transdecoder_dir/longest_orfs.pep"), 
    database = "combined_proteomes",
    out_file = "transcriptome__.blastp.outfmt6",
    other_args = c("-max_target_seqs", 1, "-evalue", 10,  "-num_threads", 2),
    wd = "01_translation",
    produces = file_out("01_translation/transcriptome__.blastp.outfmt6")
  ),
  
  # Output the final ORFs, preferentially retaining those with blast hits
  orf_predict = baitfindR::transdecoder_predict(
    transcriptome_file = file_in("data/transcriptome__"),
    blast_result = file_in("01_translation/transcriptome__.blastp.outfmt6"),
    wd = "01_translation",
    infile = file_in("01_translation/transcriptome__.transdecoder_dir/longest_orfs.pep"),
    outfile = file_out("01_translation/transcriptome__.transdecoder.cds")
  ),
  
  # Shorten names from transdecoder using R function
  # Y&S function works on whole folder, so not conducive to drake.
  # baitfindR function works on single file (transcriptome) at a time.
  seq_with_short_names = baitfindR::fix_names_from_transdecoder(
    transdecoder_output = file_in("01_translation/transcriptome__.transdecoder.cds")
  ),
  
  # Write out fasta files with shortened names
  seq_with_short_names_out = ape::write.FASTA(
    x = seq_with_short_names_transcriptome__, 
    file = file_out("01_translation/transcriptome__.cds.fa")
  ),
  
  # Reduce redundancy for cds
  # Note that output is put into clustering/
  reduced_cdhitest = baitfindR::cd_hit_est(
    input = file_in("01_translation/transcriptome__.cds.fa"), 
    output = file_out("02_clustering/transcriptome__.cds.fa.cdhitest"), 
    other_args = c("-c", "0.99", "-n", "10", "-r", "0", "-T", "2")
  )
) %>%
  evaluate_plan(rules = list(transcriptome__ = codes))

# 02_clustering, 03_clusters ----------------------------------------------

# Concatenate all the .cdhitest files into a new file all_orfs.fa
concatenate_cdhitest <- 
  drake_plan(
    read_reduced_cdhitest = readr::read_file(file_in("02_clustering/transcriptome__.cds.fa.cdhitest"))
  ) %>%
  evaluate_plan(rules = list(transcriptome__ = codes)) %>%
  bind_plans(gather_plan(., target = "cdhitest_files")) %>%
  bind_plans(
    drake_plan(
      concatenated_cdhitest = baitfindR::cat_files(
        cdhitest_files,
        output_file=file_out("02_clustering/all_orfs.fa")
      )
    )
  )

# Make database from ORFs, and run blast on one transcriptome at a time.
run_allbyall_blast <- drake_plan(
  
  # Make all-by-all database including all ORFs
  all_blast_db = baitfindR::build_blast_db(
    in_seqs = file_in("02_clustering/all_orfs.fa"),
    db_type = "nucl",
    out_name = "all_orfs",
    wd = "02_clustering",
    outfile = file_out("02_clustering/all_orfs.nhr")
  ),
  
  # Blast transcriptomes against database
  blast_transcriptome = baitfindR::blast_n(
    query = file_in("02_clustering/transcriptome__.cds.fa.cdhitest"),
    database = "all_orfs",
    out_file = "transcriptome__.allbyall.blast.outfmt6",
    other_args = c("-evalue", 10, "-num_threads", 1, "-max_target_seqs", 1000),
    outfmt = "6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore",
    wd = "02_clustering",
    depends = file_in("02_clustering/all_orfs.nhr"),
    produces = file_out("02_clustering/transcriptome__.allbyall.blast.outfmt6")
  )
) %>% 
  evaluate_plan(rules = list(transcriptome__ = codes))

# Concatenate output of single blast results
concatenate_allbyall_blast <- 
  drake_plan(
    read_blast = readr::read_file(file_in("02_clustering/transcriptome__.allbyall.blast.outfmt6"))
  ) %>%
  evaluate_plan(rules = list(transcriptome__ = codes)) %>%
  bind_plans(gather_plan(., target = "cdhitestblast_files")) %>%
  bind_plans(
    drake_plan(
      combined_cdhitest = baitfindR::cat_files(
        cdhitestblast_files, 
        output_file=file_out("02_clustering/all.rawblast")
      )
    )
  )

# Filter raw blast output by hit fraction, run mcl, 
# and write out as fasta files.
run_mcl <- drake_plan (
  
  # Prepare all by all blast results for mcl
  blast_formatted_for_mcl = baitfindR::blast_to_mcl(
    blast_results = file_in("02_clustering/all.rawblast"),
    hit_fraction_cutoff = my_hit_frac,
    outfile = file_out("02_clustering/all.rawblast.hit-frac`my_hit_frac`.minusLogEvalue")
  ),
  
  # Run mcl
  mcl_clusters = baitfindR::mcl(
    mcl_input = file_in("02_clustering/all.rawblast.hit-frac`my_hit_frac`.minusLogEvalue"),
    i_value = my_i_value,
    e_value = 5,
    other_args = c("--abc", "-te", "2"),
    mcl_output = file_out("02_clustering/hit-frac`my_hit_frac`_I`my_i_value`_e5")
  ),
  
  # Write fasta files for each cluster from mcl output to 03_clusters/
  fasta_clusters = baitfindR::write_fasta_files_from_mcl(
    all_fasta = file_in("02_clustering/all_orfs.fa"),
    mcl_outfile = file_in("02_clustering/hit-frac`my_hit_frac`_I`my_i_value`_e5"),
    minimal_taxa = 4,
    outdir = "03_clusters",
    get_hash = TRUE,
    overwrite = TRUE
  ), 
  
  # Align each cluster, trim alignment, and infer a tree
  # ("basic trees" which will be further pruned downstream)
  basic_trees = baitfindR::fasta_to_tree(
    overwrite = TRUE,
    seq_folder = "03_clusters/",
    number_cores = 1,
    seq_type = "dna",
    bootstrap = FALSE,
    get_hash = TRUE,
    infile = fasta_clusters)
) %>%
  evaluate_plan(rules = list("`my_hit_frac`" = my_hit_frac, 
                             "`my_i_value`" = my_i_value),
                rename = FALSE)

# 04_homologs, 05_orthologs -----------------------------------------------

# Identify homologs and orthologs, write out orthologs (baits)
sort_homologs_orthologs <- drake_plan(
  
  # Trip extremely long tips from tree
  trimmed_trees = baitfindR::trim_tips(
    tree_folder = "03_clusters", 
    tree_file_ending = ".tre",
    relative_cutoff = 0.2,
    absolute_cutoff = 0.4, 
    overwrite = TRUE,
    depends = basic_trees), 
  
  # Retain only one tip per taxon
  masked_trees = baitfindR::mask_tips_by_taxonID_transcripts(
    tree_folder = "03_clusters",
    aln_folder = "03_clusters",
    depends = trimmed_trees),
  
  # Break up homolog tree into orthologs by cutting 
  # at long internal branches
  homolog_trees = baitfindR::cut_long_internal_branches( 
    tree_folder = "03_clusters/",
    tree_file_ending = ".mm",
    internal_branch_length_cutoff = 0.3,
    minimal_taxa = 4,
    outdir = "04_homologs/",
    depends = masked_trees),
  
  # Write out sequences from trees
  homolog_seqs = baitfindR::write_fasta_files_from_trees(
    all_fasta = "02_clustering/all_orfs.fa",
    tree_folder = "04_homologs",
    tree_file_ending = ".subtree",
    outdir = "04_homologs",
    depends = homolog_trees),
  
  # Reassemble these into trees
  homolog_clean_trees = baitfindR::fasta_to_tree(
    seq_folder = "04_homologs",
    number_cores = 1,
    seq_type = "dna",
    bootstrap = FALSE,
    overwrite = TRUE,
    depends = homolog_seqs),
  
  # Prune orthologs using "1-to-1" method
  ortho_121_trees = baitfindR::filter_1to1_orthologs(
    tree_folder = "04_homologs", 
    tree_file_ending = ".tre",
    minimal_taxa = 3,
    outdir = "05_orthologs/tre",
    overwrite = TRUE,
    depends = homolog_clean_trees), 
  
  # Write out orthologs
  ortholog_fasta = baitfindR::write_ortholog_fasta_files( 
    all_fasta = file_in("02_clustering/all_orfs.fa"),
    tree_folder = "05_orthologs/tre",
    outdir = "05_orthologs/fasta",
    minimal_taxa = 3,
    overwrite = TRUE,
    depends = ortho_121_trees)
)

# 06_intron_masking -------------------------------------------------------

# Generate fasta files of intron-masked genes 
# from reference genomes.
mask_genes <- drake_plan (
  
  # Find introns in reference genomes
  introns = baitfindR::find_bed_regions (
    gff3_file = get_gff3("genome_list"),
    genome_file = get_genome("genome_list"),
    source_select = get_sources("genome_list"),
    prefix = "genome_list",
    out_type = "write_all",
    out_dir = "06_intron_masking",
    depends = ortholog_fasta,
    output = c(
      file_out("06_intron_masking/genome_list_genes"),
      file_out("06_intron_masking/genome_list_introns"),
      file_out("06_intron_masking/genome_list_exons")
    )
  ),
  
  # Mask introns in reference genomes
  masked_genome = baitfindR::mask_regions_in_fasta(
    bed_file = file_in("06_intron_masking/genome_list_introns"),
    fasta_file = get_genome("genome_list"),
    out_fasta_file = file_out("06_intron_masking/genome_list_masked")
  ),
  
  # Extract intron-masked genes from genomes
  masked_genes = baitfindR::extract_regions_from_fasta(
    bed_file = file_in("06_intron_masking/genome_list_genes"),
    fasta_file = file_in("06_intron_masking/genome_list_masked"),
    out_fasta_file = file_out("06_intron_masking/genome_list_masked_genes")
  )
) %>%
  evaluate_plan(
    wildcard = "genome_list",
    values = genomes)

# Concatenate all the intron-masked genes into a single file
concatenate_masked_genes <- 
  drake_plan(
    read = readr::read_file(file_in("06_intron_masking/genome_list_masked_genes"))
  ) %>%
  evaluate_plan(
    wildcard = "genome_list",
    values = genomes) %>%
  bind_plans(gather_plan(., target = "masked_gene_files")) %>%
  bind_plans(
    drake_plan(
      concatenate = baitfindR::cat_files(
        masked_gene_files, 
        output_file=file_out("06_intron_masking/all_masked_genes")
      )
    )
  )

# Make blast-db of masked genes
make_masked_genes_blast_db <- drake_plan(
  masked_db = baitfindR::build_blast_db(
    in_seqs = file_in("06_intron_masking/all_masked_genes"),
    db_type = "nucl",
    out_name = "masked_genes",
    parse_seqids = TRUE,
    wd = "06_intron_masking",
    outfile = file_out("06_intron_masking/masked_genes.nhr"))
)

# Mask introns in baits, filter final baits
mask_and_filter_baits <- drake_plan (
  
  # Filter baits by family
  # (each alignment must contain at least
  # one sample per ingroup family)
  family_filtered_baits = baitfindR::filter_fasta(
    seq_folder = "05_orthologs/fasta",
    taxonomy_data = taxonomy_data,
    filter_col = "family",
    sample_col = "code",
    depends1 = file_in("06_intron_masking/masked_genes.nhr"),
    depends2 = ortholog_fasta),
  
  # Write out unaligned, family-filtered baits
  family_filtered_baits_out = baitfindR::write_fasta_files(
    fasta_list = family_filtered_baits,
    out_dir = "06_intron_masking/taxonomy_filtered"),
  
  # Align family-filtered baits in folder
  aligned_baits = baitfindR::mafft_wrapper(
    fasta_folder = "06_intron_masking/taxonomy_filtered",
    number_cores = 2,
    overwrite = TRUE,
    get_hash = TRUE,
    depends = family_filtered_baits_out
  ),
  
  # Clean up alignments in folder
  cleaned_aligned_baits = baitfindR::phyutility_wrapper(
    fasta_folder = "06_intron_masking/taxonomy_filtered",
    min_col_occup = 0.3,
    overwrite = TRUE,
    depends = aligned_baits),
  
  # Blast clean, aligned baits against intron-masked genes.
  # Results will be written to 06_intron_masking/taxonomy_filtered.
  blast_baits = baitfindR::blast_n_list(
    fasta_folder = "06_intron_masking/taxonomy_filtered",
    fasta_pattern = "\\.aln-cln$",
    database_path = "06_intron_masking/masked_genes",
    out_ext = "outfmt6",
    outfmt = "6 qseqid qlen sseqid slen evalue bitscore",
    overwrite = TRUE,
    depends = cleaned_aligned_baits
  ),
  
  # Extract top blast hit for each bait.
  # Results will be written to 06_intron_masking/blast_filtered.
  best_hits = baitfindR::extract_blast_hits(
    blast_results_dir = "06_intron_masking/taxonomy_filtered/",
    blast_results_pattern = "\\.outfmt6$",
    blast_cols = c("qseqid","qlen","sseqid","slen","evalue","bitscore"),
    database_path = "06_intron_masking/masked_genes",
    out_dir = "06_intron_masking/blast_filtered/",
    out_ext = "bestmatch",
    depends = blast_baits
  ),
  
  # Re-align filtered baits with top blast hits
  combined_alignments = baitfindR::realign_with_best_hits(
    best_hits_dir = "06_intron_masking/blast_filtered",
    best_hits_pattern = "\\.bestmatch$",
    fasta_dir = "06_intron_masking/taxonomy_filtered",
    fasta_pattern = "\\.fa$",
    depends = best_hits
  ),
  
  # Fill-in introns and remove outgroup
  combined_alignments_filled = purrr::map (
    combined_alignments,
    ~ baitfindR::fill_introns(
      ., 
      ref_pattern = "Sacu|Azfi|1|2|3|4|5", 
      outgroup = outgroup, 
      trim_outgroup = TRUE)
  ),
  
  # Calculate summary statistics for alignments
  combined_alignments_data = map_df(
    combined_alignments_filled, 
    baitfindR::calculate_alignment_stats, 
    include_aln = TRUE,
    .id = "bait_id"),
  
  # Filter alignments to get final baits
  # Here we filter by min. number of introns > 1,
  # and choosing top five ranked by % pars. inform. characters
  final_baits_data = filter(combined_alignments_data, num_introns > 1) %>%
    arrange(desc(pars_inf)) %>%
    slice(1:5),

  # Write out final baits
  final_baits_out = baitfindR::write_fasta_files(
    fasta_list = final_baits_data$alignment,
    fasta_names = final_baits_data$bait_id,
    ext = "fasta",
    out_dir = "07_baits")
  
)

# Output report
write_report <- drake_plan(
  report = rmarkdown::render(
    knitr_in("report.Rmd"),
    output_format = "html_document",
    output_file = file_out("report.html"),
    quiet = TRUE))

main_plan <- bind_plans(
  build_blastp_db, 
  run_transdecoder, 
  concatenate_cdhitest, 
  run_allbyall_blast, 
  concatenate_allbyall_blast, 
  run_mcl,
  sort_homologs_orthologs,
  mask_genes,
  concatenate_masked_genes,
  make_masked_genes_blast_db,
  mask_and_filter_baits,
  write_report)

rm(build_blastp_db, 
   run_transdecoder, 
   concatenate_cdhitest, 
   run_allbyall_blast, 
   concatenate_allbyall_blast, 
   run_mcl,
   sort_homologs_orthologs,
   mask_genes,
   concatenate_masked_genes,
   make_masked_genes_blast_db,
   mask_and_filter_baits,
   write_report)