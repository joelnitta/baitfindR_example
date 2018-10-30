# main_plan

# Make blast protein db from combined Arabidopsis, Azolla, Lygodium, 
# and Salvinia proteomes.
# Blast databases include multiple files. Arbitrarily choose .phr as output 
# file for tracking.
build_blastp_db <- drake_plan(
  blastp_database = baitfindR::build_blast_db(
    in_seqs = file_in(here("data/combined_proteomes.fasta")), 
    out_name = here("01_translation/combined_proteomes"), 
    db_type = "prot", 
    other_args = "-parse_seqids",
    outfile = file_out(here("01_translation/combined_proteomes.phr"))
  )
)

# Run transdecoder, shorten output names, and prepare for clustering
run_transdecoder <- drake_plan(
  
  # Find candidate long ORFs (long open reading frames) in transcriptomes.
  # TransDecoder.LongOrfs outputs multiple files to .transdecoder_dir folder.
  # Use longest_orfs.pep as representative output file for tracking.
  longest_orfs = baitfindR::transdecoder_long_orfs(
    transcriptome_file = file_in(here("data/transcriptome__")),
    wd = here("01_translation"),
    other_args = "-S",
    outfile = file_out(here("01_translation/transcriptome__.transdecoder_dir/longest_orfs.pep"))
  ),
  
  # Blast candidate long ORFs against the protein database
  orf_blast = baitfindR::blast_p(
    query = file_in(here("01_translation/transcriptome__.transdecoder_dir/longest_orfs.pep")), 
    database = here("01_translation/combined_proteomes"),
    # set infile to blast db .phr file to maintain dependency on blast db
    depends = file_in(here("01_translation/combined_proteomes.phr")),
    other_args = c("-max_target_seqs 1", "-evalue 10",  "-num_threads 2"),
    out_file = file_out(here("01_translation/transcriptome__.blastp.outfmt6"))
  ),
  
  # Output the final ORFs, preferentially retaining those with blast hits
  orf_predict = baitfindR::transdecoder_predict(
    transcriptome_file = file_in(here("data/transcriptome__")),
    blast_result = file_in(here("01_translation/transcriptome__.blastp.outfmt6")),
    wd = here("01_translation"),
    infile = file_in(here("01_translation/transcriptome__.transdecoder_dir/longest_orfs.pep")),
    outfile = file_out(here("01_translation/transcriptome__.transdecoder.cds"))
  ),
  
  # Shorten names from transdecoder using R function
  # Y&S function worked on whole folder, so not conducive to drake.
  # baitfindR function works on single file (transcriptome) at a time.
  seq_with_short_names = baitfindR::fix_names_from_transdecoder(
    transdecoder_output = file_in(here("01_translation/transcriptome__.transdecoder.cds"))
  ),
  
  # Write out fasta files with shortened names
  seq_with_short_names_out = ape::write.FASTA(
    x = seq_with_short_names_transcriptome__, 
    file = file_out(here("01_translation/transcriptome__.cds.fa"))
  ),
  
  # Reduce redundancy for cds
  # Note that output is put into clustering/
  reduced_cdhitest = baitfindR::cd_hit_est(
    input = file_in(here("01_translation/transcriptome__.cds.fa")), 
    output = file_out(here("02_clustering/transcriptome__.cds.fa.cdhitest")), 
    other_args = c("-c", "0.99", "-n", "10", "-r", "0", "-T", "2")
  )
) %>%
  evaluate_plan(rules = list(transcriptome__ = codes))

# Concatenate all the .cdhitest files into a new file all_orfs.fa
concatenate_cdhitest <- 
  drake_plan(
    read_reduced_cdhitest = read_file(file_in("02_clustering/transcriptome__.cds.fa.cdhitest"))
  ) %>%
  evaluate_plan(rules = list(transcriptome__ = codes)) %>%
  bind_plans(gather_plan(., target = "cdhitest_files")) %>%
  bind_plans(
    drake_plan(
      concatenated_cdhitest = baitfindR::cat_files(
        cdhitest_files,
        output_file=file_out(here("02_clustering/all_orfs.fa"))
      )
    )
  )

# All-by-all blast: make database and run blast on one transcriptome at a time
run_allbyall_blast <- drake_plan(
  
  # Make all-by-all database including all ORFs
  all_blast_db = baitfindR::build_blast_db(
    in_seqs = file_in(here("02_clustering/all_orfs.fa") ),
    out_name = here("02_clustering/all_orfs"),
    other_args = NULL,
    db_type = "nucl",
    outfile = file_out(here("02_clustering/all_orfs.nhr"))
  ),
  
  # Blast transcriptomes against database
  blast_transcriptome = baitfindR::blast_n(
    query = file_in(here("02_clustering/transcriptome__.cds.fa.cdhitest")),
    database = here("02_clustering/all_orfs"),
    out_file = file_out(here("02_clustering/transcriptome__.allbyall.blast.outfmt6")),
    other_args = c("-evalue 10", "-num_threads 1", "-max_target_seqs 1000"),
    depends = file_in(here("02_clustering/all_orfs.nhr"))
  )
) %>% 
  evaluate_plan(rules = list(transcriptome__ = codes))

# All-by-all blast: concatenate output of single blast results
concatenate_allbyall_blast <- 
  drake_plan(
    read_blast = readr::read_file(file_in(here("02_clustering/transcriptome__.allbyall.blast.outfmt6")))
  ) %>%
  evaluate_plan(rules = list(transcriptome__ = codes)) %>%
  bind_plans(gather_plan(., target = "cdhitestblast_files")) %>%
  bind_plans(
    drake_plan(
      combined_cdhitest = baitfindR::cat_files(
        cdhitestblast_files, 
        output_file=file_out(here("02_clustering/all.rawblast"))
      )
    )
  )

# Filter raw blast output by hit fraction, run mcl, and write out as fasta files
run_mcl <- drake_plan (
  
  # Prepare all by all blast results for mcl
  blast_formatted_for_mcl = baitfindR::blast_to_mcl(
    blast_results = file_in(here("02_clustering/all.rawblast")),
    hit_fraction_cutoff = my_hit_frac,
    outfile = file_out(here(glue("02_clustering/all.rawblast.hit-frac{my_hit_frac}.minusLogEvalue")))
    ),
  
  # Run mcl
  mcl_clusters = baitfindR::mcl(
    mcl_input = file_in(here(glue("02_clustering/all.rawblast.hit-frac{my_hit_frac}.minusLogEvalue"))),
    i_value = my_i_value,
    e_value = 5,
    other_args = c("--abc", "-te", "2"),
    mcl_output = file_out(here(glue("02_clustering/hit-frac{my_hit_frac}_I{my_i_value}_e5")))
    ),
  
  # Write fasta files for each cluster from mcl output to 03_clusters/
  fasta_clusters = baitfindR::write_fasta_files_from_mcl(
    all_fasta = file_in(here("02_clustering/all_orfs.fa")),
    mcl_outfile = file_in(here(glue("02_clustering/hit-frac{my_hit_frac}_I{my_i_value}_e5"))),
    minimal_taxa = 4,
    outdir = here("03_clusters"),
    get_hash = TRUE,
    overwrite = TRUE
    ), 

  # Align each cluster, trim alignment, and infer a tree
  # ("basic trees" which will be further pruned downstream)
  basic_trees = baitfindR::fasta_to_tree(
    overwrite = TRUE,
    seq_folder = here("03_clusters/"),
    number_cores = 1,
    seq_type = "dna",
    bootstrap = FALSE,
    get_hash = TRUE,
    infile = fasta_clusters)
)

#### START PLAN02 HERE

# Make homolog trees
make_homologs <- drake_plan(
  
  trimmed_trees = trim_tips(
    tree_folder = here("03_clusters"), 
    tree_file_ending = ".tre",
    relative_cutoff = 0.2,
    absolute_cutoff = 0.4, 
    overwrite = TRUE), 
  
  masked_trees = mask_tips_by_taxonID_transcripts(
    tree_folder = here("03_clusters"),
    aln_folder = here("03_clusters"),
    depends = trimmed_trees),
  
  homolog_trees = cut_long_internal_branches( 
    tree_folder = here("03_clusters/"),
    tree_file_ending = ".mm",
    internal_branch_length_cutoff = 0.3,
    minimal_taxa = 4,
    outdir = here::here("04_homologs/"),
    depends = masked_trees),
  
  homolog_seqs = write_fasta_files_from_trees(
    all_fasta = here::here("02_clustering/all_orfs.fa"),
    tree_folder = here::here("04_homologs"),
    tree_file_ending = ".subtree",
    outdir = here::here("04_homologs"),
    depends = homolog_trees),
  
  homolog_clean_trees = fasta_to_tree(
    seq_folder = here::here("04_homologs"),
    number_cores = 1,
    seq_type = "dna",
    bootstrap = FALSE,
    overwrite = TRUE,
    depends = homolog_seqs),
  
  # Prune orthologs using "1-to-1" method
  ortho_121_trees = filter_1to1_orthologs(
    tree_folder = here::here("04_homologs"), 
    tree_file_ending = ".tre",
    minimal_taxa = 3,
    outdir = here::here("05_orthologs/tre"),
    overwrite = TRUE,
    depends = homolog_clean_trees), 
  
  # Write out orthologs
  ortholog_fasta = write_ortholog_fasta_files( 
    all_fasta = file_in(here::here("02_clustering/all_orfs.fa")),
    tree_folder = here::here("05_orthologs/tre"),
    outdir = here::here("05_orthologs/fasta"),
    minimal_taxa = 3,
    overwrite = TRUE,
    depends = ortho_121_trees),
  
  # Count number of ingroup sequences per pruning method
  filtered = filter_fasta(
    seq_folder = here::here("05_orthologs/fasta"),
    taxonomy_data = onekp_data,
    min_taxa = 3,
    sample_col = "code",
    depends = ortholog_fasta)
)

# plan_03

# In plan_03, we will account for introns in the candidate bait sequences.

# Setup drake cache -------------------------------------------------------

if (file.exists(".plan_03_cache")) {
  plan_03_cache <- this_cache(".plan_03_cache")
} else {
  plan_03_cache <- new_cache(".plan_03_cache")
}

# Define basic input values -----------------------------------------------

# Vector of reference genomes to use for making masked blast db
genomes <- c("arabidopsis", "azolla", "salvinia")

# Ortholog pruning method selected after examining output of Plan 2
my_prune_method <- "ortho_MI"

# Subfolders for the output of this plan
subfolders <- c("taxonomy_filtered", "blast_filtered")

# Setup folders ---------------------------------------------------------

plan_03_folders <- drake_plan(
  top_folder = make_dir(dir_name = "introns", 
                        outfile = file_out(here::here("introns/.name")))) %>%
  bind_plans(
    drake_plan(
      subfolder = make_dir(dir_name = "introns/subfolder__",
                           infile = file_in(here::here("introns/.name")),
                           outfile = file_out(here::here("introns/subfolder__/.name")))
    )
  ) %>%
  evaluate_plan(rules = list(subfolder__ =  subfolders))

# Make intron-masked blast DB ---------------------------------------------

# Generate fasta files of intron-masked genes from reference genomes
mask_genes_plan <- drake_plan (
  
  # Find introns in reference genomes
  introns = find_bed_regions (
    gff3_file = get_gff3("genome_list"),
    genome_file = get_genome("genome_list"),
    source_select = get_sources("genome_list"),
    prefix = "genome_list",
    out_type = "write_all",
    out_dir = here::here("introns"),
    output = c(
      file_out(here::here("introns/genome_list_genes")),
      file_out(here::here("introns/genome_list_introns")),
      file_out(here::here("introns/genome_list_exons"))
    )
  ),
  
  # Mask introns in reference genomes
  masked_genome = mask_genome(
    introns_file = file_in(here::here("introns/genome_list_introns")),
    genome_file = get_genome("genome_list"),
    masked_genome_file = file_out(here::here("introns/genome_list_masked")),
    wd = here::here("introns")
  ),
  
  # Extract intron-masked genes from genomes
  masked_genes = extract_masked_genes(
    genes_file = file_in(here::here("introns/genome_list_genes")),
    masked_genome_file = file_in(here::here("introns/genome_list_masked")),
    masked_genes_file = file_out(here::here("introns/genome_list_masked_genes")),
    wd = here::here("introns")
  )
) %>%
  evaluate_plan(
    wildcard = "genome_list",
    values = genomes)

# Concatenate all the intron-masked genes into a single file
concatenate_masked_genes_plan <- drake_plan(
  read = readr::read_file(file_in(here::here("introns/genome_list_masked_genes")))) %>%
  evaluate_plan(wildcard = "genome_list",
                values = genomes) %>%
  bind_plans(gather_plan(., target = "masked_gene_files")) %>%
  bind_plans(drake_plan(
    concatenate = baitfindR::cat_files(masked_gene_files, 
                                       output_file=file_out(here::here("introns/all_masked_genes")))))

# Make blast-db of masked genes
masked_genes_blast_db_plan <- drake_plan(
  masked_db = baitfindR::build_blast_db(
    in_seqs = file_in(here::here("introns/all_masked_genes") ),
    out_name = here::here("introns/masked_genes"),
    other_args = "-parse_seqids",
    db_type = "nucl",
    outfile = file_out(here::here("introns/masked_genes.nhr"),
                       here::here("introns/masked_genes.nin"),
                       here::here("introns/masked_genes.nsq")))
)

# Mask introns in baits ---------------------------------------------------

mask_baits_plan <- drake_plan (
  
  # Filter baits by family
  # (each alignment must contain at least
  # one sample per ingroup family)
  family_filtered_baits = filter_fasta(seq_folder = here::here(paste0(my_prune_method, "/fasta")),
                                       taxonomy_data = onekp_data,
                                       filter_col = "family",
                                       sample_col = "code",
                                       depends = masked_db),
  
  # Write out unaligned, family-filtered baits
  family_filtered_baits_out = write_fasta_files(
    fasta_list = family_filtered_baits,
    out_dir = here::here("introns/taxonomy_filtered")),
  
  # Align family-filtered baits in folder
  aligned_baits = mafft_wrapper(
    fasta_folder = here::here("introns/taxonomy_filtered"),
    number_cores = 2,
    overwrite = TRUE,
    get_hash = TRUE,
    depends = family_filtered_baits_out
  ),
  
  # Clean up alignments in folder
  cleaned_aligned_baits = phyutility_wrapper(
    fasta_folder = here::here("introns/taxonomy_filtered"),
    min_col_occup = 0.3,
    overwrite = TRUE,
    depends = aligned_baits),
  
  # Blast clean, aligned baits against intron-masked genes
  blast_baits = blastn_list(
    fasta_folder = here::here("introns/taxonomy_filtered"),
    fasta_ending = "\\.aln-cln$",
    database = here::here("introns/masked_genes"),
    overwrite = TRUE,
    depends = cleaned_aligned_baits
  ),
  
  # Extract top blast hit for each bait
  best_hits = extract_blast_hits(
    blast_db = here::here("introns/masked_genes"),
    out_dir = here::here("introns/blast_filtered/"),
    blast_results_folder = here::here("introns/taxonomy_filtered/"),
    blast_results_ending = "\\.outfmt6$",
    depends = blast_baits
  ),
  
  # Re-align filtered baits with top blast hits
  combined_alignments = realign_with_best_hits(
    blast_filtered_folder = here::here("introns/blast_filtered"),
    taxonomy_filtered_folder = here::here("introns/taxonomy_filtered"),
    depends = best_hits
  ),
  
  # Fill-in introns and remove outgroup
  combined_alignments_filled = fill_introns_loop (
    alignment_list = combined_alignments,
    outgroup = "FQGQ"
  ),
  
  # Calculate summary statistics for alignments
  combined_alignments_data = assemble_bait_data (combined_alignments_filled),
  
  # Filter alignments to get final baits
  final_baits_data = filter_alignments (combined_alignments_data),
  
  final_baits = pull(final_baits_data, alignment)
  
)

# output report
# report_plan <- drake_plan(
#   rmarkdown::render(
#     knitr_in("report_01.Rmd"),
#     output_format = "html_document",
#     output_file = file_out("report_01.html"), 
#     quiet = TRUE))

main_plan <- bind_plans(build_blastp_db, run_transdecoder, concatenate_cdhitest, run_allbyall_blast, concatenate_allbyall_blast, run_mcl)

rm(build_blastp_db, run_transdecoder, concatenate_cdhitest, run_allbyall_blast, concatenate_allbyall_blast, run_mcl)
