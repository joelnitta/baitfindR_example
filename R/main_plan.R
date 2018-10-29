# main_plan

# Define basic input values -----------------------------------------------

# Vector of transcriptome file names ("taxonID" of Y&S).
# According to Y&S, "Use short taxonIDs with 4-6 letters 
# and digits with no special characters"
# (same as those downloaded and pre-processed in data_plan)
codes <- baitfindR::onekp_data$code

# Values to try for mcl I value and Y&S hit-frac-cutoff
# These are also used in combination as names of subfolders in clusters/
hit_frac_list <- c(0.3, 0.4)
i_value_list <- c(1.4, 2)

# Make drake plan ---------------------------------------------------------

# Make blast protein db
# Blast databases include multiple files. Arbitrarily choose .phr as output 
# file for tracking.
build_blastp_db <- drake_plan(
  blastp_database = baitfindR::build_blast_db(
    in_seqs = file_in(here::here("data/arabidopsis_lygodium.fasta")), 
    out_name = here::here("01_translation/arabidopsis_lygodium"), 
    db_type = "prot", 
    other_args = "-parse_seqids",
    outfile = file_out(here::here("01_translation/arabidopsis_lygodium.phr"))
  )
)

# Run transdecoder, shorten output names, and prepare for clustering
run_transdecoder <- drake_plan(
  
  # Find candidate long ORFs (long open reading frames) in transcriptomes.
  # TransDecoder.LongOrfs outputs multiple files to .transdecoder_dir folder.
  # Use longest_orfs.pep as representative output file for tracking.
  longest_orfs = baitfindR::transdecoder_long_orfs(
    transcriptome_file = file_in(here::here("data/transcriptome__")),
    wd = here::here("01_translation"),
    other_args = "-S",
    # set infile to blast db .phr file to maintain dependency on blast db
    depends = file_in(here::here("01_translation/arabidopsis_lygodium.phr")),
    outfile = file_out(here::here("01_translation/transcriptome__.transdecoder_dir/longest_orfs.pep"))),
  
  # Blast candidate long ORFs against the protein database
  orf_blast = baitfindR::blast_p(
    query = file_in(here::here("01_translation/transcriptome__.transdecoder_dir/longest_orfs.pep")), 
    database = here::here("01_translation/arabidopsis_lygodium"),
    other_args = c("-max_target_seqs 1", "-evalue 10",  "-num_threads 2"),
    out_file = file_out(here::here("01_translation/transcriptome__.blastp.outfmt6"))),
  
  # Output the final ORFs, preferentially retaining those with blast hits
  orf_predict = baitfindR::transdecoder_predict(
    transcriptome_file = file_in(here::here("data/transcriptome__")),
    blast_result = file_in(here::here("01_translation/transcriptome__.blastp.outfmt6")),
    wd = here::here("01_translation"),
    infile = file_in(here::here("01_translation/transcriptome__.transdecoder_dir/longest_orfs.pep")),
    outfile = file_out(here::here("01_translation/transcriptome__.transdecoder.cds"))),
  
  # Shorten names from transdecoder using R function
  # Y&S function worked on whole folder, so not conducive to drake.
  # baitfindR function works on single file at a time.
  seq_with_short_names = baitfindR::fix_names_from_transdecoder(
    transdecoder_output = file_in(here::here("01_translation/transcriptome__.transdecoder.cds"))),
  
  # Write out fasta files with shortened names
  seq_with_short_names_out = ape::write.FASTA(
    seq_with_short_names_transcriptome__, 
    file = file_out(here::here("01_translation/transcriptome__.cds.fa"))),
  
  # Reduce redundancy for cds
  # Note that output is put into clustering/
  reduced = baitfindR::cd_hit_est(
    input = file_in(here::here("01_translation/transcriptome__.cds.fa")), 
    output = file_out(here::here("02_clustering/transcriptome__.cds.fa.cdhitest")), 
    other_args = c("-c", "0.99", "-n", "10", "-r", "0", "-T", "2"))) %>%
  evaluate_plan(rules = list(transcriptome__ = codes))

# Concatenate all the .cdhitest files into a new file all.fa
concatenate_cdhitest <- drake_plan(
  read = readr::read_file(file_in(here::here("02_clustering/transcriptome__.cds.fa.cdhitest")))) %>%
  evaluate_plan(rules = list(transcriptome__ = codes)) %>%
  bind_plans(gather_plan(., target = "cdhitest_files")) %>%
  bind_plans(drake_plan(
    concatenate = baitfindR::cat_files(cdhitest_files, 
                                       output_file=file_out(here::here("02_clustering/all.fa")))))

# All-by-all blast
run_allbyall_blast <- drake_plan(
  
  # Make all-by-all database
  all_blast_db = baitfindR::build_blast_db(
    in_seqs = file_in(here::here("02_clustering/all.fa") ),
    out_name = here::here("02_clustering/all.fa"),
    other_args = NULL,
    db_type = "nucl",
    outfile = file_out(here::here("02_clustering/all.fa.nhr"),
                       here::here("02_clustering/all.fa.nin"),
                       here::here("02_clustering/all.fa.nsq"))),
  
  # Blast transcriptomes
  blast_transcriptome = baitfindR::blast_n(
    query = file_in(here::here("02_clustering/transcriptome__.cds.fa.cdhitest")),
    database = here::here("02_clustering/all.fa"),
    out_file = file_out(here::here("02_clustering/transcriptome__.allbyall.blast.outfmt6")),
    other_args = c("-evalue 10", "-num_threads 1", "-max_target_seqs 1000"),
    infile = file_in(here::here("02_clustering/all.fa.nhr")))) %>% 
  evaluate_plan(rules = list(transcriptome__ = codes))

# All-by-all blast: concatenate output of single blast results
concatenate_allbyall_blast <- drake_plan(
  read_blast = readr::read_file(file_in(here::here("02_clustering/transcriptome__.allbyall.blast.outfmt6")))) %>%
  evaluate_plan(rules = list(transcriptome__ = codes)) %>%
  bind_plans(gather_plan(., target = "cdhitestblast_files")) %>%
  bind_plans(drake_plan(
    combined_cdhitest = baitfindR::cat_files(
      cdhitestblast_files, 
      output_file=file_out(here::here("02_clustering/all.rawblast")))))

# Filter raw blast output by hit fraction, run mcl, and write out as fasta files
run_mcl <- drake_plan (
  
  # prepare all by all blast results for mcl
  blast_formatted_for_mcl = baitfindR::blast_to_mcl(
    blast_results = file_in(here::here("02_clustering/all.rawblast")),
    hit_fraction_cutoff = ..MY_HIT_FRAC..,
    outfile = file_out(here::here("02_clustering/all.rawblast.hit-frac..MY_HIT_FRAC...minusLogEvalue"))),
  
  # run mcl
  mcl_clusters = baitfindR::mcl(
    mcl_input = file_in(here::here("02_clustering/all.rawblast.hit-frac..MY_HIT_FRAC...minusLogEvalue")),
    i_value = ..MY_I_VALUE..,
    e_value = 5,
    other_args = c("--abc", "-te", "2"),
    mcl_output = file_out(here::here("02_clustering/hit-frac..MY_HIT_FRAC.._I..MY_I_VALUE.._e5"))),
  
  # Write fasta files for each cluster from mcl output
  # These are named cluster1.fa, cluster2.fa, etc., and will be output to the subdirectories named after
  # the hit_frac and mcl settings used (e.g., hit-frac0.3_I1.4_e5/, etc) in clusters/
  fasta_clusters = baitfindR::write_fasta_files_from_mcl(
    all_fasta = file_in(here::here("02_clustering/all.fa")),
    mcl_outfile = file_in(here::here("02_clustering/hit-frac..MY_HIT_FRAC.._I..MY_I_VALUE.._e5")),
    minimal_taxa = 4,
    outdir = here::here("03_clusters/hit-frac..MY_HIT_FRAC.._I..MY_I_VALUE.._e5"),
    get_hash = TRUE,
    overwrite = TRUE)) %>% 
  evaluate_plan(list(
    "..MY_HIT_FRAC.." = hit_frac_list,
    "..MY_I_VALUE.." = i_value_list))

# Align each cluster, trim alignment, and infer a tree ("basic trees" which be further pruned downstream)
run_fasta_to_tree <- drake_plan(
  basic_trees = baitfindR::fasta_to_tree(
    overwrite = TRUE,
    seq_folder = here::here("03_clusters/hit-frac..MY_HIT_FRAC.._I..MY_I_VALUE.._e5"),
    number_cores = 1,
    seq_type = "dna",
    bootstrap = FALSE,
    get_hash = TRUE,
    infile = fasta_clusters_..MY_HIT_FRAC.._..MY_I_VALUE..)) %>%
  evaluate_plan(list(
    "..MY_HIT_FRAC.." = hit_frac_list,
    "..MY_I_VALUE.." = i_value_list)) %>%
  bind_plans(gather_plan(., target = "gathered_trees"))

# output report
report_plan <- drake_plan(
  rmarkdown::render(
    knitr_in("report_01.Rmd"),
    output_format = "html_document",
    output_file = file_out("report_01.html"), 
    quiet = TRUE))

plan_01 <- bind_plans(build_blastp_db, run_transdecoder, concatenate_cdhitest, run_allbyall_blast, concatenate_allbyall_blast, run_mcl, run_fasta_to_tree, report_plan)

rm(build_blastp_db, run_transdecoder, concatenate_cdhitest, run_allbyall_blast, concatenate_allbyall_blast, run_mcl, run_fasta_to_tree, report_plan)
