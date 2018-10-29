# example_data

# Download and pre-process raw data for the example.

# Setup drake cache -------------------------------------------------------

if (file.exists(".example_data_cache")) {
  example_data_cache <- this_cache(".example_data_cache")
} else {
  example_data_cache <- new_cache(".example_data_cache")
}

# Define basic input values -----------------------------------------------

# Vector of 1KP transcriptome codes
# These will be downloaded from the 1KP website. 
# Here we use a subset of eupolypod II ferns including 
# species in Aspleniaceae, Athyriaceae, and Woodsiaceae
# as the ingroup, and a single eupolypod I fern as the outgroup
codes <- baitfindR::onekp_data$code

# (Note that the proteomes used for the blastp database are
# also downloaded, but these don't change much so they are 
# hard-coded in the drake plan)

# For the example, to what fraction should transcriptomes 
# be down-sized? e.g., 0.05 = 5%
trim_frac <- 0.25
set.seed(9542) # for reproducibility

# Make drake plans --------------------------------------------------------

# Process transcriptomes
transcriptomes_plan_part_1 <- drake_plan(
  
  # Download transcriptomes
  download = download_file(
    url = "http://206.12.96.204/onekp/transcriptome__-SOAPdenovo-Trans-assembly.fa.bz2",
    destfile = file_out(here::here("data_raw/transcriptome__-SOAPdenovo-Trans-assembly.fa.bz2")),
  ),
  
  # Downsize transcriptomes
  downsize = downsize_transcriptome(file = file_in(here::here("data_raw/transcriptome__-SOAPdenovo-Trans-assembly.fa.bz2")), 
                                    keep_frac = trim_frac)
) %>%
  evaluate_plan(rules = list(transcriptome__ = codes))

# Write out transcriptomes
transcriptomes_plan_part_2 <- drake_plan(
  write = ape::write.FASTA(
    x = downsize_transcriptome__,
    file = here::here("data/transcriptome__"))) %>%
  evaluate_plan(rules = list(transcriptome__ = codes)) %>%
  bind_plans(gather_plan(., target = "final_transcriptomes"))

transcriptomes_plan <- bind_plans(transcriptomes_plan_part_1, transcriptomes_plan_part_2)

# Process proteomes
example_proteomes <- drake_plan (
  
  # Download Lygodium
  download_lygodium = download_file(
    url = "http://bioinf.mind.meiji.ac.jp/kanikusa/data/download/lygodium_predicted_potein_ver1.0RC.fasta.tar.gz",
    destfile = file_out(here::here("data_raw/lygodium_predicted_potein_ver1.0RC.fasta.tar.gz"))),
  
  # Unzip Lygodium
  unzipped_lygodium = untar_tracked(
    tarfile = file_in(here::here("data_raw/lygodium_predicted_potein_ver1.0RC.fasta.tar.gz")), 
    compressed = "gzip", 
    exdir = "data_raw",
    outfile = file_out(here::here("data_raw/lygodium_predicted_potein_ver1.0RC.fasta"))),
  
  # Load lygodium
  lygodium = read_fasta_tracked(
    file = file_in(here::here("data_raw/lygodium_predicted_potein_ver1.0RC.fasta")),
    type = "AA"),
  
  # Fix duplicate sequence names in Lygodium
  trimmed_lygodium = trim_proteome(lygodium),
  
  # Download Arabidopsis
  download_arabidopsis = download_file(
    url = "https://www.arabidopsis.org/download_files/Sequences/TAIR10_blastsets/TAIR10_pep_20110103_representative_gene_model_updated",
    destfile = file_out(here::here("data_raw/TAIR10_pep_20110103_representative_gene_model_updated.fasta")),
    depends = trimmed_lygodium),
  
  # Load Arabidopsis
  arabidopsis = read_fasta_tracked(
    file = file_in(here::here("data_raw/TAIR10_pep_20110103_representative_gene_model_updated.fasta")),
    type = "AA"),
  
  # Combine and write out
  combined_proteomes = c(arabidopsis, trimmed_lygodium),
  
  write_combined = write_fasta_tracked(
    seq = combined_proteomes, 
    file_name = file_out(here::here("data/arabidopsis_lygodium.fasta")))
)

# Get genomes and annotations
example_genomes <- drake_plan(
  arabidopsis_gff = download_file(
    url = "ftp://ftp.ensemblgenomes.org/pub/release-40/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.40.gff3.gz",
    destfile = file_out(here::here("data_raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz"))
  ),
  
  arabidopsis_genome = download_file(
    url = "ftp://ftp.ensemblgenomes.org/pub/release-40/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz",
    destfile = file_out(here::here("data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"))
  ),
  
  unzipped_arabidopsis_genome = R.utils::gunzip(
    filename = file_in(here::here("data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz")),
    destname = file_out(here::here("data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.fasta"))),
  
  renamed_arabidopsis_genome = rename_arabidopsis_genome(
    fasta_in = file_in(here::here("data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.fasta")),
    fasta_out = file_out(here::here("data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.renamed.fasta"))
  ),
  
  azolla_gff = download_file(
    url = "ftp://ftp.fernbase.org/Azolla_filiculoides/Azolla_asm_v1.1/Azolla_filiculoides.gene_models.highconfidence_v1.1.gff",
    destfile = file_out(here::here("data_raw/Azolla_filiculoides.gene_models.highconfidence_v1.1.gff"))
  ),
  
  azolla_genome = download_file(
    url = "ftp://ftp.fernbase.org/Azolla_filiculoides/Azolla_asm_v1.1/Azolla_filiculoides.genome_v1.2.fasta",
    destfile = file_out(here::here("data_raw/Azolla_filiculoides.genome_v1.2.fasta"))
  ),
  
  salvinia_gff = download_file(
    url = "ftp://ftp.fernbase.org/Salvinia_cucullata/Salvinia_asm_v1.2/Salvinia_cucullata.gene_models.highconfidence_v1.2.gff",
    destfile = file_out(here::here("data_raw/Salvinia_cucullata.gene_models.highconfidence_v1.2.gff"))
  ),
  
  salvinia_genome = download_file(
    url = "ftp://ftp.fernbase.org/Salvinia_cucullata/Salvinia_asm_v1.2/Salvinia_cucullata.genome_v1.2.fasta",
    destfile = file_out(here::here("data_raw/Salvinia_cucullata.genome_v1.2.fasta"))
  )
)
