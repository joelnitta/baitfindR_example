# example_data.R

# A drake plan to download and pre-process data for this example analysis

# Download genomes and annotations
genomes_data_plan <- drake_plan(
  
  arabidopsis_gff = download.file(
    url = "ftp://ftp.ensemblgenomes.org/pub/release-40/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.40.gff3.gz",
    destfile = file_out(here("data_raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz"))
  ),
  
  arabidopsis_genome = download.file(
    url = "ftp://ftp.ensemblgenomes.org/pub/release-40/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz",
    destfile = file_out(here("data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"))
  ),
  
  unzipped_arabidopsis_genome = R.utils::gunzip(
    filename = file_in(here("data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz")),
    destname = file_out(here("data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.fasta"))),
  
  renamed_arabidopsis_genome = rename_arabidopsis_genome(
    fasta_in = file_in(here("data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.fasta")),
    fasta_out = file_out(here("data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.renamed.fasta"))
  ),
  
  azolla_gff = download.file(
    url = "ftp://ftp.fernbase.org/Azolla_filiculoides/Azolla_asm_v1.1/Azolla_filiculoides.gene_models.highconfidence_v1.1.gff",
    destfile = file_out(here("data_raw/Azolla_filiculoides.gene_models.highconfidence_v1.1.gff"))
  ),
  
  azolla_genome = download.file(
    url = "ftp://ftp.fernbase.org/Azolla_filiculoides/Azolla_asm_v1.1/Azolla_filiculoides.genome_v1.2.fasta",
    destfile = file_out(here("data_raw/Azolla_filiculoides.genome_v1.2.fasta"))
  ),
  
  salvinia_gff = download.file(
    url = "ftp://ftp.fernbase.org/Salvinia_cucullata/Salvinia_asm_v1.2/Salvinia_cucullata.gene_models.highconfidence_v1.2.gff",
    destfile = file_out(here("data_raw/Salvinia_cucullata.gene_models.highconfidence_v1.2.gff"))
  ),
  
  salvinia_genome = download.file(
    url = "ftp://ftp.fernbase.org/Salvinia_cucullata/Salvinia_asm_v1.2/Salvinia_cucullata.genome_v1.2.fasta",
    destfile = file_out(here("data_raw/Salvinia_cucullata.genome_v1.2.fasta"))
  )
)

# Download and downsize transcriptomes
# Use wildcard transcriptomes__ to expand 
# across vector of transcriptome codes.
transcriptomes_data_plan <- drake_plan(
  
  # Download transcriptomes
  download = download.file(
    url = "http://206.12.96.204/onekp/transcriptome__-SOAPdenovo-Trans-assembly.fa.bz2",
    destfile = file_out(here("data_raw/transcriptome__-SOAPdenovo-Trans-assembly.fa.bz2"))
  ),
  
  # Downsize transcriptomes
  downsize = downsize_transcriptome(
    file = file_in(here("data_raw/transcriptome__-SOAPdenovo-Trans-assembly.fa.bz2")), 
    keep_frac = trim_frac),
  
  # Write-out transcriptomes
  write = ape::write.FASTA(
    x = downsize_transcriptome__,
    file = here("data/transcriptome__"))
) %>%
  evaluate_plan(rules = list(transcriptome__ = codes))

# Download and process proteomes
proteomes_data_plan <- drake_plan (
  
  # Download Lygodium proteome
  download_lygodium = download.file(
    url = "http://bioinf.mind.meiji.ac.jp/kanikusa/data/download/lygodium_predicted_potein_ver1.0RC.fasta.tar.gz",
    destfile = file_out(here("data_raw/lygodium_predicted_potein_ver1.0RC.fasta.tar.gz"))),
  
  # Unzip Lygodium proteome
  unzipped_lygodium = untar_tracked(
    tarfile = file_in(here("data_raw/lygodium_predicted_potein_ver1.0RC.fasta.tar.gz")), 
    compressed = "gzip", 
    exdir = "data_raw",
    outfile = file_out(here("data_raw/lygodium_predicted_potein_ver1.0RC.fasta"))),
  
  # Load Lygodium proteome
  lygodium = read.FASTA(
    file = file_in(here("data_raw/lygodium_predicted_potein_ver1.0RC.fasta")),
    type = "AA"),
  
  # Fix duplicate sequence names in Lygodium proteome
  trimmed_lygodium = trim_proteome(lygodium),
  
  # Download Arabidopsis proteome
  download_arabidopsis = download.file(
    url = "https://www.arabidopsis.org/download_files/Sequences/TAIR10_blastsets/TAIR10_pep_20110103_representative_gene_model_updated",
    destfile = file_out(here("data_raw/TAIR10_pep_20110103_representative_gene_model_updated.fasta"))),
  
  # Load Arabidopsis proteome
  arabidopsis_to_combine = read.FASTA(
    file = file_in(here("data_raw/TAIR10_pep_20110103_representative_gene_model_updated.fasta")),
    type = "AA"),
  
  # Download Azolla proteome
  azolla_proteome = download.file(
    url = "ftp://ftp.fernbase.org/Azolla_filiculoides/Azolla_asm_v1.1/Azolla_filiculoides.protein.highconfidence_v1.1.fasta",
    destfile = file_out(here("data_raw/Azolla_filiculoides.protein.highconfidence_v1.1.fasta"))
  ),
   
  # Load Azolla proteome
  azolla_to_combine = read.FASTA(
    file = file_in(here("data_raw/Azolla_filiculoides.protein.highconfidence_v1.1.fasta")),
    type = "AA"),
  
  # Download Salvinia proteome
  salvinia_proteome = download.file(
    url = "ftp://ftp.fernbase.org/Salvinia_cucullata/Salvinia_asm_v1.2/Salvinia_cucullata.protein.highconfidence_v1.2.fasta",
    destfile = file_out(here("data_raw/Salvinia_cucullata.protein.highconfidence_v1.2.fasta"))
  ),
  
  # Load Salvinia proteome
  salvinia_to_combine = read.FASTA(
    file = file_in(here("data_raw/Salvinia_cucullata.protein.highconfidence_v1.2.fasta")),
    type = "AA"),
  
  # Combine and write out proteomes
  combined_proteomes = c(arabidopsis_to_combine, 
                         trimmed_lygodium,
                         azolla_to_combine,
                         salvinia_to_combine),
  
  write_combined = write.FASTA(
    x = combined_proteomes, 
    file = file_out(here("data/combined_proteomes.fasta")))
)

# Combine into final data plan
example_data <- bind_plans(genomes_data_plan, transcriptomes_data_plan, proteomes_data_plan)
rm(genomes_data_plan, transcriptomes_data_plan, proteomes_data_plan)