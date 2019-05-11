# Make folders for storing intermediate output
# Add a .keep file to each for tracking with git.

make_dir <- function (dir_name, ...) {
  dir_name <- fs::path_abs(dir_name)
  if(!dir.exists(dir_name)) {
    dir.create(dir_name)
    sink(file = fs::path(dir_name, ".keep"), type = "output")
    cat(dir_name)
    sink()
  }
}

out_dirs <-
  c("01_translation",
    "02_clustering",
    "03_clusters",
    "04_homologs",
    "05_orthologs",
    "06_intron_masking",
    "07_baits")

purrr::walk(out_dirs, make_dir)

sub_dirs <-
  c("05_orthologs/fasta",
    "05_orthologs/tre",
    "06_intron_masking/blast_filtered",
    "06_intron_masking/taxonomy_filtered")

purrr::walk(sub_dirs, make_dir)
