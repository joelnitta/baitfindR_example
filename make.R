# baitfindR_drake

# An example project showing how baitfindR can be combined with drake to 
# automate the Yang and Smith (2014) (hereafter, Y&S) workflow to find
# candidate baits for sequence capture from a set of transcriptomes and 
# reference genomes.

### Set-up

# Load packages
source("R/packages.R") 

# Set working directory
setwd(here::here())

# Update drake settings
pkgconfig::set_config("drake::strings_in_dots" = "literals")

# Update baitfindR settigs
baitfindR::set_ys_path("/home/phylogenomic_dataset_construction")

# Load functions
source("R/functions.R")

# load plans
source("R/example_data.R")
source("R/plan_01.R")
source("R/plan_02.R")
source("R/plan_03.R")

### Run analyses (adjust number of jobs according to computer hardware)

# Download and pre-process example data
make(example_data_folders)
make(example_transcriptomes, cache = example_data_cache)
make(example_proteomes, cache = example_data_cache)
make(example_genomes, cache = example_data_cache)

# Determine optimal settings for hit_frac_cutoff and i_value
make(plan_01_folders)
make(plan_01, cache = plan_01_cache)

# Determine optimal homolog pruning method
make(plan_02_folders)
make(plan_02, cache = plan_02_cache)

# Mask introns
make(plan_03_folders)
make(plan_03, cache = plan_03_cache)

# Filter final baits
#make(plan_04_folders)
#make(plan_04, cache = plan_04_cache)
