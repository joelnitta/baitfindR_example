# baitfindR_simple

# An example project showing how baitfindR can be combined with drake to
# automate the Yang and Smith (2014) (hereafter, Y&S) workflow to find
# candidate baits for sequence capture from a set of transcriptomes and
# reference genomes.

# Set working directory
setwd(here::here())

# Load packages
source("R/packages.R")

# Update drake settings
pkgconfig::set_config("drake::strings_in_dots" = "literals")

# Update baitfindR settigs
baitfindR::set_ys_path("/home/phylogenomic_dataset_construction")

# Load functions
source("R/functions.R")

### Define basic input values

# Set seed for reproducibility
set.seed(5394)

# Vector of 4-letter transcriptome codes.
# These will be downloaded from the 1KP website.
# Here we use a subset of eupolypod II ferns including
# species in Aspleniaceae, Athyriaceae, and Woodsiaceae
# as the ingroup, and a single eupolypod I fern as the outgroup
codes <- baitfindR::onekp_data$code

# Specify outgroup
outgroup <- "FQGQ"

# Set fraction for randomly downsizing transcriptomes (10%)
keep_frac = 0.10

# Values to use for mcl I value and Y&S hit-frac-cutoff
my_hit_frac <- 0.4
my_i_value <- 2

# Vector of reference genomes to use for making masked blast db
genomes <- c("arabidopsis", "azolla", "salvinia")

### Load plans
source("R/example_data.R")
source("R/main_plan.R")

### Download and pre-process example data
make(example_data_plan)

### Run analyses
#
# `jobs` is set to 1, but may be increased up to the number of
# cores (CPUs) available.
#
# Check the number of cores available with this command:
# future::availableCores()

make(main_plan, jobs = 1)
