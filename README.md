
baitfindR\_simple
=================

This is an example of how to use [baitfindR](https://github.com/joelnitta/baitfindR) in combination with [drake](https://github.com/ropensci/drake) to run the [Yang & Smith (2014) workflow](https://bitbucket.org/yangya/phylogenomic_dataset_construction) (hereafter, 'Y&S') to find candidate baits for sequence capture from a set of transcriptomes and reference genomes.

The benefits of running the Y&S workflow this way are:

-   `drake` tracks the input and output of each step. If you change something and need to re-run the analysis, you don't have to do everything from scratch. `drake` will only re-run the parts that need updating.

-   It is very easy to shorten run times by running steps in parallel on a cluster. Simply increase `n` in `make(main_plan, jobs = n)` in `make.R`.

-   `baitfindR` provides R scripts to help with bait selection that aren't available in Y&S.

Note that on 2018-08-30 a [new version of the Y&S workflow](https://bitbucket.org/yanglab/phylogenomic_dataset_construction/) was released. I will work on updating `baitfindR` to match, but for now it continues to use the [original version](https://bitbucket.org/yangya/phylogenomic_dataset_construction).

Workflow
--------

### Example data

The `example_data` plan downloads transcriptomes from the [1000 plants, or '1KP' project](https://sites.google.com/a/ualberta.ca/onekp/) as input to the Y&S workflow. By default, it downloads eight fern transcriptomes, as well proteomes and genomes of *Arabidopsis thaliana*, *Azolla filiculoides*, *Salvinia cucullata*, and *Lygodium japonicum* (proteome only).

Because it takes a long time to run the analyses on complete transcriptomes (in particular the BLAST steps), the transcriptomes are randomly downsized to 25% of their original size.

For more information about the eight sample transcriptomes, run `baitfindR::onekp_data`.

### Main plan

The `main_plan` first analyses the eight transcriptomes using Y&S workflow scripts to obtain a set of low-copy orthologs (candidate baits). Next, `baitfindR` scripts are used to insert introns into the orthologs, filter them based on criteria such as taxon sampling and fraction of parsimony-informative characters, and produce a set of final baits.

Intermediate output files are saved to the numbered folders (`01_translation`, `02_clustering`, etc). Final baits are saved to `07_baits`.

This plan uses pre-set values for `hit_frac_cutoff` and `I-value`, and the "one-to-one" method of ortholog pruning. To test the effects of using different values and pruning methods, see the [baitfindR\_complex example project](https://github.com/joelnitta/baitfindR_complex).

Docker image
------------

A docker image is provided to run the code. This is the preferred method, as there are many dependencies which are otherwise tedious to install, and may not even work correctly if versions have changed.

To use the image, first [install Docker](https://docs.docker.com/install/).

Pull the latest build:

    docker pull joelnitta/baitfindr_tidyverse

Dependencies
------------

These are the dependencies needed if you choose to install everything locally instead of using docker.

### R Packages

You will need to install all packages (and their dependencies) listed in `R/packages.R`. Most of these are available on CRAN and can be installed using `install.packages()`.

`baitfindR` and `jntools` must be installed from github using `devtools`.

``` r
devtools::install_github("joelnitta/baitfindR")
devtools::install_github("joelnitta/jntools")
```

### Other software

You will need all the [dependencies for Y&S](https://bitbucket.org/yangya/phylogenomic_dataset_construction/src/master/tutorials/part1_dependencies.md).

Usage
-----

### Running locally (not using Docker)

1.  Clone this repository to your computer:

        git clone https://github.com/joelnitta/baitfindR_simple.git
        cd baitfindR_drake

2.  Open `make.R`, [change the default paths for Y&S scripts](https://joelnitta.github.io/baitfindR/reference/set_ys_path.html), and save.

3.  Run `make.R` by either running `Rscript make.R` from the command line, or opening `make.R` in RStudio and clicking on 'Source'.

### Running with Docker (interactive mode)

1.  Clone the repo as above.

2.  Launch the container (where `/path/to/baitfindR_drake` is the full path to the cloned repository on your machine):

        docker run -it -v /path/to/baitfindR_drake:/home/rstudio/ joelnitta/baitfindr_tidyverse bash

3.  You should now be in the container. Run `make.R`:

        Rscript make.R

### Running with Docker (detached mode)

The scripts take a few hours to finish, so it may be preferable to run docker in a detached state in the background.

1.  Clone the repo as above.

2.  Launch the container with the `-d` flag (you need to provide ROOT, USER, and PASSWORD as if you were going to use RStudio server since the docker image is based off of `rocker/tidyverse`):

        docker run -d -v /path/to/baitfindR_drake:/home/rstudio/ --name baitfindr -e ROOT=true -e USER=rstudio -e PASSWORD=clever_pw joelnitta/baitfindr_tidyverse

3.  Enter the container:

        docker exec -it baitfindr bash

4.  Run `make.R` with `nohup` (and a log) so it continues running after we exit the shell:

        nohup Rscript make.R > make.log 2>&1 &
        exit

Settings
--------

Basic input values are defined in `make.R`:

-   `codes` is a vector of 4-letter codes used to identify [1KP transcriptome samples](http://www.onekp.com/samples/list.php). You may change this to any combination of valid 1KP codes. Default = eight fern transcriptomes.

-   `trim_frac` is the fraction to which the transcriptomes will be randomly downsized prior to analysis, to keep runtimes short for this example. Default = 25%.

-   `my_hit_frac` sets the minimum percentage overlap between query and target sequence in blast results to be retained for mcl clustering in the Y&S `blast_to_mcl.py` script. Default = 0.4.

-   `my_i_value` is the inflation value to use for the mcl clustering algorithm. Default = 2.

-   `genomes` is a vector of genome names that are used to look up the locations of fasta files and GFF files during the `mask_genes` sub-plan.

Analzying your own data
-----------------------

It should be fairly straightfoward if the user wants to analyze their own transcriptomes but leave the reference proteomes and genomes as they are. In that case, comment-out or delete `transcriptomes_data_plan` in `example_data.R`. Instead, provide your own transcriptome fasta-files in `data/` and a vector of file names for `codes`. Names for the transcriptome files should follow the rules for taxonIDs in Y&S (4--6 characters, letters and digits only).

Changing the reference proteomes and genomes is possible but will require more extensive editing of plans.

Parallel computing
------------------

It may be possible to shorten runtimes by increasing the number of jobs as described above on a personal computer if it has multiple cores. My trusty macbook pro retina can handle up to 4 jobs without any problems. Be sure to check your hardware specs before changing this though!
