# baitfindR_example

This is an example of how to use [baitfindR](https://github.com/joelnitta/baitfindR) in combination with [drake](https://github.com/ropensci/drake) to run the [Yang & Smith (2014) workflow](https://bitbucket.org/yangya/phylogenomic_dataset_construction) (hereafter, 'Y&S') to find candidate baits for sequence capture from a set of transcriptomes and reference genomes.

The benefits of running the Y&S workflow this way are:

- `drake` tracks the input and output of each step. If you change something and need to re-run the analysis, you don't have to do everything from scratch. `drake` will only re-run the parts that need updating.

- It is very easy to shorten run times by running steps in parallel on a cluster. Simply increase `n` in `make(main_plan, jobs = n)` in `make.R`.

- `baitfindR` provides R scripts to help with bait selection that aren't available in Y&S.

Note that on 2018-08-30 a [new version of the Y&S workflow](https://bitbucket.org/yanglab/phylogenomic_dataset_construction/) was released. I will work on updating `baitfindR` to match, but for now it continues to use the [original version](https://bitbucket.org/yangya/phylogenomic_dataset_construction).

Also note that running this example uses a significant amount of disk space, approx. 3.5 gb total. Please make sure enough space is available before running `make.R`. It should take 1-2 hours to finish.

## Workflow

### Example data

`example_data_plan` downloads transcriptomes from the [1000 plants, or '1KP' project](https://sites.google.com/a/ualberta.ca/onekp/) as input to the Y&S workflow. By default, it downloads eight fern transcriptomes, as well proteomes and genomes of *Arabidopsis thaliana*, *Azolla filiculoides*, *Salvinia cucullata*, and *Lygodium japonicum* (proteome only).

Because it takes a long time to run the analyses on complete transcriptomes (in particular the BLAST steps), the transcriptomes are randomly downsized to 10% of their original size. This can be changed in `make.R`.

For more information about the eight sample transcriptomes, run `baitfindR::onekp_data`.

### Main plan

The `main_plan` first analyses the eight transcriptomes using Y&S workflow scripts to obtain a set of low-copy orthologs (candidate baits). Next, `baitfindR` scripts are used to insert introns into the orthologs, filter them based on criteria such as taxon sampling and fraction of parsimony-informative characters, and produce a set of final baits.

Intermediate output files are saved to the numbered folders (`01_translation`, `02_clustering`, etc). Final baits are saved to `07_baits`.

## Docker image

[A docker image](https://hub.docker.com/r/joelnitta/baitfindr) is provided to run the code. This is the preferred method, as there are many dependencies which are otherwise tedious to install, and may not work correctly if versions have changed.

To use the image, first [install Docker](https://docs.docker.com/install/).

Pull the latest build:

```
docker pull joelnitta/baitfindr
```

## Dependencies

These are the dependencies needed if you choose to install everything locally instead of using docker.

### R Packages

You will need to install all packages (and their dependencies) listed in `R/packages.R`. Most of these are available on CRAN and can be installed using `install.packages()`.

`baitfindR` must be installed from github using `devtools`.

```{r dependencies, eval=FALSE}
devtools::install_github("joelnitta/baitfindR")
```

### Other software

You will need all the [dependencies for Y&S](https://bitbucket.org/yangya/phylogenomic_dataset_construction/src/master/tutorials/part1_dependencies.md).

## Usage

### Running locally (not using Docker)

Clone this repository to your computer:

```
git clone https://github.com/joelnitta/baitfindR_example.git
cd baitfindR_example
```

Open `make.R`, [change the default paths for Y&S scripts](https://joelnitta.github.io/baitfindR/reference/set_ys_path.html), and save.

Run `make.R` by either running `Rscript make.R` from the command line, or opening `make.R` in RStudio and clicking on 'Source'.

### Running with Docker (interactive mode)

Clone the repo as above.

Launch the container (where `/path/to/repo` is the full path to the cloned repository on your machine, e.g. `/home/me/baitfindR_example`):

```
docker run -it -v /path/to/repo:/home/rstudio/baitfindR_example joelnitta/baitfindr bash
```

You should now be in the container. Run `make.R`:

```
cd baitfindR_example
Rscript make.R
```

### Running with Docker (detached mode)

The scripts take a few hours to finish, so it may be preferable to run docker in a detached state in the background.

Clone the repo as above.

Launch the container in detached mode.  There are two options to do this. 

- Use `docker run` with the `-d` flag:
```
docker run -d -v /path/to/repo:/home/rstudio/baitfindR_example --name baitfindr_example_drake_1 -e DISABLE_AUTH=true joelnitta/baitfindr
```

- Or, use `docker-compose` after navigating to the repo (this way requires less typing):
```
cd /path/to/repo
docker-compose up -d
```

Enter the container:

```
docker exec -it baitfindr_example_drake_1 bash
```

Run `make.R` with `nohup` (and a log) so it continues running after we exit the shell:

```
cd baitfindR_example
nohup Rscript make.R > make.log 2>&1 &
exit
```

You will know the plan is finished running when the report (`report.html`) is successfully compiled. This takes 1-2 hours with the example dataset.

Clean up the container when you're done:

```
docker kill baitfindr_example_drake_1
```

## Settings

Basic input values are defined in `make.R`:

- `codes` is a vector of 4-letter codes used to identify [1KP transcriptome samples](http://www.onekp.com/samples/list.php). You may change this to any combination of valid 1KP codes. Default = eight fern transcriptomes.

- `outgroup` is the code for the outgroup taxon. Default = FQGQ (*Polystichum acrostichoides*).

- `keep_frac` is the fraction to which the transcriptomes will be randomly downsized prior to analysis, to keep runtimes short for this example. Default = 10%.

- `my_hit_frac` sets the minimum percentage overlap between query and target sequence in blast results to be retained for mcl clustering in the Y&S `blast_to_mcl.py` script. Default = 0.4.

- `my_i_value` is the inflation value to use for the mcl clustering algorithm. Default = 2.

- `genomes` is a vector of genome names that are used to look up the locations of fasta files and GFF files during the `mask_genes` sub-plan.

## Analzying your own data

It should be fairly straightfoward if the user wants to analyze their own transcriptomes but leave the reference proteomes and genomes as they are. In that case, comment-out or delete `transcriptomes_plan` in `example_data.R`. Instead, provide your own transcriptome fasta-files in `data/` and a vector of file names for `codes` and `outgroup` in `make.R`. Names for the transcriptome files should follow the rules for taxonIDs in Y&S (4--6 characters, letters and digits only).

Changing the reference proteomes and genomes is possible but will require more extensive editing of plans.

## Parallel computing

It may be possible to shorten runtimes by increasing the number of jobs as described above on a personal computer if it has multiple cores. My trusty macbook pro retina can handle up to 4 jobs without any problems.  Be sure to check your hardware specs before changing this though!
