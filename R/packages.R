# Load the packages this workflow depends on

library(drake)
library(here)
library(glue)
library(knitr)
library(ape)
library(bedr)
library(tidyverse)
library(baitfindR)
library(jntools)

# Fix conflicts
expand <- tidyr::expand
gather <- tidyr::gather
here <- here::here
filter <- dplyr::filter
