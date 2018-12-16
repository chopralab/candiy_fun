#!/usr/bin/Rscript

library(readr)
library(parallel)
library(plyr)
library(dplyr)
library(tidyr)

source("utils.R")
source("ranker.R")
source("random_analysis.R")
source("extract.R")
source("graph.R")

# adjust as needed
options(mc.cores = 8)

extract.and.graph <- function(ts) {
  ts$extract$by.mesh <- extract.by.mesh(ts)
  ts$extract$by.drug <- extract.by.drug(ts)
  ts$extract$by.cdrug <- extract.by.cdrug(ts$extract$by.drug)
  ts$extract$by.catg <-extract.by.category(ts)
  
  ts$graphs$n.m.e <- ts$extract$by.mesh %>% figure.by.mesh(ts$top)
  ts$graphs$n.d.e <- ts$extract$by.drug %>% figure.by.drug(ts$top)
  ts$graphs$n.c.e <- ts$extract$by.cdrug %>% figure.by.drug(ts$top)
  ts$graphs$m.cat <- ts$extract$by.catg %>% figure.by.category(ts$top)
}

load.preset.project.files(10) -> top10
load.preset.project.files(25) -> top25
load.preset.project.files(40) -> top40
load.preset.project.files(100) -> top100

cat('Data loaded\n')

rank.interactions(top10)
rank.interactions(top25)
rank.interactions(top40)
rank.interactions(top100)

cat('Rankings calculated\n')

random.analysis(top10, 16)
random.analysis(top25, 16)
random.analysis(top40, 16)
random.analysis(top100,16)

cat('Random analysis complete\n')

categories <- read.table("categories.tsv", header = F)

extract.and.graph(top10)
extract.and.graph(top25)
extract.and.graph(top40)
extract.and.graph(top100)

cat('Graphs made\n')

arrange.all.graphs()
