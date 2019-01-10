#!/usr/bin/Rscript

library(readr)
library(parallel)
library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(cowplot)
library(grid)
library(gridExtra)

source("utils.R")
source("ranker.R")
source("random_analysis.R")
source("extract.R")
source("graph.R")
source("table.R")

# adjust as needed
options(mc.cores = 8)

extract.and.graph <- function(ts) {
  ts$extract$by.mesh <- extract.by.mesh(ts)
  ts$extract$by.drug <- extract.by.drug(ts)
  ts$extract$by.cdrug <- extract.by.cdrug(ts$extract$by.drug)
  ts$extract$by.catg <-extract.by.category(ts)
  ts$extract$by.heat <- extract.heat(ts)
  
  ts$graphs$n.m.e <- ts$extract$by.mesh %>% figure.by.mesh(ts$top)
  ts$graphs$n.d.e <- ts$extract$by.drug %>% figure.by.drug(ts$top)
  ts$graphs$n.c.e <- ts$extract$by.cdrug %>% figure.drugs.by.category(ts$top)
  ts$graphs$m.cat <- ts$extract$by.catg %>% figure.by.category(ts$top)
  ts$graphs$m.all <- ts$extract$by.heat %>% figure.by.heat(ts$top)
}

if (!file.exists("rand2.rda")) {

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

  random.analysis(top10, 2)
  random.analysis(top25, 2)
  random.analysis(top40, 2)
  random.analysis(top100,2)

  cat('Random analysis complete\n')
} else {
  load("rand2.rda")
  
  rank.interactions(top10)
  rank.interactions(top25)
  rank.interactions(top40)
  rank.interactions(top100)
}

extract.and.graph(top10)
extract.and.graph(top25)
extract.and.graph(top40)
extract.and.graph(top100)

cat('Graphs made\n')

arrange.all.graphs() -> all.graphs
all.graphs %>% save.out.all.graphs

list() -> mesh.to.disease
top100$interactions %>% filter(MESH %in% top100$diseases) %>% select(MESH, disorder) %>% distinct -> temp
mesh.to.disease[temp$MESH %>% as.character] <- temp$disorder %>% as.character
unlist(mesh.to.disease) -> mesh.to.disease

all.top(table.by.mesh) %>%
  spread(top, per.psych) %>%
  arrange(desc(`10`), desc(`25`), desc(`40`), desc(`100`)) %>%
  write_csv(path = "bymesh.csv")

all.top(table.by.drug) %>%
  spread(top, per.mental) %>%
  arrange(desc(`10`), desc(`25`), desc(`40`), desc(`100`)) %>%
  write_csv(path = "bydrug.csv")

all.top(table.by.category) %>%
  spread(top, rank.per) %>%
  arrange(desc(`10`), desc(`25`), desc(`40`), desc(`100`)) %>%
  write_csv(path = "bycatg.csv")
