#!/usr/bin/Rscript

load.project.files <- function( interactions.file,
          diseases.of.interest.file,
          drugs.of.interest.file) {

  ts <- new.env()

  read_csv(interactions.file) -> ts$interactions
  ts$interactions$MESH <- factor(ts$interactions$MESH)
  ts$interactions$compound <- factor(ts$interactions$compound)
  ts$interactions$disorder <- factor(ts$interactions$disorder)
  ts$interactions$match <- factor(ts$interactions$match)

  read_tsv(drugs.of.interest.file, col_names = c("category", "ids", "cids", "name")) -> drugs.frame
  drugs.frame$ids -> ts$drug.ids

  list() -> ts$categories
  ts$categories[drugs.frame$ids] <- drugs.frame$category

  list() -> ts$drug.names
  ts$drug.names[drugs.frame$ids] <- drugs.frame$name

  read_lines(diseases.of.interest.file) -> diseases
  as.factor(diseases) -> ts$diseases

  count.known.drugs(ts$interactions) -> ts$known.drugs

  ts
}

load.preset.project.files <- function( top ) {
  load.project.files(paste0("all_disorders_",top,".csv.gz"),
                     "cando_mental_disorder.lst",
                     "psychoactives.tsv") -> sys
  sys$top <- top
  sys
}

count.known.drugs <- function(interactions) {
  interactions %>%
    select(MESH, col) %>%
    distinct %>%
    group_by(MESH) %>%
    tally %>%
    rename(known.count = n)
}

count.predictions <- function(interactions) {
  interactions %>%
    filter(score != 0) %>%
    select(MESH, mcol) %>%
    distinct %>%
    group_by(MESH) %>%
    tally %>%
    rename(prediction.count = n)
}

all.top <- function(table_func) {
  table_func(top10) %>%
    rbind(table_func(top25)) %>%
    rbind(table_func(top40)) %>%
    rbind(table_func(top100))
}

do.stat.test <- function( ts ) {
  attach(ts, warn.conflicts = F)
  
  ts$ks.tests$by.mesh.rdrug <-  
    ks.test(extract$by.mesh$norm.checked, extract$by.mesh$rdrug, alternative = "g")
  ts$ks.tests$by.mesh.rmesh <-
    ks.test(extract$by.mesh$norm.checked, extract$by.mesh$rmesh, alternative = "g")
  
  ts$ks.tests$by.drug.rdrug <-
    ks.test(extract$by.drug$norm.checked, extract$by.drug$rdrug, alternative = "g")
  ts$ks.tests$by.drug.rmesh <-
    ks.test(extract$by.drug$norm.checked, extract$by.drug$rmesh, alternative = "g")
 
  ts$ts.tests$by.mesh.rdrug <-
    t.test(extract$by.mesh$norm.checked, extract$by.mesh$rdrug,
           alternative = "g", paired = T)
  ts$ts.tests$by.mesh.rmesh <-
    t.test(extract$by.mesh$norm.checked, extract$by.mesh$rmesh,
           alternative = "g", paired = T)

  ts$ts.tests$by.drug.rdrug <-
    t.test(extract$by.drug$norm.checked, extract$by.drug$rdrug,
           alternative = "g", paired = T)
  ts$ts.tests$by.drug.rmesh <-
    t.test(extract$by.drug$norm.checked, extract$by.drug$rmesh,
           alternative = "g", paired = T)

  detach()
}
