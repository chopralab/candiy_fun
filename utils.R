#!/usr/bin/Rscript

load.project.files <- function( interactions.file,
          diseases.of.interest.file,
          drugs.of.interest.file) {

  ts <- new.env()

  read_csv(interactions.file) -> ts$interactions
  ts$interactions$MESH <- factor(ts$interactions$MESH)
  ts$interactions$compound <- factor(ts$interactions$compound)
  ts$interactions$disorder <- factor(ts$interactions$disorder)

  read.csv(drugs.of.interest.file, header=F, sep='') -> drugs.frame

  drugs.frame$V1 -> ts$drug.ids
  drugs.frame$V3 -> ts$drugs

  read_lines(diseases.of.interest.file) -> diseases
  as.factor(diseases) -> ts$diseases

  ts
}

load.preset.project.files <- function( top ) {
  load.project.files(paste0("all_disorders_",top,".csv.gz"),
                     "cando_mental_disorder.lst",
                     "psychoactive_list_2.tsv") -> sys
  sys$top <- top
  sys
}

convert.col <- function (interactions,
                         first.col.values,
                         first.col.name,
                         second.col.name) {
  
  aaply ( as.array(first.col.values), 1,
          function(x) {
            interactions[ interactions[[first.col.name]] == x, ] -> selected
            as.character( selected[[second.col.name]][1] )
          } )
}

count.known.drugs <- function ( interactions, meshes) {
  aaply( as.character(meshes), 1,
         function(x) {
           length(unique(interactions[ interactions$MESH == x, ]$col))
         }
  )
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
    ks.test(extract$by.drug$norm.checked, extract$by.drug$rdrug, alternative = "g")
  
  detach()
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
