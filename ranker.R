#!/usr/bin/Rscript

general.ranker <- function( sub.table,
                            rank.factor, rank.factor.selection,
                            match.col, match.check.list) {
  rank.factor = enquo(rank.factor)
  match.col = enquo(match.col)
  sub.table %>%
    filter( is.element(!!rank.factor, rank.factor.selection)) %>%
    group_by(!!rank.factor, !!match.col) %>%
    summarise(rank = sum(score != 0)) %>%
    ungroup %>%
    mutate(is.checked = is.element(!!match.col, match.check.list)) %>%
    rename(bywhat = !!rank.factor, match.name = !!match.col)
}

rank.interactions <- function(ts) {
  
  hack <- 
    tibble(MESH = ts$diseases[1],
         mcol = ts$drug.ids, score = 0) %>%
    rbind(tibble(MESH = ts$diseases,
                 mcol = ts$drug.ids[1], score = 0))

  sub.table <-
    ts$interactions %>%
    select(MESH, mcol, score) %>%
    rbind(hack)

  ts$analysis$mesh.by.drug <-
    sub.table %>% general.ranker(mcol, ts$drug.ids, MESH, ts$diseases)

  ts$analysis$drug.by.mesh <-
    sub.table %>% general.ranker(MESH, ts$diseases, mcol, ts$drug.ids)
}
