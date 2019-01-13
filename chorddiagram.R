library(circlize)

related.meshes <- function(ts, consensus = 0, min.score = 0) {
  ts$interactions %>%
    filter(mcol %in% ts$drug.ids) %>%
    filter(MESH %in% ts$diseases) %>% 
    filter(score < min.score) %>%
    group_by(mcol, MESH) %>%
    tally() %>%
    filter(n >= consensus) %>%
    ungroup ->
    tally.df
  
  lapply(tally.df$MESH %>% unique, function(current.mesh) {
    tally.df %>% filter(MESH == current.mesh) %>% select(mcol) -> pred.drug.ids
    pred.drug.ids$mcol -> pred.drug.ids

    tally.df %>%
      filter(MESH != current.mesh) %>%
      filter(mcol %in% pred.drug.ids) %>%
      select(mcol, MESH) %>%
      group_by(MESH) %>%
      tally %>%
      ungroup %>%
      mutate(to = current.mesh) %>%
      select(to, from = MESH, n)
  }) %>% bind_rows() %>%
    mutate(from = mesh.to.disease[from %>% as.character()],
           to = mesh.to.disease[to %>% as.character()] )
}

make.to.from.unique <- function(df) {
  df %>%
    mutate(to.from = interaction(
      if_else(to > from, to, from),
      if_else(to > from, from ,to))) %>%
    distinct(to.from, .keep_all = T) %>%
    select(-to.from) %>%
    arrange(from, to)
}

make.all.chords <- function(seed=1) {
  layout(matrix(c(1,3,2,4), 2, 2))
  
  set.seed(seed)
  indication.col <- rand_color(length(top10$diseases))
  names(indication.col) <- mesh.to.disease
  
  related.meshes(top10, 2) %>%
    filter(n != 1) %>%
    make.to.from.unique() %>%
    chordDiagram(annotationTrack = c("name", "grid"), grid.col = indication.col)

  related.meshes(top25, 3) %>%
    filter(n != 1) %>%
    make.to.from.unique() %>%
    chordDiagram(annotationTrack = c("name", "grid"), grid.col = indication.col)
  
  related.meshes(top40, 4) %>%
    filter(n != 1) %>%
    make.to.from.unique() %>%
    chordDiagram(annotationTrack = c("name", "grid"), grid.col = indication.col)
  
  related.meshes(top100, 6) %>%
    filter(n != 1) %>%
    make.to.from.unique() %>%
    chordDiagram(annotationTrack = c("name", "grid"), grid.col = indication.col)
  
  layout(matrix(1))
}
