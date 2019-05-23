
randomize.column <- function ( interactions, col ) {
  
  attach(interactions, warn.conflicts = F)
  
  as.list(as.character( sample( unique(col ) ) ) ) -> random.mapping
  names(random.mapping) <- as.character(unique(col))
  list2env(random.mapping) -> random.mapping
  
  lapply(as.character(col),
        function(x) {
          random.mapping[[x]]
        }
  ) %>% unlist-> random.column
  
  detach()
  
  return(random.column)
}

random.by.analysis <- function( sub.table,
                                rand.col, list1,
                                norm.col, list2,
                                total.runs=1000) {
  norm.col <- enquo(norm.col)
  mclapply( 1:total.runs, function(x) {
    sub.table$rand <- randomize.column(sub.table, rand.col)

    rColRanks <-
      general.ranker(sub.table,
                     !!norm.col, list2,
                     rand,   list1
                     )

    extract.checked(rColRanks) %>% 
      select(-checked, -norm.not.checked, -not.checked) %>%
      spread(bywhat, norm.checked, fill = 0)
  } ) %>% bind_rows
}

random.mc.analysis <- function( sub.table,
                                rand.col, list1,
                                norm.col, list2,
                                total.runs=1000) {
  norm.col <- enquo(norm.col)
  mclapply( 1:total.runs, function(x) {
    sub.table$rand <- randomize.column(sub.table, rand.col)
    
    general.ranker(sub.table,
                   rand,   list1, 
                   !!norm.col, list2
                   ) -> rColRanks
    
    extract.checked(rColRanks) %>% 
      select(-checked, -norm.not.checked, -not.checked) %>%
      spread(bywhat, norm.checked, fill = 0)
  } ) %>% bind_rows
}

random.analysis <- function(ts, rand.runs = 1000) {
  hack <- 
    tibble(MESH = ts$diseases[1],
           mcol = ts$drug.ids, score = 0) %>%
    rbind(tibble(MESH = ts$diseases,
                 mcol = ts$drug.ids[1], score = 0))
  
  sub.table <-
    ts$interactions %>%
    select(MESH, mcol, score) %>%
    rbind(hack)

  ts$random.analysis$rmcol.by.mesh <-
    random.by.analysis(sub.table, 
                       mcol,   ts$drug.ids,
                       MESH, ts$diseases, 
                       total.runs = rand.runs)
  
  ts$random.analysis$rmcol.by.drug <-
    random.mc.analysis(sub.table,
                       mcol, ts$drug.ids,
                       MESH, ts$diseases,
                       total.runs = rand.runs)

  ts$random.analysis$rmesh.by.drug <-
    random.by.analysis(sub.table,
                       MESH, ts$diseases,
                       mcol, ts$drug.ids,
                       total.runs = rand.runs)

  ts$random.analysis$rmesh.by.mesh <-
    random.mc.analysis(sub.table,
                       MESH, ts$diseases,
                       mcol, ts$drug.ids,
                       total.runs = rand.runs)
}
