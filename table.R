table.by.mesh <- function(ts) {
  ts$extract$by.mesh %>%
    as.tibble %>%
    mutate(indication = mesh.to.disease[bywhat %>% as.character]) %>%
    select(MESH = bywhat, indication, per.psych = norm.checked, known.count) %>%
    arrange(desc(per.psych)) %>%
    mutate(top = ts$top) %>%
    filter(per.psych != 0)
}

table.by.drug <- function(ts) {
  ts$extract$by.cdrug %>%
    mutate(compound = ts$drug.names[bywhat] %>% unlist) %>%
    select(compound, drug.id = bywhat, category, per.mental = norm.checked) %>%
    arrange(desc(per.mental)) %>%
    mutate(top = ts$top) %>%
    filter(per.mental != 0)
}

table.by.heat <- function(ts) {
  ts$extract$by.heat %>%
    filter(rank != 0) %>%
    arrange(desc(norm.checked), desc(norm.checked2)) %>%
    mutate(compound = ts$drug.names[drug.id] %>% unlist) %>%
    mutate(indication = mesh.to.disease[as.character(bywhat)]) %>%
    select(indication, compound, category, drug.id, per.psych = norm.checked, rank) %>%
    mutate(top = ts$top)
}

table.by.category <- function(ts) {
  ts$extract$by.catg %>%
    filter(rank != 0) %>%
    arrange(desc(rank.per)) %>%
    mutate(indication = mesh.to.disease[as.character(match.name)]) %>%
    select(indication, category, known.count, rank.per) %>%
    mutate(top = ts$top)
}