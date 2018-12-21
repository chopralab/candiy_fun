
extract.checked <- function(analysis) {
  analysis %>%
    group_by(bywhat) %>%
    summarise(checked = sum(is.checked * rank),
              not.checked = sum((!is.checked) * rank)) %>%
    mutate(norm.checked = 100 * (checked) / (checked + not.checked),
           norm.not.checked = (100 - norm.checked))
}

extract.by.mesh <- function(ts) {
  extract <- extract.checked(ts$analysis$drug.by.mesh)
  extract$rdrug <- colMeans(ts$rmcol.by.mesh)
  extract$rmesh <- colMeans(ts$rmesh.by.mesh)
  
  extract[ is.na(extract) ] <- 0
  extract <- arrange(extract, desc(norm.checked), desc(rdrug), desc(rmesh))
  extract$by.what <- factor( extract$bywhat, levels = extract$bywhat)
  
  # count known drugs
  extract <- merge(extract, ts$known.drugs %>% rename(by.what = MESH))
  
  extract
}

extract.by.drug <- function(ts) {
  extract <- extract.checked(ts$analysis$mesh.by.drug)
  extract$rdrug <- colMeans(ts$rmcol.by.drug)
  extract$rmesh <- colMeans(ts$rmesh.by.drug)

  extract[ is.na(extract) ] <- 0
  extract <- arrange(extract, desc(norm.checked), desc(rdrug), desc(rmesh))
  extract$by.what <- factor(extract$bywhat,
                                    levels = extract$bywhat)
  
  extract
}

extract.by.cdrug <- function(drug.extract) {
  drug.extract$category <- convert.col(
    categories, as.character(drug.extract$by.what), 'V2', 'V1')
  
  drug.extract <- arrange(drug.extract, category)
  drug.extract$by.what <- factor(drug.extract$by.what,
                                 levels = drug.extract$by.what)

  drug.extract
}

extract.by.category <- function(ts) {
  pred <- count.predictions(ts$interactions)
  known <- ts$known.drugs

  ts$analysis$mesh.by.drug %>%
    filter(is.checked) %>%
    mutate(category = convert.col(categories, bywhat, 'V2', 'V1' )) %>%
    merge(pred %>%  rename(match.name = MESH)) %>%
    merge(known %>% rename(match.name = MESH)) %>%
    group_by(category, match.name) %>%
    summarise(rank = sum(rank),
              prediction.count = unique(prediction.count),
              known.count = unique(known.count)) %>%
    ungroup %>%
    mutate(rank.per = rank / prediction.count) %>%
    arrange(desc(rank.per)) %>%
    mutate(match.name = factor(match.name, levels = unique(match.name)))
}

extract.all <- function(ts) {
  # Top10 here on purpose, ensures X-axis is constant
  extract <- extract.checked(top40$analysis$drug.by.mesh)
  extract[ is.na(extract) ] <- 0
  
  extract2 <- extract.checked(ts$analysis$mesh.by.drug)
  extract2[ is.na(extract2) ] <- 0
  
  extract2 <-
    extract2 %>%
    select(match.name = bywhat,
           norm.checked2 = norm.checked)

  pred <- count.predictions(ts$interactions)

  ts$analysis$drug.by.mesh %>%
    filter(is.checked) %>%
    #filter(rank!=0) %>%
    merge(pred %>% rename(bywhat = MESH)) %>%
    merge(extract) %>%
    #merge(extract2) %>%
    mutate(norm.checked2 = rank / ts$top) %>%
    as.tibble() %>%
    arrange(desc(norm.checked), rank) %>%
    mutate(bywhat = factor(bywhat, levels = unique(bywhat))) %>%
    mutate(match.name = factor(match.name, levels = unique(match.name))) %>%
    mutate(category = convert.col(categories, match.name, 'V2', 'V1' ))
    
}
