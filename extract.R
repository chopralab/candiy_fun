
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
  extract$known.count <- count.known.drugs(ts$interactions, extract$by.what)
  
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
  ts$analysis$mesh.by.drug %>%
    filter(is.checked) %>%
    mutate(category = convert.col(categories, bywhat, 'V2', 'V1' )) %>%
    group_by(category, match.name) %>%
    summarise(rank = sum(rank)) %>%
    ungroup %>%
    arrange(desc(rank)) %>%
    mutate(match.name = factor(match.name, levels = unique(match.name)))
}
