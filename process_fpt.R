library(tidyverse)

all.fpt <- read_table("all_new.protein_compound_interact_real.fpt.xz", col_names=F)

all.fpt$X3 <- 0

save(all.fpt, file="all_fpt_vascunicol_zeroed.rda")

