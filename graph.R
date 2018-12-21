figure.by.mesh <- function(extract.data, top, text.size = 0.6) {
  extract.data %>% ggplot(aes(x = by.what)) +
    geom_col(aes(y = norm.checked, fill = known.count), width = 1) +
    geom_line(aes(y = rdrug, group = 1), color = "#00BA38", size = 1) +
    #geom_line(aes(y=rmesh, group = 1), color="#619CFF", size=1) +
    geom_hline(color = "darkgreen", alpha = 0.5, yintercept = mean(extract.data$rdrug)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_continuous(low = "#F8766D", high = "#619CFF", trans = "log", breaks = c(1, 10, 60)) +
    labs(fill = "Number of Known Compounds") +
    coord_cartesian(ylim = c(0, 100)) +
    xlab("Indication") + ylab("% psychoactives") + ggtitle(paste0("Top ", top)) +
    guides(color = FALSE) +
    theme(legend.position = "bottom",
          legend.key.width = unit(.12, "npc"),
          plot.title = element_text(margin = margin(b = -10, unit = "pt")),
          axis.ticks.x = element_blank(), axis.text.x = element_blank())
}

figure.by.drug <- function(data, top) {
  ggplot(data, aes(x = by.what)) +
    geom_col(aes(y = norm.checked, fill = "Normal"), width = 1) +

    # Random controls    
    geom_col(aes(y = rdrug, fill = "Randomized Compounds"), alpha = 0.4, width = 1) +
    geom_hline(color = "darkgreen", alpha = 0.5, yintercept = mean(data$rdrug)) +
    geom_col(aes(y = rmesh, fill = "Randomized Indications"), alpha = 0.4, width = 1) +
    geom_hline(color = "darkgreen", alpha = 0.5, yintercept = mean(data$rmesh)) +

    # Aesthetics
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF"), guide = guide_legend(title = "")) +
    xlab("Compound") + ylab("% mental indications") + ggtitle(paste0("Top ", top)) + guides(color = FALSE) +
    theme(legend.position = "bottom",
          plot.title = element_text(margin = margin(b = -10, unit = "pt")),
          axis.ticks.x = element_blank(), axis.text.x = element_blank())
}

figure.drugs.by.category <- function(data, top) {
  figure.by.drug(data, top) +
  facet_grid(~category, labeller = labeller(category = category_names), switch = "both", scales = "free")
}

# figure.by.category <- function(data, top) {
#   ggplot(data, aes(x=category, y=match.name, fill=rank)) +
#     geom_tile() +
#     xlab("Category of Psychoactive") +
#     ylab("Indication") +
#     scale_fill_continuous(low = "white") +
#     ggtitle(paste0("Top", top)) +
#     scale_y_discrete(breaks=NULL) +
#     scale_x_discrete(labels = category_names) +
#     labs(fill="Number of predictions") +
#     theme(legend.position = "bottom", legend.key.width = unit(.12, "npc"),
#           legend.spacing = unit(0.1, "npc"), legend.justification = "center",
#           plot.title = element_text(margin=margin(b = 10, unit = "pt")),
#           axis.text.x = element_text(angle = 10, size = 9))
# }

figure.by.category <- function(data, top) {
  ggplot(data, aes(match.name, rank.per, fill = known.count)) +
    facet_grid(category ~ ., labeller = labeller(category = category_names), switch = "both") +
    geom_col() +
    scale_y_sqrt(breaks = c(0.1, 0.5, 1.0), expand = c(0, 0), limits = c(0, 1.0)) +
    scale_fill_continuous(low = "#F8766D", high = "#619CFF", trans = "log", breaks = c(1, 10, 60)) +
    labs(fill = "Number of Known Drugs") +
    xlab("Indication") + ylab("% psychoactives") + ggtitle(paste0("Top ", top)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 9),
          legend.position = "bottom", legend.title = element_text(size = 9),
          legend.key.width = unit(.12, "npc"), legend.text = element_text(size = 9),
          legend.key.height = unit(.02, "npc"), strip.text.y = element_text(size = 9, angle = -105),
          strip.placement = "outsie", strip.background = element_blank(),
          plot.title = element_text(margin = margin(b = -10, unit = "pt")))
}

figure.by.heat <- function(data, top) {
  ggplot(data, aes(match.name, bywhat, fill = norm.checked2 * 100)) +
    facet_grid(~category, labeller = labeller(category = category_names), switch = "both", scales = "free", space="free") +
    geom_tile() +
    scale_fill_viridis_c( trans = "log1p", breaks = c(1, 5, 10, 20, 30), option = "D", limits = c(1, 30)) +
    labs(fill = "% occurance") +
    ylab("Indication") + xlab("Psychoactive") + ggtitle(paste0("Top ", top)) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          legend.position = "bottom", legend.title = element_text(size = 9, vjust=0.80), legend.text=element_text(size = 9),
          legend.key.width = unit(.12, "npc"), strip.text.x=element_text(angle=15, size=9),
          strip.placement = "outsie", strip.background = element_blank())
}

save.out.all.graphs <- function(arranged.graphs) {
  ggsave("figure2.pdf", arranged.graphs$by.mesh, height = 7.46)
  ggsave("figure3.pdf", arranged.graphs$by.drug, height = 7.46)
  ggsave("figure4.pdf", arranged.graphs$by.cate, height = 7.46)
  ggsave("figure5.pdf", arranged.graphs$by.heat, height = 7.46)
  ggsave("figure6.pdf", arranged.graphs$by.mcat, height = 7.46)
}

arrange.all.graphs <- function() {
  arranged.graphs = new.env()

  get_legend(top10$graphs$n.m.e) -> mesh.legend

  grid.arrange(top10$graphs$n.m.e + guides(fill = F), top25$graphs$n.m.e + guides(fill = F),
               top40$graphs$n.m.e + guides(fill = F), top100$graphs$n.m.e + guides(fill = F),
               mesh.legend, ncol = 2, nrow = 3, heights = c(2.5, 2.5, 0.5),
               layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5))) -> arranged.graphs$by.mesh

  get_legend(top10$graphs$n.d.e) -> drug.legend

  grid.arrange(top10$graphs$n.d.e + guides(fill = F), top25$graphs$n.d.e + guides(fill = F),
               top40$graphs$n.d.e + guides(fill = F), top100$graphs$n.d.e + guides(fill = F),
               drug.legend, ncol = 2, nrow = 3, heights = c(2.5, 2.5, 0.5),
               layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5))) -> arranged.graphs$by.drug

  get_legend(top10$graphs$n.c.e) -> cate.legend

  grid.arrange(top10$graphs$n.c.e + guides(fill = F), top25$graphs$n.c.e + guides(fill = F),
               top40$graphs$n.c.e + guides(fill = F), top100$graphs$n.c.e + guides(fill = F),
               cate.legend, ncol = 2, nrow = 3, heights = c(2.5, 2.5, 0.5),
               layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5))) -> arranged.graphs$by.cate

  get_legend(top10$graphs$m.cat) -> mcat.legend

  grid.arrange(top10$graphs$m.cat + guides(fill = F) + theme(axis.title = element_blank()),
               top25$graphs$m.cat + guides(fill = F) + theme(axis.title = element_blank()),
               top40$graphs$m.cat + guides(fill = F) + theme(axis.title = element_blank()),
               top100$graphs$m.cat + guides(fill = F) + theme(axis.title = element_blank()),
               textGrob("Indication"), ncol = 2, nrow = 3, heights = c(3, 3, 0.3),
               layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)), bottom = mcat.legend,
               left = textGrob("% psychoactive", rot = 90)) -> arranged.graphs$by.mcat

  get_legend(top10$graphs$m.all) -> mcat.legend

  grid.arrange(top10$graphs$m.all + guides(fill = F) + theme(axis.title = element_blank()),
               top25$graphs$m.all + guides(fill = F) + theme(axis.title = element_blank()),
               top40$graphs$m.all + guides(fill = F) + theme(axis.title = element_blank()),
               top100$graphs$m.all + guides(fill = F) + theme(axis.title = element_blank()),
               textGrob("Indication"), ncol = 2, nrow = 3, heights = c(3, 3, 0.5),
               layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)), bottom = mcat.legend,
               left = textGrob("% psychoactive", rot = 90)) -> arranged.graphs$by.heat

  return(arranged.graphs)
}
