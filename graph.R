library(cowplot)
library(gridExtra)

figure.by.mesh <- function ( extract.data, top, text.size=0.6 ) {
  extract.data %>% ggplot(aes(x=by.what)) +
    geom_bar(stat="identity", aes(y=norm.checked, fill=log(known.count)), width = 1) +
    geom_line(aes(y=rdrug, group = 1), color="#00BA38", size=1) +
    annotate( "segment",
              x=as.numeric(extract.data$by.what[1]),
              xend=as.numeric(tail(extract.data$by.what, 1)),
              y=mean(extract.data$rdrug),
              yend=mean(extract.data$rdrug),
              color = "darkgreen", alpha = 0.5 ) +
    scale_x_discrete(breaks=NULL) + 
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_continuous(low="#F8766D", high="#619CFF") +
    labs (fill = "Number of Known Compounds (Logarithmic)") +
    coord_cartesian(ylim=c(0,100)) +
    xlab("Indication") + ylab("% psychoactives") + ggtitle(paste0("Top ", top)) +
    guides(color=FALSE) + theme(legend.position = "top", plot.title = element_text(margin=margin(b = -10, unit = "pt")))
}

figure.by.drug <- function ( data, top ) {
  ggplot(data, aes(x=by.what)) +
    geom_bar(stat="identity", aes(y=norm.checked, fill="Normal"), width = 1) +
    geom_bar(aes(y=rdrug, fill = "Randomized Compounds"), stat = "identity", alpha=0.4, width = 1) +
    annotate( "segment",
              x=as.numeric(data$by.what[1]),
              xend=as.numeric(tail(data$by.what, 1)),
              y=mean(data$rdrug),
              yend=mean(data$rdrug),
              color="darkgreen", alpha = 0.5) +
    geom_bar(aes(y=rmesh, fill = "Randomized Indications"), stat = "identity", alpha=0.4, width = 1) +
    annotate( "segment",
              x=as.numeric(data$by.what[1]),
              xend=as.numeric(tail(data$by.what, 1)),
              y=mean(data$rmesh),
              yend=mean(data$rmesh),
              color="darkblue", alpha = 0.5) +
    scale_x_discrete(breaks=NULL) + 
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual( values = c("#F8766D", "#00BA38", "#619CFF"), guide = guide_legend(title = "") ) +
    coord_cartesian(ylim=c(0,100)) +
    xlab("Compound") + ylab("% mental indications") + ggtitle(paste0("Top ", top)) + guides(color=FALSE) +
    theme(legend.position = "top", plot.title = element_text(margin=margin(b = -10, unit = "pt")),
          axis.line.x = element_line(color="black",size = 0.1,linetype = "solid"), axis.line.y = element_line(size = 0.1))
}

figure.by.category <- function(data, top) {
  ggplot(data, aes(x=category, y=match.name, fill=rank)) +
    geom_tile() +
    xlab("Category of Psychoactive") +
    ylab("Indication") +
    scale_fill_continuous(low = "white") +
    ggtitle(paste0("Top", top)) +
    theme(legend.position = "top", plot.title = element_text(margin=margin(b = 10, unit = "pt"))) +
    scale_y_discrete(breaks=NULL) +
    labs(fill="Number of predictions")
}

save.out.all.graphs <- function(arranged.graphs) {
  ggsave("figure2.svg", arranged.graphs$by.mesh, height = 7.46)
  ggsave("figure3.svg", arranged.graphs$by.drug, height = 7.46)
  ggsave("figure4.svg", arranged.graphs$by.cate, height = 7.46)
  ggsave("figure5.svg", arranged.graphs$by.mcat, height = 7.46)
}

arrange.all.graphs <- function() {
  arranged.graphs = new.env()
  
  get_legend(top10$graphs$n.m.e) -> mesh.legend
  
  grid.arrange(top10$graphs$n.m.e + guides(fill=F), top25 $graphs$n.m.e + guides(fill=F),
               top40$graphs$n.m.e + guides(fill=F), top100$graphs$n.m.e + guides(fill=F),
               mesh.legend, ncol=2, nrow = 3, heights = c(2.5, 2.5, 0.5),
               layout_matrix = rbind( c(1,2), c(3,4), c(5,5))) -> arranged.graphs$by.mesh
  
  get_legend(top10$graphs$n.d.e) -> drug.legend
  
  grid.arrange(top10$graphs$n.d.e + guides(fill=F), top25 $graphs$n.d.e + guides(fill=F),
               top40$graphs$n.d.e + guides(fill=F), top100$graphs$n.d.e + guides(fill=F),
               drug.legend, ncol=2, nrow = 3, heights = c(2.5, 2.5, 0.5),
               layout_matrix = rbind( c(1,2), c(3,4), c(5,5))) -> arranged.graphs$by.drug
  
  get_legend(top10$graphs$n.c.e) -> cate.legend
  
  grid.arrange(top10$graphs$n.c.e + guides(fill=F), top25 $graphs$n.c.e + guides(fill=F),
               top40$graphs$n.c.e + guides(fill=F), top100$graphs$n.c.e + guides(fill=F),
               cate.legend, ncol=2, nrow = 3, heights = c(2.5, 2.5, 0.5),
               layout_matrix = rbind( c(1,2), c(3,4), c(5,5))) -> arranged.graphs$by.cate
  
  get_legend(top10$graphs$m.cat) -> mcat.legend
  
  grid.arrange(top10$graphs$m.cat + guides(fill=F), top25 $graphs$m.cat + guides(fill=F),
               top40$graphs$m.cat + guides(fill=F), top100$graphs$m.cat + guides(fill=F),
               mcat.legend, ncol=2, nrow = 3, heights = c(2.5, 2.5, 0.5),
               layout_matrix = rbind( c(1,2), c(3,4), c(5,5))) -> arranged.graphs$by.mcat
  
  
  return (arranged.graphs)
}
