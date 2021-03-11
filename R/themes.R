# Define themes
base_theme <- theme_minimal()+
  theme(
    strip.background = element_blank(),
    #strip.background = element_rect(colour = "black", fill = "lightgray"),
    strip.text = element_text(size=9, family = ""),
    axis.text.x =element_text(angle=45, hjust=1, vjust=1),
    axis.ticks = element_line(colour = "grey20"),
    #axis.title.x=element_blank(),
    plot.background = element_blank(),
    text = element_text(size=9, family = ""),
    axis.text = element_text(size=8, family = ""),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    panel.grid = element_line(size = rel(1))
    #panel.grid = element_blank()
    
  )
