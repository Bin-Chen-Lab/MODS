boxplot_multilevel= function(df,x,y, ylab, title){
  library(ggplot2)
  library(ggthemes)
  library(ggpubr)
  
  if(is.null(df)){
    print("Usage is plot_boxplot(df,x,y, ylab, title)")
  }else{
   p=ggboxplot(df, x = x, y = y,          ##change the input file and output file names
              color = "black", fill = x, palette = "jco", legend = "none")+
      stat_compare_means(method = "anova", size= 8, label = "p.signif", vjust = 0.5)+
      geom_jitter(position=position_jitter(width=.1, height=0), size = 0.5)+
      theme_bw()+theme(axis.text.x = element_text(size=12, face="bold", angle = 45, hjust = 1), 
                       axis.text.y = element_text(size=12, face="bold"),
                       axis.title.y = element_text(size=12, face="bold"),
                       axis.title.x = element_blank(),
                       legend.position = "none",
                       plot.title = element_text(size=14, face="bold",lineheight = 0.9),
                       plot.subtitle = element_text(size=14, face="bold"),
                       panel.border = element_rect(fill=NA, colour = "black", size=1),
                       axis.ticks = element_line(size = 0.5))+ylab(ylab)+
      labs(title = title)
    print(p)
  }
  
}