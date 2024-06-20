
plot_dotplot= function(df,x,y, test.method,ylab, title){
  library(ggplot2)
  library(ggthemes)
  library(ggpubr)
  
  if(is.null(df)){
    print("Usage is plot_dotplot(df,x,y, ylab, ggttle)")
  }else{
    p=ggdotplot(df, x=x, y=y, fill = x, add= "mean_sd")+#avoid plotting outliers twice
      stat_compare_means(method = test.method, paired = F, size=4, label.x = 1.4)+
      theme_bw()+theme(axis.text = element_text( size=12, face="bold"),
                       axis.title.y = element_text(size=12, face="bold"),
                       axis.title.x = element_blank(),
                       legend.position = "none",
                       plot.title = element_text(size=14, face="bold",lineheight = 0.9),
                       panel.border = element_rect(fill=NA, colour = "black", size=1),
                       axis.ticks = element_line(size = 0.5))+ylab(ylab)+
      labs(title = title)
    
    p+scale_fill_manual(values = c("royalblue1", "peru"))
  }
  
}
