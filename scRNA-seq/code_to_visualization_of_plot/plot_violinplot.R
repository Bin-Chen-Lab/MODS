library(ggplot2)
library(ggthemes)
library(ggpubr)
plot_violin= function(df,x,y, test.method,ylab, title){
  if(is.null(df)){
    print("Usage is plot_violin(df,x,y, ylab, ggttle)")
  }else{
    p=ggviolin(df, x=x, y=y, fill = x, size =0.5, alpha = 0.5, add = "boxplot")+
      #geom_jitter(position=position_jitter(width=.1, height=0), size = 0.5)+#avoid plotting outliers twice
      stat_compare_means(method = test.method, paired = F, size=20, label.x = 1.4, label = "p.signif")+
      theme_bw()+theme(axis.text = element_text( size=24, face="bold"),
                       axis.title.y = element_text(size=22, face="bold"),
                       axis.title.x = element_blank(),
                       legend.position = "none",
                       plot.title = element_text(size=24, face="bold",lineheight = 0.9),
                       panel.border = element_rect(fill=NA, colour = "black", size=1),
                       axis.ticks = element_line(size = 0.5))+ylab(ylab)+
      labs(title = title)
    
    p+scale_fill_manual(values = c("royalblue1", "peru"))
  }
  
}
