# Odds ratio calulation
#########
#function
#######
# d is a data frame with 4 columns
# d$x gives variable names
# d$y gives center point
# d$ylo gives lower limits
# d$yhi gives upper limits
forestplot <- function(d, xlab="Odds Ratio", ylab="Study"){
  require(ggplot2)
  p <- ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi, fill=x)) + 
    geom_pointrange(width=0.1, color="blue", fill="blue") + 
    coord_flip() +
    geom_hline(aes(x=1, yintercept= 1), lty=1) +
    theme_classic()+
    ylab(xlab) +
    xlab(ylab)#switch because of the coord_flip() above
  return(p)
}
#key_ECO_MODS: dataframe contains expression of DE genes in ECMO and MODS patients
# class column containes label 1 and 0 for ECMO and MODS patients, respectively.
#column: gene
#rows: Sample ids
# This needs to repeat for all the genes
load("~/Desktop/ECMO_Vs_MODS_data/key_ECO_MODS.RData")
mylogit1 <- glm(key_ECO_MODS$class ~ scale(key_ECO_MODS[,1]), family = "binomial")
test1 = exp(cbind(OR = coef(mylogit1), confint(mylogit1)))
coef(summary(mylogit1))[,4][2]

#calculate the 95% confidance interval values 
exp(summary(mylogit1)$coefficients["scale(key_ECO_MODS[,1]",1] + 
      +     qnorm(c(0.025,0.5,0.975)) * summary(mylogit1)$coefficients["scale(key_ECO_MODS[, 1])",2])


##This will append the OR to the test2 
mylogit1 <- glm(key_ECO_MODS$class ~ scale(key_ECO_MODS[,2]), family = "binomial")
test2 = exp(cbind(OR = coef(mylogit1), confint(mylogit1)))
test2 <- rbind(t(test1[2,]),test2)
row.names(test2)[2] <- colnames(key_ECO_MODS)[2]

colnames(test2) = c( "y", "ylo", "yhi","x")
png("~/Desktop/ECMO_Vs_MODS_data/OR_for_key_genes.png")
forestplot(test2)+theme(axis.text = element_text(size=14, face="bold"))
dev.off()
