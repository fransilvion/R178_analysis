#genes
#ENSG00000100926 - TMSF1
#ENSG00000092036 - HAUS4
#ENSG00000184304 - PRKD1
#ENSG00000259132 (paralogy of HAUS4)
#ENSG00000259522 (paralogy to TM9SF1)
x <- plotCounts(ddsTxi, gene=which(rownames(res) == "ENSG00000259132"), intgroup="Affected", returnData=TRUE)

x %>%
  ggplot(aes(x=Affected, y=count, color=Affected))+
  geom_point(size=3)+
  geom_boxplot(alpha=0.25)+
  theme_bw()+
  ggtitle("HAUS4_paralog - no quality control")+
  ylab("Count")+
  theme(axis.title=element_text(size=12, face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))