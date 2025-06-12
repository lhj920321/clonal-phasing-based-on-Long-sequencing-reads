rm(list = ls())

library("Biostrings")
library("dplyr")
library("ggplot2")

## read haplotype
Pers.refs <- list()
Pers.refs[["P1"]] <- "P1-E-j"
#Pers.refs[["P2"]] <- "P2-E-t"
Pers.refs[["P4"]] <- "P4-E-j"
Pers.refs[["P5"]] <- "P5-E-d"
Pers.refs[["P6"]] <- "P6-E-t"
Pers.refs[["P7"]] <- "P7-E-j"
Pers.refs[["P18"]] <- "P18-C_bR-x"
Pers.refs[["P19"]] <- "P19-C_bR-j"
Pers.refs[["P20"]] <- "P20-C_bR-t"
Pers.refs[["P21"]] <- "P21-C_sL-j"
#Pers.refs[["P22"]] <- "P22-C_sL-d"
Pers.refs[["P23"]] <- "P23-C_sL-j"
Pers.refs[["P24"]] <- "P24-C_sC-t"
Pers.refs[["P25"]] <- "P25-C_sC-t"

samp_list <- list()
samp_list[["P1"]] <- c("P1-E-d","P1-E-x","P1-E-j","P1-E-t")
samp_list[["P2"]] <- c("P2-E-t")
samp_list[["P4"]] <- c("P4-E-d","P4-E-x","P4-E-j")
samp_list[["P5"]] <- c("P5-E-d")
samp_list[["P6"]] <- c("P6-E-d","P6-E-t")
samp_list[["P7"]] <- c("P7-E-d","P7-E-j")
samp_list[["P18"]] <- c("P18-C_bR-x","P18-C_bR-j")
samp_list[["P19"]] <- c("P19-C_bR-d","P19-C_bR-j","P19-C_bR-t")
samp_list[["P20"]] <- c("P20-C_bR-d","P20-C_bR-x","P20-C_bR-j","P20-C_bR-t")
samp_list[["P21"]] <- c("P21-C_sL-d","P21-C_sL-x","P21-C_sL-j","P21-C_sL-t")
samp_list[["P22"]] <- c("P22-C_sL-d","P22-C_sL-t")
samp_list[["P23"]] <- c("P23-C_sL-j","P23-C_sL-x")
samp_list[["P24"]] <- c("P24-C_sC-x","P24-C_sC-t")
samp_list[["P25"]] <- c("P25-C_sC-j","P25-C_sC-t")

Happlot.path <- "/home/try/Documents/projects/liuhj-Other/phasing-2024/blocks-phasing/Phased-Haps-plot/"

##revombi
phased.Hap.f <- paste0(Happlot.path,person,".phased.Hap.forPlot.xlsx")
phased.Hapspltdf <- openxlsx::read.xlsx(phased.Hap.f)
phased.Hapspltdf[,c("Start","End","Start.Posi","End.Posi","line.whith.Freq")] <- lapply(phased.Hapspltdf[,c("Start","End","Start.Posi","End.Posi","line.whith.Freq")],as.numeric)
phased.Hapspltdf <- phased.Hapspltdf %>% mutate(Length = End.Posi-Start.Posi) %>% 
  mutate(Locus.n = End-Start + 1)

Chuck.long <- phased.Hapspltdf %>% filter(Locus.n >= 3,
                                          Comp_to_RefMajorHap == "NS") 

median()
p.hist.Len <-  ggplot(Chuck.long, aes(x = Length)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", color = "black", alpha = 0.7) +  # 绘制直方图
  labs(
    title = "",
    x = "Length of chruck in haplotype with recombination",  
    y = "No. of chruck"  
  ) +
  #ylim(5,12) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        #panel.grid.major.y = element_line(size=0.25, color="#D2D3D3"),
        #axis.text.x = element_text(angle=90, size=5),
        axis.text.x = element_text(angle=90, size=7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        legend.position="none")
ggsave(p.hist.Len,filename = paste0(Happlot.path,"phased.Hap-recomb.jpg"),width = 4 ,height = 2.3)
ggsave(p.hist.Len,filename = paste0(Happlot.path,"phased.Hap-recomb.pdf"),width = 4 ,height = 2.3)




all.snv.c <- c()
chruk.snv.c <- c()
samp.c <- c()
pers.c <- c()

for (person in names(Pers.refs)) {
#for (person in c("P1")) {
  print(person)
  ref.samp <- Pers.refs[[person]]
  sort.Samp <- c(ref.samp,setdiff(samp_list[[person]],ref.samp) )
  relaibLocus.P <- "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/"
  block.table <- read.table(paste0( relaibLocus.P,ref.samp,".blocks.txt"),
                            sep = "\t",quote = "",header = T)
  
 
  
  ## recombi region, SNV percent 
  SNV.table <- read.table(paste0( "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/",person,".reliable.locus.CCS.table.txt"),
                          sep = "\t",quote = "",header = T)
  for (samp in sort.Samp) {
    SNVs.Freq.smp <- SNV.table[,c("X.Posi", gsub("-",".",paste0(samp,"_2_",ref.samp)))]
    colnames(SNVs.Freq.smp) <- c("X.Posi","Freq")
    SNV.locus <- SNVs.Freq.smp %>% filter(Freq != "NO")
    SNV.locus[,colnames(SNV.locus)] <- lapply( SNV.locus[,colnames(SNV.locus)],as.numeric)
    Tot.SNV <- SNV.locus %>% filter(Freq >= 0.02)
    
    all.Chuck.SNVsites <- c()
    for (Chuck.n in seq(1,nrow(Chuck.long)) ) {
      Chuck.S <- Chuck.long[Chuck.n,"Start.Posi"]
      Chuck.E <- Chuck.long[Chuck.n,"End.Posi"]
      Chuck.SNV.sites <- Tot.SNV %>% filter(X.Posi >= Chuck.S,
                                                  X.Posi <= Chuck.E)
      all.Chuck.SNVsites <- c(all.Chuck.SNVsites,Chuck.SNV.sites$X.Posi)
    }
    
    
    n.all.SNV <-  length(Tot.SNV$X.Posi)
    n.Chuck.SNV <-  length(unique(all.Chuck.SNVsites))
      
    
    all.snv.c <- c(all.snv.c, n.all.SNV)
    chruk.snv.c <- c(chruk.snv.c,n.Chuck.SNV)
    samp.c <- c(samp.c,samp)
    pers.c <- c(pers.c,person)

  }
}


chruck.snv.df <- data.frame("person" = pers.c,
                            "sample" = samp.c,
                            "all.snv.No" = all.snv.c,
                            "chruk.snv.No" = chruk.snv.c)
chruck.snv.df <- chruck.snv.df %>% mutate(perent.chruk.snv = round(chruk.snv.No*100/all.snv.No,2 )) %>% 
  mutate(X.kilo = round(all.snv.No/1000,4))


p.point <- ggplot(chruck.snv.df, aes(x = X.kilo, y = perent.chruk.snv)) + # 
  geom_point(alpha=0.8, size= 1,shape = 1) +  
  geom_smooth(method = "loess",color = "red",size = 0.5) +
  scale_color_manual(values = c("blue","black")) + 
  geom_vline(xintercept=c(0.5),col="lightgrey",lwd=0.5) + 
  scale_x_continuous(breaks = seq(0,10) ) +
  #geom_hline(yintercept = Neg.beta.cutoff, lty="dotdash",col="black",lwd=0.8) +  
  #geom_text_repel(aes(label = label)) +
  labs(x="SNV no. per sample(kilo)", y="SNVs within chrunks(%)") +  #x、y轴标签
  #ggtitle(sampFlag) + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank(),
        panel.grid.major = element_blank()) +
  ylim(0,100)

ggsave(p.point,filename = paste0(Happlot.path,"phased.Hap-chrunk.SNV.jpg"),width = 4 ,height = 2.3)
ggsave(p.point,filename = paste0(Happlot.path,"phased.Hap-chrunk.SNV.pdf"),width = 4 ,height = 2.3)



