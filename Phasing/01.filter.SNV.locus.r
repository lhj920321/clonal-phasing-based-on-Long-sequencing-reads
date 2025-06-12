
library("dplyr")

Persons <- c(1,2,4,5,6,7,seq(18,25))
for (persid in Persons) {
  person <- paste0("P",persid)
  print(person)
  #arg <- commandArgs(T)
  ## CCS
  CCSF <- paste0("/disk/lhj_4T/lhj/HP/process/CCS_map2PersonRespectSamp/table/",person,"/all.iSNV_with_SNP.pyResults.txt") 
  CCS.df <-  read.table(CCSF,sep = "\t",comment.char = "",header =  T)
  samps.CCS <- colnames(CCS.df)[3:ncol(CCS.df)]
  
  ## NGS
  NGSF <- paste0("/disk/lhj_4T/lhj/HP/process/NGS_map2PersonResptGnm/table/",person,"/all.iSNV_with_SNP.pyResults.txt")
  NGS.df <-  read.table(NGSF,sep = "\t",comment.char = "",header =  T)
  #samps.NGS <- colnames(NGS.df)[3:ncol(CCS.df)]
  
  Ref.samp <- strsplit(samps.CCS[1],"_2_")[[1]][2]
  Ref.samp <- gsub("[.]","-",Ref.samp)
  
  
  # ##repeat region positive
  # RepeatRGs.F <-  paste0("/disk/lhj_4T/lhj/HP/process/Genome_repeatRegion/RepeatRegions/",Ref.samp,".repeat-regions-Stat.txt")
  # RepeatRGs <-  read.table(RepeatRGs.F,sep = "\t",header =  T)
  
  ## reliable locus of CCS and NGS
  kaopu.locus <- list()
  share.locus <- c()
  for (CCS.smp in samps.CCS) {
    CCS.subdf <-  CCS.df[,c("X.Posi",CCS.smp)] 
    NGS.subdf <-  NGS.df[,c("X.Posi",CCS.smp)] 
    colnames(CCS.subdf) <- c("X.Posi","CCS")
    colnames(NGS.subdf) <- c("X.Posi","NGS")
    mg.df <-  merge(CCS.subdf,NGS.subdf,all = T) %>% na.omit()
    mg.df[mg.df=="NO"] <- "0"
    mg.df$CCS <- as.numeric(mg.df$CCS)
    mg.df$NGS <- as.numeric(mg.df$NGS)
    mg.df <- mg.df %>% mutate(diff = CCS - NGS) %>% filter(abs(diff) < 0.2)
    if(length(share.locus) == 0){
      share.locus <- mg.df$X.Posi
    }else{
      share.locus <- intersect(share.locus,mg.df$X.Posi)
    }
    kaopu.locus[[CCS.smp]] <- mg.df
  }
  
  filter.table <- CCS.df %>% filter(X.Posi %in% share.locus)
  filter.table$new.locus <- rownames(filter.table)
  
  write.table(filter.table,paste0( "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/",person,".reliable.locus.CCS.table.txt"),
              sep = "\t",quote = F,row.names = F)
}



