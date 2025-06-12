rm(list = ls())

library("Biostrings")
library("dplyr")
library("ggplot2")
# Function to read a FASTA file
read_fasta <- function(file_path) {
  lines <- readLines(file_path)  # Read all lines from the file
  sequences <- list()            # Initialize a list to store sequences
  header <- NULL                 # Initialize a variable to store the current header
  for (line in lines) {
    if (startsWith(line, ">")) {  # Check if the line is a header
      header <- sub("^>", "", line)  # Remove the ">" from the header
      sequences[[header]] <- ""      # Initialize an empty sequence for this header
    } else {
      sequences[[header]] <- paste0(sequences[[header]], line)  # Append the sequence
    }
  }
  return(sequences)  # Return the list of sequences
}


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


O.path <- "/home/try/Documents/projects/liuhj-Other/phasing-2024/blocks-phasing/Phased-Haps-plot/"

for (person in names(Pers.refs)) {
#for (person in c("P24")) {
  print(person)
  block.color <- data.frame()
  
  ref.samp <- Pers.refs[[person]]
  SNV.table <- read.table(paste0( "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/",person,".reliable.locus.CCS.table.txt"),
                          sep = "\t",quote = "",header = T)
  
  sort.Samp <- c(ref.samp,setdiff(samp_list[[person]],ref.samp) )
  

  ##ref sample major Haps 
  Hap.OP <- paste0("/home/try/Documents/projects/liuhj-Other/phasing-2024/blocks-phasing/",ref.samp,"/") 
  relaibLocus.P <- "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/"
  block.table <- read.table(paste0( relaibLocus.P,ref.samp,".blocks.txt"),
                            sep = "\t",quote = "",header = T)
  ##### read hap files, Ref.majorHap
  Ref.majorHap <- data.frame()
  for (n.block in seq(1,nrow(block.table)) ) {
    #print( paste0("block : ",n.block) )
    B <- block.table[n.block,"block"]
    B.S <- block.table[n.block,"Block.Start"]
    B.E <- block.table[n.block,"Block.End"]
    SNVs.Freq.smp <- SNV.table[,c("X.Posi", gsub("-",".",paste0(ref.samp,"_2_",ref.samp)))]
    block.SNV.sites <- SNVs.Freq.smp %>% filter(X.Posi >= B.S,
                                                X.Posi <= B.E)
    block.sites <- block.SNV.sites$X.Posi
    
    # for (block.site in block.sites) {
    #   Ref.majorHap[[as.character(block.site) ]] <- list()
    # }
    ##phased
    PhasedHaps.F1 <- paste0( Hap.OP,ref.samp,".",B,".",B.S,"-",B.E,".fasta")
    PhasedHaps.F2 <- paste0(Hap.OP ,ref.samp,".",B,".",B.S,"-",B.E,".ReadHaps.fasta")
    
    if (file.exists(PhasedHaps.F1)) {
      PhasedHaps.F <- PhasedHaps.F1
      fasta_data <- read_fasta(PhasedHaps.F)
      majorHap <-  names(fasta_data)[1]
      majorHap.Freq <- strsplit( strsplit(majorHap,":")[[1]][2],"_" )[[1]][1] %>% as.numeric()
      if( length(names(fasta_data)) >=2  ){
        for (OtheHap in seq(2,length(names(fasta_data)) ) ) {
          otherHap <-  names(fasta_data)[OtheHap]
          otherHap.Freq <- strsplit( strsplit(otherHap,":")[[1]][2],"_" )[[1]][1] %>% as.numeric()
          if(otherHap.Freq > majorHap.Freq){
            majorHap.Freq <- otherHap.Freq
            majorHap <- otherHap
          }
        }
      }
      MajorBase <- fasta_data[[majorHap]]
    } else{
      PhasedHaps.F <- PhasedHaps.F2
      fasta_data <- read_fasta(PhasedHaps.F)
      majorHap <-  names(fasta_data)[1]
      MajorBase <- fasta_data[[majorHap]]
    }
    
    ## Ref samp major Hap
    MajorSeq <- unlist( strsplit(MajorBase,"")  )
    Majordf <- data.frame(
      posi = block.sites,  
      base = MajorSeq
    )
    colnames(Majordf)[2] <- "Major.Hap"
    if( nrow(Ref.majorHap) == 0 ){ 
      Ref.majorHap <-  Majordf 
    }else{
      Ref.majorHap <- rbind(Ref.majorHap,Majordf)
    }
  }
  

  
  ##
  Samp.y.Start <- 1
  plt.Haps <- data.frame()
  for (samp in sort.Samp) {
    Hap.OP <- paste0("/home/try/Documents/projects/liuhj-Other/phasing-2024/blocks-phasing/",samp,"/") 
    relaibLocus.P <- "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/"
    block.table <- read.table(paste0( relaibLocus.P,samp,".blocks.txt"),
                              sep = "\t",quote = "",header = T)
    ## read hap files
    Phased.site.Base <- list()
    y.val.min <- Samp.y.Start
    block.HapOlst <- list()
    xuni.breks <- c()
    #for (n.block in seq(1,1) ) {
    for (n.block in seq(1,nrow(block.table)) ) {
      #print( paste0("block : ",n.block) )
      B <- block.table[n.block,"block"]
      B.S <- block.table[n.block,"Block.Start"]
      B.E <- block.table[n.block,"Block.End"]
      
      SNVs.Freq.smp <- SNV.table[,c("X.Posi", gsub("-",".",paste0(samp,"_2_",ref.samp)))]
      SNVs.Freq.smp$xuniLocus <- rownames(SNVs.Freq.smp)
      block.SNV.sites <- SNVs.Freq.smp %>% filter(X.Posi >= B.S,
                                                  X.Posi <= B.E)
      block.sites <- block.SNV.sites$X.Posi
      
      for (block.site in block.sites) {
        Phased.site.Base[[as.character(block.site) ]] <- list()
      }

      ##phased
      PhasedHaps.F1 <- paste0( Hap.OP,samp,".",B,".",B.S,"-",B.E,".fasta")
      PhasedHaps.F2 <- paste0(Hap.OP ,samp,".",B,".",B.S,"-",B.E,".ReadHaps.fasta")
      
      ##major hap  & base stat for each cite
      Hap.df <- data.frame()
      Hap.Nmlst <- c()
      Hap.Freq <- list()
      #FreqSortLst <- list()
      if (file.exists(PhasedHaps.F1)) {
        PhasedHaps.F <- PhasedHaps.F1
        fasta_data <- read_fasta(PhasedHaps.F)
        majorHap <-  names(fasta_data)[1]
        majorHap.Freq <- strsplit( strsplit(majorHap,":")[[1]][2],"_" )[[1]][1] %>% as.numeric()
        Hap.Nmlst <- c(Hap.Nmlst,majorHap)
        Hap.Freq[[majorHap]] <- majorHap.Freq
        #FreqSortLst[[as.character(majorHap.Freq) ]] <- c(majorHap)
        if( length(names(fasta_data)) >=2  ){
          for (OtheHap in seq(2,length(names(fasta_data)) ) ) {
            otherHap <-  names(fasta_data)[OtheHap]
            otherHap.Freq <- strsplit( strsplit(otherHap,":")[[1]][2],"_" )[[1]][1] %>% as.numeric()
            Hap.Nmlst <- c(Hap.Nmlst,otherHap)
            Hap.Freq[[otherHap]] <- otherHap.Freq
            #FreqSortLst[[as.character(majorHap.Freq)]] <- c(FreqSortLst[[as.character(majorHap.Freq)]],otherHap)
            if(otherHap.Freq > majorHap.Freq){
              majorHap.Freq <- otherHap.Freq
              majorHap <- otherHap
            }
          }
        }
      } else{
        PhasedHaps.F <- PhasedHaps.F2
        fasta_data <- read_fasta(PhasedHaps.F)
        majorHap <-  names(fasta_data)[1]
        majorHap.Freq <- strsplit( strsplit(majorHap,":")[[1]][2],"_" )[[1]][1] %>% as.numeric()
        Hap.Nmlst <- c(Hap.Nmlst,majorHap)
        Hap.Freq[[majorHap]] <- majorHap.Freq
        #FreqSortLst[[as.character(majorHap.Freq)]] <- c(majorHap)
      }
      
      Freqs <- list()
      Freq.c <- c()
      for (Hhp in names(Hap.Freq)) {
        fq <- Hap.Freq[[Hhp]]
        Freq.c <- c(Freq.c,fq) %>% unique()
        if(!fq %in% Freqs){
          Freqs[[ as.character(fq) ]] <- c(Hhp)
        } else{
          Freqs[[ as.character(fq) ]] <- c(Freqs[[fq]],Hhp)
        }
      }
      Freq.c <- sort(Freq.c,decreasing = T)
      
      
      ##major hap, bases 
      for (Hap in names(fasta_data)) {
        baseSeq <- unlist( strsplit(fasta_data[[Hap]],"")  )
        df <- data.frame(
          posi = block.sites,  # seq_along(baseSeq),
          base = baseSeq
        )
        colnames(df) <- c("citeIdx",Hap)
        if(nrow(Hap.df) == 0){
          Hap.df <- df
        }else{
          Hap.df <- merge(Hap.df,df)
        }
      }
      
      ##count for A G C T
      Hap.df.s <- Hap.df
      Hap.df.s$A_count <- rowSums(data.frame(Hap.df.s[, Hap.Nmlst])  == "A")
      Hap.df.s$G_count <- rowSums(data.frame(Hap.df.s[, Hap.Nmlst]) == "G")
      Hap.df.s$C_count <- rowSums(data.frame(Hap.df.s[, Hap.Nmlst]) == "C")
      Hap.df.s$T_count <- rowSums(data.frame(Hap.df.s[, Hap.Nmlst]) == "T")
      
      stat.cols <- c("A_count","G_count","C_count","T_count")
      Hap.df.s$zero.Count <- rowSums(Hap.df.s[, stat.cols] == 0)
      
      ##mut or not
      Hap.df.s <- Hap.df.s %>% mutate(Mut_or_not = case_when(zero.Count == 3 ~ "No",
                                                             TRUE ~ "Mut"))
      Hap.df.s <-  Hap.df.s %>% dplyr::select(-c("A_count","G_count","C_count","T_count"))
      #MjbaseSeq <- unlist( strsplit(fasta_data[[majorHap]],"")  )
      ##add major Hap bases
      #Hap.df.s$Major.Hap <- MjbaseSeq
      Hap.df.s <- merge(Hap.df.s,Ref.majorHap,by.x = "citeIdx",by.y = "posi",all.x = T)
      
      
      # site. base         
      for (block.site in block.sites) {
        MutNot <-  Hap.df.s[Hap.df.s$citeIdx == as.numeric(block.site),"Mut_or_not"] 
        BASE <-  Hap.df.s[Hap.df.s$citeIdx == as.numeric(block.site),"Major.Hap"] 
        if(MutNot == "No"){
          Phased.site.Base[[ as.character(block.site) ]][[BASE]] <- 1
        }
      }
      
      ## sort Hap by Freq
      Hap.df.s <- merge(Hap.df.s,SNVs.Freq.smp,by.x = "citeIdx",by.y = "X.Posi",all.x = T)
      xuni.brek <- Hap.df.s[1,"xuniLocus"] %>% as.numeric()
      xuni.breks <- c(xuni.breks,xuni.brek)
      Phased.block.plt <- data.frame()
      y.val <- Samp.y.Start
      
      for (Freq.idx in seq(1,length(Freq.c)) ) {
        Freq <- Freq.c[Freq.idx]
        Freq.haps <- Freqs[[as.character(Freq) ]]
        for (Freq.hap in Freq.haps) {
          y.val <- y.val - 1
          if(y.val < y.val.min){y.val.min <- y.val}
          Hap.df.s <- Hap.df.s %>% mutate(Base = case_when(.data[[Freq.hap]] == Major.Hap ~ "S",
                                                           TRUE ~ "NS") )
          Hap.df.s$Freq.f <- Freq
           # Hap.df.s <- Hap.df.s %>% mutate(Freq.f = case_when(Mut_or_not == "No" ~ 1,
          #                                                    TRUE ~ Freq))
          bas.s <- Hap.df.s$Base
          rle_result <- rle(bas.s)
          result <- character(length(bas.s))
          start <- 1
          for (i in seq_along(rle_result$lengths)) {
            end <- start + rle_result$lengths[i] - 1
            label <- paste0( Hap.df.s$xuniLocus[start],"-",Hap.df.s$xuniLocus[end], ":", rle_result$values[i],"-",Hap.df.s$Freq.f[end],":",Hap.df.s$citeIdx[start],"-",Hap.df.s$citeIdx[end])
            result[start:end] <- label
            start <- end + 1
          }
          Rg.result <- result %>% unique()
          split_result <- strsplit(Rg.result, "[:-]")
          plt.df <- data.frame(
            Start = sapply(split_result, function(x) x[1]),  
            End = sapply(split_result, function(x) x[2]), 
            Start.Posi = sapply(split_result, function(x) x[5]),  
            End.Posi = sapply(split_result, function(x) x[6]), 
            Comp_to_RefMajorHap = sapply(split_result, function(x) x[3]),
            line.whith.Freq = sapply(split_result, function(x) x[4]) ) 
          plt.df$Block <- B
          plt.df$Y <- y.val
          plt.df$sample <- samp
          block.color <- rbind(block.color,plt.df)
          
          colnames(Hap.df.s)[ncol(Hap.df.s)-1] <- paste0("Base__",Freq.hap)
          colnames(Hap.df.s)[ncol(Hap.df.s)] <- paste0("Freq__",Freq.hap)

        }
      }
      block.HapOlst[[paste0("block.",n.block)]] <- Hap.df.s
    } #block
    Samp.y.Start <- y.val.min - 1
    openxlsx::write.xlsx(block.HapOlst,paste0(O.path,person,".",samp,".phased.Hap.Info.xlsx"))
    
  } # sample
  
  
  
  
  openxlsx::write.xlsx(block.color,paste0(O.path,person,".phased.Hap.forPlot.xlsx"))
  
  block.color[,c("Start","End","Start.Posi","End.Posi","line.whith.Freq")] <- lapply(block.color[,c("Start","End","Start.Posi","End.Posi","line.whith.Freq")],as.numeric)
  SampBlockplot <- ggplot(data=block.color,aes(x=Start,y=Y,fill=factor(Comp_to_RefMajorHap))) +   #
    geom_rect(aes(xmin=Start, xmax=End,ymin=Y-0.5*line.whith.Freq,ymax=Y+ 0.5*line.whith.Freq),
              colour="NA") +
    #scale_x_continuous(breaks = xticksData$xticks, labels = xticksData$absPos ) +
    scale_fill_manual(name="compare with ref major haplotype",
                      labels = c("NS" = "Not same","S"="Same" ),
                      values = c("NS" = "#0000FF","S"="#FFCC00")) +
    xlab("block region") +
    ylab("") +
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
          legend.position="none")+
    geom_vline(xintercept = unique(xuni.breks),color = "lightgrey",linewidth = 0.1,linetype = "dashed")

  ggsave(SampBlockplot,filename = paste0(O.path,person,".phased.Hap-Plot.jpg"),width = 6 ,height = 1 + 0.06 * -y.val.min)
  ggsave(SampBlockplot,filename = paste0(O.path,person,".phased.Hap-Plot.pdf"),width = 6 ,height = 1 + 0.06 * -y.val.min)
  
  
  
  SampBlockplot2 <- ggplot(data=block.color,aes(x=Start.Posi,y=Y,fill=factor(Comp_to_RefMajorHap))) +   #
    geom_rect(aes(xmin=Start.Posi, xmax=End.Posi,ymin=Y-0.5*line.whith.Freq,ymax=Y+ 0.5*line.whith.Freq),
              colour="NA") +
    #scale_x_continuous(breaks = xticksData$xticks, labels = xticksData$absPos ) +
    scale_fill_manual(name="compare with ref major haplotype",
                      labels = c("NS" = "Not same","S"="Same" ),
                      values = c("NS" = "#0000FF","S"="#FFCC00")) +
    xlab("block region") +
    ylab("") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),  #element_line(linewidth=0.25, color="grey"),
          #panel.grid.major.y = element_line(size=0.25, color="#D2D3D3"),
          #axis.text.x = element_text(angle=90, size=5),
          axis.text.x = element_text(angle=90, size=7),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          legend.position="none") +
    geom_vline(xintercept = unique(xuni.breks),color = "lightgrey",linewidth = 0.1,linetype = "dashed")
  
  ggsave(SampBlockplot2,filename= paste0(O.path,person,".phased.Hap-Plot.realCite.jpg") ,width = 6,height = 1 + 0.06 * -y.val.min)
  ggsave(SampBlockplot2,filename= paste0(O.path,person,".phased.Hap-Plot.realCite.pdf") ,width = 6,height = 1 + 0.06 * -y.val.min)
  
}  #person






