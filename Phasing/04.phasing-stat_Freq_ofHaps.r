

library("Biostrings")
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

samps.c <- c()
R.c <- c()
Pvalue <- c()

for (person in names(Pers.refs)) {
#for (person in c("P24")) {
  print(person)
  ref.samp <- Pers.refs[[person]]
  SNV.table <- read.table(paste0( "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/",person,".reliable.locus.CCS.table.txt"),
                          sep = "\t",quote = "",header = T)
  # ref genenome
  ref.gnmP <- "/disk/lhj_4T/lhj/HP/process/assembly/all_genomes/"
  ref.gnmF <- paste0(ref.gnmP,ref.samp,".fasta")
  ref.gnm <-  readDNAStringSet(ref.gnmF)[[1]]  #  read_fasta(ref.gnmF)[1]
  Bases <- list()
  for (SNVs.site in SNV.table$X.Posi) {
    ref.base <-  subseq(ref.gnm , start = as.numeric(SNVs.site),end = as.numeric(SNVs.site))
    Bases[[ paste0("C_",SNVs.site) ]] <- ref.base
  }
  Bases.df <- Bases %>% as.data.frame() %>% t()
  colnames(Bases.df) <- "Ref.base" 
  Bases.df <- Bases.df %>% as.data.frame()
  Bases.df$X.Posi.1 <- rownames(Bases.df)
  Bases.df$X.Posi <- str_split_fixed(Bases.df$X.Posi.1,"_",2)[,2]
  Bases.df$X.Posi.1 <- NULL
  SNV.table <- merge(SNV.table,Bases.df)
  
  
  #for (samp in c("P1-E-d")) {
  for (samp in samp_list[[person]]) {
    Hap.OP <- paste0("/home/try/Documents/projects/liuhj-Other/phasing-2024/blocks-phasing/",samp,"/") 
    
    relaibLocus.P <- "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/"
    block.table <- read.table(paste0( relaibLocus.P,samp,".blocks.txt"),
                              sep = "\t",quote = "",header = T)
    SNVs.Freq.smp <- SNV.table[,c("X.Posi", gsub("-",".",paste0(samp,"_2_",ref.samp)),"Ref.base")]
    
    ## read hap files
    Phased.site.Freqs <- list()
    #for (n.block in seq(3,3) ) {
    for (n.block in seq(1,nrow(block.table)) ) {
        #for (n.block in seq(5,5) ) {
      print( paste0("block : ",n.block) )
      B <- block.table[n.block,"block"]
      B.S <- block.table[n.block,"Block.Start"]
      B.E <- block.table[n.block,"Block.End"]
      ##phased
      PhasedHaps.F1 <- paste0( Hap.OP,samp,".",B,".",B.S,"-",B.E,".fasta")
      PhasedHaps.F2 <- paste0(Hap.OP ,samp,".",B,".",B.S,"-",B.E,".ReadHaps.fasta")
      block.SNV.sites <- SNVs.Freq.smp %>% filter(X.Posi >= B.S,
                                                  X.Posi <= B.E)
      block.sites <- block.SNV.sites$X.Posi
      if (file.exists(PhasedHaps.F1)) {
        PhasedHaps.F <- PhasedHaps.F1
        fasta_data <- read_fasta(PhasedHaps.F)
        #Block.Freq <- list()
        for (Hpid in names(fasta_data)) {
          HapFreq <- strsplit( strsplit(Hpid,":")[[1]][2],"_")[[1]][1] %>% as.numeric()
          Hp.seq <- fasta_data[[Hpid]]
          for (Hp.seq.idx in seq(1,nchar(Hp.seq))) {
            idx.Locus <- block.sites[Hp.seq.idx] %>% as.character()
            idx.Base <- substr(Hp.seq,Hp.seq.idx,Hp.seq.idx)
            if(!idx.Locus %in% names(Phased.site.Freqs)){
              Phased.site.Freqs[[idx.Locus]] <- list()
              Phased.site.Freqs[[idx.Locus]][[idx.Base]] <- HapFreq
            }
            else{
              if(!idx.Base %in% names(Phased.site.Freqs[[idx.Locus]]) ){
                Phased.site.Freqs[[idx.Locus]][[idx.Base]] <- HapFreq
              }
              else{
                Phased.site.Freqs[[idx.Locus]][[idx.Base]] <- Phased.site.Freqs[[idx.Locus]][[idx.Base]] + HapFreq 
              }
            }
          }
        }
        
      } 
      else{
        PhasedHaps.F <- PhasedHaps.F2
        fasta_data <- read_fasta(PhasedHaps.F)
        HapFreq <- 1
        Hp.seq <- fasta_data[[names(fasta_data)[1]]]
        for (Hp.seq.idx in seq(1,nchar(Hp.seq))) {
          idx.Locus <- block.sites[Hp.seq.idx] %>% as.character()
          idx.Base <- substr(Hp.seq,Hp.seq.idx,Hp.seq.idx)
          if(!idx.Locus %in% names(Phased.site.Freqs)){
            Phased.site.Freqs[[idx.Locus]] <- list()
            Phased.site.Freqs[[idx.Locus]][[idx.Base]] <- HapFreq
          }
      }
    }
    }
    
    Mut.list <- list()
    for (ps.cit in names(Phased.site.Freqs)) {
      Bss <- names( Phased.site.Freqs[[ps.cit]] ) 
      Ref.bas <- SNVs.Freq.smp[SNVs.Freq.smp$X.Posi == as.numeric(ps.cit),"Ref.base"]
      for (Bs in Bss) {
        if(Bs != Ref.bas & Bs != "Z"){
          if(!ps.cit %in% Mut.list){
            Mut.list[[ps.cit]] <- Phased.site.Freqs[[ps.cit]][[Bs]]
          } else{
            Mut.list[[ps.cit]] <- Mut.list[[ps.cit]] + Phased.site.Freqs[[ps.cit]][[Bs]]
          }
        }
      }
    }
    for (ps.cit in names(Phased.site.Freqs)) {
      if(!ps.cit %in% names(Mut.list)){
        Mut.list[[ps.cit]] <- 0
      }
    }
    
    Hap.mut <- Mut.list %>% as.data.frame() %>% t() %>% as.data.frame()
    colnames(Hap.mut) <- "Freq.FromHaps"
    Hap.mut$X.Posi <- rownames(Hap.mut) 
    Hap.mut$X.Posi <- str_split_fixed(Hap.mut$X.Posi,"X",2)[,2]
    
    SNVs.Freq.smp.2 <-  merge(SNVs.Freq.smp,Hap.mut)
    colnames(SNVs.Freq.smp.2) <- c("X.Posi","Freq.raw","base.ref","Freq.phased")
    SNVs.Freq.smp.2[SNVs.Freq.smp.2 == "NO"] <- 0
    SNVs.Freq.smp.2$Freq.raw <- SNVs.Freq.smp.2$Freq.raw %>% as.numeric()
    SNVs.Freq.smp.2 <- SNVs.Freq.smp.2 %>% mutate(diff = Freq.phased - Freq.raw)
    openxlsx::write.xlsx(SNVs.Freq.smp.2,paste0("/home/try/Documents/projects/liuhj-Other/phasing-2024/blocks-phasing/",samp,".phased.Freq.table.xlsx"))
    
    SNVs.Freq.smp.2$X.Posi <- NULL
    SNVs.Freq.smp.2$base.ref <- NULL
    library("ggplot2")
    library(ggpubr)
    p <- ggplot(SNVs.Freq.smp.2, aes(x = Freq.raw, y = Freq.phased)) +
      #geom_line(stat = "summary", fun = mean, linewidth = 1) +  # Line for mean values
      geom_point(stat = "summary", fun = mean, size = 1,alpha = 0.3) +  # Points for mean values
      theme_minimal() +
      labs(title = "", x = "Freq.raw", y = "Freq.phased") +
      xlim(0,1) +
      ylim(0,1) + 
      theme(panel.border = element_rect(color = "black",fill = NA,linewidth = 1)) 
    p <- p + stat_cor(method = "pearson", label.x = 0.1, label.y = 0.8)
    
    corr.R <- cor.test(SNVs.Freq.smp.2$Freq.raw, SNVs.Freq.smp.2$Freq.phased)
    
    samps.c <- c(samps.c,samp)
    R.c <- c(R.c,corr.R$estimate)
    Pvalue <- c(Pvalue,corr.R$p.value)
    

    ggsave(p,filename = paste0(paste0("/home/try/Documents/projects/liuhj-Other/phasing-2024/blocks-phasing/",samp,".phased.Freq.jpg") ),width = 3,height = 3)
  }
}


corr.df <-  data.frame("Sample" = samps.c,
           "Pearson.r" = R.c,
           "P.value" = Pvalue)
openxlsx::write.xlsx(corr.df,paste0("/home/try/Documents/projects/liuhj-Other/phasing-2024/blocks-phasing/AllSamples.phased.Freq.table.xlsx"))

