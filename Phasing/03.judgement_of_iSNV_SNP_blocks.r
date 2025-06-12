library("dplyr")



readHap.P <- "/home/try/Documents/projects/liuhj-Other/phasing-2024/haplotype-read/"
# samp <- "P1-E-t"
# ref.samp <- "P1-E-j"
# person <- "P1"


Pers.refs <- list()
Pers.refs[["P1"]] <- "P1-E-j"
Pers.refs[["P2"]] <- "P2-E-t"
Pers.refs[["P4"]] <- "P4-E-j"
Pers.refs[["P5"]] <- "P5-E-d"
Pers.refs[["P6"]] <- "P6-E-t"
Pers.refs[["P7"]] <- "P7-E-j"
Pers.refs[["P18"]] <- "P18-C_bR-x"
Pers.refs[["P19"]] <- "P19-C_bR-j"
Pers.refs[["P20"]] <- "P20-C_bR-t"
Pers.refs[["P21"]] <- "P21-C_sL-j"
Pers.refs[["P22"]] <- "P22-C_sL-d"
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


#for (person in names(Pers.refs)) {
for (person in c("P24")) {
    print(person)
  ref.samp <- Pers.refs[[person]]
  for (samp in samp_list[[person]]) {
    print(samp)
    readHap.F <- paste0(readHap.P,"readHap.",samp,"_2_",ref.samp,"/2.read-haplotype.Stat.txt") 
    readHaps <- read.table(readHap.F,sep = "\t",quote = "",header = F)
    colnames(readHaps) <- c("count","S.E","haplotype")
    readHaps$id <- rownames(readHaps)
    readHaps <- readHaps[ !is.na(readHaps$haplotype),]
    
    ##SNV locus
    SNV.table <- read.table(paste0( "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/",person,".reliable.locus.CCS.table.txt"),
                            sep = "\t",quote = "",header = T)
    library("stringr")
    ## read haplotype
    
    readHaps$id <- rownames(readHaps)
    readHaps$S <- str_split_fixed(readHaps$S.E, '-', 3)[,2]
    readHaps$E <- str_split_fixed(readHaps$S.E, '-', 3)[,3]
    readHaps$count <- as.numeric(readHaps$count)
    readHaps$S <- as.numeric(readHaps$S)
    readHaps$E <- as.numeric(readHaps$E)
    
    locus.list <- list()
    locus.depth <- list()
    for (Locus in SNV.table$X.Posi) {
      #for (Locus in c("11431")) {
      Locus <- as.numeric(Locus)
      right.df <- readHaps %>% filter(S <=  Locus  & E >=  Locus) 
      locus.list[[paste0("P.",Locus)]] <- c(right.df$id)
      locus.depth[[paste0("P.",Locus)]] <- sum(right.df$count)
    }
    
    
    
    #### blocks phasing
    minSharedHapDeep <- 10
    blockLst <- list()
    for (n in seq(1,(nrow(SNV.table)-1) )) {
      #for (n in seq(93,93)) {
      locus <- SNV.table[n,"X.Posi"]
      Next.locus <- SNV.table[n+1,"X.Posi"]
      Flg.locus <- paste0("P.",locus)
      Flg.Nextlocus <- paste0("P.",Next.locus)
      readHaps.locus <- locus.list[[Flg.locus]]
      readHaps.Nextlocus <- locus.list[[Flg.Nextlocus]]
      shared.Read <-  intersect(readHaps.locus,readHaps.Nextlocus)
      
      shared.df <- readHaps %>% filter(id %in% shared.Read) 
      shared.dep <- sum(shared.df$count)
      depth.locus <- locus.depth[[Flg.locus]]
      depth.Nextlocus <- locus.depth[[Flg.Nextlocus]]
      deapMean <-  mean(c( depth.locus,depth.Nextlocus) )
      
      if(shared.dep >= deapMean/2){
        Flag <- "ke"
        if(length(names(blockLst)) == 0){
          BlockFlag <- locus  }
      }else{
        Flag <- "No"
        BlockFlag <- Next.locus
      }
      blockLst[[paste0("P.",BlockFlag) ]] <- Next.locus
    }
    
    
    ##output 
    Block.Start <- c()
    Block.End <- c()
    for (block in names(blockLst)) {
      Block.Start <- c(Block.Start,block)
      Block.End <- c(Block.End,blockLst[[block]])
    }
    blocks <-  data.frame("Sample" = samp,
                          "Block.Start.1"= Block.Start,
                          "Block.End" = Block.End)
    blocks$Block.Start <- str_split_fixed(blocks$Block.Start.1, '[.]', 2)[,2]
    blocks$block.No <- rownames(blocks)
    blocks <- blocks %>% mutate(block = paste0("block.",block.No)) %>% mutate(Block.length =  as.numeric(Block.End)  - as.numeric(Block.Start) + 1 )
    blocks$Block.Start.1 <- NULL
    blocks$block.No <- NULL
    
    O.df <- blocks[,c("Sample","block","Block.Start","Block.End","Block.length")]
    write.table(O.df,paste0( "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/",samp,".blocks.txt"),
                sep = "\t",quote = F,row.names = F)
    
  }
}







