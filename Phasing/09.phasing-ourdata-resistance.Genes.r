##resistans 

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

resist.Genes <- list()
resist.Genes[["Cla"]] <- c("23S ribosomal RNA")
resist.Genes[["Lef"]] <- c("DNA gyrase subunit A")

Ref.Anno.P <- "/disk/lhj_4T/lhj/HP/process/assembly/all_gff/"
relaibLocus.P <- "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/"
readHap.P <- "/home/try/Documents/projects/liuhj-Other/phasing-2024/haplotype-read/"
O.path <- "/home/try/Documents/projects/liuhj-Other/phasing-2024/resistance.Gene.phasedHaps/"
dir.create(O.path)

write_fasta <- function(sequences, file_name) {
  # 使用 writeLines 将每条序列按行写入文件
  lines <- unlist(lapply(sequences, function(seq) {
    c(seq$header, seq$sequence)
  }))
  # 写入 FASTA 文件
  writeLines(lines, con = file_name)
}

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


panduan.phase  <-  function(person,samp,SNV.table,drug,durg.Gene,O.path){
  print(person)
  ref.samp <- Pers.refs[[person]]
  #for (samp in samp_list[[person]]) {
    print(samp)
    readHap.F <- paste0(readHap.P,"readHap.",samp,"_2_",ref.samp,"/2.read-haplotype.Stat.txt") 
    readHaps <- read.table(readHap.F,sep = "\t",quote = "",header = F)
    colnames(readHaps) <- c("count","S.E","haplotype")
    readHaps$id <- rownames(readHaps)
    readHaps <- readHaps[ !is.na(readHaps$haplotype),]
    
    ##SNV locus
    # SNV.table <- read.table(paste0( "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/",person,".reliable.locus.CCS.table.txt"),
    #                         sep = "\t",quote = "",header = T)
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
    O.df <- na.omit(O.df)
    write.table(O.df,paste0( O.path,samp,"_",drug,"_",durg.Gene,".blocks.txt"),
                sep = "\t",quote = F,row.names = F)
  #}
}


phase <- function(person,samp,SNV.table,drug,durg.Gene,O.path){
    print(person)
    ref.samp <- Pers.refs[[person]]
    #for (samp in "P1-E-d") {
    #for (samp in samp_list[[person]]) {
      print(samp)
      Hap.OP.1 <- paste0("/home/try/Documents/projects/liuhj-Other/phasing-2024/blocks-phasing-resistance/") 
      dir.create(Hap.OP.1)
      Hap.OP <- paste0(Hap.OP.1,samp,"/") 
      dir.create(Hap.OP)
      readHap.P <- "/home/try/Documents/projects/liuhj-Other/phasing-2024/haplotype-read/" 
      relaibLocus.P <- "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/"
      readHap.F <- paste0(readHap.P,"readHap.",samp,"_2_",ref.samp,"/2.read-haplotype.Stat.txt") 
      readHaps <- read.table(readHap.F,sep = "\t",quote = "",header = F)
      colnames(readHaps) <- c("count","S.E","haplotype")
      readHaps$id <- rownames(readHaps)
      readHaps <- readHaps[ !is.na(readHaps$haplotype),]
      
      readHaps$id <- rownames(readHaps)
      readHaps$S <- str_split_fixed(readHaps$S.E, '-', 3)[,2]
      readHaps$E <- str_split_fixed(readHaps$S.E, '-', 3)[,3]
      readHaps$count <- as.numeric(readHaps$count)
      readHaps$S <- as.numeric(readHaps$S)
      readHaps$E <- as.numeric(readHaps$E)
      
      # ##SNV locus
      # SNV.table <- read.table(paste0(relaibLocus.P ,person,".reliable.locus.CCS.table.txt"),
      #                         sep = "\t",quote = "",header = T)
  
      drug <- Res.durg
      durg.Gene <- Res.Gene
      ## blocks
      block.table <- read.table(paste0(O.path, samp,"_",drug,"_",durg.Gene,".blocks.txt"),
                                sep = "\t",quote = "",header = T)
      
      ## phasing
      blocks.phasing.df <- list()
      #block.table <- block.table %>% filter(Block.Start != Block.End)
      for (n.block in seq(1,nrow(block.table)) ) {
        #for (n.block in seq(44,44) ) {
        print( paste0("block : ",n.block) )
        B <- block.table[n.block,"block"]
        B.S <- block.table[n.block,"Block.Start"]
        B.E <- block.table[n.block,"Block.End"]
        SNV.table.Cites <- SNV.table$X.Posi[SNV.table$X.Posi >= B.S & SNV.table$X.Posi <= B.E]
        SNV.table.Cites.df <- SNV.table.Cites %>% as.data.frame()
        colnames(SNV.table.Cites.df) <- "locus"
        Fit.df <-  readHaps %>% filter(S >= B.S & B.E >= S | S <= B.S & B.E <= E | B.S <= E & B.E >= E,
                                       count >= 3)
        Fit.df <- Fit.df[order(Fit.df$count,decreasing = T),]
        
        if(nrow(Fit.df) >= 1){
          ##filter haplotype
          n <- 0
          for (hap.n in seq(1,nrow(Fit.df))) {
            n <- n + 1
            hap <- Fit.df[hap.n,"haplotype"]
            Ct <- Fit.df[hap.n,"count"]
            hap.Nm <- paste0("Hap_",n,".",Ct) 
            bases <-  str_split(hap,"_")[[1]]
            pos.c <- c()
            base.c <- c()
            for (Pbase in bases) {
              base <- gsub('[0-9]','',Pbase[1])
              Pos <- gsub('[A-Z]','',Pbase[1])
              pos.c <- c(pos.c,Pos)
              base.c <- c(base.c,base)
            }
            hap.df <- data.frame("locus" = pos.c,"base" = base.c)
            colnames(hap.df)[2] <- paste0(hap.Nm)
            SNV.table.Cites.df <- merge(SNV.table.Cites.df,hap.df,all.x = T)
          }
          rownames(SNV.table.Cites.df) <- SNV.table.Cites.df$locus
          SNV.table.Cites.df$locus <- NULL
          
          
          ## core phasing script
          ### 1. find ref seq
          dfFormg <- SNV.table.Cites.df
          dfFormg.t <- t(dfFormg) %>% as.data.frame()
          dfFormg.t[is.na(dfFormg.t)] <- "Z"
          dfFormg.t$combined_column <- apply(dfFormg.t, 1, function(row) paste(row, collapse = ""))
          dfFormg.t <- dfFormg.t[order(dfFormg.t$combined_column),]
          
          
          ##output sort dataframe
          Out.list <- list()
          O.index <- 0
          for (blockHaps in rownames(dfFormg.t) ) {
            O.index <- O.index + 1
            name <- paste0(">",blockHaps)
            seq <- paste(dfFormg.t[blockHaps,"combined_column"], collapse = "")
            Out.list[[O.index]] <- list(header = name, sequence = seq)
          }
          write_fasta(Out.list, paste0(Hap.OP ,samp,"_",drug,"_",durg.Gene,".",B,".",B.S,"-",B.E,".ReadHaps.fasta") )
          
          ## 
          dfFormg.t.2 <- dfFormg.t
          dfFormg.t.2$combined_column <- NULL
          dfFormg.t.2 <- t(dfFormg.t.2) %>% as.data.frame()
          dfFormg.t.2$locus <- rownames(dfFormg.t.2)
          
          REFs <- dfFormg.t.2[,c("locus",colnames(dfFormg.t.2)[1])] 
          n.otherHaps <- length(colnames(dfFormg.t.2))-1
          
          if(n.otherHaps >1){
            #print(n.otherHaps)
            for (comCol in colnames(dfFormg.t.2)[2:n.otherHaps] ) {
              #for (comCol in colnames(dfFormg.t.2)[2:3] ) {
              # print("comp")
              # print(comCol)
              match <- c()
              for (refCol in colnames(REFs)[2:length(colnames(REFs))] ) {
                # print("refCol")
                # print(refCol)
                DF <- merge(REFs[,c("locus",refCol)], dfFormg.t.2[,c("locus",comCol)])
                #print(colnames(DF))
                C1 <- colnames(DF)[2]
                C2 <- colnames(DF)[3]
                DF <- DF %>% mutate(comp = case_when(.data[[C1]] == "Z" | .data[[C2]] == "Z" ~ "Z",
                                                     .data[[C1]] == .data[[C2]] ~ "Same",
                                                     TRUE ~ "Not-same") )
                Diff.Ct <- nrow(DF[DF$comp == "Not-same",])
                C1.readN <- strsplit(C1,"[.]")[[1]][2] %>% as.numeric()
                C2.readN <- strsplit(C2,"[.]")[[1]][2] %>% as.numeric()
                if( comCol == "Hap_377.3"){
                  print("match")
                  print(refCol)
                  print("comp")
                  print(comCol)
                  print(Diff.Ct)
                }
                if( nrow( DF[DF$comp == "N",] ) != nrow(DF) &  Diff.Ct == 0 ){
                  match <- c(match,refCol)
                }
              }
              if( comCol == "Hap_377.3"){
                print("match")
                print(match)
              }
              
              # 
              # print(match)
              # print(length(match) )
              if(length(match) == 0){ REFs <- merge(REFs,dfFormg.t.2[,c("locus",comCol)],all.x = T)  }
              if( length(match) == 1 ){
                DF <- merge(REFs[,c("locus",match[1])], dfFormg.t.2[,c("locus",comCol)])
                C1 <- colnames(DF)[2]
                C2 <- colnames(DF)[3]
                DF <- DF %>% mutate(comp = case_when(.data[[C1]] =="Z" | .data[[C2]] =="Z" ~ "N",
                                                     .data[[C1]] == .data[[C2]] ~ "Same",
                                                     TRUE ~ "Not-same") )
                DF <- DF[order(DF$locus), ]
                Diff.Ct <- nrow(DF[DF$comp == "Not-same",])
                C1.readN <- strsplit(C1,"[.]")[[1]][2] %>% as.numeric()
                C2.readN <- strsplit(C2,"[.]")[[1]][2] %>% as.numeric()
                
                if(C1.readN >= C2.readN){maxC = C1} else{maxC = C2}
                DF <- DF %>% mutate(mergeBase = case_when(.data[[C1]] =="Z" & .data[[C2]] == "Z" ~ "Z",
                                                          !.data[[C1]] =="Z" & .data[[C2]] == "Z" ~ .data[[C1]],
                                                          .data[[C1]] =="Z" & !.data[[C2]] == "Z" ~ .data[[C2]],
                                                          comp == "Same" ~ .data[[C1]],
                                                          comp == "Not-same" ~ .data[[maxC]],
                                                          TRUE ~ "") )
                mg.DF <- DF[,c("locus","mergeBase")]
                colnames(mg.DF)[2] <- paste0(strsplit(C1,"[.]")[[1]][1],".",C1.readN + C2.readN)
                REFs <- merge(REFs,mg.DF,all.x = T)
                REFs <- REFs %>% dplyr::select(-c(C1)) }
              
              if(length(match) > 1){
                sum <- 0
                for (match.Ref in match) {
                  readN <- strsplit(match.Ref,"[.]")[[1]][2] %>% as.numeric()
                  sum <- sum + readN  }
                for (match.Ref in match) {
                  #print("she")
                  #print(head(REFs[,c("locus",match.Ref)]))
                  DF <- merge(REFs[,c("locus",match.Ref)], dfFormg.t.2[,c("locus",comCol)])
                  C1 <- colnames(DF)[2]
                  C2 <- colnames(DF)[3]
                  C1.readN <- strsplit(C1,"[.]")[[1]][2] %>% as.numeric()
                  C2.readN <- strsplit(C2,"[.]")[[1]][2] %>% as.numeric()
                  
                  DF <- DF %>% mutate(mergeBase = case_when(.data[[C1]] =="Z" & .data[[C2]] == "Z" ~ "Z",
                                                            !.data[[C1]] =="Z" & .data[[C2]] == "Z" ~ .data[[C1]],
                                                            .data[[C1]] =="Z" & !.data[[C2]] == "Z" ~ .data[[C2]],
                                                            .data[[C1]] == .data[[C2]] ~ .data[[C1]],
                                                            TRUE ~ "") )
                  mg.DF <- DF[,c("locus","mergeBase")]
                  colnames(mg.DF)[2] <- paste0(strsplit(C1,"[.]")[[1]][1],".",C1.readN + ceiling(C2.readN/sum) )
                  #print("222222")
                  #print(head(mg.DF))
                  REFs <- merge(REFs,mg.DF,all.x = T)
                  #print(head(REFs))
                  REFs <- REFs %>% dplyr::select(-c(C1)) }
              }
              #print(colnames(REFs) )  
            }
            #REFs <- REFs[, !apply(REFs, 2, function(col) any(grepl("Z", col)))]
          }
          
          
          ##filter Phased Haps
          all.Haps <-  colnames(REFs)[2:length(colnames(REFs))]
          Len.c <- c()
          for (Hap in all.Haps) {
            freq <- strsplit(Hap,"[.]")[[1]][length( strsplit(Hap,"[.]")[[1]] )] %>% as.numeric()
            Len.c <- c(Len.c,freq)
          }
          Len.sumNum <- sum(Len.c)
          Len.cPerct <-  Len.c/Len.sumNum
          l.idx <-  which(Len.cPerct<0.02)
          OK.index <-  setdiff(seq(1,length(Len.cPerct)),l.idx)
          Ok.Haps <- all.Haps[OK.index]
          OK.lensSum <- Len.c[OK.index] %>% sum()
          
          
          new.Hap.ID <- c()
          III <- 0 
          for (Hap in Ok.Haps) {
            III <- III + 1
            Num <- strsplit(Hap,"[.]")[[1]][length( strsplit(Hap,"[.]")[[1]] )] %>% as.numeric()
            freqt <- round(Num/OK.lensSum,2)
            newID <- paste0("Phased.Hap_",III,"_Freq:",freqt,"_",Num)
            new.Hap.ID <- c(new.Hap.ID,newID)
          }
          
          
          ##out phasing Haps
          REFs.O <- REFs[,c("locus",Ok.Haps)]
          colnames(REFs.O) <- c("locus",new.Hap.ID) 
          rownames(REFs.O) <- REFs.O$locus
          REFs.O$locus <- NULL
          REFs.O <- REFs.O[as.character(SNV.table.Cites) ,]
          REFs.t <- t(REFs.O) %>% as.data.frame()
          REFs.t$combined_column <- apply(REFs.t, 1, function(row) paste(row, collapse = ""))
          
          if( length(Ok.Haps) == 1){
            rownames(REFs.t) <- new.Hap.ID
          }
          
          Out.Haps <- list()
          O.index <- 0
          for (blockHaps in rownames(REFs.t) ) {
            O.index <- O.index + 1
            name <- paste0(">",samp,".",blockHaps)
            seq <- paste(REFs.t[blockHaps,"combined_column"], collapse = "")
            Out.Haps[[O.index]] <- list(header = name, sequence = seq)
          }
          write_fasta(Out.Haps, paste0( Hap.OP,samp,"_",drug,"_",durg.Gene,".",B,".",B.S,"-",B.E,".fasta") )
        }
      }
    }




#######------------------------------------------------------------------------
Samp.cc <- c()
resGene <- c()
resGene.iSNV <- c()

for (person in names(Pers.refs)) {
#for (person in c("P24")) {
  print(person)
  ref.samp <- Pers.refs[[person]]
  sort.Samp <- c(ref.samp,setdiff(samp_list[[person]],ref.samp) )
   
  ## recombi region, SNV percent 
  SNV.table <- read.table(paste0(relaibLocus.P ,person,".reliable.locus.CCS.table.txt"),
                          sep = "\t",quote = "",header = T)
  
  RefSamp.annoF <- paste0(Ref.Anno.P,ref.samp,".gff/",ref.samp,".gff")

  RefSamp.anno <- read.delim(RefSamp.annoF, header = FALSE, comment.char = "#", sep = "\t", 
                         col.names = c("seqid", "source", "type", "start", "end", 
                                       "score", "strand", "phase", "attributes"))
  
  ref.gnmP <- "/disk/lhj_4T/lhj/HP/process/assembly/all_genomes/"
  ref.gnmF <- paste0(ref.gnmP,ref.samp,".fasta")
  ref.gnm <-  readDNAStringSet(ref.gnmF)[[1]]  #  read_fasta(ref.gnmF)[1]
  
  
  ##ref resistance genes phasing
  for (Res.durg in names(resist.Genes)) {
    print(Res.durg)
    Res.Genes <- resist.Genes[[Res.durg]]
    drug.O.df <- data.frame()
    for (Res.Gene in Res.Genes) {
      print(Res.Gene)
      Res.Gene.l <-  RefSamp.anno %>% filter( grepl( Res.Gene, attributes) )
      drug.O.df <- rbind(drug.O.df,Res.Gene.l)
      for (cpy in seq(1,nrow(Res.Gene.l))) {
        Res.S <- Res.Gene.l$start[cpy]
        Res.E <- Res.Gene.l$end[cpy]
        
        SNV.table.Res <- SNV.table %>% filter(X.Posi >= Res.S,
                                              X.Posi <= Res.E)
        
        Samp.cc <- c(Samp.cc,person)
        resGene <- c(resGene,Res.durg)
        resGene.iSNV <- c(resGene.iSNV,nrow(SNV.table.Res))
        
        Bases <- list()
        for (SNVs.site in seq(Res.S,Res.E)) {
          ref.base <-  subseq(ref.gnm , start = as.numeric(SNVs.site),end = as.numeric(SNVs.site))
          Bases[[ paste0("C_",SNVs.site) ]] <- ref.base
        }
        
        if(nrow(SNV.table.Res) >= 1){
          Res.durg <- paste0(Res.durg,"-",cpy)
          print("hhhh")
          Res.blocks <- c()
          print(dim(SNV.table.Res))
          openxlsx::write.xlsx(SNV.table.Res,paste0(O.path,person,".",ref.samp,"_",Res.durg,"_",Res.Gene,".iSNV.table.xlsx"))
          
          
          for (samp in sort.Samp) {
            panduan.phase(person,samp,SNV.table.Res,Res.durg,Res.Gene,O.path)
            phase(person,samp,SNV.table.Res,Res.durg,Res.Gene,O.path)
            
            
            block.table <- read.table(paste0(O.path, samp,"_",Res.durg,"_",Res.Gene,".blocks.txt"),
                                      sep = "\t",quote = "",header = T)
            ## phasing
            Out.list <- list()
            for (n.block in seq(1,nrow(block.table)) ) {
              print( paste0("block : ",n.block) )
              B <- block.table[n.block,"block"]
              B.S <- block.table[n.block,"Block.Start"]
              B.E <- block.table[n.block,"Block.End"]
              
              Hap.OP.1 <- paste0("/home/try/Documents/projects/liuhj-Other/phasing-2024/blocks-phasing-resistance/") 
              Hap.OP <- paste0(Hap.OP.1,samp,"/") 
              Phased.HapsF <-  paste0(Hap.OP, samp,"_",Res.durg,"_",Res.Gene,".",B,".",B.S,"-",B.E,".fasta")
              Phased.Haps <- read_fasta(Phased.HapsF)
              idex <- 0
              for (Phased.Hap in names(Phased.Haps)) {
                idex <- idex + 1
                seq.Hap <- Phased.Haps[[Phased.Hap]]
                seq.Hap.Cites <- SNV.table.Res$X.Posi
                seq.Hapseq <- c()
                for (Cite in seq(Res.S,Res.E)) {
                  if(Cite %in% seq.Hap.Cites){
                   BBBs <-  strsplit(seq.Hap[which(seq.Hap.Cites == Cite)], split = "")[[1]]  
                  }else{
                  BBBs <- Bases[[paste0("C_",Cite)]] %>% as.character()
                  }
                  seq.Hapseq <- c(seq.Hapseq,BBBs)
                }
                seq.Hapseq.O <- paste0(seq.Hapseq,collapse = "")
                Out.list[[idex]] <- list(header = paste0(">",Phased.Hap), sequence = seq.Hapseq.O)
              }
            }
            write_fasta(Out.list, paste0(Hap.OP ,samp,"_",Res.durg,"_",Res.Gene,".",B,".",B.S,"-",B.E,".AllSeqOfGene.Haps.fasta") )
 
            
           } 
         }
      }
    }
    openxlsx::write.xlsx(drug.O.df,paste0(O.path,person,".",ref.samp,"_",Res.durg,".xlsx"))
  }
}



O.resiSNV <- data.frame("Person" = Samp.cc,
                        "Resistance.Gene" = resGene,
                        "iSNV.locus" = resGene.iSNV)


openxlsx::write.xlsx(O.resiSNV,paste0(O.path,"all.Person.res.Genes.iSNV.xlsx"))


## Lev
O.resiSNV.Lev <- O.resiSNV %>% filter(Resistance.Gene == "Lef",
                                      iSNV.locus != 0)
O.resiSNV.Lev <- O.resiSNV.Lev[order(O.resiSNV.Lev$iSNV.locus,decreasing = T),]

O.resiSNV.Lev$Person <- factor(O.resiSNV.Lev$Person,levels = O.resiSNV.Lev$Person)
p.bar  <- ggplot(O.resiSNV.Lev, aes(x = Person,y = iSNV.locus  )) +
  geom_bar(stat = "identity") +
  labs( #title = "Histogram of haplotype number",
       #x = "Person",
       y = "No. of iSNV locus in gyrA") +
  geom_text(aes(label = iSNV.locus), vjust = -0.5, color = "black", size = 3) +  
  theme_minimal() +
  ylim(0,45) +
  theme(
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),   # 移除次网格线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75)  # 添加黑色边框
  )


ggsave(p.bar,filename = paste0(O.path,"RestGene.Lev.iSNV.Num.barplot.pdf"),width = 4,height = 2.5)


## 23S rRNA
O.resiSNV.cla <- O.resiSNV %>% filter(Resistance.Gene == "Cla",
                                      iSNV.locus != 0)
O.resiSNV.cla <- O.resiSNV.cla[order(O.resiSNV.cla$iSNV.locus,decreasing = T),]
O.resiSNV.cla$Person <- factor(O.resiSNV.cla$Person,levels = O.resiSNV.cla$Person)
p.bar  <- ggplot(O.resiSNV.cla, aes(x = Person,y = iSNV.locus  )) +
  geom_bar(stat = "identity") +
  labs( #title = "Histogram of haplotype number",
    #x = "Person",
    y = "No. of iSNV locus in 23S rRNA") +
  geom_text(aes(label = iSNV.locus), vjust = -0.5, color = "black", size = 3) +  
  theme_minimal() +
  ylim(0,45) +
  theme(
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),   # 移除次网格线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75)  # 添加黑色边框
  )
ggsave(p.bar,filename = paste0(O.path,"RestGene.Cla.iSNV.Num.barplot.pdf"),width = 4,height = 2.5)








