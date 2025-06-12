library("dplyr")
library("stringr")

##trans to hapSeq
# 将序列写入文件
write_fasta <- function(sequences, file_name) {
  # 使用 writeLines 将每条序列按行写入文件
  lines <- unlist(lapply(sequences, function(seq) {
    c(seq$header, seq$sequence)
  }))
  # 写入 FASTA 文件
  writeLines(lines, con = file_name)
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
Pers.refs[["P22"]] <- "P22-C_sL-d"
Pers.refs[["P23"]] <- "P23-C_sL-j"
Pers.refs[["P24"]] <- "P24-C_sC-t"
Pers.refs[["P25"]] <- "P25-C_sC-t"

samp_list <- list()
samp_list[["P1"]] <- c("P1-E-d","P1-E-x","P1-E-j","P1-E-t")
#samp_list[["P2"]] <- c("P2-E-t")
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


for (person in names(Pers.refs)) {
#for (person in c("P24")) {
  print(person)
  ref.samp <- Pers.refs[[person]]
  #for (samp in "P1-E-d") {
  for (samp in samp_list[[person]]) {
    print(samp)
    Hap.OP <- paste0("/home/try/Documents/projects/liuhj-Other/phasing-2024/blocks-phasing/",samp,"/") 
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
    
    ##SNV locus
    SNV.table <- read.table(paste0(relaibLocus.P ,person,".reliable.locus.CCS.table.txt"),
                            sep = "\t",quote = "",header = T)
    
    ## blocks
    block.table <- read.table(paste0( relaibLocus.P,samp,".blocks.txt"),
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
      write_fasta(Out.list, paste0(Hap.OP ,samp,".",B,".",B.S,"-",B.E,".ReadHaps.fasta") )
      
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
        write_fasta(Out.Haps, paste0( Hap.OP,samp,".",B,".",B.S,"-",B.E,".fasta") )
        
      
    }


}}






