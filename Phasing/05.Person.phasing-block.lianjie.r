library(Biostrings)

person <- "P1"
samp <- "P1-E-d"
ref.samp <- "P1-E-j"



Hap.OP <- paste0("/home/try/Documents/projects/liuhj-Other/phasing-2024/blocks-phasing/",samp,"/") 
dir.create(Hap.OP)

##SNV locus
relaibLocus.P <- "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/"
SNV.table <- read.table(paste0(relaibLocus.P ,person,".reliable.locus.CCS.table.txt"),
                        sep = "\t",quote = "",header = T)
## blocks
block.table <- read.table(paste0( relaibLocus.P,samp,".blocks.txt"),
                          sep = "\t",quote = "",header = T)


## phasing
blocks.phasing.df <- list()
Haps.lst <- list()
peidui.list <- list()
#block.table <- block.table %>% filter(Block.Start != Block.End)
max.break <- 0
sampBreaks <- c()
#for (n.block in seq(1,2) ) {
for (n.block in seq(1,nrow(block.table)-1) ) {
  print( paste0("block : ",n.block) )
  B <- block.table[n.block,"block"]
  B.S <- block.table[n.block,"Block.Start"]
  B.E <- block.table[n.block,"Block.End"]
  print(paste0( B,".",B.S,"-",B.E))
  SNV.table.Cites <- SNV.table$X.Posi[SNV.table$X.Posi >= B.S & SNV.table$X.Posi <= B.E]
  SNV.table.Cites.df <- SNV.table.Cites %>% as.data.frame()
  colnames(SNV.table.Cites.df) <- "locus"
  Phased.HapsF <- paste0( Hap.OP,samp,".",B,".",B.S,"-",B.E,".fasta")
  readHaps.F <- paste0(Hap.OP ,samp,".",B,".",B.S,"-",B.E,".ReadHaps.fasta")
  if (file.exists(Phased.HapsF)) {
    dna_sequences <- readDNAStringSet(Phased.HapsF)
  } else {
    dna_sequences <- readDNAStringSet(readHaps.F)
  }
  seq_letters <- sapply(dna_sequences, as.character)
##next block
  B.next <- block.table[n.block+1,"block"]
  B.S.next <- block.table[n.block+1,"Block.Start"]
  B.E.next <- block.table[n.block+1,"Block.End"]
  Phased.HapsF.next <- paste0( Hap.OP,samp,".",B.next,".",B.S.next,"-",B.E.next,".fasta")
  readHaps.F.next <- paste0(Hap.OP ,samp,".",B.next,".",B.S.next,"-",B.E.next,".ReadHaps.fasta")
  if (file.exists(Phased.HapsF.next)) {
    dna_seq.next <- readDNAStringSet(Phased.HapsF.next)
  } else {
    dna_seq.next <- readDNAStringSet(readHaps.F.next)
  }
  seq.letters.next <- sapply(dna_seq.next, as.character)
  print(seq_letters)
  print(seq.letters.next)
  
  
###-This haps freqs
  This.Fs <- c()
  for (Hap.NM in names(seq_letters)) {
    if( grepl("[:]",Hap.NM) ){
      Freq <-  strsplit( strsplit(Hap.NM,"[:]")[[1]][2],"[_]" )[[1]][1] %>% as.numeric()
    } else{
      Freq <- 1
    }
    This.Fs <- c(This.Fs,Freq)
  }
  This.Fs <- This.Fs %>% sort(decreasing = T)
  print(This.Fs)
  
  This.NMs <- c()
  for (This.F in This.Fs) {
    for (Hap.NM in names(seq_letters)) {
      F1.idex <- F1.idex + 1
      if( grepl("[:]",Hap.NM) ){
        Freq <-  strsplit( strsplit(Hap.NM,"[:]")[[1]][2],"[_]" )[[1]][1] %>% as.numeric()
      } else{
        Freq <- 1
      }
      if (This.F == Freq & !Hap.NM %in% This.NMs ){
        This.NMs <- c(This.NMs,Hap.NM)
      }
    }
  }
  print(This.NMs)
  
  This.Names <- list()
  This.Freqs <- list()
  This.Haps <- list()
  This.sumFs <- c()
  F1.idex <- 0
  F1.s <- 0
  for (Hap.NM in This.NMs) {
    F1.idex <- F1.idex + 1
    if( grepl("[:]",Hap.NM) ){
      Freq <-  strsplit( strsplit(Hap.NM,"[:]")[[1]][2],"[_]" )[[1]][1] %>% as.numeric()
    } else{
      Freq <- 1
    }
    F1.s <- F1.s + Freq
    This.sumFs <- c(This.sumFs,F1.s)
    This.Names[[as.character(F1.idex)]] <- Hap.NM
    This.Freqs[[as.character(F1.idex)]] <- Freq
    This.Haps[[as.character(F1.idex)]] <- seq_letters[[Hap.NM]]
  }
  
  
  
  ###Next haps freqs
  Next.Fs <- c()
  for (Hap.NM in names(seq.letters.next)) {
    if( grepl("[:]",Hap.NM) ){
      Freq <-  strsplit( strsplit(Hap.NM,"[:]")[[1]][2],"[_]" )[[1]][1] %>% as.numeric()
    } else{
      Freq <- 1
    }
    Next.Fs <- c(Next.Fs,Freq)
  }
  Next.Fs <- Next.Fs %>% sort(decreasing = T)
  print(Next.Fs)
  
  
  Next.NMs <- c()
  for (Next.F in Next.Fs) {
    for (Hap.NM in names(seq.letters.next)) {
      F1.idex <- F1.idex + 1
      if( grepl("[:]",Hap.NM) ){
        Freq <-  strsplit( strsplit(Hap.NM,"[:]")[[1]][2],"[_]" )[[1]][1] %>% as.numeric()
      } else{
        Freq <- 1
      }
      if (Next.F == Freq & !Hap.NM %in% Next.NMs ){
        Next.NMs <- c(Next.NMs,Hap.NM)
      }
    }
  }
  print(Next.NMs)
  
  Next.Names <- list()
  Next.Freqs <- list()
  Next.Haps <- list()
  Next.sumFs <- c()
  F1.idex <- 0
  F1.s <- 0
  for (Hap.NM in Next.NMs) {
    F1.idex <- F1.idex + 1
    if( grepl("[:]",Hap.NM) ){
      Freq <-  strsplit( strsplit(Hap.NM,"[:]")[[1]][2],"[_]" )[[1]][1] %>% as.numeric()
    } else{
      Freq <- 1
    }
    F1.s <- F1.s + Freq
    Next.sumFs <- c(Next.sumFs,F1.s)
    Next.Names[[as.character(F1.idex)]] <- Hap.NM
    Next.Freqs[[as.character(F1.idex)]] <- Freq
    Next.Haps[[as.character(F1.idex)]] <- seq.letters.next[[Hap.NM]]
  }
  

  ### find all breaks
  UniqBreks <- c()
  allBreks <- c(This.sumFs,Next.sumFs) %>% unique() %>% sort()
  if(length(allBreks) != 1){
    for (iii in seq(1,length(allBreks)-1)) {
      nxi <- iii + 1
      print(allBreks[iii])
      print(allBreks[nxi])
      if( abs( allBreks[iii] - allBreks[nxi]) <= 0.03 ){
        UniqBreks <- c(UniqBreks,allBreks[nxi])
        UniqBreks <- UniqBreks[UniqBreks != allBreks[iii]]
        #iii <- iii + 2
      }
      if( abs( allBreks[iii] - allBreks[nxi]) > 0.03 ){
        UniqBreks <- c(UniqBreks,allBreks[iii])
        UniqBreks <- c(UniqBreks,allBreks[nxi])
      }
    }
    UniqBreks <- unique(UniqBreks)
    print(UniqBreks)
  } else{
    UniqBreks <- allBreks
  }
  
  print(B.S)
  max.break <- c(max.break,length(UniqBreks))
  sampBreaks <- c(sampBreaks,UniqBreks)
  ### duiying of all breaks
  all.BKs <-  c(0,UniqBreks)
  peidui <- list()
  for (BKid in seq(1,length(all.BKs)-1) ) {
    print("")
    print("")
    print(BKid)
    BK.S <- all.BKs[[BKid]]
    BK.E <- all.BKs[[BKid + 1]]
    print(BK.S)
    print(BK.E)
    print("--")
    This.BKs <-  c(0,This.sumFs)
    for (T.BKid in seq(1,length(This.BKs)-1) ) {
      T.BK <- This.BKs[[T.BKid]]
      T.N.BK <- This.BKs[[T.BKid + 1]]
      if (  BK.S >= T.BK-0.03 & BK.E <= T.N.BK ){
        print("This pipei")
        print(T.BK)
        print(T.N.BK)
        T.Hao.ID <- T.BKid
        break 
        }
    }
    Next.BKs <-  c(0,Next.sumFs)
    for (N.BKid in seq(1,length(Next.BKs)-1) ) {
      TN.BK <- Next.BKs[[N.BKid]]
      TN.N.BK <- Next.BKs[[N.BKid + 1]]
      if (BK.S >= TN.BK-0.05 & BK.E <= TN.N.BK ){
        print("Next pipei")
        print(TN.BK)
        print(TN.N.BK)
        TN.Hao.ID <- N.BKid
        break
      }
    } 
    print("duiying ID")
    print(T.Hao.ID)
    print(TN.Hao.ID)
    peidui[[paste0(BK.S,"-",BK.E)]] <- paste0(T.Hao.ID,":",TN.Hao.ID)
  }
  
}


  
  
#     
#   
#   for (This.Freq in This.sumFs) {
#     print("--")
#     print(This.Freq)
#     yijing.Next <- c()
#     for (Next.Freq in Next.sumFs) {
#       print(Next.Freq)
#       yijing.Next <- c(yijing.Next,Next.Freq)
#       if(abs(Next.Freq-This.Freq) > 0.05 & !Next.Freq %in% allBreks & Next.Freq %in% yijing.Next ) {
#         allBreks <- c(allBreks,Next.Freq)
#         # print(Next.Freq)
#         # print(This.Freq)
#         print("hh")
#         break
#   }}}
#   
#   
#   
#   
#   This.Fsum <- 0
#   for (This.i in names(This.Freqs)) {
#     Next.Fsum <- 0
#     for (Next.i in names(Next.Freqs)) {
#       
#     }
#   }
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   Next.Names <- list()
#   Next.Freqs <- list()
#   Next.Haps <- list()
#   Next.sumFs <- c()
#   F1.idex <- 0
#   for (Hap.NM in names(seq.letters.next)) {
#     F1.idex <- F1.idex + 1
#     Freq <-  strsplit( strsplit(Hap.NM,"[:]")[[1]][2],"[_]" )[[1]][1] %>% as.numeric()
#     Next.sumFs <- c(Next.sumFs,F1.s)
#     Next.Names[[as.character(F1.idex)]] <- Hap.NM
#     Next.Freqs[[as.character(F1.idex)]] <- Freq
#     Next.Haps[[as.character(F1.idex)]] <- seq.letters.next[[Hap.NM]]
#   }
#   
#   
#   
#   
#   
#   
#   ThisFreqs.brk  <- list()
#   ll <- 0
#   for ( ii in seq(1,length(This.Freqs))) {
#     ThisFreqs.brk[[ as.character( ll )  ]] <- This.Haps[ ii]
#     ll <- ll + This.Freqs[ii]
#   }
# 
#   
#   ## next
#   Next.Freqs <- c()
#   Next.Haps <- list()
#   #Next.Freq.Flag <- list()
#   freqSm <- 0
#   for (Hap.NM in names(seq.letters.next)) {
#     Freq <-  strsplit( strsplit(Hap.NM,"[:]")[[1]][2],"[_]" )[[1]][1] %>% as.numeric()
#     Next.Freqs <- c(Next.Freqs,Freq)
#     Next.Haps[[Hap.NM]] <- seq.letters.next[[Hap.NM]]
#   }
#   # NextFreqs.brk  <- list()
#   # ll <- 0
#   # for ( ii in seq(1,length(Next.Freqs))) {
#   #   NextFreqs.brk[[ as.character( ll )  ]] <- Next.Haps[ ii]
#   #   ll <- ll + Next.Freqs[ii]
#   # }
#   # 
#   
#   ### this freq breaks
#   Freq.brk <- c() 
#   Freq.sum <- 0
#   for (Freq in sort(This.Freqs,decreasing = T)   ) {
#     Freq.sum <- Freq.sum + Freq
#     Freq.brk <- c(Freq.brk,Freq.sum) 
#   }
#   All.break <- Freq.brk[seq(1,length(Freq.brk)-1)]
#   ##next freq breaks
#   Freq.nxt.brk <- c() 
#   Freq.sum <- 0
#   for (Freq.nxt in sort(Next.Freqs,decreasing = T) ) {
#     Freq.sum <- Freq.sum + Freq.nxt
#     print(Freq.sum)
#     if(abs( Freq.sum - max(All.break) )  > 0.05 ){
#       print(">")
#       All.break <- c(All.break,Freq.sum)
#     }
#   }
#   print( paste0("break.Num : ", length(All.break)-1))
# 
#   
#   
#   
#   for (variable in vector) {
#     
#   }
#   
#   
#   
#     
#   ### This and next peidui 
#   peidui <- list()
#   C <- 1
#   comp.L <- 0
#   qian.Freq <- 0
#   for (comp.R in breaks) {
#     print("oooo")
#     print(comp.R)
#     pipei.Freq <- round(comp.R - qian.Freq,3) 
#     print("iii")
#     print(pipei.Freq)
#     ##this
#     peidui[[ paste0(as.character(C),"_",pipei.Freq )  ]] <- list()
#     print( paste0(as.character(C),"_",pipei.Freq ) )
#     
#     this.L <- 0
#     for (HapPert in HapreadPerct) {
#       this.R <- this.L + HapPert
#       if( this.L <= comp.L & this.R >= comp.R){
#         thisHap <- Hap.Num.lst[[HapreadNums[which(HapreadPerct == HapPert)]  %>% as.character()]]
#         peidui[[paste0(as.character(C),"_",pipei.Freq) ]][[thisHap]] <- '' 
#         break
#       }
#       this.L <- this.R
#     }
#     ##next
#     next.L <- 0
#     for (HapPert.nxt in next.HapreadPerct) {
#       if( next.L <= comp.L & next.R >= comp.R){
#         NxtHap <- Hap.nxt.Num.lst[[next.HapreadNums[which(next.HapreadPerct == HapPert.nxt)]  %>% as.character()]]
#         peidui[[ paste0(as.character(C),"_",pipei.Freq)  ]][[ thisHap ]] <-  NxtHap
#         break
#       }
#       next.L <- next.R
#     }
#     qian.Freq <- comp.R
#     comp.L <-  comp.R
#     C <- C + 1
#   }
#   peidui.list[[B]] <- peidui
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###lian jie 
# Hap.maxNum <- max(max.break)
# 
# Final.long.Haps <- list()
# First.Haps <- peidui.list[[ names(peidui.list)[1] ]]
# print(names(First.Haps))
# for (First.Hap  in names(First.Haps) ) {
#   print( paste0("1st: ",First.Hap) )
#   First.Hap.ID <- names(First.Haps[[First.Hap]] )[[1]]
#   print( paste0("1st HapID: ",First.Hap.ID) )
#   Match.Hap.ID <- First.Haps[[First.Hap]][[First.Hap.ID]]
#   Final.long.Haps[[paste(First.Hap.ID,Match.Hap.ID,sep = "|") ]] <- c("a")
#   
#   print(Match.Hap.ID)
#   print(paste(First.Hap.ID,Match.Hap.ID,sep = "|"))
#   print("match")
#   #for ( Other.Bs   in names(peidui.list)[2: length(names(peidui.list)) ]   ) {
#   for ( Other.Bs in names(peidui.list)[2: 4]   ) {
#     currrtB.Haps <- peidui.list[[Other.Bs]]
#     print(paste0("Block:",Other.Bs)  )
#     sameCount <- 0
#     for (currrtB.Hap  in names(currrtB.Haps) ) {
#       print(currrtB.Hap)
#       currrtB.Hap.ID <- names(currrtB.Haps[[currrtB.Hap]] ) [[1]]
#       matchNext.idx <- currrtB.Haps[[currrtB.Hap]][[currrtB.Hap.ID]]
#       print( paste0("current hapID:",currrtB.Hap.ID) )
#       print( paste0("pipei next hapID:",matchNext.idx) )
#       
#       for (Final.longH in names(Final.long.Haps)) {
#         print("--")
#         print(Final.longH)
#         Final.longH.last <- strsplit(Final.longH, "[|]")[[1]][ length(strsplit(Final.longH,"[|]")[[1]]) ]
#         print(Final.longH.last)
#         print(currrtB.Hap.ID)
#         if(currrtB.Hap.ID == Final.longH.last ){
#           print("same")
#           sameCount <- sameCount + 1
#           new.ID <- paste(Final.longH,matchNext.idx,sep = "|")
#           if(sameCount == 1){ Final.long.Haps[[Final.longH]]  <- NULL }
#           Final.long.Haps[[new.ID]] <- c(Final.long.Haps[[new.ID]],c("b"))
#         }
#       }
#         
#     }
#   }
# }
# 








