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


P <- c()
S <- c()
Bc <- c()
N <- c()
for (person in names(Pers.refs)) {
  #for (person in c("P23")) {
  print(person)
  
  for (samp in samp_list[[person]]) {
    
    Hap.OP <- paste0("/home/try/Documents/projects/liuhj-Other/phasing-2024/blocks-phasing/",samp,"/") 
    relaibLocus.P <- "/home/try/Documents/projects/liuhj-Other/phasing-2024/reliable_locus/"
    block.table <- read.table(paste0( relaibLocus.P,samp,".blocks.txt"),
                              sep = "\t",quote = "",header = T)
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
      if (file.exists(PhasedHaps.F1)) {
        PhasedHaps.F <- PhasedHaps.F1
      } 
      else{
        PhasedHaps.F <- PhasedHaps.F2
      }
      fasta_data <- read_fasta(PhasedHaps.F)
      N.Haps <- length(names(fasta_data))
      
      P <- c(P,person)
      S <- c(S,samp)
      Bc <- c(Bc,B) 
      N <- c(N,N.Haps)
    }
  }
}


Stat.df <- data.frame("Person" = P,
                      "Sample" = S,
                      "Block" = Bc,
                      "Hap.n" = N)
Stat.df.plt <- Stat.df %>% group_by(Sample,Hap.n) %>% summarise(count = n())

p.bar  <- ggplot(Stat.df.plt, aes(x = Hap.n,y = count)) +
  geom_bar(stat = "identity",position = "dodge") +
  labs(title = "Histogram of haplotype number",
       x = "Hap number",
       y = "Frequency") +
  theme_minimal() + 
  facet_wrap( ~ Sample)


ggsave(p.bar,filename = paste0(paste0("/home/try/Documents/projects/liuhj-Other/phasing-2024/blocks-phasing/phased.Haps.Num.jpg") ),width = 10,height = 10)






# 
# 
# 
# library(ggplot2)
# 
# 
# main_line <- data.frame(
#   x = c(1, 2, 3, 4, 5),
#   y = c(1,1,1,1,1)
# )
# 
# 
# branch_lines <- data.frame(
#   x_start = c(1, 1),       # 分叉起点的 x 坐标
#   y_start = c(1, 1),       # 分叉起点的 y 坐标
#   x_end = c(2, 2),         # 分叉终点的 x 坐标
#   y_end = c(2, 2),         # 分叉终点的 y 坐标
#   weight = c(2, 5)         # 分叉线的权重
# )
# 
# 
# 
# 
# pp <- ggplot() +
# geom_path(data = main_line, aes(x = x, y = y),
#           color = "black", size = 1.5) +
# geom_segment(data = branch_lines,
#              aes(x = x_start, y = y_start,
#                  xend = x_end, yend = y_end,
#                  color = weight, size = weight)) +
# scale_color_gradient(low = "blue", high = "red") +  # 颜色根据权重变化
#   scale_size_continuous(range = c(1, 5)) +           # 粗细根据权重变化
# labs(title = "Line with Weighted Branches",
#      x = "X Axis",
#      y = "Y Axis",
#      color = "Weight",
#      size = "Weight") +
# theme_minimal()

 

