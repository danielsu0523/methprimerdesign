Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/home/user/Downloads/ncbi-blast-2.15.0+/bin", sep = ":"))
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/home/user/Downloads/primer3/src", sep = ":"))


BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
BiocManager::install("TAPseq", force = TRUE)
BiocManager::install("Biostrings", force = TRUE)

devtools::install_github("jensenlab/primer3")
devtools::install_github("ruthkr/rnafolding")

#library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(primer3)
library(TAPseq)
# An R-package to design PCR primers for TAP-seq 
library(GenomicRanges)
library(BiocParallel)
library(rnafolding)
library(shiny)
library(BSgenome)


# chromosome 11 truncated transcript sequences and annotations
data("chr11_truncated_txs_seq")

# create TsIOList object for the first two sequence templates
tapseq_io <- TAPseqInput(chr11_truncated_txs_seq[1:2], product_size_range = c(350, 500))

# design primers
## Not run: 
#tapseq_io <- designPrimers(tapseq_io, primer3_core = getOption("TAPseq.primer3_core"))
tapseq_io <- designPrimers(tapseq_io, primer3_core = "/home/user/Downloads/primer3/src/primer3_core")

## End(Not run)

# designed primers are stored in the tapseq_primers slot
tapseq_primers(tapseq_io)


# IRanges object with 10 ranges and 8 metadata columns:
#   start       end     width |   penalty             sequence        tm gc_percent self_any_th
# <integer> <integer> <integer> | <numeric>          <character> <numeric>  <numeric>   <numeric>
#   AKIP1.primer_left_0       626       645        20 |  0.211770 ATCCCTTGGTCCTGGAGGCA    62.788         60        0.00
# AKIP1.primer_left_1       627       646        20 |  0.405167 TCCCTTGGTCCTGGAGGCAG    63.405         65        0.00
# AKIP1.primer_left_2       515       534        20 |  0.405167 TGTGCCAGAGGAAGGAGGGG    63.405         65        0.00
# AKIP1.primer_left_3       561       580        20 |  0.442332 GGCGAGTCGAAGCTGCACAT    63.442         60       18.61
# AKIP1.primer_left_4       557       576        20 |  0.473386 CAGAGGCGAGTCGAAGCTGC    63.473         65       18.61
# ARFIP2.primer_left_0      1594      1613        20 |  0.296325 CCCAGAAGTTGCTGCCCTGT    62.704         60        0.00
# ARFIP2.primer_left_1      1588      1607        20 |  0.583893 TGGAGGCCCAGAAGTTGCTG    62.416         60        0.00
# ARFIP2.primer_left_2      1538      1557        20 |  0.642407 CTGGGGCCTGACACCAGTTT    62.358         60       14.98
# ARFIP2.primer_left_3      1568      1587        20 |  0.655122 GCTATGGTGGGAAGAGGGCC    62.345         65        0.00
# ARFIP2.primer_left_4      1590      1609        20 |  0.721365 GAGGCCCAGAAGTTGCTGCC    63.721         65        0.00
# self_end_th hairpin_th end_stability
# <numeric>  <numeric>     <numeric>
#   AKIP1.primer_left_0           0      34.44          4.75
# AKIP1.primer_left_1           0      34.44          4.85
# AKIP1.primer_left_2           0       0.00          4.79
# AKIP1.primer_left_3           0       0.00          3.21
# AKIP1.primer_left_4           0       0.00          5.25
# ARFIP2.primer_left_0           0      34.56          4.00
# ARFIP2.primer_left_1           0      37.58          4.41
# ARFIP2.primer_left_2           0      45.21          2.66
# ARFIP2.primer_left_3           0       0.00          5.80
# ARFIP2.primer_left_4           0      41.17          4.85
####

#######################
######## Primer3 ######
#######################
# github.com/jensenlab/primer3 
# primer3 provides a direct interface to the thermodynamic calculations in Primer3 (melting temperatures, primer dimers,
# secondary structures, etc.) The package is self-contained and does not require a separate installation of Primer3

calculate_tm(c("AAGCCGCGTACGA","AAGAGCGATGACG"))
# [1] 44.8188 38.743
calculate_hairpin("CCCCATCCGATCAGGGGGG")
# $structure_found
# [1] TRUE
# 
# $temp
# [1] 50.76991
# 
# $ds
# [1] -99.40729
# 
# $dh
# [1] -32200
# 
# $dg
# [1] -1368.829
# 
# $align_end_1
# [1] -32200
# 
# $align_end_2
# [1] -96
calculate_homodimer("CCCCATCCGATCAGGGGGG")
# $structure_found
# [1] TRUE
# 
# $temp
# [1] 7.101923
# 
# $ds
# [1] -240.7316
# 
# $dh
# [1] -77600
# 
# $dg
# [1] -2937.099
# 
# $align_end_1
# [1] 18
# 
# $align_end_2
# [1] 19
calculate_dimer("CCCCATCCGATCAGGGGGG", "TTTTCCTTTAAAATTGGGGTAATAGCCCATATCCTCCACATCTCAGAAGGTTCATTTT")
# $structure_found
# [1] TRUE
# 
# $temp
# [1] -12.13625
# 
# $ds
# [1] -105.2097
# 
# $dh
# [1] -36900
# 
# $dg
# [1] -4269.206
# 
# $align_end_1
# [1] 5
# 
# $align_end_2
# [1] 44



### Find Target Sequence from UCSC genome browser in R, 300 bp around the probe of interest

# beware of the strandness
seq_chr4 <- getSeq(Hsapiens, "4", start = 47000 - 3000, end = 47000 - 1, strand = "-")
seq_chr4
# 3000-letter DNAString object
# seq: TTTTCCTTTAAAATTGGGGTAATAGCCCATATCCTCCACATCTCAGAAGGTTCATTTTTATTGGTC...GTGAGTATTTTAATCAAATTCTCAAACTTGAGAAAGGGGTTATGGGAGTCCTCACTTCTTAGTAGC
toString(seq_chr4)
# [1] "TTTTCCTTTAAAATTGGGGTAATAGCCCATATCCTCCACATCTCAGAAGGTTCATTTTTATTGGTCCAAAGGAAGTATATGAAAACACTTCCAAAAAAAGTGTAGGTGTTATCAGTATATCATACCATGTGAGAGAAGGTGTAAAGCAACAAAAGGAGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTTGAGATGGAGTTTCACTCTTGTTGCCCAGGCTGGAGTGCAATGGCGCCATCTCAGCTCACTGCAACCTCTGCCTCCCAGGTTCAAGCGATTCTCCTGCTTCAGCCTCCCAAGCAGCTGGGATTACAGGCATGCACCACCATGACCAGCCAATTTTGTGTTGTTAGTAGAGATGGGGTTTCTCCATGTTGGTCAGGTTGGTCTCGAACGCCCGACCTCAGGTGATCCGCCTGCCTTGGCTTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACTGCGCCTGGACAGGAGTTTGGTTTTTAAAAATTTAACACAAGAACATGTATATGCATTTGAATGTACACCAGCCTTTAAAGCAGTCACCATTGGAAGCCATATATACCAAATATCATGAGGGAGAAACTGGAACAAAAAGATGAGTTGGCTATAACTTACGGAAAAGTGAAGAAAATTTCAAAAACCTCTTTATGTGAATGCTCAACTATATTTAAAACACATCAATTCTATGAGGTAGGTATTAACCCTTTCTTACAGAAGAAGGAACTGAACCCCATCAAGGCTACACATAATAAGTGGCAGAGGTGAAATTTGAACCCAGATCTACCAACTCCAAGTTTCTTTTCAATTTGACTATAACTGCATTTGTGATTGCTTGTGACATAAGACATACGAGAGAGTAGATTTGAATTATGTATAAGCAGCTGAGAAAAATGCAATTGAGGGACTTTTCTGAATTGAGAAAAGGATGACATTCTTTGGAGAGCAGACTTTCTACCCTCGCCATCTATCATTGAATGTTATATTCATTGCTGCCACTAGTTAAATATTAGCTAGACATTGACTGGTAGTTACAATCAATGACACAGTCACTGGCCTAATTACTTTGGTGATTTTTATACCATGGTGTCTGACATTCCCAGGTGTTTTGTATAAGGCACACGTATTTCAAGAGCTCTTCAGCCAGGCCTAATCTCTATGCTGATGCCTCAGTCTTGCTTGCTGTTGGAAAGCTTGGCATATTCTCCAGACATAGTCAGAGGACTATATTTCAAAGCTATGTGTGGGTGCACTGAATTCTATGCCTCAGATGATTCCACTGTGGAAAATGATGTGGAAACCCAGTTGAAACGCTCTCACTTCAGTTAACTTTTATATCACAAATCTTAAGTCTGGAAAAGGAAAACAAACTCTTGTGTCCATACTTATTCTTTCTTAAGACGTGAAAACAATAGTACATTTTAACATGGCTGAAGAAAAATTAATTTATGTGACTGTTTATAGAATATAATATAGTTTTATAAGATAAAAATTAAGGACAATATTTTACTTTTTTTTTTTTTTGAGACAAAGTTTTGTTCTTGTTGCCCAGGCTGGAATGCAATGGCATGATCTCAGCTCACTGCAACCTCTGCCTCCCAGGTTTAAGCAATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGATTACCAGAGTACACCACCAGGCCCGGCTAATTTTGTACTTCTTTAGTAGAGACAGGGTTTCACCATGTTGGTCAGGCTGGTCTTGAATCCCTGAACTCAAGTGATCCACCCCTGCACACTCGGCCTCCCAAAATGCTGGAATATTTTACTTTTGTATACATCTAAAAAGTAATTATTAATTTTGAAGATTGAGTCTTCTAAGTTAAAATAGACAGTTGTAACTCATAGACAGAACGTTCACAACTCACAGATACCTAAATAGTGCAGTTTTTGCTTTTAGTTCTGAAATGAAGATTGTCATGAAGCAGCTTATCACTCTCTGTGCCTTCACCATTCTCCTGGGTTATGATTTGGTGAAAAGGCTGCGGGAAAAAATGTGGAACATGTTCAACAGATAGGGAGAAGTGAGAAAGCCGATGGTACCTACAGAAATGAAATCCTTGGCTGGGAGTGGTGGCTCACGCCTGTAATCCCAACACTTTGGGAGGCCGAGGCAGGTGGATCACCTGAGGTCAGGAGTTTGAGACTAGTGTGGCCAACATGGTGAAACCCCTTCTCTACTAAAAATACCAAAAATTAGCCGGGTGTAGTGGCGGACACCTGTAATCCCAGCTACTCAGGAGGCTGAGACAGGAGAATCACTTGAACCTGGGAGGGAGAGGTTGCAGTGAGCTGAGATTGCGCCATTACACTCTGGCCTGGGCAACAAGAGTGAAACTCTGTCTCAAAAAAAAAAAAAAAAAAAAGAAATCCTCTCCAATCAACAAATTAGGTGACTGAATTCTTTCACAATTCTTGCATTTCCATTCAGGACTTGGTGGCTATTCTTTTCAAAGTAGAGGCATGTTGTGCCTTATAGTACCACATACAACTTAACAAAAATATGTAGTTATAAAACTTGCCTTCCTGGAGGGTATTTGTAGCTAAAAGTAAACTAGATTTCAACTTAGATTGCCAGTAAAACACATTATAAACTATAAAGTATTACAAAAATAGGTATATAATATTAAGGTTACATGTGTATCTAAATCTTCCCTAAATCTGAATTTGTTTTATTTGACAAAGTCTATAAACAACCCCAGGACAATAAAAATTGTATCCTCAGAAAGTTTGTCCCAATCCCCTTTCATCCTATCCCCAGCTGCATCTGCCTACAAGTCCCCAGCCTGCCCAGGCTCTGTAGCTTCTATACACCTGTTCCATTTCCAACTGTTCCTGTGTTGTATATTTTATAACAAACTGTTATACACAAGTACAGTGTTTGGCTGAGTTCTGTGAGTATTTTAATCAAATTCTCAAACTTGAGAAAGGGGTTATGGGAGTCCTCACTTCTTAGTAGC"
seq_chr4@shared
seq_chr4@offset
seq_chr4@length
seq_chr4@elementMetadata
seq_chr4@metadata

#### perform in silico bisulfite conversion
## have done the Python version, re-write to R code

#######

## IDT OligoAnalyzer for secondary structures analysis
# Input ? Output ?
# Calculate for GC content, melting temperature (Tm), molecular weight, extinction coefficient, ug/OD, nmol/OD, identify secondary
# structure potential, minimize dimerization, use NCBI BLAST
# Analysis, Hairpin, Self-Dimer, Hetero-Dimer, NCBI blast, TM Mismatch
# parameter sets, Target type: DNA/RNA, Oligo Conc, Na conc, Mg conc, dNTP conc,

# UNAFold www.unafold.org

# Homo-Dimer Analysis
# delta G is calculated by taking into account the longest stretch of complementary bases.
# Dotted lines represent additional complementary bases for that dimer structure, but their presence does not impact calculated delta G values.
# Actual delta G values may vary based on presence of additional complementary bases. The Maximum Delta G value
# refers to the free energy of the oligo sequence binding to its perfect complement.

# Similar tools: Mfold, RNAstructure, DNASIS, ViennaRNA (RNAfold)

# RNAfold

BiocManager::install("LncFinder")
library(LncFinder)
?run_RNAfold


 
g_m <- file(paste0(strsplit(fname, "\\.")[][], "_meth.txt"), "w")
g_u <- file(paste0(strsplit(fname, "\\.")[][], "_unmeth.txt"), "w")

data <- file("Primer3_data/primer3.txt")
data
# primer3_me <- file(paste0(strsplit(fname, "\\.")[][], "_me_primer3.txt"), "w")
# primer3_un <- file(paste0(strsplit(fname, "\\.")[][], "_un_primer3.txt"), "w")

m_index <- c()
i_2 <- ""
i <- "CCCTTCGAATCCGCCG"

for(i in readLines(data)){
  index <- gregexpr("CG", i)[[1]]
  print(index)
  # for(m in index){
  #   m_index <- c(m_index, m)
  # }
  cat("original:", "\n", i, "\n")
  #i_2 <- i
  # for(k in seq_along(m_index)){
  #   substring(i_1, m_index[k], m_index[k]) <- "T"
  # }
  i_1 <- gsub("C", "T", i)
  #cat("bisulfite converted:", "\t", i_1, "\n")
  i_2 = i_1
  #for(j in index){
    #cat("@",j, i_2, "\n")
  #i_2[index] <- "C"
  stringsplit <- strsplit(i_2, "")[[1]]
  stringsplit[index] <- "C"
  result_i_2 <- paste(stringsplit, collapse = "")
  #print(result_i_2)
    #i_2 <- cat(substring(i_2,1,j), "C", substring(i_2,j+1,nchar(i_2)))
  #}
  #cat("Unmethylated converted:", "\n", i_1, "\n")
  cat("Unmethylated converted:", "\n", i_1, "\n")
  cat("methylated converted:", "\n", result_i_2, "\n")
}

result_i_2
close(f)


data <- rnorm(100)

# Design a interactive primer design utility
library(shiny)
ui <- fluidPage(
  titlePanel("Primer Design"),
  sidebarLayout(
    sidebarPanel(
      selectInput("plotType", "Choose a plot type:",
                  choices = c("Histogram" = "hist",
                              "Boxplot" = "boxplot"))
    ),
    mainPanel(
      plotOutput("plot")
    )
  )
)

server <- function(input, output){
  output$plot <- renderPlot({
   if(input$plotType == "hist") {
     hist(data, main = "Histogram")
   }else if (input$plotType == "boxplot"){
     boxplot(data, main = "Boxplot")
   }
  })
}

shinyApp(ui = ui, server = server)






###########################
# Load the package
library(rnafolding)

# Define the FASTA path from the sample 5S rRNA sequence
seq_5S <- system.file("extdata", "5S.fasta", package = "rnafolding")
seq_5S

# Fold 5S with Sliding Windows
windows_5S <- rnafolding::fold(
  filename = seq_5S,
  winsize = 50,
  stepsize = 5
)
windows_5S

pdf("windows_5S.pdf")
plot_structure_map(
  windows_5S,
  plot_list = c("entropy", "bpp")
)
dev.off()
?plot_structure_map

############################