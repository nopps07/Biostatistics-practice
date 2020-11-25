#Git link

#Shall we check if push and pull is working

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
a

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
a

BiocManager::available()
BiocManager::install()
a


library(BiocManager)
packageVersion("BiocManager")


#The Role of S4 in Bioconductor

#S4
#Positive
#-Formal definition of classes
#-Bioconductor reusability
#-Has validation of types
#-Naming conventions
#Example: mydescriptor <- new("GenomeDescription)

#Negative
#Complex structure compared to S3

#How do we know whether an object is S4 or not
#isS4(mydescriptor)
#TRUE

#str(mydescriptor)
#Formal class ....

#S4 class describes a representation
#-name. slots(methods/fields). contains(inheritance definition)

#example
MyEpicProject <- setClass(# Define class name with UpperCamelCase
                          "MyEpicProject",
                          # Define slots, helpful for validation
                          slots = c(ini = "Date",
                                    end = "Date",
                                    milestone = "character"),
                          # Define inheritance
                          contains = "MyProject")


#S4 accesors (methods)
#.S4methods(class = "genomeDescription")
#showMethods(classes = "GenomeDescription), where = search())

#object summary
#show(myDescriptor)
BiocManager::install("BSgenome")

library("BSgenome")
#
showClass("BSgenome")


# What is a_genome's main class?
class(a_genome)  # "BSgenome"

# What is a_genome's other classes?
is(a_genome)  # "BSgenome", "GenomeDescription"

# Is a_genome an S4 representation?
isS4(a_genome)  # TRUE

#Investigate the a_genome
show(a_genome)

#Investigate some other accesors
organism(a_genome)
provider(a_genome)
seqinfo(a_genome)

#Introducing biology of genomic datasets

install("BSgenome.Scerevisiae.UCSC.sacCer3")
installed.genomes()
available.genomes()

library(BSgenome.Scerevisiae.UCSC.sacCer3)
yeast <- BSgenome.Scerevisiae.UCSC.sacCer3


str(yeast)
summary(yeast)
length(yeast)
names(yeast)
seqlengths(yeast)

#S4 method getSeq() requires a BSgenome object
getSeq(yeast)
#Select chromosome sequence by name, one or many
getSeq(yeast, "chrM")
#Select start, end and or width
#end = 10, selects first 10 base pairs of each chromosome
getSeq(yeast, end = 10)


yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3
# Get the head of seqnames and tail of seqlengths for yeastGenome
head(seqnames(yeastGenome))
tail(seqlengths(yeastGenome))
# Print chromosome M, alias chrM
yeastGenome$chrM
# Count characters of the chrM sequence
nchar(yeastGenome$chrM)

# Get the first 30 bases of each chromosome
getSeq(yeastGenome, end = 30)

# other options
getSeq(yeastGenome, names = "chrM", start = 10, end = 50)


#THe University of California, Santa Cruze (UCSC)
#Genome Browser has made available the most genomes for BSgenome
#totaling 74 of various species!!


#Introduction to Biostrings

#biological string containers
#memory efficient to store and manipulate sequence of characters
#containers that can be inherited
#The BString comes from big string
showClass("XString")
showClass("BString")
showClass("BStringSet")

#Biostring alphabets
DNA_BASES #DNA 4 bases
RNA_BASES #RNA 4 bases
AA_STANDARD # 20 Amino acids

DNA_ALPHABET # contains IUPAC_CODE_MAP
RNA_ALPHABET # contains IUPAC_CODE_MAP
AA_ALPHABET  # contains AMINO_ACID_CODE

#Transcription DNA to RNA
#DNA single string
3dna_seq <- DNAString("ATGATCTCGTAA")
dna_seq

#T's to U's
rna_seq <- RNAString(dna_seq)
rna_seq

#Translation RNA to amino acids
aa_seq <- translate(rna_seq)
aa_seq
#Three RNA bases form one AA: AUG = M, AUC = I, UCG = S, UAA = *

#shortcut
aa_seq2 <- translate(dna_seq)
aa_seq2

#The Zika virus
install.packages("read.gb")
library(read.gb)
zikaVirus <- read.gb("sequence.gb", DNA = TRUE, Type = "full")

getwd()
setwd("C:/Users/grate/R")

summary(zikaVirus)
head(zikaVirus)

install.packages("read.fasta")

zikaVirus2 <- Biostrings::readDNAStringSet("sequence.fasta")
getSeq(zikaVirus2)
alphabet(zikaVirus)
summary(zikaVirus2)


#alphabet()
#alphabetFrequency()
#alphabet(file, baseOnly = TRUE)

zv_dna_seq <- subseq(unlist(zikaVirus2), end = 21)
zv_dna_seq
#ook it works
zv_rna_seq <- RNAString(zv_dna_seq)
zv_aa_seq <- translate(zv_rna_seq)
zv_aa_seq <- translate(zv_dna_seq)
#both are same
zv_aa_seq


#Sequence handling
#XString to store a single sequence
#BString for any string
#DNAString for DNA
#RNAString for RNA
#AAString for amino acids

#Create a string Set and collate it
length(zikaVirus2)
width(zikaVirus2)
zikaVirus_seq <- unlist(zikaVirus2)
length(zikaVirus_seq)

#From a single sequence to a set
zikaSet <- DNAStringSet(zikaVirus_seq,
                        start = c(1, 101, 201),
                        end = c(100, 200, 300))
zikaSet
length(zikaSet)
width(zikaSet)

#Double DNA?
#complement()
#give you a paired strand

#Rev a sequence
#rev(): reverse the order of your sequences at the same time

#reverseComplement()
#reverse(complement())


#exerciese
subZikv <- subseq(zikaVirus_seq, end= 30)
subZikv

reverse(subZikv)
complement(subZikv)
reverseComplement(subZikv)
translate(zikv)

#Why are we interested in Patterns??
#Zebra, fingerprint, sunflower head
#Patterns in Biology are outstanding and we can learn more about them using sequencing.
#sequence repeats, protein and codons, poly-A tails, conserved sequences, binding sites, and more
#our goal in analysing sequence patterns is to discover their occurence, frequency, periodicity, and length

#What can we find with patterns?
#Gene strat
#Protein end
#Regions that enhance or silence gene expression
#Conserved regions between organisms
#Genetic variation

#matchPattern(pattern, subject)
#1 string to 1 string

#vmatchPattern(pattern, subject)
#1 set of strings to 1 string
#1 string to a set of strings

#Palindromes
#findPalindromes() #Find palindromic regions in a single sequence

#Translation has six possibilities
#From a single DNA string, there are 6 possible string frames.
#Three +, Three -
#A negative strand is the reverse complement of a positive sequence strand.

#single base sliding window?

#exercise
#For Sets
vmatchPattern(pattern = "ACATGGGCCTACCATGGGAG",
              subject = zikaVirus2, max.mismatch = 1)
#For single sequences
matchPattern(pattern = "ACATGGGCCTACCATGGGAG",
             subject = zikaVirus_seq, max.mismatch =  1)

length(zikaVirus2)
length(zikaVirus_seq)

#Finding Palindromes
findPalindromes(zikaVirus_seq)
findPalindromes(subZikv)

#NS5
NS50 <- Biostrings::readDNAStringSet("https://www.uniprot.org/uniprot/A0A0B4ZYS0.fasta")
NS5 <- unlist(NS50) 

vcountPattern(pattern = NS5,
              subject = zikaSet, max.mismatch =  15)

#select the frame that contains the match
#zikaSet is not appropriate here. So just use it as an example
selectedSet <- zikaSet[3]

#Convert this frame into a single sequence
selectedSeq <- unlist(selectedSet)

#a set
vmatchPattern(pattern = NS5, subject = selectedSet,
              max.mismatch = 15)

#a single sequence
matchPattern(pattern = NS5, subject = selectedSeq,
             max.mismatch = 15)


#IRanges and Genomic Structures
library(IRanges)

myIRanges <- IRanges(start = 20, end = 30)
myIRanges

#More IRanges examples
myIRanges_width <- IRanges(start = c(1,20), width = c(30,11))
myIRanges_end <- IRanges(start = c(1, 20), end = 30)

myIRanges_width
myIRanges_end

#Equation: width  = end - start + 1

#Rle - run length encoding
#Rle stands for Run length encoding
#Computes and stores the lengths and values of a vector or factor
#Rle is general S4 container used to save long repetitive vectors efficiently

some_numbers <- c(3, 2, 2, 2, 3, 3, 4, 2)
Rle(some_numbers)
#This is useful to represent sequence rages as very commonly they will have repetitions

#IRanges with logical vector
IRanges(start = c(FALSE, FALSE, TRUE, TRUE))
#IRanges with logical Rle
gi <- c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE)
myRle <- Rle(gi)
myRle

IRanges(start = myRle)

#IRanges are hierarchical data structures can contain metadata
#useful to store genes, transcripts, polymorphisms, GC content, and more.
#start, end or width as numeric vectors (or NULL)
#start argument as a logical vertor or logical Rle object
#Rle stands for Run length encoding and is storage efficient
#IRanges arguments get recycled (fill in the blanks)
#equation for sequence range: width = end - start + 1

#example
# start vector 1 through 5, end 100 
IRnum1 <- IRanges(start = 1:5, end = 100)

# end 100, width 89 and 10
IRnum2 <- IRanges(end = 100, width = c(89, 10))

# logical argument start = Rle(c(F, T, T, T, F, T, T, T))
IRlog1 <- IRanges(start = Rle(c(F, T, T, T, F, T, T, T)))

# Printing objects in a list
print(list(IRnum1 = IRnum1, IRnum2 = IRnum2, IRlog1 = IRlog1))

# Create the first sequence seq_1
seq_1 <- IRanges(start = 10, end = 37)

# Create the second sequence seq_2
seq_2 <- IRanges(start = c(5, 35, 50),
                 end = c(12, 39, 61),
                 names = LETTERS[1:3])
# Check the width and length



#Gene of interest

#examples of genomic intervals
#1. Reads aligned to a reference
#2. Genes of interest
#3. Exonic regions
#4. Single nucleotide polymorphisms (SNPs)
#5. Regions of transcription or binding sites, RNA-seq or CHIP-seq

#Genomic Ranges
library(GenomicRanges)
myGR <- GRanges("chr1:200-300")
#GRanges class is a container to save genomic intervals by chromosome
#seqnames() & seqinfo()
#From data to GRanges

#Genomic Ranges accessors
methods(class = "GRanges")

#used for chromosome names
seqnames()

#returns an IRanges object for ranges
ranges()

#stores metadata columns
mcols()

#generic function to store sequence information
seqinfo()

#stores the genome name
genome()

#Accessors are both setter and getter function
#Accessors can be inherited thanks to S4 definitions

#Chromosome X GRanges
install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg <- TxDb.Hsapiens.UCSC.hg38.knownGene


#select genes from chromosome X
hg_chrXg <- genes(hg, filter = list(tx_chrom = c("chrx")))
hg_chrXg

hg_chrXgp <- genes(hg, filter = list(tx_chrom = c("chrx"), tx_strand = "+"))
sort(hg_chrXgp)
hg_chrXgp

#exercise
myGR <- as(seq_intervals, "GRanges")
myGR
#Q. seq_intervals <- what function is used to make this?

seqnames = c("chrI", "chrI", "chrII", "chrII")
start = c(11, 12, 13, 14)
end = c(36, 37, 38, 39)
seq_intervals <- data.frame(seqnames, start, end)
seq_intervals
#tibble?
?tibble

seqinfo(myGR)
mcols(myGR)


#Manipulating collections of GRanges
#The GrangesList-class is a container for storing a collection of GRanges
#Efficient for storing a large number of elements

#To construct a GRangeList
#as(mylist, "GRangesList")
#GRangesList(myGranges1, myGranges2, ...)

#To convert back to GRanges
#unlist(myGRangesList)
#Accessors method(class = "GRangeslist")

#Examples of GRangesLists:
#transcripts by gene
#exons by transcripts
#read alignments
#sliding windows

#GRanges object with 595 genes
hg_chrXg
slidingWindows(hg_chrXg, width = 20000, step = 10000)


#Genomic features and TxDb
#GenomicFeatures uses transcript database (TxDb) objects to store metadata,
#manage genomic locations and relationships between features and its identifiers

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
(hg <- TxDb.Hsapiens.UCSC.hg38.knownGene)

hg <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(hg) <- c("chrX") #prefilter results to chrX

#transcripts
transcripts(hg, columns = c("tx_id", "tx_name", filter = NULL))

#exons
exons(hg, columns("tx_id", "exon_id"), filter = list(tx_id = "179161"))

#examples?
exonsBytx <- exonsBy(hg, by = "tx") # exons by transcript

abcd1_179161 <- exonsBytx[["179161"]] # transcript id

width(abcd1_179161) #width of each exon, the purple regions of the figure


#countOverlaps results in an integer vector of counts
countOverlaps(query, subject)
#findOverlaps results in a Hits object
findOverlaps(query, subject)
#subsetByOverlaps returns a GRangesList object
subsetByOverlaps(query, subject)

#The length of the gene -1 will return at least 2 windows.
#The length of the last window can be partial.

#Exercise
#there is an overlap between chromosome X and the gene ABCD1
#Let's find its gene id and its location, also called 'locus'

rangefound <- subsetByOverlaps(hg_chrX, ABCD1)
names(rangefound)
ABCD1
rangefound

hg <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(hg) <- c("chrX")
hg_chrXt <- transcriptsBy(hg, by = "gene")
#select gene '215' from the transcripts
hg_chrXt$'215'


# Unlist hg_ChrX and save result as myGR
myGR <- unlist(hg_ChrX)

# Compare classes of hg_ChrX and myGR
class(hg_ChrX)
class(myGR)

# Compare length of hg_ChrX and myGR
length(hg_ChrX)
length(myGR)