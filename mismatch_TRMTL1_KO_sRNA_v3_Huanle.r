####################################################
### Missmatch profiles for mouse TRMT1L WTvsKO sRNA
### Eva Maria Novoa
### Jan 2018
####################################################


############################
## 1. LOAD MISSMATCH DATA
############################

setwd("~/Dropbox/Eva-Huanle/trmt1l_sRNASeq/HAMR_Huanle_m22G/try2/mismatch_sorted_hamr_files")

from_file_to_clean_dat<-function(myfile) {
    
    dat<-read.table(myfile, header=T)
    dat<-dat[,1:9]
    colnames(dat)<-c("chr","pos","strand","ref_nuc","A","C","G","T","num_miss")
    print(head(dat))
    dat$coverage<-rowSums(dat[,5:8])
    dat$missmatch_freq<-dat$num_miss/dat$coverage
    row.names(dat)<-paste(paste(dat$chr,dat$pos,sep="."),dat$strand,sep=".")
    return(dat)
}

#WT1<-from_file_to_clean_dat("WT1_mismatches_sorted.txt")
#WT2<-from_file_to_clean_dat("WT2_mismatches_sorted.txt")
#WT3<-from_file_to_clean_dat("WT3_mismatches_sorted.txt")
#HO1<-from_file_to_clean_dat("HO1_mismatches_sorted.txt")
#HO2<-from_file_to_clean_dat("HO2_mismatches_sorted.txt")
#HO3<-from_file_to_clean_dat("HO3_mismatches_sorted.txt")

#save(WT1,file="WT1.Rdata")
#save(WT2,file="WT2.Rdata")
#save(WT3,file="WT3.Rdata")
#save(HO1,file="HO1.Rdata")
#save(HO2,file="HO2.Rdata")
#save(HO3,file="HO3.Rdata")

load("HO1.Rdata")
load("HO2.Rdata")
load("HO3.Rdata")
load("WT1.Rdata")
load("WT2.Rdata")
load("WT3.Rdata")

## Stats

dim(HO1) # 498603  // try2: 3427542
dim(HO2) # 270043  // try2: 1641089
dim(HO3) # 272220  // try2: 1730202
dim(WT1) # 256686  // try2: 1617063
dim(WT2) # 2594034 // try2: 43660041
dim(WT3) # 232160  // try2: 1433599


############################
## 2. MERGE BY REPLICABILITY
############################


## 2.1 Reproducible positions
#############################

HO.reps.IDs<-intersect(intersect(row.names(HO1),row.names(HO2)),row.names(HO3))
WT.reps.IDs<-intersect(intersect(row.names(WT1),row.names(WT2)),row.names(WT3))

length(HO.reps.IDs) # 162962 // try2: 781419
length(WT.reps.IDs) # 152779 // try2: 777775

colnames(HO2)<-paste(colnames(HO2),".rep2",sep="")
colnames(HO3)<-paste(colnames(HO3),".rep3",sep="")

colnames(WT2)<-paste(colnames(WT2),".rep2",sep="")
colnames(WT3)<-paste(colnames(WT3),".rep3",sep="")

# Merge with all the other datasets too
########################################

HO1_HO2<-merge(HO1,HO2, by=0, all=FALSE)
row.names(HO1_HO2)<-HO1_HO2$Row.names
HO1_HO2<-HO1_HO2[,2:dim(HO1_HO2)[2]]

HO.reps<-merge(HO1_HO2,HO3, by=0, all=FALSE)
row.names(HO.reps)<-HO.reps$Row.names
HO.reps<-HO.reps[,2:dim(HO.reps)[2]]


WT1_WT2<-merge(WT1,WT2, by=0, all=FALSE)
row.names(WT1_WT2)<-WT1_WT2$Row.names
WT1_WT2<-WT1_WT2[,2:dim(WT1_WT2)[2]]

WT.reps<-merge(WT1_WT2,WT3, by=0, all=FALSE)
row.names(WT.reps)<-WT.reps$Row.names
WT.reps<-WT.reps[,2:dim(WT.reps)[2]]


# SANITY CHECK: Check that we are capturing tRNAPhe
HO.reps[row.names(HO.reps)=="chr10.80249010.-",]
WT.reps[row.names(WT.reps)=="chr10.80249010.-",]

HO.reps[row.names(HO.reps)=="chr13.21161429.+",]
WT.reps[row.names(WT.reps)=="chr13.21161429.+",]

HO.reps[row.names(HO.reps)=="chr13.21882931.-",]
WT.reps[row.names(WT.reps)=="chr13.21882931.-",]

dim(HO.reps) # 781419
dim(WT.reps) # 777775

#save(HO.reps,file="HO.reps.Rdata")
#save(WT.reps,file="WT.reps.Rdata")

load("HO.reps.Rdata")
load("WT.reps.Rdata")

###################################################
## 3. FILTER DATA (min COVERAGE, not being a SNP )
###################################################

# a) For replicable positions

from_alldat_to_missmatched_noSNPS<- function(dat) {
	
    print ("Original dataset:")
    print(dim(dat)[1])
    
    # 1. Filter based on minimum coverage
	dat_missmatched<-dat[dat$coverage >= 50,]
	dat_missmatched<-dat[dat$coverage.rep2 >= 50,]
	dat_missmatched<-dat[dat$coverage.rep3 >= 50,]
    
    print ("Filtering by coverage:")
    print(dim(dat_missmatched)[1])
	
	# 2. Filter based on  min_num_missmatch
	dat_missmatched<-dat_missmatched[dat_missmatched$num_miss >=5 ,]#  Minimum of 5 SNPs
	dat_missmatched<-dat_missmatched[dat_missmatched$num_miss.rep2 >=5 ,]#  Minimum of 5 SNPs
	dat_missmatched<-dat_missmatched[dat_missmatched$num_miss.rep3 >=5 ,]#  Minimum of 5 SNPs

    print ("Filtering by min_num_missmatch:")
    print(dim(dat_missmatched)[1])

	# 3. Filter based on  top mismatch_frequency
	dat_missmatched<-dat_missmatched[dat_missmatched$missmatch_freq <= 0.4 ,]
	dat_missmatched<-dat_missmatched[dat_missmatched$missmatch_freq.rep2 <= 0.4 ,]
	dat_missmatched<-dat_missmatched[dat_missmatched$missmatch_freq.rep3 <= 0.4 ,]

    print ("Filtering by top_mismatch_frequency:")
    print(dim(dat_missmatched)[1])

	# 4. Filter based on  bottom mismatch_frequency
	dat_missmatched<-dat_missmatched[dat_missmatched$missmatch_freq >= 0.05 ,] 
	dat_missmatched<-dat_missmatched[dat_missmatched$missmatch_freq.rep2 >= 0.05 ,] 
	dat_missmatched<-dat_missmatched[dat_missmatched$missmatch_freq.rep3 >= 0.05 ,] 

    print ("Filtering by bottom_mismatch_frequency:")
    print(dim(dat_missmatched)[1])

	# 5. Add ratio info
	dat_missmatched$TGratio<-dat_missmatched$T/dat_missmatched$G
	dat_missmatched$TGratio.rep2<-dat_missmatched$T.rep2/dat_missmatched$G.rep2
	dat_missmatched$GTratio.rep3<-dat_missmatched$T.rep3/dat_missmatched$G.rep3
	
	# Print stats
    print ("Mismatch distribution based on reference nucleotide:")
    print(table(dat_missmatched$ref_nuc))

	print ("Final filtered dataset:")
	print(dim(dat_missmatched)[1])
	
	return (dat_missmatched)
}


HO.reps.filtered<-from_alldat_to_missmatched_noSNPS(HO.reps)
WT.reps.filtered<-from_alldat_to_missmatched_noSNPS(WT.reps)

dim(HO.reps.filtered) # 32074
dim(WT.reps.filtered) # 29705

## Check that I didn't throw away tRNAPhe:
HO.reps.filtered[row.names(HO.reps.filtered)=="chr10.80249010.-",]
WT.reps.filtered[row.names(WT.reps.filtered)=="chr10.80249010.-",]

HO.reps.filtered[row.names(HO.reps.filtered)=="chr13.21161429.+",]
WT.reps.filtered[row.names(WT.reps.filtered)=="chr13.21161429.+",]

HO.reps.filtered[row.names(HO.reps.filtered)=="chr13.21882931.-",]
WT.reps.filtered[row.names(WT.reps.filtered)=="chr13.21882931.-",]


# b) For individual replicates


from_alldat_to_missmatched_noSNPS_single_rep<- function(dat) {
    
    print ("Original dataset:")
    print(dim(dat)[1])
    
    # 1. Filter based on minimum coverage
    dat_missmatched<-dat[dat$coverage >= 50,]
    
    print ("Filtering by coverage:")
    print(dim(dat_missmatched)[1])
    
    # 2. Filter based on  min_num_missmatch
    dat_missmatched<-dat_missmatched[dat_missmatched$num_miss >=5 ,]#  Minimum of 5 SNPs
    
    print ("Filtering by min_num_missmatch:")
    print(dim(dat_missmatched)[1])
    
    # 3. Filter based on  top mismatch_frequency
    dat_missmatched<-dat_missmatched[dat_missmatched$missmatch_freq <= 0.4 ,]
    
    print ("Filtering by top_mismatch_frequency:")
    print(dim(dat_missmatched)[1])
    
    # 4. Filter based on  bottom mismatch_frequency
    dat_missmatched<-dat_missmatched[dat_missmatched$missmatch_freq >= 0.05 ,]
    
    print ("Filtering by bottom_mismatch_frequency:")
    print(dim(dat_missmatched)[1])
    
    # 5. Add ratio info
    dat_missmatched$TGratio<-dat_missmatched$T/dat_missmatched$G
    
    # Print stats
    print ("Mismatch distribution based on reference nucleotide:")
    print(table(dat_missmatched$ref_nuc))
    
    print ("Final filtered dataset:")
    print(dim(dat_missmatched)[1])
    
    return (dat_missmatched)
}

HO1.filtered<-from_alldat_to_missmatched_noSNPS_single_rep(HO1)
HO2.filtered<-from_alldat_to_missmatched_noSNPS_single_rep(HO2)
HO3.filtered<-from_alldat_to_missmatched_noSNPS_single_rep(HO3)

WT1.filtered<-from_alldat_to_missmatched_noSNPS_single_rep(WT1)
WT2.filtered<-from_alldat_to_missmatched_noSNPS_single_rep(WT2)
WT3.filtered<-from_alldat_to_missmatched_noSNPS_single_rep(WT3)

dim(HO1.filtered) # 72363
dim(HO2.filtered) # 54812
dim(HO3.filtered) # 57852
dim(WT1.filtered) # 55443
dim(WT2.filtered) # 271506
dim(WT3.filtered) # 49694


########################################################
### 4. INTERSECTION OF WT AND KO --> TRMT1L CANDIDATES
########################################################

# Intersection--> how many in WT that are lost in KO?

## a) For replicable, non-filtered
####################################

# Common
WTvsKO.reps.IDs<-intersect(row.names(HO.reps), row.names(WT.reps))

# Different
WTnotinKO.reps<-WT.reps[!(row.names(WT.reps) %in% row.names(HO1)),]
WTnotinKO.reps<-WTnotinKO.reps[!(row.names(WTnotinKO.reps) %in% row.names(HO2)),]
WTnotinKO.reps<-WTnotinKO.reps[!(row.names(WTnotinKO.reps) %in% row.names(HO3)),]

# Stats
length(WTvsKO.reps.IDs)
dim(WTnotinKO.reps) # 18485
prop.table(table(WTnotinKO.reps$ref_nuc))
#A         C         G         T
#0.2463619 0.2385177 0.3197187 0.1954017

## b) For replicable, filtered
####################################

# Common
WTvsKO.reps.filtered.IDs<-intersect(row.names(HO.reps.filtered), row.names(WT.reps.filtered))

# Different
WTnotinKO.reps.filtered<-WT.reps.filtered[!(row.names(WT.reps.filtered) %in% row.names(HO1.filtered)),]
WTnotinKO.reps.filtered<-WTnotinKO.reps.filtered[!(row.names(WTnotinKO.reps.filtered) %in% row.names(HO2.filtered)),]
WTnotinKO.reps.filtered<-WTnotinKO.reps.filtered[!(row.names(WTnotinKO.reps.filtered) %in% row.names(HO3.filtered)),]

# Stats
length(WTvsKO.reps.filtered.IDs)
dim(WTnotinKO.reps.filtered) # 1159
prop.table(table(WTnotinKO.reps.filtered$ref_nuc))

#        A         C         G         T
#0.3347714 0.1518550 0.3149267 0.1984469

############################################################
### 4. MISSMATCH SIGNATURES
#############################################################


library("ggtern")

myggtern_plot_A<-function(mydf,C,G,T,mytitle) {
    plot<-ggtern(mydf,aes(C,G,T))+ geom_point(aes(color = log(coverage)), mydf) + ggtitle(mytitle)
    return(plot)
    #plot+ theme_bw() + theme_nogrid()
}

myggtern_plot_C<-function(mydf,A,G,T,mytitle) {
    plot<-ggtern(mydf,aes(A,G,T))+ geom_point(aes(color = log(coverage)), mydf) + ggtitle(mytitle)
    return(plot)
    #plot+ theme_bw() + theme_nogrid()
}

myggtern_plot_G<-function(mydf,A,C,T,mytitle) {
    plot<-ggtern(mydf,aes(A,C,T))+ geom_point(aes(color = log(coverage)), mydf) + ggtitle(mytitle)
    return(plot)
    #plot+ theme_bw() + theme_nogrid()
}

myggtern_plot_T<-function(mydf,A,C,G,mytitle) {
    plot<-ggtern(mydf,aes(A,C,G))+ geom_point(aes(color = log(coverage)), mydf) + ggtitle(mytitle)
    return(plot)
    #plot+ theme_bw() + theme_nogrid()
}

myggtern_plot_for_4bases<-function(dat_freq,sampleName){
    dat.A<-dat_freq[dat_freq$ref_nuc=="A",]
    dat.C<-dat_freq[dat_freq$ref_nuc=="C",]
    dat.G<-dat_freq[dat_freq$ref_nuc=="G",]
    dat.T<-dat_freq[dat_freq$ref_nuc=="T",]
    pdf(paste(sampleName,"TERN_PLOT.pdf",sep="."),height=7,width=7)
    plot(myggtern_plot_A(dat.A,C,G,T,paste("A missmatches",sampleName,sep=", ")))
    plot(myggtern_plot_C(dat.C,A,G,T,paste("C missmatches",sampleName,sep=", ")))
    plot(myggtern_plot_G(dat.G,A,C,T,paste("G missmatches",sampleName,sep=", ")))
    plot(myggtern_plot_T(dat.T,A,C,G,paste("T missmatches",sampleName,sep=", ")))
    dev.off()
}


myggtern_plot_for_4bases(WTnotinKO.reps.filtered,"WTnotinKO.reps.filtered")
myggtern_plot_for_4bases(WTnotinKO.reps,"WTnotinKO.reps")



##################################################################
### 5. WHERE ARE IN THESE WTnoinKO.reps.filtered positions???
##################################################################

# a) Overlap with reference annotations (previously prepared)
#############################################################
Gmismatches<-WTnotinKO.reps.filtered[WTnotinKO.reps.filtered$ref_nuc=="G",] # Candidate positions to study: 365

# Convert to GRanges
library(GenomicRanges)
Gmismatches.gr<-makeGRangesFromDataFrame(Gmismatches,
    keep.extra.columns=TRUE,
    starts.in.df.are.0based=FALSE,
    start.field="pos",
    end.field="pos"
)

# Load TX annotations
library(GenomicFeatures)
txdb<-loadDb(file="~/reference_data/mm10/txdb_makeTxDbFromGFF.mm10.gencodeM14.sqlite")
snoRNA_txdb<-loadDb(file="~/reference_data/mm10/txdb_makeTxDbFromGFF.mm10.gencodeM14.snoRNAs.sqlite")

mm10_transcripts <- transcripts(txdb)
snoRNA_mm10_transcripts <- transcripts(snoRNA_txdb)

# Intersect with TX annotations
table(Gmismatches.gr %over% mm10_transcripts) # 231 overlaps
#FALSE  TRUE
#231   134
table(Gmismatches.gr %over% snoRNA_mm10_transcripts) # 11 overlaps
#FALSE  TRUE
#364     1

# b) Closest genes analysis
############################
library(ChIPpeakAnno)
library(EnsDb.Mmusculus.v79)

## Annotate my GR dataset (Gmismatches)
annoData <- toGRanges(EnsDb.Mmusculus.v79) # USING V79 ENSEMBL !!!!
seqlevelsStyle(Gmismatches.gr) <- seqlevelsStyle(annoData)# keep the seqnames in the same style
anno <- annotatePeakInBatch(Gmismatches.gr, AnnotationData=annoData)

# Add gene symbol to GR
library(org.Mm.eg.db)
anno$symbol <- xget(anno$feature, org.Mm.egSYMBOL)
anno <- addGeneIDs(anno, orgAnn="org.Mm.eg.db",feature_id_type="ensembl_gene_id",IDs2Add=c("symbol"))
head(anno)

# Stats of overlapping results
table(anno[anno$insideFeature=="inside"]$symbol)





