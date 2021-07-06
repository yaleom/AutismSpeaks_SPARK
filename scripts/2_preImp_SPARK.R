
####==========================================================================================####
####========================Step2: Pre-imputation                 ============================####
####==========================================================================================####
library(data.table)
###Delta2 version
genodir="/afm01/UQ/Q3457/SPARK/genotypes/"
outdir="/afm01/UQ/Q3490/SPARK/QCd/"
bfile="SPARK.WES1.release.2021_03.genotype"
scriptdir="/scratch/90days/uqywan67/auti_proj/AutismSpeaks_SPARK/scripts/"
plotdir="/scratch/90days/uqywan67/auti_proj/AutismSpeaks_SPARK/outputs/"
rangefile="/afm01/UQ/Q3490/SSC/QCd/high-LD-regions-hg19-GRCh37.txt"

strandfile="/afm01/UQ/Q3490/SPARK/GSA-24v1-0_A1-b38.strand"
strandbuild37="/afm01/UQ/Q3490/SPARK/GSA-24v1-0_A1-b37.strand"
ref="/afm01/UQ/Q3490/SSC/QCd/human_g1k_v37.fasta"
OUTDIR <- outdir

###2.1 generate QCd file~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
system(paste0("plink --bfile ", genodir, bfile, " --chr 1-23, 25 --geno 0.05 --maf 0.01 --mind 0.1 --hwe 1e-6 --make-bed --out ", OUTDIR, bfile, "_allchrs"))

#Update chr23 and chr25 to chrX
bim <- read.table(paste0(OUTDIR, bfile, "_allchrs.bim"))
bim[bim$V1 == 23 | bim$V1 == 25,]$V1 <- "X"
chrs <- bim[,c(2,1)]
write.table(chrs, file = paste0(OUTDIR, bfile, "_allchrs.chr"), col.names = F, row.names = F, quote = F, sep = "\t")
system(paste0("plink --bfile ", OUTDIR, bfile, "_allchrs --output-chr M --update-chr ", OUTDIR, bfile, "_allchrs.chr --make-bed --out ", OUTDIR, bfile, "_allchrs_update"))



###2.2 update SNP IDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
##some SNP ids including GSA, removing those IDs
##but the code for indels should be updated/removed, here we removed with only 29
pos <- read.table(paste0(strandfile), stringsAsFactors = F)
bim3 <- fread(paste0(OUTDIR, bfile, "_allchrs_update.bim"), stringsAsFactors = F)
als <- c("A", "G", "C", "T")
bim3_2 <- bim3[bim3$V5 %in% als,]
#bim3_2[,new := gsub("GSA-", "", V2)]

# snps <- bim3_2[,.(V2, new)]
# fwrite(snps, file =  paste0(OUTDIR, bfile, "_updateName.snplist"), sep = "\t", col.names = F)

write.table(bim3_2$V2, file = paste0(OUTDIR, bfile, "_allchrs_update.snplist"), col.names = F, row.names = F, quote = F, sep = "\t")
system(paste0("plink --bfile ", OUTDIR, bfile, "_allchrs_update --extract ", OUTDIR, bfile, "_allchrs_update.snplist --output-chr M --update-chr ", OUTDIR, bfile, "_allchrs.chr --make-bed --out ", OUTDIR, bfile, "_allchrs_update"))

## generate plink bfile based on update-name and remove duplicates
#system(paste0("plink --bfile ", OUTDIR, bfile, "_allchrs_update --update-name ", OUTDIR, bfile, "_updateName.snplist --make-bed --out ",  OUTDIR, bfile, "_allchrs_update" ))


###2.3 update build ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
##update from build 36 to build 37 and flip minus strand
system(paste0("/afm01/UQ/Q3490/SSC/QCd/update_build.sh ", OUTDIR, bfile, "_allchrs_update  ", strandbuild37, " ",  OUTDIR, bfile, "_allchrs_flip"))

#remove chr24 and update chr23 to chrX
system(paste0("plink --bfile ", OUTDIR, bfile, "_allchrs_flip --output-chr M --chr 1-23 --make-bed --out ", OUTDIR, bfile, "_allchrs_flip2"))

#==============================================================================
print("step 2.3 finished: update build")





###2.4 update kgp IDs to rsids ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

library(plyr)
#BiocManager::install("SNPlocs.Hsapiens.dbSNP142.GRCh37")
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
bimname <- paste0(OUTDIR, bfile, "_allchrs_flip2.bim")
bimname2 <- paste0(bfile, "_allchrs_flip2")
bim <- read.table(bimname)
bim$index <- 1:nrow(bim)

# For each chromosome download the SNP locations and match to .bim file
a <- ddply(bim, .(V1), .progress="text", function(x)
{
  x <- mutate(x)
  chr <- paste("ch", x$V1[1], sep="")
  snps <- getSNPlocs(chr)
  snps <- subset(snps, loc %in% x$V4, select=-c(alleles_as_ambig))
  snps$RefSNP_id <- paste("rs", snps$RefSNP_id, sep="")
  snps <- subset(snps, !duplicated(loc))
  snps <- subset(snps, !duplicated(RefSNP_id))

  x <- merge(x, snps, by.x="V4", by.y="loc", all.x=T)

 x <- x[order(x$index), ]
  index <- !is.na(x$RefSNP_id)
  x$V2[index] <- x$RefSNP_id[index]
  x <- subset(x, select=c(V1, V2, V3, V4, V5, V6))
  return(x)
})

# If there are duplicated SNPs for any reason then label them as _2, _3 etc
temp <- rle(a$V2)
temp2 <- paste0(rep(temp$values, times = temp$lengths), "_", unlist(lapply(temp$lengths, seq_len)))
temp2 <- gsub("_1", "", temp2)
a$V2 <- temp2

##CAUTION: the SNPs in the bim file are sorted in 1, 11, 12, .....19, 2, 21, 22, 3, 4, 5, 6, 7,..
##reorder the SNPs, note chrX
sorted.a=a[with(a, order(as.numeric(V1),V4)),]

###remove GSA- in SNP ids
sorted.a$V2 <- gsub("GSA-", "", sorted.a$V2)
write.table(a, file = paste0(bimname, "2"), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(sorted.a, file = paste0(OUTDIR, bimname2, "_updateIDs.bim"), col.names = F, row.names = F, quote = F, sep = "\t")

## rename .fam .bed as .bim
system(paste0("mv ", OUTDIR, bimname2, ".fam ", OUTDIR, bimname2, "_updateIDs.fam"))
system(paste0("mv ", OUTDIR, bimname2, ".bed ", OUTDIR, bimname2, "_updateIDs.bed"))

#==============================================================================
print("step 2.4 finished: update SNP IDs")





###2.5 Genotype QCs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
system(paste0("plink --bfile ", OUTDIR, bfile, "_allchrs_flip2_updateIDs --autosome --make-bed --out ", OUTDIR, bfile, "_qcd_autosome"))
system(paste0("plink --bfile ", OUTDIR, bfile, "_allchrs_flip2_updateIDs --chr 23, 25 --make-bed --out ", OUTDIR, bfile, "_qcd_chrX"))

#==============================================================================
print("step 5 finished: genotype QC")





###2.6  transform to .vcf format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
##(should separated to chromosome if use michigen imputation server)
#system(paste0("plink --bfile ",  OUTDIR, bfile, "_allchrs_flip2 --set-hh-missing --recode vcf --out ", OUTDIR, bfile, "_allchrs_flip2_recode"))
#system(paste0("awk '{if($0 ~ /#/){if($0 ~ /ID=23/){gsub(\"ID=23\", \"ID=X\", $0); print;next} else{print}} else {gsub(23, \"X\", $1); print}}' ",  OUTDIR, bfile, "_allchrs_flip2_recode.vcf > ",  OUTDIR, bfile, "_allchrs_flip2_recodeII.vcf
#"))
#system(paste0("vcf-sort ", OUTDIR, bfile, "_allchrs_flip2_recodeII.vcf | bgzip -c > ", OUTDIR, bfile, "_allchrs_flip2_recode.vcf.gz"))

#Note error occurs when combined autosome and chrX in fixed reference step, therefore, autosome and chrX will be separated.

######autosome only
system(paste0("plink --bfile ",  OUTDIR, bfile, "_allchrs_flip2_updateIDs --autosome --recode vcf --out ", OUTDIR, bfile, "_allchrs_flip2_recode"))
system(paste0("~/vcftools/src/perl/vcf-sort ", OUTDIR, bfile, "_allchrs_flip2_recode.vcf | ~/.local/bin/bin/bgzip -c > ", OUTDIR, bfile, "_allchrs_flip2_recode.vcf.gz"))
##using vcftools/HTSlib (https://gist.github.com/adefelicibus/f6fd06df1b4bb104ceeaccdd7325b856)

## remove the SNPs in hh file for X chromosome if a hh file was generated. otherwise the imputation will fail.
system(paste0("plink --bfile ",  OUTDIR, bfile, "_allchrs_flip2_updateIDs --chr 23  --set-hh-missing --recode vcf --out ", OUTDIR, bfile, "_allchrs_flip2_recode_chrX"))
system(paste0("sed 's/23/X/' ", OUTDIR, bfile, "_allchrs_flip2_recode_chrX.vcf > ", OUTDIR, bfile, "_allchrs_flip2_recode_chrXII.vcf"))
system(paste0("~/vcftools/src/perl/vcf-sort ", OUTDIR, bfile, "_allchrs_flip2_recode_chrXII.vcf | ~/.local/bin/bin/bgzip -c > ", OUTDIR, bfile, "_allchrs_flip2_recode_chrX.vcf.gz"))

system(paste0("rm ",OUTDIR, bfile, "_allchrs_flip2_recode.vcf"))
system(paste0("rm ",OUTDIR, bfile, "_allchrs_flip2_recode.log"))
system(paste0("rm ",OUTDIR, bfile, "_allchrs_flip2_recodeII.vcf"))

system(paste0("rm ",OUTDIR, bfile, "_allchrs_flip2_recode_chrX.vcf"))
system(paste0("rm ",OUTDIR, bfile, "_allchrs_flip2_recode_chrXII.vcf"))
system(paste0("rm ",OUTDIR, bfile, "_allchrs_flip2_recode_chrX.log"))


#==============================================================================
print("step 6 finished: transformed to vcf")
#Note remember to change 23 into X in the vcf file if not.



###2.7  update REF allele based on HRC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
input <- paste0(OUTDIR, bfile, "_allchrs_flip2_recode.vcf.gz")
output <- paste0(OUTDIR, bfile, "_allchrs_flip2_recode_fixed.vcf.gz")
system(paste0("bcftools +fixref ", input, " -- -f ", ref))
system(paste0("bcftools +fixref ", input, " -Oz -o ", output," -- -d -f ", ref, " -m flip"))
system(paste0("bcftools +fixref ", output, " -- -f ", ref))
##using bcftools from https://samtools.github.io/bcftools/howtos/install.html

##for chrX
input2 <- paste0(bfile, "_allchrs_flip2_recode_chrX.vcf.gz")
output2 <- paste0(bfile, "_allchrs_flip2_recode_chrX_fixed.vcf.gz")
system(paste0("bcftools +fixref ", input2, " -- -f ", ref))
system(paste0("bcftools +fixref ", input2, " -Oz -o ", output2," -- -d -f ", ref, " -m flip"))
system(paste0("bcftools +fixref ", output2, " -- -f ", ref))

#==============================================================================
print("step 2.7 finished: update ref based on HRC ")

#remove temp files
system(paste0("rm ",OUTDIR, bfile, "_allchrs.*"))
system(paste0("rm ",OUTDIR, bfile, "_allchrs_update.*"))
system(paste0("rm ",OUTDIR, bfile, "_allchrs_updateName.*"))
system(paste0("rm ",OUTDIR, bfile, "_allchrs_flip.*"))
system(paste0("rm ",OUTDIR, bfile, "_allchrs_updateName_exDups.*"))
system(paste0("rm ",OUTDIR, bfile, "_allchrs_flip2.*"))




