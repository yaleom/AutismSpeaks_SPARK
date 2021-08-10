####=================================================================####
#
# Autism Speaks: Data Analysis Pipeline
#                                  
####=================================================================####


###Delta2 version
genodir="/afm01/UQ/Q3457/SPARK/genotypes"
outdir="/afm01/UQ/Q3490/SPARK/QCd"
bfile="SPARK.WES1.release.2021_03.genotype"
scriptdir="/scratch/90days/uqywan67/auti_proj/AutismSpeaks_SPARK/scripts"
plotdir="/scratch/90days/uqywan67/auti_proj/AutismSpeaks_SPARK/outputs"
rangefile="/afm01/UQ/Q3490/SSC/QCd/high-LD-regions-hg19-GRCh37.txt"
####==========================================================================================####
####========================Step1: Pre-QC using raw genotype files============================####
####==========================================================================================####
###1.1 Check sex~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#get a list of QCd SNPs
plink --bfile $genodir/${bfile} \
--geno 0.05 --hwe 1e-6 --mind 0.1 --maf 0.01 \
--make-just-bim --out $outdir/${bfile}_qcd

#calculate F-statistic based on QCd chr23 SNPs only
plink --bfile $genodir/${bfile} \
 --check-sex  \
 --extract $outdir/${bfile}_qcd.bim \
 --out $outdir/${bfile}_checksex

#F > 0.8 coded as male, F < 0.2 coded as female; Note those individuals with extreme values

##update sex using SNPSEX if necessary [Usually check with phenotype files first]
awk '{print $1, $2, $4}' $outdir/${bfile}_checksex.sexcheck  | awk 'NR>1' > $outdir/${bfile}_toUpdateSex.txt

plink --bfile $genodir/${bfile} \
--update-sex $outdir/${bfile}_toUpdateSex.txt \
--must-have-sex \
--make-bed \
--out $genodir/${bfile}

######################################skip this for now######################################

##summary checksex information 
Rscript $scriptdir/plot_checkSex.R $outdir/${bfile}_checksex.sexcheck $plotdir/sexD.png


###1.2 IBD estimation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#prune SNPs and using only ~50K SNPs, exclude high LD regions
r2=0.02 #might try a range of different r2 values to get ~50K SNPs for further analysis
plink --bfile $genodir/${bfile} \
 --geno 0.05 --maf 0.01 --mind 0.2 --hwe 1e-6 \
 --indep-pairwise 50 10 $r2 \
 --chr 1-22 \
 --exclude range $rangefile \
 --out $outdir/${bfile}_pruned

#estimate based on same family
plink --bfile $genodir/${bfile} \
--extract $outdir/${bfile}_pruned.prune.in \
--genome rel-check \
--out $outdir/${bfile}_pruned_IBD

#plot IBD distribution
Rscript $scriptdir/plot_IBD.R $outdir/${bfile}_pruned_IBD.genome $plotdir/IBD.png

#Note if there is any inbreeding, compare the chrX-F VS IBD estimates

###1.3 Check genome-wide heterozygosity VS missing rates~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#using pruned SNPs
plink --bfile $genodir/${bfile} \
--extract $outdir/${bfile}_pruned.prune.in \
--het \
--out $outdir/${bfile}_pruned_het

plink --bfile $genodir/${bfile} \
--missing \
--autosome \
--out $outdir/${bfile}_missing 

#plot genome-wide heterozygosity/F VS missing rates
Rscript $scriptdir/plot_genomewideF_VS_missing.R $outdir/${bfile}_pruned_het.het $outdir/${bfile}_missing.imiss $plotdir/Het_VS_missing.png $plotdir/inbreeding_VS_missing.png




