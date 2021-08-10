####==========================================================================================####
####========================Step3: After Imputation                ===========================####
####==========================================================================================####

#Note see the scripts in /QRISdata/Q3490/SSC/scripts

#===============3.1 convert to plink format [Note the delimiters of FID_IID "_", if two "_" then use "--const-fid 0]
##for X-chromosome
wkdir="/afm01/UQ/Q3457/SPARK/hrc_imputed"
imp_folder="CHRX.vcfs"
cohort="SPARK"
qsubshcom "Rscript vcfToplink_format.R --wkdir ${wkdir} --imp_folder ${imp_folder} --cohort ${cohort} --chr X" 1 60g vcftoplink-X 24:00:00 ""

##for autosomal chromosomes
imp_folder="SPARK-ALL.vcfs"
cohort="SPARK"
wkdir="/afm01/UQ/Q3457/SPARK/hrc_imputed"
qsubshcom "Rscript vcfToplink_format.R --wkdir ${wkdir} --imp_folder ${imp_folder} --cohort ${cohort} --chr {TASK_ID}" 1 60g vcftoplink 24:00:00 "-array=1-22"


#===============3.2 merge bfiles &.info 
cohort="SPARK"
wkdir="/afm01/UQ/Q3457/SPARK/hrc_imputed"
qsubshcom "Rscript merge_plink.R --wkdir ${wkdir} --threads 1 --cohort ${cohort} " 1 120G merge-plink 24:00:00 ""

#===============3.3 re-check sex information based on imputed chr-X
#calculate F-statistic based on QCd chr23 SNPs only
genodir="/afm01/UQ/Q3457/SPARK/hrc_imputed/plink_bfile"
bfile="SPARK_imputed_chrX"
snps="/afm01/UQ/Q3457/SPARK/hrc_imputed/QCd/SPARK_imputed_chrX_geno.05_maf.01_hwe1e-6_info.3.bim"
outdir="/afm01/UQ/Q3490/SPARK/QCd"

plink --bfile $genodir/${bfile} \
 --check-sex  \
 --extract ${snps} \
 --out $outdir/${bfile}_checksex_imputed
 
##only 71 SNPs left....something must be wrong

#===============3.4  Frequency check for imputed data
#for common SNPs with MAF >1% intersecting with HRC reference panel
root="/afm01/UQ/Q3457/SPARK/hrc_imputed/QCd"
frqfile=${root}/SPARK_imputed_autosome_geno.05_maf.01_hwe1e-6_info.3_freq.frq
frqCheckOutfile=${root}/SPARK_imputed_autosome_qcd_hrc_frqcheck
QcdBimfile=${root}/SPARK_imputed_autosome_geno.05_maf.01_hwe1e-6_info.3.bim
FigureFile="/scratch/90days/uqywan67/auti_proj/AutismSpeaks_SPARK/outputs/frqcheck.png"
scriptdir="/scratch/90days/uqywan67/auti_proj/AutismSpeaks_SPARK/scripts"
Rscript ${scriptdir}/plot_freqCheck.R $frqfile $frqCheckOutfile $QcdBimfile $FigureFile


#===============3.5 PCA analysis
#based on HM3 SNPs using the scripts in 4_PCA_SPARK.R













