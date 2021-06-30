####==========================================================================================####
####========================Step3: After Imputation                ===========================####
####==========================================================================================####

#Note see the scripts in /QRISdata/Q3490/SSC/scripts


#===============3.1 convert to plink format [Note the delimiters of FID_IID "_"]
#===============3.2 extract .info information
#===============3.3 fix missing IDs in the .bim file
#===============3.4 fix missing IDs in the .info file

#for i in {1..3};do
# qsubshcom "Rscript afterImp_QC_aut.R ${i} {TASK_ID}" 1 30G aftImp-aut-${i} 24:00:00 "-array=1-22"
#done
#qsubshcom "Rscript afterImp_QC_chrX.R {TASK_ID}" 1 30G aftImp-chrX 24:00:00 "-array=1-3"

#===============3.5 merge bfiles &.info for each cohort
#-need to update .fam if there are multiple "_" in FID_IID
#-extract a list of .info > 0.3
#qsubshcom "Rscript afterImp_QC2_aut.R {TASK_ID}" 1 120G aftImp2-aut 24:00:00 "-array=1-3"
#qsubshcom "Rscript afterImp_QC2_chrX.R {TASK_ID}" 1 60G aftImp-chrX 24:00:00 "-array=1-3"

#===============3.6 merge bfiles for all cohorts
#-based on .info > 0.3 for each cohort
#qsubshcom "Rscript afterImp_QC3_aut.R" 1 300G aftImp3-aut 24:00:00 ""


#===============3.7 Frequency check for imputed data
#for common SNPs with MAF >1% intersecting with HRC reference panel
Rscript /30days/uqywan67/SSC/scripts/plot_freqCheck.R frqfile frqCheckOutfile QcdBimfile FigureFile


#===============3.8 PCA analysis
#based on HM3 SNPs or genotype QCd SNPs











