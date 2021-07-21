####==========================================================================================####
####========================Step4: PCA analysis                   ============================####
####==========================================================================================####

genodir="/afm01/UQ/Q3457/SPARK/hrc_imputed/QCd/" ##genotype dir
outdir="/afm01/UQ/Q3457/SPARK/hrc_imputed/QCd/grm/"  ## grm output dir
scriptdir="/scratch/90days/uqywan67/auti_proj/AutismSpeaks_SPARK/scripts/" ##script dir for weight_frq
bfile="SPARK_imputed_autosome_geno.05_maf.01_hwe1e-6_info.3" ##bfile prefix
projc="SPARK" ##cohort name for cohort to be projected
r2="0.3"##LD pruned r2 threshold
data=paste0(genodir, bfile)
kgdir="/afm01/UQ/Q3117/1kgAll-phase3/"
ref="1000G_phase3_hapmap3_autosome"
refdata=paste0(kgdir, ref)
refc="1000G"

###4.1 pruend SNPs first with ~100K SNPs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
system(paste0("plink2 --bfile ",  data, "  --extract /afm01/UQ/Q3457/SSC/hrc_imputed/QCd/hapmap3_autosome.snplist  --indep-pairwise 100 50 ", r2, " --threads 4 --out ", outdir, bfile, "_qcd_pruned", r2))


###4.2 common SNPs with 1KG~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
bim <- read.table(paste0(refdata, ".bim"), stringsAsFactors = F)
bim2 <- read.table(paste0(data, ".bim"), stringsAsFactors = F)
bim$comb1 <- paste0(bim$V2, bim$V5,bim$V6)
bim$comb2 <- paste0(bim$V2, bim$V6,bim$V5)
bim2$comb <- paste0(bim2$V2, bim2$V5, bim2$V6)

coms <- bim[bim$comb1 %in% bim2$comb | bim$comb2 %in% bim2$comb,]
prune1 <- read.table(paste0(outdir, bfile, "_qcd_pruned", r2, ".prune.in"), stringsAsFactors = F)
snps <- coms[coms$V2 %in% prune1$V1,]$V2

write.table(snps, file = paste0(outdir, refc, "_com_", projc, "_qcd_pruned", r2, ".snplist"), col.names = F, row.names = F, quote = F)

common_snp <- paste0(outdir, refc, "_com_", projc, "_qcd_pruned", r2, ".snplist")

###4.3 calculate and update frq~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#calculate frq
system(paste0("plink --bfile ", data, " --extract ", common_snp, " --freq --out ", outdir, projc, "_pruned", r2))
system(paste0("plink --bfile ", refdata, " --extract ", common_snp, " --freq --out ", outdir, refc, "_pruned", r2))

##update freq
frq_projc <- paste0(outdir, projc, "_pruned", r2)
frq_refc <- paste0(outdir, refc, "_pruned", r2)

system(paste0("Rscript ", scriptdir, "weight_freq.R ", frq_projc, " ", frq_refc, " ", outdir, projc, "_com_", refc, ".frq"))

###4.4 calculate GRM/PCs in 1KG~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
system(paste0("gcta64 --bfile ", refdata, " --extract ", common_snp, "  --update-freq ",  outdir, projc, "_com_", refc, ".frq --make-grm-bin --out ", outdir, refc, "_pruned", r2, " --threads 4"))
system(paste0("gcta64 --grm  ", outdir, refc, "_pruned", r2, " --pca 20 --out ", outdir, refc, "_pruned", r2, ".pca20 --threads 4"))


###4.5 get the pc loading~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
system(paste0("gcta64 --bfile ", refdata, " --extract ", common_snp, "  --update-freq ",  outdir,  projc, "_com_", refc, ".frq --pc-loading ", outdir, refc, "_pruned", r2, ".pca20 --out ", outdir, refc, "_", projc, "_loading.pruned", r2, " --threads 4"))

###4.6 project the pcs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
system(paste0("gcta64 --bfile ", data, " --extract ", common_snp, " --project-loading ", outdir, refc, "_", projc, "_loading.pruned", r2, " 20 --out ", outdir, refc, "_", projc, "_proj_pruned", r2, " --threads 4"))




