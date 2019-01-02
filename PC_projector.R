#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--target_plink_chr", action="store", default=NA, type='character',
		help="Path to per chromosome target PLINK files [required]"),
make_option("--ref_plink_chr", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK files [required]"),
make_option("--plink", action="store", default='plink', type='character',
		help="Path PLINK software binary [optional]"),
make_option("--output", action="store", default='./PC_projector_output/Output', type='character',
		help="Path for output files [optional]"),
make_option("--batch_size", action="store", default=10000, type='numeric',
		help="Number of individuals in each batch [optional]"),
make_option("--n_pcs", action="store", default=20, type='numeric',
		help="Number of PCs to extract [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

tmp<-sub('.*/','',opt$output)
opt$output_dir<-sub(paste0(tmp,'*.'),'',opt$output)
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# PC_projector.R V1.0
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

###
# Merge the per chromosome reference genetic data and update IDs to be clearly distinguishable from the target
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Merging per chromosome reference data...')
sink()

# Create merge list
ref_merge_list<-paste0(opt$ref_plink_chr,1:22)
write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)

# Create file to update IDs
ref_fam<-data.frame(fread(paste0(opt$ref_plink_chr,'22.fam')))
ref_ID_update<-data.frame(OLD_FID=ref_fam$V1,
                          OLD_IID=ref_fam$V2,
                          NEW_FID=paste0('REF_',ref_fam$V1),
                          NEW_IID=paste0('REF_',ref_fam$V2))
                          
write.table(ref_ID_update, paste0(opt$output_dir,'ref_id_update.txt'), row.names=F, col.names=F, quote=F)

system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --make-bed --update-ids ',opt$output_dir,'ref_id_update.txt --out ',opt$output_dir,'ref_merge'))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Create SNP list for LD pruning
# Remove regions of long range LD which can confound estimates of ancestry estimates (REF: PMC2443852)
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Identifying LD independent SNPs based on reference data...')
sink()

# Read in the bim file
ref_bim<-data.frame(fread(paste0(opt$output_dir,'ref_merge.bim')))

# Create file removing these regions.
long_ld_exclude<-ref_bim$V2[ (ref_bim$V1 == 1 & ref_bim$V4 >= 48e6 & ref_bim$V4 <= 52e6) |
                                  (ref_bim$V1 == 2 & ref_bim$V4 >= 86e6 & ref_bim$V4 <= 100.5e6) |
                                  (ref_bim$V1 == 2 & ref_bim$V4 >= 134.5e6 & ref_bim$V4 <= 138e6) |
                                  (ref_bim$V1 == 2 & ref_bim$V4 >= 183e6 & ref_bim$V4 <= 190e6) |
                                  (ref_bim$V1 == 3 & ref_bim$V4 >= 47.5e6 & ref_bim$V4 <= 50e6) |
                                  (ref_bim$V1 == 3 & ref_bim$V4 >= 83.5e6 & ref_bim$V4 <= 87e6) |
                                  (ref_bim$V1 == 3 & ref_bim$V4 >= 89e6 & ref_bim$V4 <= 97.5e6) |
                                  (ref_bim$V1 == 5 & ref_bim$V4 >= 44.5e6 & ref_bim$V4 <= 50.5e6) |
                                  (ref_bim$V1 == 5 & ref_bim$V4 >= 98e6 & ref_bim$V4 <= 100.5e6) |
                                  (ref_bim$V1 == 5 & ref_bim$V4 >= 129e6 & ref_bim$V4 <= 132e6) |
                                  (ref_bim$V1 == 5 & ref_bim$V4 >= 135.5e6 & ref_bim$V4 <= 138.5e6) |
                                  (ref_bim$V1 == 6 & ref_bim$V4 >= 25.5e6 & ref_bim$V4 <= 33.5e6) |
                                  (ref_bim$V1 == 6 & ref_bim$V4 >= 57e6 & ref_bim$V4 <= 64e6) |
                                  (ref_bim$V1 == 6 & ref_bim$V4 >= 140e6 & ref_bim$V4 <= 142.5e6) |
                                  (ref_bim$V1 == 7 & ref_bim$V4 >= 55e6 & ref_bim$V4 <= 66e6) |
                                  (ref_bim$V1 == 8 & ref_bim$V4 >= 8e6 & ref_bim$V4 <= 12e6) |
                                  (ref_bim$V1 == 8 & ref_bim$V4 >= 43e6 & ref_bim$V4 <= 50e6) |
                                  (ref_bim$V1 == 8 & ref_bim$V4 >= 112e6 & ref_bim$V4 <= 115e6) |
                                  (ref_bim$V1 == 10 & ref_bim$V4 >= 37e6 & ref_bim$V4 <= 43e6) |
                                  (ref_bim$V1 == 11 & ref_bim$V4 >= 46e6 & ref_bim$V4 <= 57e6) |
                                  (ref_bim$V1 == 11 & ref_bim$V4 >= 87.5e6 & ref_bim$V4 <= 90.5e6) |
                                  (ref_bim$V1 == 12 & ref_bim$V4 >= 33e6 & ref_bim$V4 <= 40e6) |
                                  (ref_bim$V1 == 12 & ref_bim$V4 >= 109.5e6 & ref_bim$V4 <= 112e6) |
                                  (ref_bim$V1 == 20 & ref_bim$V4 >= 32e6 & ref_bim$V4 <= 34.5e6)]

write.table(long_ld_exclude, paste0(opt$output_dir,'long_ld.exclude'), col.names=F, row.names=F, quote=F)
  
# Identify LD independent SNPs.
system(paste0(opt$plink,' --bfile ',opt$output_dir,'ref_merge --exclude ',opt$output_dir,'long_ld.exclude --indep-pairwise 1000 5 0.2 --out ',opt$output_dir,'ref_merge'))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Calculate allele frequencies in the reference data
###
# This ensures the PCs are not dependent on the target sample.

system(paste0(opt$plink,' --bfile ',opt$output_dir,'ref_merge --freqx --out ',opt$output_dir,'ref_merge'))

###
# Calculate PCs in the reference sample for scaling the target sample factor scores.
###

# Extract LD independent SNPs
system(paste0(opt$plink,' --bfile ',opt$output_dir,'ref_merge --extract ',opt$output_dir,'ref_merge.prune.in --make-bed --out ',opt$output_dir,'ref_merge_pruned'))

# Calculate PCs
system(paste0(opt$plink,' --bfile ',opt$output_dir,'ref_merge_pruned --read-freq ',opt$output_dir,'ref_merge.frqx --pca ',opt$n_pcs,' --out ',opt$output_dir,'ref_merge_pruned'))

# Read in reference PC scores and calculate the mean and SD
PCs_ref<-data.frame(fread(paste0(opt$output_dir,'ref_merge_pruned.eigenvec')))
names(PCs_ref)<-c('FID','IID',paste0('PC',1:100))
PCs_ref_centre_scale<-data.frame(PC=names(PCs_ref[-1:-2]),
								  Mean=sapply(PCs_ref[,-1:-2], function(x) mean(x)),
								  SD=sapply(PCs_ref[,-1:-2], function(x) sd(x)),
								  row.names=seq(1:100))

###
# Define splits for the target sample batches
###

targ_fam<-data.frame(fread(paste0(opt$target_plink_chr,'22.fam')))

batches<-split(1:dim(targ_fam)[1], ceiling(seq_along(1:dim(targ_fam)[1])/opt$batch_size))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Target data will be split into',length(batches),'batches.\n')
sink()

# Create an object to store all PCs for target sample
PCs_targ<-NULL

for(batch in 1:length(batches)){

  ###
  # Extract batch from the target data, update IDs to be clearly distinguishable, and merge per chromosome files
  ###
  targ_fam_batch<-targ_fam[batches[[batch]],]
  targ_fam_batch_ID_update<-data.frame( OLD_FID=targ_fam_batch$V1,
                                        OLD_IID=targ_fam_batch$V2,
                                        NEW_FID=paste0('TARG_',targ_fam_batch$V1),
                                        NEW_IID=paste0('TARG_',targ_fam_batch$V2))
                                        
  write.table(targ_fam_batch_ID_update[3:4], paste0(opt$output_dir,'targ_batch',batch,'_keep.txt'), row.names=F, col.names=F, quote=F)
  write.table(targ_fam_batch_ID_update, paste0(opt$output_dir,'targ_batch',batch,'_id_update.txt'), row.names=F, col.names=F, quote=F)
  
  for(chr in 1:22){
    system(paste0(opt$plink,' --bfile ',opt$target_plink_chr,chr,' --make-bed --keep ',opt$output_dir,'targ_batch',batch,'_keep.txt --update-ids ',opt$output_dir,'targ_batch',batch,'_id_update.txt --out ',opt$output_dir,'targ_batch',batch,'_chr',chr))
  }
  
  targ_batch_merge_list<-paste0(opt$output_dir,'targ_batch',batch,'_chr',1:22)
  write.table(targ_batch_merge_list, paste0(opt$output_dir,'targ_batch',batch,'_mergelist.txt'), row.names=F, col.names=F, quote=F)
  
  system(paste0(opt$plink,' --merge-list ',opt$output_dir,'targ_batch',batch,'_mergelist.txt --out ',opt$output_dir,'targ_batch',batch,'_merge'))

  ###
  # Merge the target_batch and ref data
  ###

  system(paste0(opt$plink,' --bfile ',opt$output_dir,'targ_batch',batch,'_merge --bmerge ',opt$output_dir,'ref_merge --make-bed --out ',opt$output_dir,'ref_targ_batch',batch,'_merge'))

  # Delete per chromosome target data
  system(paste0('rm ', opt$output_dir,'targ_batch',batch,'_chr*'))
  
  ###
  # Extract SNPs identified as LD independent.
  ###

  system(paste0(opt$plink,' --bfile ',opt$output_dir,'ref_targ_batch',batch,'_merge --extract ',opt$output_dir,'ref_merge.prune.in --make-bed --out ',opt$output_dir,'ref_targ_batch',batch,'_merge_pruned'))
  
  ###
  # Perform PCA using reference individuals only
  ###
  
  # Create cluster membership file
  targ_ref_fam<-data.frame(fread(paste0(opt$output_dir,'ref_targ_batch',batch,'_merge_pruned.fam')))
  targ_fam_batch_clust<-data.frame(targ_ref_fam[,1:2])
  targ_fam_batch_clust$Cluster<-'TARG'
  targ_fam_batch_clust$Cluster[grepl('REF_',targ_fam_batch_clust$V1)]<-'REF'
  write.table(targ_fam_batch_clust, paste0(opt$output_dir,'ref_targ_batch',batch,'_merge_pruned.clusters'), col.names=F, row.names=F, quote=F)
  
  system(paste0(opt$plink,' --bfile ',opt$output_dir,'ref_targ_batch',batch,'_merge_pruned --read-freq ',opt$output_dir,'ref_merge.frqx --pca ',opt$n_pcs,' --pca-cluster-names REF --within ',opt$output_dir,'ref_targ_batch',batch,'_merge_pruned.clusters --out ',opt$output_dir,'ref_targ_batch',batch,'_merge_pruned'))

  ###
  # Read in the PCs to be scaled to the reference and combined with the other batches.
  ###
  
  PCs_ref_targ_batch<-data.frame(fread(paste0(opt$output_dir,'ref_targ_batch',batch,'_merge_pruned.eigenvec')))
  PCs_targ_batch<-PCs_ref_targ_batch[grepl('TARG_', PCs_ref_targ_batch$V1),]
  names(PCs_targ_batch)<-c('FID','IID',paste0('PC',1:opt$n_pcs))

  # Scale based on mean and sd in the reference
  PCs_targ_batch_scale<-PCs_targ_batch
  for(i in names(PCs_targ_batch)[-1:-2]){
	PCs_targ_batch_scale[i]<-PCs_targ_batch[i]-PCs_ref_centre_scale$Mean[PCs_ref_centre_scale$PC == i]
	PCs_targ_batch_scale[i]<-PCs_targ_batch[i]/PCs_ref_centre_scale$SD[PCs_ref_centre_scale$PC == i]
  }

  # Convert the IDs back to match the target data
  PCs_targ_batch$FID<-sub('TARG_','',PCs_targ_batch$FID)
  PCs_targ_batch$IID<-sub('TARG_','',PCs_targ_batch$IID)

  PCs_targ<-rbind(PCs_targ,PCs_targ_batch)

  ###
  # Delete files that are not required for the next batch
  ###
  
  system(paste0('rm ',opt$output_dir,'*batch',batch,'*'))
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Batch',batch,'complete.\n')
  sink()

}

write.table(PCs_targ, paste0(opt$output,'.eigenvec'), col.names=T, row.names=F, quote=F)

###
# Clean up temporary files
###

system(paste0('rm ',opt$output_dir,'long_ld.exclude'))
system(paste0('rm ',opt$output_dir,'ref_id_update.txt'))
system(paste0('rm ',opt$output_dir,'ref_merge*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
cat('Project PCs are here: ',opt$output,'.eigenvec.\n',sep='')
sink()
