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
		help="Number of PCs to extract [optional]"),
make_option("--extract", action="store", default='NA', type='character',
		help="Number of PCs to extract [optional]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]"),
make_option("--pop_data", action="store", default=NA, type='character',
		help="Population data on reference individuals [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(caret)
library(ggplot2)
library(cowplot)

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
# Extract SNPs from the reference
###

if(!is.na(opt$extract)){
	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat(paste0('Extracting SNPs in ',opt$extract,'...'))
	sink()
	for(i in 1:22){
		system(paste0(opt$plink,' --bfile ',opt$ref_plink_chr,i,' --make-bed --extract ',opt$extract,' --out ',opt$output_dir,'ref.chr',i,' --memory ',floor(opt$memory/0.5)))
	}
	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('Done!\n')
	sink()
}

###
# Label duplicates from the reference
###

for(i in 1:22){
	success <- FALSE
	bim<-fread(paste0(opt$output_dir,'ref.chr',i,'.bim'))
	while (!success) {
		print(sum(duplicated(bim$V2)))	
		bim$V2[duplicated(bim$V2)]<-paste0(bim$V2[duplicated(bim$V2)],'_dup')
		success <- sum(duplicated(bim$V2)) == 0
	}	
	write.table(bim, paste0(opt$output_dir,'ref.chr',i,'.bim'), col.names=F, row.names=F, quote=F)
}

###
# Merge the per chromosome reference genetic data and update IDs to be clearly distinguishable from the target
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Merging per chromosome reference data...')
sink()

# Create merge list
if(!is.na(opt$extract)){
	ref_merge_list<-paste0(opt$output_dir,'ref.chr',1:22)
} else {
	ref_merge_list<-paste0(opt$ref_plink_chr,1:22)
}

write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)

# Create file to update IDs
ref_fam<-data.frame(fread(paste0(opt$ref_plink_chr,'22.fam')))
ref_ID_update<-data.frame(OLD_FID=ref_fam$V1,
                          OLD_IID=ref_fam$V2,
                          NEW_FID=paste0('REF_',ref_fam$V1),
                          NEW_IID=paste0('REF_',ref_fam$V2))
                          
write.table(ref_ID_update, paste0(opt$output_dir,'ref_id_update.txt'), row.names=F, col.names=F, quote=F)

system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --make-bed --update-ids ',opt$output_dir,'ref_id_update.txt --out ',opt$output_dir,'ref_merge --memory ',floor(opt$memory/0.5)))

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
system(paste0(opt$plink,' --bfile ',opt$output_dir,'ref_merge --exclude ',opt$output_dir,'long_ld.exclude --indep-pairwise 1000 5 0.2 --out ',opt$output_dir,'ref_merge --memory ',floor(opt$memory/0.5)))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Calculate allele frequencies in the reference data
###
# This ensures the PCs are not dependent on the target sample.

system(paste0(opt$plink,' --bfile ',opt$output_dir,'ref_merge --freqx --out ',opt$output_dir,'ref_merge --memory ',floor(opt$memory/0.5)))

###
# Calculate PCs in the reference sample for scaling the target sample factor scores.
###

# Extract LD independent SNPs
system(paste0(opt$plink,' --bfile ',opt$output_dir,'ref_merge --extract ',opt$output_dir,'ref_merge.prune.in --make-bed --out ',opt$output_dir,'ref_merge_pruned --memory ',floor(opt$memory/0.5)))

# Calculate PCs
system(paste0(opt$plink,' --bfile ',opt$output_dir,'ref_merge_pruned --read-freq ',opt$output_dir,'ref_merge.frqx --pca ',opt$n_pcs,' --out ',opt$output_dir,'ref_merge_pruned --memory ',floor(opt$memory/0.5)))

# Read in reference PC scores and calculate the mean and SD
PCs_ref<-data.frame(fread(paste0(opt$output_dir,'ref_merge_pruned.eigenvec')))
names(PCs_ref)<-c('FID','IID',paste0('PC',1:100))
PCs_ref_centre_scale<-data.frame(PC=names(PCs_ref[-1:-2]),
								  Mean=sapply(PCs_ref[,-1:-2], function(x) mean(x)),
								  SD=sapply(PCs_ref[,-1:-2], function(x) sd(x)),
								  row.names=seq(1:100))

PCs_ref_scale<-PCs_ref
for(i in names(PCs_ref)[-1:-2]){
PCs_ref_scale[i]<-PCs_ref[i]-PCs_ref_centre_scale$Mean[PCs_ref_centre_scale$PC == i]
PCs_ref_scale[i]<-PCs_ref[i]/PCs_ref_centre_scale$SD[PCs_ref_centre_scale$PC == i]
}

PCs_ref_scale$FID<-sub('REF_','',PCs_ref_scale$FID)
PCs_ref_scale$IID<-sub('REF_','',PCs_ref_scale$IID)

write.table(PCs_ref_scale, paste0(opt$output,'.reference.eigenvec'), col.names=T, row.names=F, quote=F)

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
  
	if(!is.na(opt$extract)){
	  for(chr in 1:22){
		system(paste0(opt$plink,' --bfile ',opt$target_plink_chr,chr,' --make-bed --extract ',opt$extract,' --keep ',opt$output_dir,'targ_batch',batch,'_keep.txt --update-ids ',opt$output_dir,'targ_batch',batch,'_id_update.txt --out ',opt$output_dir,'targ_batch',batch,'_chr',chr,' --memory ',floor(opt$memory/0.5)))
	  }
	} else{
 	  for(chr in 1:22){
		system(paste0(opt$plink,' --bfile ',opt$target_plink_chr,chr,' --make-bed --keep ',opt$output_dir,'targ_batch',batch,'_keep.txt --update-ids ',opt$output_dir,'targ_batch',batch,'_id_update.txt --out ',opt$output_dir,'targ_batch',batch,'_chr',chr,' --memory ',floor(opt$memory/0.5)))
	  }
	}
	
  targ_batch_merge_list<-paste0(opt$output_dir,'targ_batch',batch,'_chr',1:22)
  write.table(targ_batch_merge_list, paste0(opt$output_dir,'targ_batch',batch,'_mergelist.txt'), row.names=F, col.names=F, quote=F)
  
	###
	# Label duplicates in the target
	###

	for(i in 1:22){
			success <- FALSE
			bim<-fread(paste0(opt$output_dir,'targ_batch',batch,'_chr',i,'.bim'))
			bim$V1<-i
  	while (!success) {
			print(sum(duplicated(bim$V2)))	
			bim$V2[duplicated(bim$V2)]<-paste0(bim$V2[duplicated(bim$V2)],'_dup2')
    	success <- sum(duplicated(bim$V2)) == 0
  	}
			write.table(bim, paste0(opt$output_dir,'targ_batch',batch,'_chr',i,'.bim'), col.names=F, row.names=F, quote=F)
	}

  system(paste0(opt$plink,' --merge-list ',opt$output_dir,'targ_batch',batch,'_mergelist.txt --out ',opt$output_dir,'targ_batch',batch,'_merge --memory ',floor(opt$memory/0.5)))

  ###
  # Merge the target_batch and ref data
  ###

  log<-system(paste0(opt$plink,' --bfile ',opt$output_dir,'targ_batch',batch,'_merge --bmerge ',opt$output_dir,'ref_merge --make-bed --out ',opt$output_dir,'ref_targ_batch',batch,'_merge --memory ',floor(opt$memory/0.5)))

	if(log ==3){  
  		system(paste0(opt$plink,' --bfile ',opt$output_dir,'targ_batch',batch,'_merge --exclude ',opt$output_dir,'ref_targ_batch1_merge-merge.missnp --make-bed --out ',opt$output_dir,'targ_batch',batch,'_merge --memory ',floor(opt$memory/0.5)))
		log<-system(paste0(opt$plink,' --bfile ',opt$output_dir,'targ_batch',batch,'_merge --bmerge ',opt$output_dir,'ref_merge --make-bed --out ',opt$output_dir,'ref_targ_batch',batch,'_merge --memory ',floor(opt$memory/0.5)))
	}
  
  # Delete per chromosome target data
  system(paste0('rm ', opt$output_dir,'targ_batch',batch,'_chr*'))
  
  ###
  # Extract SNPs identified as LD independent.
  ###

  system(paste0(opt$plink,' --bfile ',opt$output_dir,'ref_targ_batch',batch,'_merge --extract ',opt$output_dir,'ref_merge.prune.in --make-bed --out ',opt$output_dir,'ref_targ_batch',batch,'_merge_pruned --memory ',floor(opt$memory/0.5)))
  
  ###
  # Perform PCA using reference individuals only
  ###
  
  # Create cluster membership file
  targ_ref_fam<-data.frame(fread(paste0(opt$output_dir,'ref_targ_batch',batch,'_merge_pruned.fam')))
  targ_fam_batch_clust<-data.frame(targ_ref_fam[,1:2])
  targ_fam_batch_clust$Cluster<-'TARG'
  targ_fam_batch_clust$Cluster[grepl('REF_',targ_fam_batch_clust$V1)]<-'REF'
  write.table(targ_fam_batch_clust, paste0(opt$output_dir,'ref_targ_batch',batch,'_merge_pruned.clusters'), col.names=F, row.names=F, quote=F)
  
  system(paste0(opt$plink,' --bfile ',opt$output_dir,'ref_targ_batch',batch,'_merge_pruned --read-freq ',opt$output_dir,'ref_merge.frqx --pca ',opt$n_pcs,' --pca-cluster-names REF --within ',opt$output_dir,'ref_targ_batch',batch,'_merge_pruned.clusters --out ',opt$output_dir,'ref_targ_batch',batch,'_merge_pruned --memory ',floor(opt$memory/0.5)))

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
  PCs_targ_batch_scale$FID<-sub('TARG_','',PCs_targ_batch_scale$FID)
  PCs_targ_batch_scale$IID<-sub('TARG_','',PCs_targ_batch_scale$IID)

	PCs_targ<-rbind(PCs_targ,PCs_targ_batch_scale)

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

if(!is.na(opt$extract)){
	for(i in 1:22){
		system(paste0('rm ',opt$output_dir,'ref.chr',i,'.bed'))
		system(paste0('rm ',opt$output_dir,'ref.chr',i,'.bim'))
		system(paste0('rm ',opt$output_dir,'ref.chr',i,'.fam'))
		system(paste0('rm ',opt$output_dir,'ref.chr',i,'.nosex'))
		system(paste0('rm ',opt$output_dir,'ref.chr',i,'.log'))
	}
}

###
# Create plot PC scores of target sample compared to the reference
###
if(!is.na(opt$pop_data)){
	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('Plotting target sample PCs on reference...')
	sink()

	# Read in population data
	pop_data<-data.frame(fread(opt$pop_data))
	names(pop_data)[1]<-'IID'

	# Plot the reference sample PCs
	ref_PCs<-data.frame(fread(paste0(opt$output,'.reference.eigenvec')))
	ref_PCs<-merge(ref_PCs, pop_data, by='IID')
	ref_PCs$FID<-NULL

	# Read in the target sample PCs
	targ_PCs<-data.frame(fread(paste0(opt$output,'.eigenvec')))
	targ_PCs$FID<-NULL

	new_cols<-names(ref_PCs[!grepl('PC|ID', names(ref_PCs))])
	new_cols_2<-data.frame(matrix(rep('Target',length(new_cols)),ncol=length(new_cols)))
	names(new_cols_2)<-names(ref_PCs[!grepl('PC|ID', names(ref_PCs))])
	targ_PCs<-cbind(targ_PCs,new_cols_2)

	# Combine the two sets
	ref_PCs_targ_PCs<-rbind(ref_PCs,targ_PCs)

	Label_groups<-names(ref_PCs_targ_PCs[!grepl('PC|IID',names(ref_PCs_targ_PCs))])

	for(i in Label_groups){
		PC_1_2<-ggplot(ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] != 'Target',], aes(x=PC1,y=PC2, colour=get(i))) + 
		  geom_point() + 
			geom_point(data=ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] == 'Target',], aes(x=PC1,y=PC2), size=3, colour='black') + 
		  ggtitle("PCs 1 and 2") +
			labs(colour="")
		PC_3_4<-ggplot(ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] != 'Target',], aes(x=PC3,y=PC4, colour=get(i))) + 
		  geom_point() + 
			geom_point(data=ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] == 'Target',], aes(x=PC3,y=PC4), size=3, colour='black') + 
		  ggtitle("PCs 3 and 4") +
			labs(colour="")
		PC_5_6<-ggplot(ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] != 'Target',], aes(x=PC5,y=PC6, colour=get(i))) + 
		  geom_point() + 
			geom_point(data=ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] == 'Target',], aes(x=PC5,y=PC6), size=3, colour='black') + 
		  ggtitle("PCs 5 and 6") +
			labs(colour="")
		PC_7_8<-ggplot(ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] != 'Target',], aes(x=PC7,y=PC8, colour=get(i))) + 
		  geom_point() + 
			geom_point(data=ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] == 'Target',], aes(x=PC7,y=PC8), size=3, colour='black') + 
		  ggtitle("PCs 7 and 8") +
			labs(colour="")
		PC_9_10<-ggplot(ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] != 'Target',], aes(x=PC9,y=PC10, colour=get(i))) + 
		  geom_point() + 
			geom_point(data=ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] == 'Target',], aes(x=PC9,y=PC10), size=3, colour='black') + 
		  ggtitle("PCs 9 and 10") +
			labs(colour="")

		png(paste0(opt$output,'.PCs_plot_',i,'.png'), units='px', res=300, width=4000, height=2500)
		print(plot_grid(PC_1_2,PC_3_4,PC_5_6,PC_7_8,PC_9_10))
		dev.off()

		print(i)
	}

	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('Done!\n')
	sink()

	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('Building C5.0 tree to predicting most likely population for target samples.\n')
	sink()

	for(i in Label_groups){
		model <- train(y=as.factor(ref_PCs[[i]]), x=ref_PCs[grepl('PC',names(ref_PCs))], method="C5.0Tree", metric='logLoss', trControl=trainControl(method="cv", number=10, classProbs= TRUE, savePredictions = 'final', summaryFunction = multiClassSummary))
		pred<-predict(object = model$finalModel, newdata = data.matrix(targ_PCs[grepl('PC',names(targ_PCs))]), type = "prob")

		write.table(data.frame(IID=targ_PCs$IID,round(pred,3)),paste0(opt$output,'.',i,'_prediction.txt'), col.names=T, row.names=F, quote=F)
		
		sink(file = paste(opt$output,'.log',sep=''), append = T)
		cat('Model predicting',i,'has an AUC of ',model$results$AUC,'.\n')
		sink()

		print(i)
	}
} 

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
cat('Project PCs are here: ',opt$output,'.eigenvec.\n',sep='')
sink()
