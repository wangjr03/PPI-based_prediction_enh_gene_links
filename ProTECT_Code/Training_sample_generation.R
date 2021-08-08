# 4 arguments
# path to significant chromatin interactions
# path to contact domain annotations
# path of the output file
# index of cell-type
library(bedtoolsr)
args <- commandArgs(T)
gs_path <- args[1]
domain_path <- args[2]
train_pairs_out <- args[3] # output path of the processed training pairs
ct <- as.numeric(args[4])

enh_path = "../data/enhancer_coords.bed" # path to the enhancer coordinates
gene_act_mat_path = "../data/RPKM_gene_exp_mat.bed" # path to gene activity matrix
gene_name_path = "../data/RPKM_gene_name.bed" # path to the gene names
gene_annotation_path = "../data/Gene_annotation.bed" # path to the gene annotation files

enhancer <- read.table(enh_path)
gene_act_mat <- read.table(gene_act_mat_path)
gene_name <- read.table(gene_name_path)
gene_anno <- read.table(gene_annotation_path,colClasses = "character")
gs <- read.table(gs_path)


#for GM12878 check 48 column
act_gene <- subset(gene_name,gene_act_mat[,ct]>0)

gene_anno <- subset(gene_anno,gene_anno[,1]%in%act_gene[,1])

get_pro <- function(x){
  
  if(x[5]=='1'){
    
    return(x[c(2,3,3,1)])
    
  }else{
    
    return(x[c(2,4,4,1)])
    
  }
  
  
  
}

pro <- t(apply(gene_anno,1,get_pro))


pro[,1] <- paste0("chr",pro[,1])

pro <- as.data.frame(pro)

pro[,2] <- as.numeric(as.character(pro[,2]))-1000

pro[,3] <- as.numeric(as.character(pro[,3]))+1000

pro[,1] <- factor(pro[,1],levels = paste0("chr",c(1:22,"X")))

pro <- pro[order(pro[,1],pro[,2]),]


# overlaps enhancers and genes with gold-standard pairs to generate positive samples
gs$id <- 1:(dim(gs)[1])
overlapped_enh.1 <- bt.intersect(enhancer,gs[,c(1:3,6)],wa=T,wb=T)
overlapped_gene.1 <- bt.intersect(pro[,1:3],gs[,c(1:3,6)],wa=T,wb=T)
overlapped_enh.2 <- bt.intersect(enhancer,gs[,c(c(1,4,5,6))],wa=T,wb=T)
overlapped_gene.2 <- bt.intersect(pro[,1:3],gs[,c(1,4,5,6)],wa=T,wb=T)

pos_pair.1 <- merge(overlapped_enh.1,overlapped_gene.2,by='V7')
pos_pair.2 <- merge(overlapped_enh.2,overlapped_gene.1,by='V7')

pos_pair <- rbind(pos_pair.1,pos_pair.2)[,c(2:4,9:10)]


N = 3000000

enh_chr = split(1:dim(enhancer)[1],enhancer[,1])
pro_chr = split(1:dim(pro)[1],pro[,1])

sel_pair = list()
for(i in 1:23){
  
  tmp_enh_id= sample(enh_chr[[i]],N,replace = T)
  tmp_pro_id = sample(pro_chr[[i]],N,replace = T)
  tmp_id= which(abs(enhancer[tmp_enh_id,2] - pro[tmp_pro_id,2])<2e6)
  sel_pair[[i]] = data.frame(enhancer[tmp_enh_id[tmp_id],],pro[tmp_pro_id[tmp_id],])
  
}
neg_pair <- do.call(rbind,sel_pair)
neg_pair <- neg_pair[,c(1,2,3,5,6)]


# enh_name <- apply(neg_pair,1,function(x){
  
#   tmp = paste0("K562|",x[1],":",x[2],"-",x[3])
  
#   return(gsub(" ","",tmp))
  
# })

# promoter_name <- apply(neg_pair,1,function(x){
  
#   tmp = paste0("K562|",x[1],":",x[4],"-",x[5])
  
#   return(gsub(" ","",tmp))
  
# })

# neg_pair <- data.frame(neg_pair[,1],neg_pair[,2],neg_pair[,3],enh_name,neg_pair[,1],neg_pair[,4],neg_pair[,5],promoter_name)


# names(neg_pair) <- c('enhancer_chrom', 'enhancer_start', 'enhancer_end', 'enhancer_name',
#                      'promoter_chrom', 'promoter_start', 'promoter_end', 'promoter_name')


#control distance
d_neg = abs(neg_pair[,4]-neg_pair[,2])
d_pos = abs(pos_pair[,4]-pos_pair[,2])

bin_count_neg = floor(d_neg/1e4)

dist_control_neg = list()

for(i in 0:199){
  
  pos_count = sum(floor(d_pos/1e4) == i)
  candidate = which(bin_count_neg == i)
  sel_idx = sample(candidate,pos_count*3)
  dist_control_neg[[i+1]] = neg_pair[sel_idx,]
  
}

dist_neg_pair = do.call(rbind,dist_control_neg)

match <- function(x){
  
  which(data[,1]==as.character(x[1])&data[,2]<as.numeric(as.character(x[3]))&data[,3]>as.numeric(as.character(x[2])))
  
}



domain <- read.table(domain_path)
colnames(pos_pair) <- colnames(dist_neg_pair) <- 1:5
pairs <- rbind(pos_pair,dist_neg_pair)
pairs$id <- c(rep(1,dim(pos_pair)[1]),rep(-1,dim(dist_neg_pair)[1]) )
get_domain <- function(x){
  
  
  id <- which(domain[,1]==as.character(x[1])&domain[,2] < as.numeric(x[3]) & domain[,3] > as.numeric(x[2]))
  
  if(length(id)>0){
    
    return(id)
    
  }else{
    
    return(0)
    
  }
  
}

enh_domain <- apply(pairs[,1:3],1,get_domain)

pro_domain <- apply(pairs[,c(1,4,5)],1,get_domain)

a <- subset(pairs,enh_domain==pro_domain&enh_domain!=0)

write.table(subset(a,a$id==1),paste0(train_pairs_out,"/ProTECT_training_pairs_pos.txt"),row.names = F,col.names = F,sep="\t",quote=F)
write.table(subset(a,a$id==-1),paste0(train_pairs_out,"/ProTECT_training_pairs_neg.txt"),row.names = F,col.names = F,sep="\t",quote=F)


# reformat for TargetFinder
enh_name <- apply(a,1,function(x){
  
  tmp = paste0("K562|",x[1],":",x[2],"-",x[3])
  
  return(gsub(" ","",tmp))
  
})

promoter_name <- apply(a,1,function(x){
  
  tmp = paste0("K562|",x[1],":",x[4],"-",x[5])
  
  return(gsub(" ","",tmp))
  
})

tf_pairs <- data.frame(a[,1],a[,2],a[,3],enh_name,a[,1],a[,4],a[,5],promoter_name)


names(tf_pairs) <- c('enhancer_chrom', 'enhancer_start', 'enhancer_end', 'enhancer_name',
                     'promoter_chrom', 'promoter_start', 'promoter_end', 'promoter_name')


write.table(a,paste0(train_pairs_out,"/TargetFinder_training_pairs.txt"),row.names = F,col.names = F,sep="\t",quote=F)
write.table(a$id,paste0(train_pairs_out,"/TargetFinder_training_pairs_labels.txt"),row.names = F,col.names = F,sep="\t",quote=F)
