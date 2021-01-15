enh_path = "" # path to the enhancer coordinates
gene_act_mat_path = "" # path to gene activity matrix
gene_name_path = "" # path to the gene names
gene_annotation_path = "" # path to the gene annotation files
pos_path = "" # path to the positive enhancer-gene pairs
domain_path = "" # path to the TAD annotation
train_pairs_out = "" # output path of the processed training pairs


enhancer <- read.table(enh_path)

gene_act_mat <- read.table(gene_act_mat_path)
gene_name <- read.table(gene_name_path)

gene_anno <- read.table(gene_annotation_path,colClasses = "character")

#for GM12878 check 48 column
act_gene <- subset(gene_name,gene_act_mat[,54]>0)

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


# enh_ind <- sample(1:dim(enhancer)[1],55000000,replace = T)
# pro_ind <- sample(1:dim(pro)[1],55000000,replace=T)
# 
# id = enhancer[enh_ind,1] == pro[pro_ind,1]
# 
# neg_pair <- data.frame(enhancer[enh_ind[which(a==1)],],pro[pro_ind[which(a==1)],])
# 
# neg_pair <- subset(neg_pair,neg_pair[,1]==neg_pair[,4])
# 
# neg_pair <- neg_pair[,c(1,2,3,5,6)]
# 
# d <- abs(neg_pair[,4]-neg_pair[,2])
# 
# neg_pair <- subset(neg_pair,d<2.1e6)

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


enh_name <- apply(neg_pair,1,function(x){
  
  tmp = paste0("K562|",x[1],":",x[2],"-",x[3])
  
  return(gsub(" ","",tmp))
  
})

promoter_name <- apply(neg_pair,1,function(x){
  
  tmp = paste0("K562|",x[1],":",x[4],"-",x[5])
  
  return(gsub(" ","",tmp))
  
})

neg_pair <- data.frame(neg_pair[,1],neg_pair[,2],neg_pair[,3],enh_name,neg_pair[,1],neg_pair[,4],neg_pair[,5],promoter_name)


names(neg_pair) <- c('enhancer_chrom', 'enhancer_start', 'enhancer_end', 'enhancer_name',
                     'promoter_chrom', 'promoter_start', 'promoter_end', 'promoter_name')

#control distance
pos_pair <- read.csv(pos_path)

d_neg = abs(neg_pair[,6]-neg_pair[,2])
d_pos = abs(pos_pair[,7]-pos_pair[,3])

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

pairs <- rbind(pos_pair,dist_neg_pair)
get_domain <- function(x){
  
  
  id <- which(domain[,1]==as.character(x[1])&domain[,2] < as.numeric(x[3]) & domain[,3] > as.numeric(x[2]))
  
  if(length(id)>0){
    
    return(id)
    
  }else{
    
    return(0)
    
  }
  
}

enh_domain <- apply(pairs[,1:3],1,get_domain)

pro_domain <- apply(pairs[,c(5,6,7)],1,get_domain)

a <- subset(pairs,enh_domain==pro_domain&enh_domain!=0)

write.table(a,train+pairs_out,row.names = F,col.names = F,row.names=F,sep="\t",quote=F)
