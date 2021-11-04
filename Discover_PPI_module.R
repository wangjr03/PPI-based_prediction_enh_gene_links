
mycluster <- function(hc,n,cluster_max){
  
  rmr_step <- c()
  
  membership <- list()
  
  merge_step_1 <- max(c(which(hc$merge[,1] < 0),which(hc$merge[,2] < 0)))
  
  count <- 0
  re_c <- c()
  for(i in 1:dim(hc$merge)[1]){
    
    tmp = hc$merge[i,]
    
    if(tmp[1] <0 & tmp[2] <0){
      
      membership[[i]] <- c(abs(tmp[1]),abs(tmp[2]))
      count = count + 1
    }else{
      
      if(tmp[1] >0 & tmp[2] >0){
        
        check_re <- union(membership[[tmp[1]]],membership[[tmp[2]]])
        
        if(length(check_re) < cluster_max){
          
          membership[[i]] <- check_re
          membership[[tmp[1]]] <- membership[[tmp[2]]] <- logical(0)
          count = count - 1
        }else{
          rmr_step <- c(rmr_step,i)
          if(length(membership[[tmp[1]]]) < length(membership[[tmp[2]]])){
            
            membership[[i]] <- membership[[tmp[1]]]
            
            membership[[tmp[1]]] <- logical(0)
            
          }else{
            
            membership[[i]] <- membership[[tmp[2]]]
            
            membership[[tmp[2]]] <- logical(0)
            
            
          }
          
        }
        
        
        
      }else{
        
        if(tmp[1] >0 & tmp[2] < 0){
          
          membership[[i]] <- union(membership[[tmp[1]]],abs(tmp[2]))
          membership[[tmp[1]]] <-logical(0)
          
          
        }else{
          
          membership[[i]] <- union(membership[[tmp[2]]],abs(tmp[1]))
          membership[[tmp[2]]] <- logical(0)
          
        }
        
        
        
      }
      
      
      
      
    }
    re_c <- c(re_c,count)
    if(count < max(re_c)/2 & count == n){
      
      break
      
    }
    #print(count)
  }
  
  #for remaining merge, never merge two clusters but continue with smaller one
  
  
  for(i in (i+1):dim(hc$merge)[1]){
    
    tmp = hc$merge[i,]
    
    if(tmp[1] <0 & tmp[2] <0){
      
      membership[[i]] <- c(abs(tmp[1]),abs(tmp[2]))
      count = count + 1
    }else{
      
      if(tmp[1] >0 & tmp[2] >0){
        
        check_re <- union(membership[[tmp[1]]],membership[[tmp[2]]])
        
        rmr_step <- c(rmr_step,i)
        if(length(membership[[tmp[1]]]) < length(membership[[tmp[2]]])){
          
          membership[[i]] <- membership[[tmp[1]]]
          
          membership[[tmp[1]]] <- logical(0)
          
        }else{
          
          membership[[i]] <- membership[[tmp[2]]]
          
          membership[[tmp[2]]] <- logical(0)
          
          
        }
        
        
        
        
      }else{
        
        if(tmp[1] >0 & tmp[2] < 0){
          
          membership[[i]] <- union(membership[[tmp[1]]],abs(tmp[2]))
          membership[[tmp[1]]] <-logical(0)
          
          
        }else{
          
          membership[[i]] <- union(membership[[tmp[2]]],abs(tmp[1]))
          membership[[tmp[2]]] <- logical(0)
          
        }
        
        
        
      }
      
      
      
      
    }
    re_c <- c(re_c,count)
    #print(count)
  }
  
  
  #generate label
  
  label <- rep(0,length(unlist(membership)))
  
  k <- 1
  
  for(i in 1:length(membership)){
    
    if(length(membership[[i]]) == 0 ){
      next
      
    }else{
      
      label[membership[[i]]] <- k
      k <- k+1
      
    }
    
    
  }
  return(label)
  
  
}

mymodule <- function(hc,G,n){
  
  m_c <- c()
  
  for(i in seq(30,n)){
    
    m_c <- c(m_c,modularity(G,mycluster(hc,i,500)))
    
  }
  
  return(m_c)
  
}


library(igraph)
args <- commandArgs(T)
ppi_file_path = args[1]
ppi_th = as.numeric(args[2])
ppi_data = read.table(ppi_file_path)
ppi_data = subset(ppi_data,ppi_data[,3]>ppi_th)
G = graph_from_edgelist(as.matrix(ppi_data)[,1:2])
E(G)$weights <- ppi_data[,3]
G=as.undirected(G,edge.attr.comb = 'min')

A = as_adjacency_matrix(G,sparse = T)
A = as.matrix(A)
diag(A) <- 1
A = as(A,'sparseMatrix')
d = apply(A,1,sum)

n = dim(A)[1]

D = diag(n)
diag(D) <- 1/d
D = as(D,'sparseMatrix')
P = D%*%A

#t step transition prob
t = 20
library(expm)
P = as.matrix(P)
Pt = P%^%t

# #calculate dist
# D_h = D^(1/2)
# dp = D_h%*%Pt
# library(parallelDist)
# input_mat <- as.matrix(t(dp))
# 
# 
# 
# input_mat <- t(input_mat)

Gram = Pt%*%t(Pt)
g = diag(Gram)
O = rep(1,dim(Gram)[1])
R = g%*%t(O) + O%*%t(g) -2*Gram
hc = hclust(as.dist(R),method='ward.D')
print('Searching for the optimal module size ...')
mc <- mymodule(hc,G,1000)

optimal_s_size = 29+which.max(mc)
optimal_l_size = which(mc>max(mc)*0.99)
optimal_l_size = optimal_l_size[length(optimal_l_size)]+29
print(paste('Optimal S-module size:',optimal_s_size))
print(paste('Optimal L-module size:',optimal_l_size))

pdf('Network_modularity.pdf')
par(mar=c(5,5,5,2))
plot(x=seq(30,29+length(mc)),y=mc,type="l",ylab='modularity',xlab='number of clusters',cex.lab=2,cex.axis=2)
points(x=optimal_s_size,y=max(mc),col='red',pch=19)
text(x = optimal_s_size, y=0.99*max(mc),paste('optimal s-module size:',which.max(mc)),cex=1,col="red")
points(x=optimal_l_size,y=mc[optimal_l_size-29],col='blue',pch=19)
text(x = optimal_l_size, y=0.99*mc[optimal_l_size-29],paste('optimal l-module size:',optimal_l_size),cex=1,col="blue")
dev.off()

s_label = mycluster(hc,optimal_s_size,500)
l_label = mycluster(hc,optimal_l_size,500)

pdf("Network_module_detection.pdf",width=100,height=100)
plot(G,vertex.size=0.1,vertex.label.cex=0.1,vertex.color = s_label,vertex.frame.color=NA)
dev.off()

hc_module <- data.frame(TF_name = V(G)$name, s_model=s_label,l_model=l_label)
write.table(hc_module,'hc_TF_community.txt',col.names=F,row.names=F,sep='\t',quote=F)