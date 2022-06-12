library(matrixStats)
require(matrixcalc)
library(pracma)
#      This Code done by : 

#   Amen A. Khabeer(1), Raffaele Vitale(2), David Morales-Jimenez(3),  Carolina Gómez-Llorente(3), José Camacho (c) (3)
#(1) University of Technology, Iraq,( 2) University of Lille, France 
#(3) University of Granada, Spain
#Amen.a.khabeer@uotechnology.edu.iq, rvitale86@gmail.com, dmorales@ugr.es, gomezll@ugr.es, 
#josecamacho@ugr.es

# Exampl of Running this code
# after import your data (  X as input and f as class)
# call the function 
#paranovaVS (X, f)
# and run the code 




# create the struct VASCA that has store all vlaues insde it

setClass("paranovao", representation(factors = "list", interactions = "list",data="matrix", design="matrix" ,n_factors="numeric",n_levels="numeric",n_interactions="numeric",effects="matrix",residuals="matrix",ord_factors="matrix",p="matrix"))
Vasca<-new("paranovao")
X=X1
f=F1
f1=f
paranovaVS <- function(X,f)
{
  N <- nrow(X) 
  M <- ncol(X)
  
  interactions <-vector()
  center <- 2
  n_perm <-1000
  ts <- 0
  
  
  size_data<-dim(X)                # size of the data matrix
  n_levels<-max(f1)                    # number of levels / factor
  n_interactions<-as.integer(nrow(interactions))      # number of interactions
  n_factors<-ncol(f1)                # number of factors
  factors<-c(1:n_factors)
  X_raw<-X1
  X_level_means<- as.list(nrow(n_factors))# level_means per factor
  X_interaction_means <-as.list((n_interactions))# cell_means for interactions
  SSQ_factors <- matrix(0,as.numeric(n_factors),size_data[2])# sum of squares for factors
  f_factors =array(0,dim=c(n_perm+1,as.numeric(n_factors),as.numeric(size_data[2]) ))         # F-value to make significance tests
  
  if(length(n_interactions)==0)
  {
    SSQ_interactions <-matrix(0,1,size_data[2])# sum of squares for interactions
    f_interactions = array(0,dim=c(n_perm+1,1,size_data[2])) # F-value to make significance tests
    
  }else
  {
    SSQ_interactions <-matrix(0,n_interactions,size_data[2])# sum of squares for interactions
    f_interactions = array(0,dim=c(n_perm+1,n_interactions,size_data[2])) # F-value to make significance tests
  }
  
  # center/standardize the data
  if (center == 1){
    Mean = (matrix(1,size_data[1],1) %*% t(as.matrix(colMeans(X_raw))))        # Overall mean
    X = (X_raw - Mean)                              # Center
    SSQ_mean = colSums(Mean^2)                    # SSQ overall means
    SSQ_X = colSums(X_raw^2)                      # Sum of squares data matrix
  }else if ( center == 2){
    stdm <- colSds(as.matrix(X_raw)) 
    stdm[stdm == 0] <- 1  
    
    Mean_std <- (matrix(1,size_data[1],1) %*% t(colMeans(as.matrix(X_raw))))/(matrix(1,size_data[1],1) %*% t(matrix(stdm)))
    X_std <- X_raw/(matrix(1,size_data[1],1) %*% t(as.matrix(stdm)))
    X <- (X_std - Mean_std)                          # Standardize
    SSQ_mean <- colSums(Mean_std^2)                # SSQ overall means
    SSQ_X <- colSums(X_std^2)
    
  }
  X_residuals<-X # initial matrix with residuals
  
  Vasa<-new("paranovao")
  
  
  
  Vasa@interactions=list(as.numeric(interactions),1)
  Vasa@data=matrix(X,nrow=nrow(X),ncol=ncol(X))
  Vasa@design=as.matrix(f)
  Vasa@n_factors=n_factors
  Vasa@n_levels=n_levels
  Vasa@n_interactions=as.numeric(n_interactions)
  
  # Collect level means for the factors indicated in the model
  
  factorx = c(1:as.integer(factors))
  
  for(factor in 1:as.integer(factors))
  {
    X_level_means[factor] <-list(level_means(as.matrix(X), Vasa, factor) )   
    
    
    
    
    SSQ_factors[factor,] = colSums(matrix(unlist(X_level_means[factor]),ncol=size_data[2])^2)
    X_residuals =X_residuals -(matrix(unlist(X_level_means[factor]),ncol=size_data[2]))
    Vasa@factors[factor] = SSQ_factors[factor] 
  }
  
  X_residuals_afterF = X_residuals 
  
  # Interactions
  print(as.numeric(n_interactions))
  for (i in 1 :0 )
  {
    
    
    
    
    
  }
  
  
  SSQ_residuals = colSums(X_residuals^2)
  
  
  if(as.integer(ts)>0){
    
    for(factor in factorx)
    {
      f_factors[1,factor,] = as.numeric(SSQ_factors[factor,]/SSQ_X)
    }
    for (i in 1 : n_interactions){
      f_interactions[1,i, ]= as.numeric(SSQ_interactions[i, ]/SSQ_X)
    }
    
    
    
  }else
  {
    for(factor in factorx)
    {
      f_factors[1,factor,] = SSQ_factors[factor,]
      
    }
    
    
    
    
  }
  
  
  perc_effects = effect_explains(SSQ_X, SSQ_mean, SSQ_factors, SSQ_interactions, SSQ_residuals)
  
  
  Vasa@effects = perc_effects
  Vasa@residuals = as.matrix(X_residuals)
  print("Factors: variables average")
  print(colMeans(as.matrix(perc_effects[,1 + (1 : n_factors)])))
  if (length(n_interactions)>0){
    print("Interactions: variables average")
    print(colMeans(as.matrix(perc_effects[,1 + (1 : n_factors)+(1 : n_interactions)])))
    
  }
  
  print("Percentage each effect contributes to the total sum of squares")
  print("Overall means: variables average")
  print(colMeans(as.matrix(perc_effects[,1])))
  print("Residuals: variables average")
  print(colMeans(as.matrix(perc_effects[,ncol(perc_effects)])))
  
  
  
  
  
  for (j in 1 : n_perm )
  {
    
    
    perms =randperm(size_data[1],size_data[1]) # permuted data (permute whole data matrix)
    
    X_residuals = X[perms, ]
    
    for (factor in 1 : n_factors)
    {# Level means for each factor
      X_level_means[factor] = list(level_means(as.matrix(X[perms, ]), Vasa, factor))
      
      SSQ_factors[factor,] = colSums(matrix(unlist(X_level_means[factor]),ncol=size_data[2])^2)
      X_residuals =X_residuals -(matrix(unlist(X_level_means[factor]),ncol=size_data[2]))
    }
    # Permutations for interactions.
    if(length( n_interactions)>0){
      for (interaction in 1 : n_interactions){
        factor_1 = interactions[interaction,1]
        factor_2 = interactions[interaction,2]
        X_interaction_means[interaction] = matrix(0,size_data[1],size_data[2])
        for (level_factor_1 in 1 : n_levels(factor_1)){                   # levels for first factor
          for (level_factor_2 in 1 : n_levels(factor_2)){               # levels for second factor
            
            tmp = matrix(0,size_data[1],1)
            found = which((f[,factor_2] == level_factor_2) && (f[,factor_1] == level_factor_1),arr.ind = TRUE)   # find rows
            
            tmp[found] = 1
            
            m = colMeans(X_residuals_afterF[perms[found],])
            tmp[found] = 1                                     
            X_interaction_means[interaction] = X_interaction_means[interaction] + tmp*m
          }
        }
        SSQ_interactions[interaction,] = colSums( (X_interaction_means[interaction])^2)
        X_residuals = X_residuals - X_interaction_means[i]
      }
      
    }
    
    if (as.integer(ts)>0){
      for (factor in 1: factors){
        f_factors[1 + j,factor,] = SSQ_factors[factor,]/SSQ_X
      }
      for (i in 1 : n_interactions){
        f_interactions[1 + j,i,] = SSQ_interactions[i,]/SSQ_X
      }
    }
    else
    {
      for (factor in 1: factors){
        f_factors[1 + j,factor,] = SSQ_factors[factor,]
      }
      if( length(n_interactions)>0){{for (i in 1 : n_interactions){
        f_interactions[1 + j,i,]= SSQ_interactions[i,]}
      }}
      
      
      
      
    }
    
  }
  
  ord_factors =  matrix(0,nrow=n_factors,ncol=size_data[2])
  #ord_interactions =   matrix(0,nrow=n_interactions,ncol=size_data[2])
  
  
  for (factor in 1 :as.integer( n_factors))
  {
    ord_factors = f_factors[order(f_factors[1 ,factor,],decreasing = TRUE)]
    
    for (j in 1 : n_perm){
      f_factors[1 + j,factor,] = f_factors[order(f_factors[1+j ,factor,],decreasing = TRUE)]
    }
  }
  
  
  Vasa@ord_factors = ord_factors
  
  
  
  
  
  
  
  
  
  
  return(Vasa)
  
  
  
  
}




#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&7

effect_explains <- function(ssq_X, ssq_mean, ssq_factors,ssq_interactions, ssq_residuals)                                                 
{
  #% Function to calculate the percentage of variance explained by
  #% each effect.
  
  #% Input
  #% sums of squares of the data matrix, the factors and the
  #% interactions
  
  #% Output
  #% vector with percentages explained variance for each effect.
  ssq_factors1=t(c(ssq_factors))
  ssq_interactions1=t(c(ssq_interactions))
  ssq1=matrix(0,length(ssq_mean),4)
  for(i in 1:length(ssq_mean))
  {
    ssq1[i,1]=(ssq_mean[i])
    ssq1[i,2]=(ssq_factors1[i])
    ssq1[i,3]=(ssq_interactions1[i])
    ssq1[i,4]=(ssq_residuals[i])
    
  }
  ssq =ssq1
  
  perc_explained_effect = 100*(ssq/( (( as.matrix(ssq_X)) %*% matrix(1,1,ncol(ssq)))  ))
  
  return(perc_explained_effect)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
level_means <- function(Y,Input,factor)
{
  
  # Calculate level means for factor in matrix Y
  size_Y = dim(Y)     
  M= matrix(0,size_Y[1],size_Y[2]) 
  
  for(level in c(1:2))
  {
    tmp = matrix(0,size_Y[1],1)
    
    
    
    found<-which(Input@design[,factor]== level,arr.ind = TRUE)
    
    m = colMeans(Y[found,])
    tmp[found] = 1                                      
    M = M+ tmp%*%m
    
  }
  
  return(M)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



































