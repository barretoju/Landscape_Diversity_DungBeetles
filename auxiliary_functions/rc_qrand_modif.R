## Beta Raup-Crick function with abundance data
## Modified from Puttker et al. 2014, which modified from Chase et al. 2011

# spXsite -> site (rows) by species (columns) abundance matrix
# reps: number of repetitions of random draws for each of the pairs of sites
# classic_metric: if TRUE, RC will be calculated as a value from 0-1, if FALSE from -1 to 1

# JB & ML (set, 2019) we added "pool" argument to include the regional pool

rc_qrand<- function(spXsite, reps=999, classic_metric=FALSE, pool= NULL) {
  
  # remove species not found in the focal landscape
  spXsite <- spXsite[, colSums(spXsite)>0]
  # to ensure spXsite is a data frame not a tibble!!
  spXsite <- as.data.frame(spXsite) 
  
  # define number of sites and species in the lanscape
  num_sites<-dim(spXsite)[1]
  num_sps <- dim(spXsite)[2]
  
  # pool:
  if(is.null(pool)){   ## this path if not regional pool (JB & ML (set, 2019))
    # number of species
    gamma<-ncol(spXsite)
    # get overall abundance for each species (colsums) to create pool of species
    abund<-colSums(spXsite)
  } else{  ## this path if regional - ## JB & ML (set, 2019)
    gamma<-nrow(pool)
    abund<-pool$abund
  }
  
  
  #################################################################
  ## transform matrix in relative value of site total abundance ##
  #################################################################
  
  # get row sums
  abund_row<-apply(spXsite, MARGIN=1, FUN=sum)
  
  # create new matrix site*species
  spXsite_trans<-matrix(NA, nrow=num_sites, ncol=num_sps)
  
  # loop over sites and species division by rowsums
  for (i in 1:num_sites){
    for (j in 1:num_sps){
      spXsite_trans[i,j]<-spXsite[i,j]/abund_row[i]
    }
  }
  
  
  #####################################################################################
  # create a pool following abundance distribution of real total pool of species #####
  ####################################################################################
  pool<-c()
  for (i in 1:gamma) {
    spec<-rep(i, abund[i])
    
    pool<-append(pool, spec)
  }
  
  ## create triangular matrix for storing RC values in the end
  RC_matrix<-matrix(data=NA, nrow=num_sites, ncol=num_sites)
  
  
  #############################################################
  ## Calculation of obs and exp values abundance shared #######
  #############################################################
  
  
  ## loop over sites
  for(i in 1:(num_sites-1))
  {
    
    ## calculate overall abundance of site i for drawing randomly afterwards
    abund1<-sum(spXsite[i,])
    
    
    ## loop again over sites
    for(j in (i+1):num_sites)
    {
      
      ## calculate overall abudance of site j for drawing randomly afterwards
      abund2<-sum(spXsite[j,])
      
      ######################################################################################
      ## calculate transformed abundance shared of original matrix between site i and j:
      ## create matrix with abund_min values
      ab_min <- matrix(data=NA, nrow=1, ncol=num_sps)
      
      ## loop over species 
      for(k in 1:num_sps) {
        abund_min<-min(c(spXsite_trans[i,k], spXsite_trans[j,k]))
        ab_min[,k]<-abund_min
        
      }
      
      ## sum and save value in matrix first value
      obs<-sum(ab_min)
      
      #################################
      ## calculate expected values
      
      ## create matrix with 1 lines and reps columns for storing the expected values 
      ab_sh<-matrix(data=NA, nrow=1, ncol=reps)
      
      ## loop over reps
      for(l in 1:reps)
      {
        
        ## create table for both random communities
        com<-matrix(NA, 2, gamma)
        
        
        ## sample abund1 individuals randomly from pool
        samp1<-sample(pool, abund1, replace=FALSE)
        
        ## Define numbers of individuals per species in samp1
        ## loop over gamma
        for(m in 1:gamma) {
          ## detect frequency of numbers (i.e. species drawn randomly) 
          num<-sum(samp1==m)
          
          ## put in table
          com[1,m]<-num
          
        }
        
        ## sample abund2 individuals randomly from pool
        samp2<-sample(pool, abund2, replace=FALSE)
        
        ## Define numbers of individuals per species in samp2
        ## loop over gamma
        for(m in 1:gamma) {
          ## detect frequency of numbers (i.e. species drawn randomly) 
          num<-sum(samp2==m)
          
          ## put in table
          com[2,m]<-num
          
        }
        
        ###################################################
        ## tranform expected abundances in relative values
        
        ## create com_trans
        com_trans<-matrix(NA, 2, gamma)
        
        ## transfomr com by deviding with site total abundances
        for (o in 1:gamma)
        {
          com_trans[1,o]<-com[1,o]/abund1
          com_trans[2,o]<-com[2,o]/abund2
        }
        
        ## define transformed abundance shared between two random communities
        ## create matrix with rand_abund_min values
        rand_ab_min<-matrix(data=NA, nrow=1, ncol=gamma)
        
        ## loop over species 
        for(n in 1:gamma) {
          rand_abund_min<-min(c(com_trans[1,n], com_trans[2,n]))
          rand_ab_min[,n]<-rand_abund_min
          
        }
        
        ## sum and save value in table
        ab_sh[,l]<-sum(rand_ab_min)
      }
      
      ##################################
      #### calculate RC
      
      ## count the times the value is larger or same than observed
      ## number of exp equal to obs
      num_exact_matching_in_null<-sum(obs==ab_sh[,1:reps])
      
      ## number of exp greater than obs
      num_greater_in_null<-sum(obs<ab_sh[,1:reps])
      
      ## calculate probability of observing original transformed shared abundance or larger shared abundance given random draws from the species pool    
      rc<-((num_greater_in_null+(num_exact_matching_in_null)/2)/reps)
      
      
      ## option of classic metric
      if(!classic_metric)
      {
        
        ## standardization of metric to range from -1 to 1 instead of 0 to 1
        rc<-(rc-.5)*2
      }
      
      ## store the metric in the results matrix:
      RC_matrix[j,i]<-rc
    } 
  }
  # Take out NA?s 
  RC_matrix=as.dist(RC_matrix)
  
  # return result
  RC_matrix
}
