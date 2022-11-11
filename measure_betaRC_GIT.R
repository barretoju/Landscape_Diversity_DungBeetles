## Script dedicated to calculate Beta Raup-Crick metrics
# Updated during BIOLCONS revision, nov 2022

library(tidyverse)
library(tidyselect)
library(purrr)
library(here)

# Load functions (BELOW)
source("/Volumes/GoogleDrive-102171931027616483972/Meu Drive/PhD/Project/analysis/CH1 - Diversity and forest loss/before_may2020/explora/beta_raupcrick/raup_crick_modif.R")

# função rc_qrand modificado de Puttker et al. 2014
source("/Volumes/GoogleDrive-102171931027616483972/Meu Drive/PhD/Project/analysis/CH1 - Diversity and forest loss/before_may2020/explora/beta_raupcrick/rc_qrand_modif.R")

# dados:
load(here("datasets", "beetles.Rdata"))

# siteXspp data (all species)
b.com

## vector for each set of species
fs.spp <- hab.data %>% filter(Classification=="FS") %>% pull(spp)
nfs.spp <-  hab.data %>% filter(Classification=="NFS") %>% pull(spp)

# sitXsp datasets for each
b.com.fs <- b.com %>% dplyr::select(Landscape, Point, all_of(fs.spp))
b.com.nfs <- b.com %>% dplyr::select(Landscape, Point, all_of(nfs.spp))


##### RAUP-CRICK ocorrência modif Chase #####

##### ALL species ####

# Pool regional - SAD todas as paisagens
reg.pool <- colSums(b.com[,3:dim(b.com)[2]]) ## making regional pool to set as occur argument
## loop to run raup crick per landscape
Landscape <- unique(b.com$Landscape)

rc.df <- data.frame(Landscape=NA, beta.rc=NA)
beta.lands <- list()

for (i in 1:length(Landscape)){
  b.com2 <- b.com[b.com$Landscape == Landscape[i],]

  beta <- raup_crick(b.com2[,3:dim(b.com2)[2]], plot_names_in_col1 = F, occur = reg.pool)
  rc.df[i,1] <- Landscape[i]
  rc.df[i,2] <- mean(beta)
  beta.lands[[i]] <- beta
}

(rc.df <- rc.df %>% rename(brc.all.raw= beta.rc))

# rmNA: remove 0 from calculations
beta.lands3 <- map(beta.lands, as.matrix)
fufa2 <- function(x)  ifelse(x==0, NA, x)
beta.lands3 <- modify_depth(beta.lands3, 1, fufa2)  # funcao que sub 0 por NA
beta.lands3 <- modify_depth(beta.lands3, 1, as.dist)

rc.df$brc.allNA <- map_dbl(beta.lands3, mean, na.rm=T)

##### FS #####

reg.pool.fs <- colSums(b.com.fs[,3:dim(b.com.fs)[2]]) ## regional pool of FS spp
rc.df.fs <- data.frame(Landscape=NA, beta.rc=NA)
beta.lands.fs <- list()

for (i in 1:length(Landscape)){
  b.com.fs2 <- b.com.fs[b.com.fs$Landscape == Landscape[i],]
  beta <- raup_crick(b.com.fs2[,3:dim(b.com.fs2)[2]], plot_names_in_col1 = F, occur = reg.pool.fs)
  rc.df.fs[i,1] <- Landscape[i]
  rc.df.fs[i,2] <- mean(beta)
  beta.lands.fs[[i]] <- beta
}

(rc.df.fs <- rc.df.fs %>% rename(brc.fs.raw = beta.rc))

# rmNA: remove 0 from calculations
ttes <- map(beta.lands.fs, as.matrix)
ttes <- modify_depth(ttes, 1, fufa2)  # funcao que sub 0 por NA
ttes <- modify_depth(ttes, 1, as.dist)

rc.df.fs$brc.fsNA <- map_dbl(ttes, mean, na.rm=T)

##### NFS #####

reg.pool.nfs <- colSums(b.com.nfs[,3:dim(b.com.nfs)[2]]) ## regional pool of nfs spp
rc.df.nfs <- data.frame(Landscape=NA, beta.rc=NA)
beta.lands.nfs <- list()

for (i in 1:length(Landscape)){
  b.com.nfs2 <- b.com.nfs[b.com.nfs$Landscape == Landscape[i],]
  beta <- raup_crick(b.com.nfs2[,3:dim(b.com.nfs2)[2]], plot_names_in_col1 = F, occur = reg.pool.nfs)
  rc.df.nfs[i,1] <- Landscape[i]
  rc.df.nfs[i,2] <- mean(beta)
  beta.lands.nfs[[i]] <- beta
}


(rc.df.nfs <- rc.df.nfs %>% rename(brc.nfs.raw = beta.rc))

# rmNA: remove 0 from calculations
ttes <- map(beta.lands.nfs, as.matrix)
ttes <- modify_depth(ttes, 1, fufa2) # funcao que sub 0 por NA
ttes <- modify_depth(ttes, 1, as.dist)
rc.df.nfs$brc.nfsNA <- map_dbl(ttes, mean, na.rm=T)

##### RAUP-CRICK ABUNDANCE modif Puttker #####

## set regional spp pool
reg.pool <- b.com %>% ungroup() %>% 
  dplyr::select(-Landscape, -Point) %>% 
  summarise_all(sum) %>% 
  gather("spp","abund",1:51)

reg.pool <- as.data.frame(reg.pool)

rc.df.abund.all <- data.frame(Landscape=NA, brc.abund.all.raw=NA)
beta.lands.abund.all <- list()

for (i in 1:length(Landscape)){
  b.com2 <- b.com[b.com$Landscape == Landscape[i],]
  beta <- rc_qrand(b.com2[,3:dim(b.com2)[2]], pool = reg.pool)
  rc.df.abund.all[i,1] <- Landscape[i]
  rc.df.abund.all[i,2] <- mean(beta, na.rm = T)
  beta.lands.abund.all[[i]] <- beta
}

# ## if sites have zero occurrence, results in NA - in which we therefore remove
# rmNA: remove 0 from calculations
ttes <- map(beta.lands.abund.all, as.matrix)
ttes <- modify_depth(ttes, 1, fufa2) # funcao que sub 0 por NA
ttes <- modify_depth(ttes, 1, as.dist)
rc.df.abund.all$brc.abund.allNA <- map_dbl(ttes, mean, na.rm= T)

#### FS ####

## set regional spp pool
reg.pool.fs <- b.com.fs %>% ungroup() %>% dplyr::select(-Landscape, -Point) %>% 
  summarise_all(sum) %>% gather("spp","abund",1:28)
reg.pool.fs <- as.data.frame(reg.pool.fs)

rc.df.abund.fs <- data.frame(Landscape=NA, brc.abund.fs.raw=NA)
beta.lands.abund.fs <- list()

for (i in 1:length(Landscape)){
  b.com.fs2 <- b.com.fs[b.com.fs$Landscape == Landscape[i],]
  b.com.fs2 <- as.matrix(b.com.fs2)
  beta <- rc_qrand(b.com.fs2[,3:dim(b.com.fs2)[2]], pool = reg.pool.fs)
  rc.df.abund.fs[i,1] <- Landscape[i]
  rc.df.abund.fs[i,2] <- mean(beta, na.rm = T)
  beta.lands.abund.fs[[i]] <- beta
}

# rmNA: remove 0 from calculations
ttes <- map(beta.lands.abund.fs, as.matrix)
ttes <- modify_depth(ttes, 1, fufa2) # funcao que sub 0 por NA
ttes <- modify_depth(ttes, 1, as.dist)
rc.df.abund.fs$brc.abund.fsNA <- map_dbl(ttes, mean, na.rm= T)

#### NFS ####

## set regional spp pool
reg.pool.nfs <- b.com.nfs %>% ungroup() %>% dplyr::select(-Landscape, -Point) %>% summarise_all(sum) %>% gather("spp","abund",1:23)
reg.pool.nfs <- as.data.frame(reg.pool.nfs)

rc.df.abund.nfs <- data.frame(Landscape=NA, brc.abund.nfs.raw=NA)
beta.lands.abund.nfs <- list()

for (i in 1:length(Landscape)){
  b.com.nfs2 <- b.com.nfs[b.com.nfs$Landscape == Landscape[i],]
  beta <- rc_qrand(b.com.nfs2[,3:dim(b.com.nfs2)[2]], pool = reg.pool.nfs)
  rc.df.abund.nfs[i,1] <- Landscape[i]
  rc.df.abund.nfs[i,2] <- mean(beta, na.rm = T)
  beta.lands.abund.nfs[[i]] <- beta
}

# rmNA: remove 0 from calculations
ttes <- map(beta.lands.abund.nfs, as.matrix)
ttes <- modify_depth(ttes, 1, fufa2) # funcao que sub 0 por NA
ttes <- modify_depth(ttes, 1, as.dist)
rc.df.abund.nfs$brc.abund.nfsNA <- map_dbl(ttes, mean, na.rm= T)

##### save data: ####
#### gather all data in one set to then save it
(rc.df <- rc.df %>% 
   left_join(rc.df.fs, by= "Landscape") %>%  # rich-based rc fs spp
   left_join(rc.df.nfs, by= "Landscape") %>% # rich-based rc nfs
   left_join(rc.df.abund.all, by= "Landscape") %>% # abund-based RC for all
   left_join(rc.df.abund.fs, by= "Landscape") %>% # abund-based fs
   left_join(rc.df.abund.nfs, by= "Landscape") %>% # abund-based nfs
   left_join(env.data, by= "Landscape"))

# write.csv(rc.df, file = here("datasets", "betaRC.csv"), row.names=FALSE)


############## FUNCTIONS FROM PUTTKER ET AL 2015 ########

# Raup-Crick modified by CHASE ET AL 2011 and BARRETO ET AL. 2020

raup_crick=function(spXsite, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=9999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE, occur= NULL){
  
  ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended.  The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).  Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  If ties are split (as we recommend) the dissimilarity (default) and similarity (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work. The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.  
  
  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model	
  
  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }
  
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  
  ##make the spXsite matrix into a pres/abs. (overwrites initial spXsite matrix):
  ceiling(spXsite/max(spXsite))->spXsite
  
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  if(is.null(occur)){   ### to set a local species pool (per landscape) ## JB & ML (set, 2019)
    occur<-apply(spXsite, MARGIN=2, FUN=sum)
  }
  
  ##NOT recommended- this is a non-trivial change to the metric:
  ##sets all species to occur with equal frequency in the null model
  ##e.g.- discards any occupancy frequency information
  if(set_all_species_equal){
    occur<-rep(1,gamma)
  }
  
  ## determine how many unique species richness values are in the dataset
  ##this is used to limit the number of null communities that have to be calculated
  alpha_levels<-sort(unique(apply(spXsite, MARGIN=1, FUN=sum))) ## ordena pela riqueza por local
  
  ##make_null:
  
  ##alpha_table is used as a lookup to help identify which null distribution to use for the tests later.  It contains one row for each combination of alpha richness levels. 
  
  alpha_table<-data.frame(c(NA), c(NA))
  names(alpha_table)<-c("smaller_alpha", "bigger_alpha")
  col_count<-1
  
  ##null_array will hold the actual null distribution values. Each element of the array corresponds to a null distribution for each combination of alpha values.  The alpha_table is used to point to the correct null distribution- the row numbers of alpha_table correspond to the [[x]] indices of the null_array.  Later the function will find the row of alpha_table with the right combination of alpha values.  That row number is used to identify the element of null_array that contains the correct null distribution for that combination of alpha levels. 
  
  null_array<-list()
  
  ##looping over each combination of alpha levels:
  
  for(a1 in 1:length(alpha_levels)){
    for(a2 in a1:length(alpha_levels)){
      
      ##build a null distribution of the number of shared species for a pair of alpha values:
      null_shared_spp<-NULL
      for(i in 1:reps){
        
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        
        ##add alpha1 number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, alpha_levels[a1], replace=FALSE, prob=occur)]<-1
        
        
        ##same for com2:
        com2[sample(1:gamma, alpha_levels[a2], replace=FALSE, prob=occur)]<-1
        
        ##how many species are shared in common?
        null_shared_spp[i]<-sum((com1+com2)>1)
        
      }
      
      ##store null distribution, record values for alpha 1 and 2 in the alpha_table to help find the correct null distribution later:
      null_array[[col_count]]<-null_shared_spp
      
      alpha_table[col_count, which(names(alpha_table)=="smaller_alpha")]<-alpha_levels[a1]
      alpha_table[col_count, which(names(alpha_table)=="bigger_alpha")]<-alpha_levels[a2]
      
      #increment the counter for the columns of the alpha table/ elements of the null array
      col_count<-col_count+1
      
    }
    
  }
  
  ##create a new column with both alpha levels to match on:
  alpha_table$matching<-paste(alpha_table[,1], alpha_table[,2], sep="_")
  
  #####################
  ##do the test:
  
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  
  
  ##for each pair of sites (duplicates effort now to make a full matrix instead of a half one- but this part should be minimal time as compared to the null model building)
  for(i in 1:n_sites){
    for(j in 1:n_sites){
      
      ##how many species are shared between the two sites:
      n_shared_obs<-sum((spXsite[i,]+spXsite[j,])>1)
      
      ## what was the observed richness of each site?
      obs_a1<-sum(spXsite[i,])
      obs_a2<-sum(spXsite[j,])
      
      ##place these alphas into an object to match against alpha_table (sort so smaller alpha is first)
      obs_a_pair<-sort(c(obs_a1, obs_a2))
      
      ##match against the alpha table- row index identifies which element of the null array contains the correct null distribution for the observed combination of alpha values:
      null_index<-which(alpha_table$matching==paste(obs_a_pair[1], obs_a_pair[2], sep="_"))
      
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null<-sum(null_array[[null_index]]==n_shared_obs)
      
      ##how many null values are bigger than the observed value?
      num_greater_in_null<-sum(null_array[[null_index]]>n_shared_obs)
      
      rc<-(num_greater_in_null)/reps
      
      if(split_ties){
        
        rc<-((num_greater_in_null+(num_exact_matching_in_null)/2)/reps)
      }
      
      if(!classic_metric){
        
        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
        
        rc<-(rc-.5)*2
      }
      
      ## at this point rc represents an index of dissimilarity- multiply by -1 to convert to a similarity as specified in the original 1979 Raup Crick paper
      if(report_similarity & !classic_metric){
        rc<- rc*-1
      }
      
      ## the switch to similarity is done differently if the original 0 to 1 range of the metric is used:
      if(report_similarity & classic_metric){
        rc<- 1-rc
      }
      
      ##store the metric in the results matrix:
      results[i,j]<-round(rc, digits=2)
      
    }
  }
  
  if(as.distance.matrix){
    results<-as.dist(results)
  }	
  
  return(results)
}

 ######### ABUNDANCE-BASED:

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
  } else{  ## this path if regional - THAT'S OUR!!  ## JB & ML (set, 2019)
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
