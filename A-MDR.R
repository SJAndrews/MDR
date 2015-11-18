##-----------------------------------------##
# Functions
##-----------------------------------------##

#classify genotypic combinations as high or low risk
#pi0, equation pg 3
pi0.func <- function(x){
  if (x[1]/(x[1] + x[2]) >=  sum(y)/ length(y)){x[3] <- 1}
  else {x[3] <- 0}
}

#counts the number of matches of an individual data vector with the rows of a target matrix, and is called by the function mdr
#from package MDR
compare <- function (mat, vec, k) 
{
  b <- 1
  match <- 1:dim(mat)[1]
  while (b <= (k + 1)) {
    match <- match[mat[match, b] == as.numeric(vec[b])]
    b <- b + 1
  }
  return(length(match))
}


##-----------------------------------------##

dat <- as.data.frame(rbind(c(0,2,1,1),
                           c(1,1,1,0),
                           c(2,0,1,1),
                           c(2,0,2,0)), header = TRUE)
colnames(dat) <- c("x1", "x2", "x3", "disease")

#rm(PermData, n_interaction, group, RefOR, RefRR, Refstat, pOR, pRR, pstat, level, maxlevel, interaction, tj1, tj2, tRF)
#<-read.table("C:/SampleData.txt", header = TRUE)
data <- subset(mdr1, select = c(SNP.1, SNP.2, SNP.3, SNP.4, SNP.5, Response))
n.repeat <- 10
missing <- 999
loci <- 1
alpha <- 0.05
i2 <- 1
CMDR<-function(data=DATA,n.repeat=100,missing=999,loci=2,alpha=0.05, genotype = c(0,1,2)){
  
  PermData<-rbind(t(apply(matrix(data[,dim(data)[2]], n.repeat, dim(data)[1],byrow=TRUE),1,sample)),data[,dim(data)[2]]) #permuted data
  interaction <- t(combn(dim(data)[2]-1, loci))                                      #empty matrix for interactions 
  n_interaction <- dim(interaction)[1]                                               #setting number of snpxsnp interactions 
  group<-data[,dim(data)[2]]                                                         #identifing disease 
  RefOR<-RefRR<-Refstat<-matrix(missing, nrow=n_interaction,ncol=n.repeat+1)         #empty matrix for permuted results 
  pOR<-pRR<-pstat<-c(missing,n_interaction)                                          #empty matrix for final results 
  g <- length(genotype)     								                                         #number of genotypes
  geno <- list(genotype)	
  geno.comb <- expand.grid(rep(geno, loci))                                          #possible genotypic combinations
  hr.lr <- matrix(missing, nrow = n_interaction, ncol = (1 + g^loci)) #hr and lr genotypic combinations for interactions
  hr.lr[,1] <- apply(interaction, 1, paste, collapse = "-")                          #potentional interactions 
  colnames(hr.lr) <- c("interaction", apply(geno.comb, 1, paste, collapse = "-"))       #names of genotypic combinations 

  #loop for running through permuted datasets
  for (i2 in 1:(n.repeat+1))
  {
    y<-PermData[i2,]                                                     #Data for permutation (cases vs controls)
    
    stat <-OR<-RR<-pstat<-pOR<-pRR<-rep(NULL,n_interaction)              #empty vector
    test<-control<-n_size<-matrix(0, nrow=n_interaction,ncol=maxlevel)   #empty matrix for n cases/controls in genotypic combination
    lowrisk<-highrisk<-matrix(999, nrow=n_interaction,ncol=2)            #empty matrix for predisposing risk table
    colnames(lowrisk) = c('case', 'control'); colnames(highrisk) = c('case', 'control') 
    RF<-matrix(0,nrow=dim(data)[1],ncol=n_interaction)                   #empty matrix
    
    case <- cbind(rep(1, g^loci), expand.grid(rep(geno, loci)))  	#genotype cases 
    ctrl <- cbind(rep(0, g^loci), expand.grid(rep(geno, loci)))		#genotype controls
    counts <- matrix(0, dim(case)[1], 3); colnames(counts) <- c('case', 'ctrl', 'ratio')	#k-way genotype combinations - ratiol of cases/controls
    
    for (j1 in 1:n_interaction))                                    
    y<-PermData[i2,]                                                     #Data for permutation (cases vs controls)
    
    stat <-OR<-RR<-pstat<-pOR<-pRR<-rep(NULL,n_interaction)              #empty vector
    lowrisk<-highrisk<-matrix(999, nrow=n_interaction,ncol=2)            #empty matrix for predisposing risk table
    colnames(lowrisk) = c('case', 'control'); colnames(highrisk) = c('case', 'control') 
    
    case <- cbind(rep(1, g^loci), expand.grid(rep(geno, loci)))  	#genotype cases 
    ctrl <- cbind(rep(0, g^loci), expand.grid(rep(geno, loci)))		#genotype controls
    counts <- matrix(0, dim(case)[1], 3); colnames(counts) <- c('case', 'ctrl', 'ratio')	#k-way genotype combinations - ratiol of cases/controls
    
    #loop for running through GxG interactions to determing high/low risk genotypes
    for (j1 in 1:n_interaction)                                    
    {
      model <- interaction[j1, ]                                          #spcifing the interaction between SNPs
      
      part <- data[,c(model)]
      part <- cbind(y, part)
      
      counts[, 1] <- apply(case, 1, compare, mat = part, k = loci) #number of cases for k-kway interaction
      counts[, 2] <- apply(ctrl, 1, compare, mat = part, k = loci) #number of controls for k-way interactions 
      counts[, 3] <- apply(counts, 1, pi0.func)        		           #ratio of cases to controls for combination
            
      hr.lr[j1, 2:ncol(hr.lr)] <- counts[,3]
      
      highrisk[j1,] <- c(sum(counts[counts[,3] == 1, 1]), sum(counts[counts[,3] == 1, 2])) #sum highrisk cases/controls
      lowrisk[j1,] <- c(sum(counts[counts[,3] == 0, 1]), sum(counts[counts[,3] == 0, 2]))  #sum lowrisk cases/controls
      
      OR[j1]<-highrisk[j1,1]*lowrisk[j1,2]/highrisk[j1,2]/lowrisk[j1,1]                    #pOR, pg3 eq1 
      RR[j1]<-highrisk[j1,1]/sum(highrisk[j1,])/lowrisk[j1,1]*sum(lowrisk[j1,])            #pRR, pg4 eq2
      stat[j1]<-chisq.test(rbind(lowrisk[j1,],highrisk[j1,]),correct=FALSE)$statistic    #pChi pg5 eq3    
    }

  RefOR[,i2]<-OR                                                                  #All and permuted and final pOR 
  RefRR[,i2]<-RR                                                                  #All and permuted and final pRR
  Refstat[,i2]<-stat                                                              #All and permuted and final pChi 
  }

  
  for (i3 in 1:length(stat))
  {
    pstat[i3]<-sum(Refstat[i3,1:n.repeat]>stat[i3])/n.repeat                        #pvalue for pChi
    pOR[i3]<-sum(RefOR[i3,1:n.repeat]>OR[i3])/n.repeat                              #pvalue for pOR 
    pRR[i3]<-sum(RefRR[i3,1:n.repeat]>RR[i3])/n.repeat                              #pvalue for PRR
  }
  #return(list(data,interaction,OR,pOR,RR,pRR,stat,pstat,tj1,tj2,tRF))
  out <- list(data,interaction,OR,pOR,RR,pRR,stat,pstat,tj1,tj2,tRF)
  names(out) <- c("data","interaction","OR","pOR","RR","pRR","stat","pstat","tj1","tj2","tRF")
  out
  #return(stat)
}

CMDR_output<-CMDR(data=dat)
CMDR_output

































