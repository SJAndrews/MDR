dat <- as.data.frame(rbind(c(0,2,1,1),
                      c(1,1,1,0),
                      c(2,0,1,1),
                      c(2,0,2,0)), header = TRUE)
colnames(dat) <- c("x1", "x2", "x3", "disease")


#<-read.table("C:/SampleData.txt", header = TRUE)


CMDR<-function(data=DATA,n.repeat=100,missing=999,loci=2,alpha=0.05){

PermData<-rbind(t(apply(matrix(data[,dim(data)[2]], n.repeat, dim(data)[1],byrow=TRUE),1,sample)),data[,dim(data)[2]]) #permuted data
n_interaction<-(choose(dim(data)[2]-1,loci))                                    #setting number of interactions 
group<-data[,dim(data)[2]]                                                      #identifing disease 
RefOR<-RefRR<-Refstat<-matrix(missing, nrow=n_interaction,ncol=n.repeat+1)      #empty matrix for permuted results 
pOR<-pRR<-pstat<-c(missing,n_interaction)                                       #empty matrix for final results 
level<-rep(missing,dim(data)[2]-1)                                              #empty vecotr of number of genotypes in SNPs
for (i in 1: (dim(data)[2]-1)){level[i]<-length(unique(data[,i]))}              #number of genotpyes in each snp into level
maxlevel<-sort(level,decreasing=TRUE)[loci-1]*sort(level,decreasing=TRUE)[loci] #maximum genotypic combinations 
interaction<-matrix(missing, nrow=n_interaction,ncol=loci)                      #empty matrix for interactions 
tj1<-tj2<-tRF<-matrix(missing,n_interaction, maxlevel)                          #tj1 empty matrix for genotype of SNP1 in genotopic combination
                                                                                #tj2 empty matrix for genotype of SNP1 in genotopic combination
                                                                                #tRF empty matrix for specifing HR/LR genotypes
for (i2 in 1:(n.repeat+1))
{
y<-PermData[i2,]                                                     #Data for permutation (cases vs controls)

stat <-OR<-RR<-pstat<-pOR<-pRR<-rep(NULL,n_interaction)              #empty vector
test<-control<-n_size<-matrix(0, nrow=n_interaction,ncol=maxlevel)   #empty matrix for n cases/controls in genotypic combination
lowrisk<-highrisk<-matrix(999, nrow=n_interaction,ncol=2)            #empty matrix for predisposing risk table
colnames(lowrisk) = c('case', 'control'); colnames(highrisk) = c('case', 'control') 
RF<-matrix(0,nrow=dim(data)[1],ncol=n_interaction)                   #empty matrix

i=0                                                                  #??????

for (j1 in 1:(dim(data)[2]-loci))                                    #loop for SNP 1 in interaction
{
for (j2 in (j1+1):(dim(data)[2]-loci+1))                             #loop for SNP 2 in interaction
{

i=i+1                                                                #??????


interaction[i,]<-c(j1,j2)                                            #spcifing the interaction between SNPs

for (k in 1:dim(unique(data[,c(j1,j2)]))[1])                         #loop for genotypic combinations in SNPxSNP matrix  
{
#number of controls in genotypic combination
control[i,k] <-sum(y[data[,j1]==unique(data[,c(j1,j2)])[k,1] & data[,j2]==unique(data[,c(j1,j2)])[k,2]]==0) 
#number of cases in genotypic combination
test[i,k]    <-sum(y[data[,j1]==unique(data[,c(j1,j2)])[k,1] & data[,j2]==unique(data[,c(j1,j2)])[k,2]]==1)
tj1[i,k]<-unique(data[,c(j1,j2)])[k,1] #SNP 1: genotype in genotypic combination
tj2[i,k]<-unique(data[,c(j1,j2)])[k,2] #SNP 2: genotype in genotypic combination
# pi0 = equation, pg 3
#if percentage of cases in genotypic combination is greater then then percentage in total cohort tRF = 1
if (test[i,k]/max((control[i,k]+test[i,k]),1)>=  sum(y)/ length(y))   
  {RF[,i][data[,j1]==unique(data[,c(j1,j2)])[k,1] & data[,j2]==unique(data[,c(j1,j2)])[k,2]]<-1
                                                                     tRF[i,k]<-1
                                                                     }
#if percentage of cases in genotypic combination is greater then then percentage in total cohort tRF = 0
if (test[i,k]/max((control[i,k]+test[i,k]),1)< sum(y)/ length(y))   
  {RF[,i][data[,j1]==unique(data[,c(j1,j2)])[k,1] & data[,j2]==unique(data[,c(j1,j2)])[k,2]]<-0
                                                                     tRF[i,k]<-0
                                                                     }

}

n_size[i,]<-control[i,]+test[i,]                                     #n individules in genotypic combination

#predisposing risk table, pg 4
lowrisk[i,]<-c(sum(test[i,][(test[i,]/n_size[i,])<sum(test[i,])/sum(n_size[i,])&n_size[i,]>0]), sum(control[i,][(test[i,]/n_size[i,])<sum(test[i,])/sum(n_size[i,])&n_size[i,]>0]))
highrisk[i,]<-c(sum(test[i,][(test[i,]/n_size[i,])>=sum(test[i,])/sum(n_size[i,])&n_size[i,]>0]), sum(control[i,][(test[i,]/n_size[i,])>=sum(test[i,])/sum(n_size[i,])&n_size[i,]>0]))

OR[i]<-highrisk[i,1]*lowrisk[i,2]/highrisk[i,2]/lowrisk[i,1]                    #pOR, pg3 eq1 
RR[i]<-highrisk[i,1]/sum(highrisk[i,])/lowrisk[i,1]*sum(lowrisk[i,])            #pRR, pg4 eq2
stat[i]<-chisq.test(rbind(lowrisk[i,],highrisk[i,]),correct=FALSE)$statistic    #pChi pg5 eq3
}
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































