

library(gtools)
library(parallel)
library(snow)

###########################################################################################################
##############################  STEP 0: loading the phenotype file ########################################
###########################################################################################################

Pheno.Sim = readRDS(file="Path:/Data/Pheno.Batch.Effect.Rda")


Age = Pheno.Sim$Age
HBA1C = Pheno.Sim$HBA1C
Group.1 = Pheno.Sim$Group.1
Group.2 = Pheno.Sim$Group.2


###########################################################################################################
###########################  STEP A: Select all Blocks by Group ###########################################
###########################################################################################################


#STEP 1: Selection all possible within batch groups (creating dataset of all possible blocks)
#block size = #of case samples within each batch
Data.1.G1 = data.frame(combinations(30, 10, Pheno.Sim$ID[Pheno.Sim$Group.1==1]))
Data.1.G1$G.1_YN = 1
Data.1.G1$G.2_YN = 0
Data.1.G1$Block.ID = 1:nrow(Data.1.G1)

saveRDS(Data.1.G1, file="Path:/Data/Data.1.G1.Rda")


Data.1.G2 = data.frame(combinations(30, 10, Pheno.Sim$ID[Group.1==0]))
Data.1.G2$G.1_YN = 0
Data.1.G2$G.2_YN = 1
Data.1.G2$Block.ID = 1:nrow(Data.1.G2)

saveRDS(Data.1.G2, file="Path:/Data/Data.1.G2.Rda")


#######################################################################################################
#################################### STEP B.1: Group #1 ###############################################
#######################################################################################################


Data.1.G1 = readRDS(file="Path:/Data/Data.1.G1.Rda")

Pheno.G1 = Pheno.Sim[Pheno.Sim$Group.1==1, ]
Pheno.G2 = Pheno.Sim[Pheno.Sim$Group.1==0, ]


#Balancing group and batch selection

a = 1:nrow(Data.1.G1)

B.Propensity.G1 = function(a){
  
#  #Loop step 1 - pull in the phenotype data
  
  #Selected into batch
  S.ID = c(Data.1.G1[a, 1:10])
  S.Data = Pheno.G1[which(Pheno.G1$ID %in% S.ID), ]
  S.Data$Selected = 1
  
  #Not selected into the batch
  NS.Data = Pheno.G1[Pheno.G1$ID %in% c(setdiff(Pheno.G1$ID, S.ID)), ]
  NS.Data$Selected = 0
  
  B.Data = rbind(S.Data, NS.Data)
  
  LR.1 = glm(Selected~Age+HBA1C , family = "binomial", data = B.Data)
  
  P.1 = predict(LR.1,type="response")
  
  #Selected vs Not Selected
  Av.P.Selected = mean(P.1[B.Data$Selected==1])
  Av.P.NSelected = mean(P.1[B.Data$Selected==0])
  Block.Diff = abs(Av.P.Selected-Av.P.NSelected)
  
  output =  c(Av.P.Selected, Av.P.NSelected, Block.Diff)
  
  names(output) = c("Av.P.Selected", "Av.P.NSelected", "Block.Diff")
  
  return(output)
  
}

#Start of parallel processing
numCores = detectCores()-1
cl=makeCluster(numCores,type="SOCK")
clusterExport(cl, c("Data.1.G1", "Pheno.G1", "Pheno.G2"))
Data.2.G1.0 = clusterApply(cl, a, B.Propensity.G1)
stopCluster(cl)


#Converting from list to dataframe
Data.2.G1  =  as.data.frame(t(matrix(unlist(Data.2.G1.0), nrow=length(unlist(Data.2.G1.0[1])))))
colnames(Data.2.G1) = c("Av.P.Selected", "Av.P.NSelected", "Block.Diff")
Data.2.G1$Block.ID = 1:nrow(Data.2.G1)

saveRDS(Data.2.G1, file="Path:/Data/Data.2.G1.Rda")


#######################################################################################################
#################################### STEP B.1: Group #2 ###############################################
#######################################################################################################


Data.1.G2 = readRDS(file="Path:/Data/Data.1.G2.Rda")

Pheno.G1 = Pheno.Sim[Pheno.Sim$Group.1==1, ]
Pheno.G2 = Pheno.Sim[Pheno.Sim$Group.1==0, ]

a = 1:nrow(Data.1.G2)

B.Propensity.G2 = function(a){
  
  #Loop step 1 - pull in the phenotype data
  
  #Selected into batch
  S.ID = c(Data.1.G2[a, 1:10])
  S.Data = Pheno.G2[which(Pheno.G2$ID %in% S.ID), ]
  S.Data$Selected = 1
  
  #Not selected into the batch
  NS.Data = Pheno.G2[Pheno.G2$ID %in% c(setdiff(Pheno.G2$ID, S.ID)), ]
  NS.Data$Selected = 0
  
  B.Data = rbind(S.Data, NS.Data)
  
  LR.1 = glm(Selected~Age+HBA1C, family = "binomial", data = B.Data)
  
  P.1 = predict(LR.1,type="response")
  
  #Selected vs Not Selected
  Av.P.Selected = mean(P.1[B.Data$Selected==1])
  Av.P.NSelected = mean(P.1[B.Data$Selected==0])
  Block.Diff = abs(Av.P.Selected-Av.P.NSelected)
  
  output =  c(Av.P.Selected, Av.P.NSelected, Block.Diff)
  
  names(output) = c("Av.P.Selected", "Av.P.NSelected", "Block.Diff")
  
  return(output)
  
}

#Start of parallel processing
#numCores = detectCores()-1
cl=makeCluster(numCores,type="SOCK")
clusterExport(cl, c("Data.1.G2", "Pheno.G1", "Pheno.G2"))
Data.2.G2.0 = clusterApply(cl, a, B.Propensity.G2)
stopCluster(cl)


#Converting from list to dataframe
Data.2.G2  =  as.data.frame(t(matrix(unlist(Data.2.G2.0), nrow=length(unlist(Data.2.G2.0[1])))))
colnames(Data.2.G2) = c("Av.P.Selected", "Av.P.NSelected", "Block.Diff")
Data.2.G2$Block.ID = 1:nrow(Data.2.G2)

saveRDS(Data.2.G2, file="Path:/Data/Data.2.G2.Rda")


###########################################################################################################
#####################################  STEP C: Order Blocks ###############################################
###########################################################################################################


#Group 1
Data.1.G1 = readRDS(file="C:/Users/carryp/R/Simulation Batch Effect/Data/Data.1.G1.Rda")
Data.2.G1 = readRDS(file="C:/Users/carryp/R/Simulation Batch Effect/Data/Data.2.G1.Rda")

Data.2.G1.A = cbind(Data.1.G1, Data.2.G1)
Data.2.G1.S = Data.2.G1.A[ order(Data.2.G1.A$Block.Diff), ]


#Group 2
Data.1.G2 = readRDS(file="C:/Users/carryp/R/Simulation Batch Effect/Data/Data.1.G2.Rda")
Data.2.G2 = readRDS(file="C:/Users/carryp/R/Simulation Batch Effect/Data/Data.2.G2.Rda")

Data.2.G2.A = cbind(Data.1.G2, Data.2.G2)
Data.2.G2.S = Data.2.G2.A[ order(Data.2.G2.A$Block.Diff), ]



###########################################################################################################
###########################  STEP D: Pick Best Unique Blocks ##############################################
###########################################################################################################


#For group 1
Block.1.1 = unlist(c(Data.2.G1.S[1, 1:10]))
Block.ID.1.1 = Data.2.G1.S$Block.ID[1]
Row.1 = 1
i = 2
length.1 = 0
length.2 = 0
print("PICKING 2nd BLOCK")
print(i)
while (length.1 <20) {
  Two = c(Block.1.1, unlist(c(Data.2.G1.S[i, 1:10])) )
  length.1 = length(unique(Two))
  Block.ID.1.2 = Data.2.G1.S$Block.ID[i]
  Block.1.2 = unlist(c(Data.2.G1.S[i, 1:10]))
  Row.2 = i
  i = i+1
}

print("PICKING 3rd BLOCK")
print(i)
while (length.2 <30) {
  Three = c(Two, 
            unlist(c(Data.2.G1.S[i, 1:10])))
  length.2 = length(unique(Three))
  Block.ID.1.3 = Data.2.G1.S$Block.ID[i]
  Block.1.3 = unlist(c(Data.2.G1.S[i, 1:10]))
  Row.3 = i
  i = i+1
}
print(i)

R.1.1 = cbind(Data.2.G1.S[which(Data.2.G1.S$Block.ID %in% Block.ID.1.1), 1:10],Block.ID.1.1) 
R.1.2 = cbind(Data.2.G1.S[which(Data.2.G1.S$Block.ID %in% Block.ID.1.2), 1:10],Block.ID.1.2) 
R.1.3 = cbind(Data.2.G1.S[which(Data.2.G1.S$Block.ID %in% Block.ID.1.3), 1:10],Block.ID.1.3) 


colnames(R.1.1) = c("ID.1", "ID.2", "ID.3", "ID.4", "ID.5", "ID.6","ID.7", "ID.8", "ID.9", "ID.10", "Block.ID")
colnames(R.1.2) = c("ID.1", "ID.2", "ID.3", "ID.4", "ID.5", "ID.6","ID.7", "ID.8", "ID.9", "ID.10", "Block.ID")
colnames(R.1.3) = c("ID.1", "ID.2", "ID.3", "ID.4", "ID.5", "ID.6","ID.7", "ID.8", "ID.9", "ID.10", "Block.ID")


Group.1.Data = rbind(R.1.1, R.1.2, R.1.3)


  #Check
  length(unique( c(R.1.1[1,1:10], R.1.2[1,1:10], R.1.3[1,1:10])  ))

saveRDS(Group.1.Data, file = "C:/Users/carryp/R/Simulation Batch Effect/Data/Optimal.Block.Group.1.Rda")



#For group 2
Block.2.1 = unlist(c(Data.2.G2.S[1, 1:10]))
Block.ID.2.1 = Data.2.G2.S$Block.ID[1]
Row.1 = 1
i = 2
length.1 = 0
length.2 = 0
print("PICKING 2nd BLOCK")
print(i)
while (length.1 <20) {
  Two = c(Block.2.1, unlist(c(Data.2.G2.S[i, 1:10])) )
  length.1 = length(unique(Two))
  Block.ID.2.2 = Data.2.G2.S$Block.ID[i]
  Block.2.2 = unlist(c(Data.2.G2.S[i, 1:10]))
  Row.2 = i
  i = i+1
}

print("PICKING 3rd BLOCK")
print(i)
while (length.2 <30) {
  Three = c(Two, 
            unlist(c(Data.2.G2.S[i, 1:10])))
  length.2 = length(unique(Three))
  Block.ID.2.3 = Data.2.G2.S$Block.ID[i]
  Block.2.3 = unlist(c(Data.2.G2.S[i, 1:10]))
  Row.3 = i
  i = i+1
}
print(i)

R.2.1 = cbind(Data.2.G2.S[which(Data.2.G2.S$Block.ID %in% Block.ID.2.1), 1:10],Block.ID.2.1) 
R.2.2 = cbind(Data.2.G2.S[which(Data.2.G2.S$Block.ID %in% Block.ID.2.2), 1:10],Block.ID.2.2) 
R.2.3 = cbind(Data.2.G2.S[which(Data.2.G2.S$Block.ID %in% Block.ID.2.3), 1:10],Block.ID.2.3) 


colnames(R.2.1) = c("ID.1", "ID.2", "ID.3", "ID.4", "ID.5", "ID.6","ID.7", "ID.8", "ID.9", "ID.10", "Block.ID")
colnames(R.2.2) = c("ID.1", "ID.2", "ID.3", "ID.4", "ID.5", "ID.6","ID.7", "ID.8", "ID.9", "ID.10", "Block.ID")
colnames(R.2.3) = c("ID.1", "ID.2", "ID.3", "ID.4", "ID.5", "ID.6","ID.7", "ID.8", "ID.9", "ID.10", "Block.ID")


Group.2.Data = rbind(R.2.1, R.2.2, R.2.3)


#Check
length(unique( c(R.2.1[1,1:10], R.2.2[1,1:10], R.2.3[1,1:10])  ))

saveRDS(Group.2.Data, file = "C:/Users/carryp/R/Simulation Batch Effect/Data/Optimal.Block.Group.2.Rda")


###########################################################################################################
###########################################  STEP E: Selecting Batches  ###################################
###########################################################################################################

#STEP #2: Selecting batches based on unique combinations of the unique blocks identified in step 1

#NOTE: Be sure to load the phenotype file

Group.1.Data = readRDS(file = "C:/Users/carryp/R/Simulation Batch Effect/Data/Optimal.Block.Group.1.Rda")
Group.2.Data = readRDS(file = "C:/Users/carryp/R/Simulation Batch Effect/Data/Optimal.Block.Group.2.Rda")

Data.3 = data.frame(expand.grid(Group.1.Data$Block.ID, 
                                Group.2.Data$Block.ID))
Data.3$Batch.ID = 1:nrow(Data.3)

#STEP #3:Estimate balance between the different groups (within batch balance betweeen groups) and between batches
Propensity.Batch = function(a){
  
  ID.Include = c( c(  Group.1.Data[which(Group.1.Data$Block.ID %in% a[1]), 1:10]  ), 
                 c(  Group.2.Data[which(Group.2.Data$Block.ID %in%  a[2]), 1:10]  ))
            
  Data.Selected = Pheno.Sim[which(Pheno.Sim$ID %in% ID.Include), ]
  Data.Selected$Selected = 1

  Data.Not.Selected = Pheno.Sim[which(Pheno.Sim$ID %in% setdiff(Pheno.Sim$ID, unlist(ID.Include))), ]
  Data.Not.Selected$Selected = 0
  
  Batch.Data = data.frame(rbind(Data.Selected, Data.Not.Selected))
  
  #Between group comparison
  LR.G = glm(Group.1~Age + HBA1C, family = "binomial", data = Data.Selected)
  P.G = predict(LR.G,type="response")

  #Estimates
  Av.G1 = mean(P.G[Data.Selected$Group.1==1])
  Av.G2 = mean(P.G[Data.Selected$Group.1==0])
  Group.Diff = abs(Av.G1- Av.G2)
  
  #Between batch comparison
  LR.B = glm(Selected~Age + HBA1C, family = "binomial", data = Batch.Data)
  P.B = predict(LR.B,type="response")
  
  #Estimates
  Av.S = mean(P.B[Batch.Data$Selected==1])
  Av.nS = mean(P.B[Batch.Data$Selected==0])
  Batch.Diff = abs(Av.S- Av.nS)

  output =  c(Av.G1, Av.G2, Group.Diff, Av.S, Av.nS, Batch.Diff)
  
  names(output) = c("Av.G1", "Av.G2", "Group.Diff", "Av.S", "Av.nS", "Batch.Diff")
  
  return(output)
  
}

Output.Batch.Comp = data.frame(t(apply(Data.3, 1, Propensity.Batch)))

Data.4 = data.frame(cbind(Data.3, Output.Batch.Comp))

Data.4$G.Diff.Rank = rank(Data.4$Group.Diff)
Data.4$B.Diff.Rank = rank(Data.4$Batch.Diff)
Data.4$Av.Rank = apply(Data.4[, c("G.Diff.Rank", "B.Diff.Rank")], 1, mean)

Data.4.S = Data.4[order(Data.4$Av.Rank),]

###########################################################################################################
##################################################  STEP F:  ##############################################
###########################################################################################################


Batch.1 = unlist(c( Group.1.Data[which(Group.1.Data$Block.ID %in% Data.4.S[1, 1]), 1:10], 
             Group.2.Data[which(Group.2.Data$Block.ID %in% Data.4.S[1, 2]), 1:10]))
Row.1 = 1
i = 2
B.length = 0
while (B.length <40) {
  B.Two = c(Batch.1, 
          c( unlist( c( Group.1.Data[which(Group.1.Data$Block.ID %in% Data.4.S[i, 1]), 1:10], 
                        Group.2.Data[which(Group.2.Data$Block.ID %in% Data.4.S[i, 2]), 1:10])) ))
  
  B.length = length(unique(B.Two))
  print(i)
  Batch.2 = c( c( Group.1.Data[which(Group.1.Data$Block.ID %in% Data.4.S[i, 1]), 1:10], 
                  Group.2.Data[which(Group.2.Data$Block.ID %in% Data.4.S[i, 2]), 1:10]))
  B.Row.2 = i
  i = i+1
}

while (B.length <60) {
  B.Three = c(B.Two, 
              c( Group.1.Data[which(Group.1.Data$Block.ID %in% Data.4.S[i, 1]), 1:10], 
                 Group.2.Data[which(Group.2.Data$Block.ID %in% Data.4.S[i, 2]), 1:10]))
              
  B.length = length(unique(B.Three))
  Batch.3 = c( c( Group.1.Data[which(Group.1.Data$Block.ID %in% Data.4.S[i, 1]), 1:10], 
                  Group.2.Data[which(Group.2.Data$Block.ID %in% Data.4.S[i, 2]), 1:10]))
  B.Row.3 = i
  print(i)
  i = i+1
}


#Creating phenotype file
Batch.1.Final = Pheno.Sim[Pheno.Sim$ID %in% Batch.1, ]
Batch.1.Final$Batch = "B1"
Batch.2.Final = Pheno.Sim[Pheno.Sim$ID %in% Batch.2, ]
Batch.2.Final$Batch = "B2"
Batch.3.Final = Pheno.Sim[Pheno.Sim$ID %in% Batch.3, ]
Batch.3.Final$Batch = "B3"

Pheno.Final.Output = rbind(Batch.1.Final, Batch.2.Final, Batch.3.Final)

Pheno.Optimal = merge(Pheno.Sim, Pheno.Final.Output, by = "ID")

saveRDS(Pheno.Final.Output, "Path:/Data/Pheno.Optimal.Rda")

