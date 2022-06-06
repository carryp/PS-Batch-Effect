
#Updated: 01-17-2022

library(sva)
library(experiment)

Pheno.Sim = readRDS(file="Path:/Data/Pheno.Batch.Effect.Rda")

Age = Pheno.Sim$Age
HBA1C = Pheno.Sim$HBA1C
Group.1 = Pheno.Sim$Group.1
Group.2 = Pheno.Sim$Group.2
ID = Pheno.Sim$ID


#Optimal permutation
Pheno.Optimal = readRDS(file="Path:/Data/Pheno.Optimal.Rda")

Pheno.Optimal.Sort = Pheno.Optimal[order(Pheno.Optimal$ID), ]

Opt.Batch = Pheno.Optimal.Sort$Batch

  #Data check to confirm balance
  table(Pheno.Optimal.Sort$Batch, Pheno.Optimal.Sort$Group.1)
  
  #Data check to confirm ordering of IDs matches, must sum to 0
  Check.1 = Pheno.Optimal.Sort$ID
  Check.2 = Pheno.Sim$ID
  sum(ifelse(Check.1==Check.2, 0, 1))
  
  
E.All = readRDS(file="Path:/Data/Expression.Data.Rda")  

Output = readRDS(file="Path:/Output/Age.HBA1C.Output.Data.Rda") 

#***NOTE*** both conditions 'optimal permutation' vs randomization use the same dataset

#Batch assignments (to reproduce results in manuscript)
Batch.Assignment =readRDS(file="Path:/Output/Batch.Assignment.Revised.Rda")

################################################################################################################
########################################### SIMULATION #########################################################
################################################################################################################

start <- Sys.time()

replicates = 1000

Simulation.Output = matrix(NA, nrow = replicates, ncol = 89)

for (n in 1:replicates){
  
  #Setting to null to prevent carry over from previous iteration
  Pheno.Sim$Batch.R = NA
  Pheno.Sim$Batch.SR = NA

  Y.True = NA 
  
  Noise = NA 
  
  Expression.Data = NA 
  
  R.I.1 = NA 
  R.I.2 = NA 
  R.I.3 = NA 
  
  SR.I.1 = NA 
  SR.I.2 = NA 
  SR.I.3 = NA 
  
  Y.Obs.R = NA 
  Output.Data.R = NA 
  
  Y.Obs.SR = NA 
  Output.Data.SR = NA 
  
  Opt.I.1 = NA 
  Opt.I.2 = NA 
  Opt.I.3 = NA 
  
  Y.Obs.Opt = NA 
  Output.Data.Obs.Optimal = NA 
  
  Output.Data.Corrected.R = NA 
  Output.Data.Corrected.SR = NA 
  Output.Data.Corrected.Optimal = NA 
  
  Diff1.1 = NA 
  Diff1.2 = NA 
  Diff1.3 = NA 
  Diff1.4 = NA 
  Diff1.5 = NA 
  Diff1.6 = NA 
  
  Diff1.7 = NA 
  Diff1.8 = NA 
  Diff1.9 = NA 
  
  Diff2.1 = NA 
  Diff2.2 = NA 
  Diff2.3 = NA 
  Diff2.4 = NA 
  Diff2.5 = NA 
  Diff2.6 = NA 
  Diff2.7 = NA 
  Diff2.8 = NA 
  Diff2.9 = NA 
  
  Diff3.1 = NA 
  Diff3.2 = NA 
  Diff3.3 = NA 
  Diff3.4 = NA 
  Diff3.5 = NA 
  Diff3.6 = NA 
  Diff3.7 = NA 
  Diff3.8 = NA 
  Diff3.9 = NA 

  mod.0 = NA 
  batch.0 = NA 
  
  mod.1 = NA 
  batch.1 = NA 
  
  mod.2 = NA 
  batch.2 = NA 
  
  LM.Output.True = NA 

  
  LM.Output.R = NA 

  
  LM.Output.SR = NA 

  
  LM.Output.Optimal = NA 

  
  LM.Output.Corrected.R = NA 

  
  LM.Output.Corrected.SR = NA 

  
  LM.Output.Corrected.Optimal = NA 
  
  LM.Output.FE.R = NA 
  
  LM.Output.FE.SR = NA 
  LM.Output.FE.Optimal = NA 

  
  #Creating the expression set (Y - truth or biological variation)
  
  N.Features = 10000
  N.Subjects = nrow(Pheno.Sim)
  
  #Expression Dataset
  
  #Setting seed so that it is reproducible
  set.seed(0430227+n)
  
  #Matrix algebra to save computation time
  Y.True = t(E.All)
  
  
  #Batch effects (2 x median biological variation)
  E = rnorm(3, 0, 2*0.2828492)
  
  Noise = matrix(rnorm(N.Subjects*N.Features, 0, 0.1), N.Subjects, N.Features)
  
  Expression.Data = data.frame(t(Y.True))

  
  print(1)
  
  #Randomization Condition
  Pheno.Sim$Batch.R = unlist(t(Batch.Assignment[n, 2:61]))
  
  #Stratified Randomization condition
  Pheno.Sim$Batch.SR = unlist(t(Batch.Assignment[n, 62:121]))
  
  #Batch effect - Randomization Condition
  R.I.1 = as.numeric(ifelse(Pheno.Sim$Batch.R=="B1",1, 0))
  R.I.2 = as.numeric(ifelse(Pheno.Sim$Batch.R=="B2",1, 0))
  R.I.3 = as.numeric(ifelse(Pheno.Sim$Batch.R=="B3",1, 0))
  
  #Adding in the batch effect - Randomization Condition
  Y.Obs.R = Y.True +  R.I.1*E[1] +  R.I.2*E[2] + R.I.3*E[3] + Noise
  
  #Creating the 'observed' data - with the batch effect error - Randomization condition
  Output.Data.R = data.frame(t(Y.Obs.R))
  colnames(Output.Data.R) = ID
  
  #Batch effect - Stratified Randomization Condition
  SR.I.1 = as.numeric(ifelse(Pheno.Sim$Batch.SR=="B1",1, 0))
  SR.I.2 = as.numeric(ifelse(Pheno.Sim$Batch.SR=="B2",1, 0))
  SR.I.3 = as.numeric(ifelse(Pheno.Sim$Batch.SR=="B3",1, 0))
  
  #Adding in the batch effect - Stratified Randomization Condition
  Y.Obs.SR = Y.True +  SR.I.1*E[1] +  SR.I.2*E[2] + SR.I.3*E[3] + Noise
  
  #Creating the 'observed' data - with the batch effect error - Randomization condition
  Output.Data.SR = data.frame(t(Y.Obs.SR))

  #Batch effect - Optimal Condition
  Opt.I.1 = ifelse(Pheno.Optimal.Sort$Batch=="B1",1, 0)
  Opt.I.2 = ifelse(Pheno.Optimal.Sort$Batch=="B2",1, 0)
  Opt.I.3 = ifelse(Pheno.Optimal.Sort$Batch=="B3",1, 0)
  
  #Adding in the batch effect - Optimal Condition
  Y.Obs.Opt = Y.True + Opt.I.1*E[1] +  Opt.I.2*E[2] + Opt.I.3*E[3] + Noise
  
  #Creating the 'observed' data - with the batch effect error - Randomization condition
  Output.Data.Optimal = data.frame( t(Y.Obs.Opt))
  
  print(2)
  
  #Correcting for batch effects  
  
 
  #batch correction using combat - Randomization condition
  mod.0 = model.matrix(~as.factor(Group.1)+Age+HBA1C, data=Pheno.Sim)
  batch.0 = as.factor(Pheno.Sim$Batch.R)
  Output.Data.Corrected.R = ComBat(as.matrix(Output.Data.R), batch.0, mod.0)
  
  #batch correction using combat - Stratified Randomization condition
  mod.1 = model.matrix(~as.factor(Group.1)+Age+HBA1C, data=Pheno.Sim)
  batch.1 = as.factor(Pheno.Sim$Batch.SR)
  Output.Data.Corrected.SR = ComBat(as.matrix(Output.Data.SR), batch.1, mod.1)
  
  #batch correction using combat - optimal condition
  mod.2 = model.matrix(~as.factor(Group.1)+Age+HBA1C, data=Pheno.Optimal.Sort)
  batch.2 = as.factor(Pheno.Optimal.Sort$Batch)
  Output.Data.Corrected.Optimal = ComBat(as.matrix(Output.Data.Optimal), batch.2, mod.2)
  
  print(3)
  
  #Running the model
  Linear.Model = function(a){
    Model <- try(  lm (a ~  Group.1 + Age + HBA1C)      )
    want = rep(NA, length=3)
    if("try-error" %in% class(Model)){want = want}else{
      
      #Getting the parameter estimates
      Beta.1 = summary(Model)$coefficients[2, 1]
      Beta.2 = summary(Model)$coefficients[3, 1]
      Beta.3 = summary(Model)$coefficient[4, 1]
      
      SE.1 = summary(Model)$coefficients[2, 2]
      SE.2 = summary(Model)$coefficients[3, 2]
      SE.3 = summary(Model)$coefficients[4, 2]
      
      P.1 = summary(Model)$coefficients[2, 4]
      P.2 = summary(Model)$coefficients[3, 4]
      P.3 = summary(Model)$coefficients[4, 4]
      
      want = c(Beta.1, Beta.2, Beta.3, SE.1, SE.2, SE.3, P.1, P.2, P.3)
    }
    
    names(want) = c("Beta.1", "Beta.2", "Beta.3", "SE.1", "SE.2", "SE.3", "P.1", "P.2", "P.3")
    
    return(want)
  }
  
  print(4)
  LM.Output.True = data.frame(t(apply(Expression.Data, 1, Linear.Model)))
  
  
  print(5)
  LM.Output.R = data.frame(t(apply(Output.Data.R, 1, Linear.Model)))
  
  
  print(6)
  LM.Output.SR = data.frame(t(apply(Output.Data.SR, 1, Linear.Model)))
  
  
  print(7)
  LM.Output.Optimal = data.frame(t(apply(Output.Data.Optimal, 1, Linear.Model)))
  
  
  print(8)
  LM.Output.Corrected.R = data.frame(t(apply(Output.Data.Corrected.R, 1, Linear.Model)))
  
  
  print(9)
  LM.Output.Corrected.SR = data.frame(t(apply(Output.Data.Corrected.SR, 1, Linear.Model)))
  

  print(10)
  LM.Output.Corrected.Optimal = data.frame(t(apply(Output.Data.Corrected.Optimal, 1, Linear.Model)))
  
  
  print(11)
  Linear.Model.2 = function(a){
    Model <- try(  lm (a ~  Group.1 +Age+HBA1C + R.I.1 + R.I.2)      )
    want = rep(NA, length=3)
    if("try-error" %in% class(Model)){want = want}else{
      
      #Getting the parameter estimates
      Beta.1 = summary(Model)$coefficients[2, 1]
      Beta.2 = summary(Model)$coefficients[3, 1]
      Beta.3 = summary(Model)$coefficient[4, 1]
      
      SE.1 = summary(Model)$coefficients[2, 2]
      SE.2 = summary(Model)$coefficients[3, 2]
      SE.3 = summary(Model)$coefficients[4, 2]
      
      P.1 = summary(Model)$coefficients[2, 4]
      P.2 = summary(Model)$coefficients[3, 4]
      P.3 = summary(Model)$coefficients[4, 4]
      
      want = c(Beta.1, Beta.2, Beta.3, SE.1, SE.2, SE.3, P.1, P.2, P.3)
    }
    
    names(want) = c("Beta.1", "Beta.2", "Beta.3", "SE.1", "SE.2", "SE.3", "P.1", "P.2", "P.3")
    
    return(want)
  }
  

  #Including batch as a fixed effect: Randomization condition
  LM.Output.FE.R = data.frame(t(apply(Output.Data.R, 1, Linear.Model.2)))

  print(12)
  
  Linear.Model.3 = function(a){
    Model <- try(  lm (a ~  Group.1+Age + HBA1C + SR.I.1 + SR.I.2)      )
    want = rep(NA, length=3)
    if("try-error" %in% class(Model)){want = want}else{
      
      #Getting the parameter estimates
      Beta.1 = summary(Model)$coefficients[2, 1]
      Beta.2 = summary(Model)$coefficients[3, 1]
      Beta.3 = summary(Model)$coefficient[4, 1]
      
      SE.1 = summary(Model)$coefficients[2, 2]
      SE.2 = summary(Model)$coefficients[3, 2]
      SE.3 = summary(Model)$coefficients[4, 2]
      
      P.1 = summary(Model)$coefficients[2, 4]
      P.2 = summary(Model)$coefficients[3, 4]
      P.3 = summary(Model)$coefficients[4, 4]
      
      want = c(Beta.1, Beta.2, Beta.3, SE.1, SE.2, SE.3, P.1, P.2, P.3)
    }
    
    names(want) = c("Beta.1", "Beta.2", "Beta.3", "SE.1", "SE.2", "SE.3", "P.1", "P.2", "P.3")
    
    return(want)
  }
  
  
  #Including batch as a fixed effect: Randomization condition
  LM.Output.FE.SR = data.frame(t(apply(Output.Data.SR, 1, Linear.Model.3)))
  
  

  print(13)
  Linear.Model.4 = function(a){
    Model <- try(  lm (a ~  Group.1 +Age+HBA1C + Opt.I.1 + Opt.I.2)      )
    want = rep(NA, length=3)
    if("try-error" %in% class(Model)){want = want}else{
      
      #Getting the parameter estimates
      Beta.1 = summary(Model)$coefficients[2, 1]
      Beta.2 = summary(Model)$coefficients[3, 1]
      Beta.3 = summary(Model)$coefficient[4, 1]
      
      SE.1 = summary(Model)$coefficients[2, 2]
      SE.2 = summary(Model)$coefficients[3, 2]
      SE.3 = summary(Model)$coefficients[4, 2]
      
      P.1 = summary(Model)$coefficients[2, 4]
      P.2 = summary(Model)$coefficients[3, 4]
      P.3 = summary(Model)$coefficients[4, 4]
      
      want = c(Beta.1, Beta.2, Beta.3, SE.1, SE.2, SE.3, P.1, P.2, P.3)
    }
    
    names(want) = c("Beta.1", "Beta.2", "Beta.3", "SE.1", "SE.2", "SE.3", "P.1", "P.2", "P.3")
    
    return(want)
  }

  
  #Including batch as a fixed effect  
  LM.Output.FE.Optimal = data.frame(t(apply(Output.Data.Optimal, 1, Linear.Model.4)))
  
 
  print(14)
  Diff1.1 = LM.Output.R$Beta.1 - LM.Output.True$Beta.1
  Diff1.2 = LM.Output.SR$Beta.1 - LM.Output.True$Beta.1
  Diff1.3 = LM.Output.Optimal$Beta.1 - LM.Output.True$Beta.1
  Diff1.4 = LM.Output.Corrected.R$Beta.1 - LM.Output.True$Beta.1
  Diff1.5 = LM.Output.Corrected.SR$Beta.1 - LM.Output.True$Beta.1
  Diff1.6 = LM.Output.Corrected.Optimal$Beta.1 - LM.Output.True$Beta.1
  
  Diff2.1 = LM.Output.R$Beta.2 - LM.Output.True$Beta.2
  Diff2.2 = LM.Output.SR$Beta.2 - LM.Output.True$Beta.2
  Diff2.3 = LM.Output.Optimal$Beta.2 - LM.Output.True$Beta.2
  Diff2.4 = LM.Output.Corrected.R$Beta.2 - LM.Output.True$Beta.2
  Diff2.5 = LM.Output.Corrected.SR$Beta.2- LM.Output.True$Beta.2
  Diff2.6 = LM.Output.Corrected.Optimal$Beta.2 - LM.Output.True$Beta.2
  
  Diff3.1 = LM.Output.R$Beta.3 - LM.Output.True$Beta.3
  Diff3.2 = LM.Output.SR$Beta.3 - LM.Output.True$Beta.3
  Diff3.3 = LM.Output.Optimal$Beta.3 - LM.Output.True$Beta.3
  Diff3.4 = LM.Output.Corrected.R$Beta.3 - LM.Output.True$Beta.3
  Diff3.5 = LM.Output.Corrected.SR$Beta.3- LM.Output.True$Beta.3
  Diff3.6 = LM.Output.Corrected.Optimal$Beta.3 - LM.Output.True$Beta.3
  
  Diff1.7 = LM.Output.FE.R$Beta.1 - LM.Output.True$Beta.1 
  Diff1.8 = LM.Output.FE.SR$Beta.1 - LM.Output.True$Beta.1
  Diff1.9 = LM.Output.FE.Optimal$Beta.1 - LM.Output.True$Beta.1

  Diff2.7 = LM.Output.FE.R$Beta.2 - LM.Output.True$Beta.2 
  Diff2.8 = LM.Output.FE.SR$Beta.2 - LM.Output.True$Beta.2
  Diff2.9 = LM.Output.FE.Optimal$Beta.2 - LM.Output.True$Beta.2

  Diff3.7 = LM.Output.FE.R$Beta.3 - LM.Output.True$Beta.3 
  Diff3.8 = LM.Output.FE.SR$Beta.3 - LM.Output.True$Beta.3
  Diff3.9 = LM.Output.FE.Optimal$Beta.3 - LM.Output.True$Beta.3
  
  print(15)
  
  #Outputting the results
  Simulation.Output[n, 1] =   n
  Simulation.Output[n, 2] =   mean(LM.Output.True$SE.1)
  Simulation.Output[n, 3] =   mean(LM.Output.R$SE.1)
  Simulation.Output[n, 4] =   mean(LM.Output.SR$SE.1)
  Simulation.Output[n, 5] =   mean(LM.Output.Optimal$SE.1)
  Simulation.Output[n, 6] =   mean(LM.Output.Corrected.R$SE.1)
  Simulation.Output[n, 7] =   mean(LM.Output.Corrected.SR$SE.1)
  Simulation.Output[n, 8] =   mean(LM.Output.Corrected.Optimal$SE.1)
  Simulation.Output[n, 9] =   mean(LM.Output.FE.R$SE.1)
  Simulation.Output[n, 10] =   mean(LM.Output.FE.SR$SE.1)
  Simulation.Output[n, 11] =   mean(LM.Output.FE.Optimal$SE.1)
  #Platform level summary: SE (X1)
  Simulation.Output[n, 12] =   mean(LM.Output.True$SE.2)
  Simulation.Output[n, 13] =   mean(LM.Output.R$SE.2)
  Simulation.Output[n, 14] =   mean(LM.Output.SR$SE.2)
  Simulation.Output[n, 15] =   mean(LM.Output.Optimal$SE.2)
  Simulation.Output[n, 16] =   mean(LM.Output.Corrected.R$SE.2)
  Simulation.Output[n, 17] =   mean(LM.Output.Corrected.SR$SE.2)
  Simulation.Output[n, 18] =   mean(LM.Output.Corrected.Optimal$SE.2)
  Simulation.Output[n, 19] =   mean(LM.Output.FE.R$SE.2)
  Simulation.Output[n, 20] =   mean(LM.Output.FE.SR$SE.2)
  Simulation.Output[n, 21] =   mean(LM.Output.FE.Optimal$SE.2)
  #Platform level summary: SE (X2)
  Simulation.Output[n, 22] =   mean(LM.Output.True$SE.3)
  Simulation.Output[n, 23] =   mean(LM.Output.R$SE.3)
  Simulation.Output[n, 24] =   mean(LM.Output.SR$SE.3)
  Simulation.Output[n, 25] =   mean(LM.Output.Optimal$SE.3)
  Simulation.Output[n, 26] =   mean(LM.Output.Corrected.R$SE.3)
  Simulation.Output[n, 27] =   mean(LM.Output.Corrected.SR$SE.3)
  Simulation.Output[n, 28] =   mean(LM.Output.Corrected.Optimal$SE.3)
  Simulation.Output[n, 29] =   mean(LM.Output.FE.R$SE.3)
  Simulation.Output[n, 30] =   mean(LM.Output.FE.SR$SE.3)
  Simulation.Output[n, 31] =   mean(LM.Output.FE.Optimal$SE.3)
  Simulation.Output[n, 32] =  max(abs(Diff1.1))
  Simulation.Output[n, 33] =  max(abs(Diff1.2))
  Simulation.Output[n, 34] =  max(abs(Diff1.3))
  Simulation.Output[n, 35] =  max(abs(Diff1.4))
  Simulation.Output[n, 36] =  max(abs(Diff1.5))
  Simulation.Output[n, 37] =  max(abs(Diff1.6))
  Simulation.Output[n, 38] =  max(abs(Diff1.7))
  Simulation.Output[n, 39] =  max(abs(Diff1.8))
  Simulation.Output[n, 40] =  max(abs(Diff1.9))
  #Bias (X1)
  Simulation.Output[n, 41] =  max(abs(Diff2.1))
  Simulation.Output[n, 42] =  max(abs(Diff2.2))
  Simulation.Output[n, 43] =  max(abs(Diff2.3))
  Simulation.Output[n, 44] =  max(abs(Diff2.4))
  Simulation.Output[n, 45] =  max(abs(Diff2.5))
  Simulation.Output[n, 46] =  max(abs(Diff2.6))
  Simulation.Output[n, 47] =  max(abs(Diff2.7))
  Simulation.Output[n, 48] =  max(abs(Diff2.8))
  Simulation.Output[n, 49] =  max(abs(Diff2.9))
  #Bias (X2)
  Simulation.Output[n, 50] =  max(abs(Diff3.1))
  Simulation.Output[n, 51] =  max(abs(Diff3.2))
  Simulation.Output[n, 52] =  max(abs(Diff3.3))
  Simulation.Output[n, 53] =  max(abs(Diff3.4))
  Simulation.Output[n, 54] =  max(abs(Diff3.5))
  Simulation.Output[n, 55] =  max(abs(Diff3.6))
  Simulation.Output[n, 56] =  max(abs(Diff3.7))
  Simulation.Output[n, 57] =  max(abs(Diff3.8))
  Simulation.Output[n, 58] =  max(abs(Diff3.9))
  #Alternative hypothesis, gene associated with age and HBA1c levels
  Simulation.Output[n, 59] = LM.Output.R["8051275", ]$Beta.1 - LM.Output.True["8051275", ]$Beta.1
  Simulation.Output[n, 60] = LM.Output.SR["8051275", ]$Beta.1 - LM.Output.True["8051275", ]$Beta.1
  Simulation.Output[n, 61] = LM.Output.Optimal["8051275", ]$Beta.1 - LM.Output.True["8051275", ]$Beta.1
  Simulation.Output[n, 62] = LM.Output.Corrected.R["8051275", ]$Beta.1 - LM.Output.True["8051275", ]$Beta.1
  Simulation.Output[n, 63] = LM.Output.Corrected.SR["8051275", ]$Beta.1 - LM.Output.True["8051275", ]$Beta.1
  Simulation.Output[n, 64] = LM.Output.Corrected.Optimal["8051275", ]$Beta.1 - LM.Output.True["8051275", ]$Beta.1
  
  Simulation.Output[n, 65] = LM.Output.R["8051275", ]$Beta.2 - LM.Output.True["8051275", ]$Beta.2   
  Simulation.Output[n, 66] = LM.Output.SR["8051275", ]$Beta.2 - LM.Output.True["8051275", ]$Beta.2      
  Simulation.Output[n, 67] = LM.Output.Optimal["8051275", ]$Beta.2 - LM.Output.True["8051275", ]$Beta.2      
  Simulation.Output[n, 68] = LM.Output.Corrected.R["8051275", ]$Beta.2 - LM.Output.True["8051275", ]$Beta.2
  Simulation.Output[n, 69] = LM.Output.Corrected.SR["8051275", ]$Beta.2- LM.Output.True["8051275", ]$Beta.2  
  Simulation.Output[n, 70] = LM.Output.Corrected.Optimal["8051275", ]$Beta.2 - LM.Output.True["8051275", ]$Beta.2 
  
  Simulation.Output[n, 71] = LM.Output.R["8051275", ]$Beta.3 - LM.Output.True["8051275", ]$Beta.3
  Simulation.Output[n, 72] =  LM.Output.SR["8051275", ]$Beta.3 - LM.Output.True["8051275", ]$Beta.3
  Simulation.Output[n, 73] = LM.Output.Optimal["8051275", ]$Beta.3 - LM.Output.True["8051275", ]$Beta.3
  Simulation.Output[n, 74] = LM.Output.Corrected.R["8051275", ]$Beta.3 - LM.Output.True["8051275", ]$Beta.3
  Simulation.Output[n, 75] = LM.Output.Corrected.SR["8051275", ]$Beta.3- LM.Output.True["8051275", ]$Beta.3
  Simulation.Output[n, 76] =  LM.Output.Corrected.Optimal["8051275", ]$Beta.3 - LM.Output.True["8051275", ]$Beta.3
  
  Simulation.Output[n, 77] = LM.Output.FE.R["8051275", ]$Beta.1 - LM.Output.True["8051275", ]$Beta.1 
  Simulation.Output[n, 78] = LM.Output.FE.SR["8051275", ]$Beta.1 - LM.Output.True["8051275", ]$Beta.1
  Simulation.Output[n, 79] = LM.Output.FE.Optimal["8051275", ]$Beta.1 - LM.Output.True["8051275", ]$Beta.1
  
  Simulation.Output[n, 80] = LM.Output.FE.R["8051275", ]$Beta.2 - LM.Output.True["8051275", ]$Beta.2 
  Simulation.Output[n, 81] = LM.Output.FE.SR["8051275", ]$Beta.2 - LM.Output.True["8051275", ]$Beta.2
  Simulation.Output[n, 82] = LM.Output.FE.Optimal["8051275", ]$Beta.2 - LM.Output.True["8051275", ]$Beta.2
  
  Simulation.Output[n, 83] = LM.Output.FE.R["8051275", ]$Beta.3 - LM.Output.True["8051275", ]$Beta.3 
  Simulation.Output[n, 84] = LM.Output.FE.SR["8051275", ]$Beta.3 - LM.Output.True["8051275", ]$Beta.3
  Simulation.Output[n, 85] =  LM.Output.FE.Optimal["8051275", ]$Beta.3 - LM.Output.True["8051275", ]$Beta.3
  
  #Output for confirming reproducibility of code
  Simulation.Output[n, 86] = E[1]
  Simulation.Output[n, 87] = E[2]
  Simulation.Output[n, 88] = E[3]
  
  Simulation.Output[n, 89] = Noise[30, 100]
  
  print("Iteration #:")
  print(n)
  
}


Simulation.Output.Data = data.frame(Simulation.Output)

colnames(Simulation.Output.Data) = 
c(
  'Iteration',
  'Av_SE_1_True',
  'Av_SE_1_R',
  'Av_SE_1_SR',
  'Av_SE_1_Optimal',
  'Av_SE_1_Corrected_R',
  'Av_SE_1_Corrected_SR',
  'Av_SE_1_Corrected_Optimal',
  'Av_SE_1_FE_R',
  'Av_SE_1_FE_SR',
  'Av_SE_1_FE_Optimal',
  
  'Av_SE_2_True',
  'Av_SE_2_R',
  'Av_SE_2_SR',
  'Av_SE_2_Optimal',
  'Av_SE_2_Corrected_R',
  'Av_SE_2_Corrected_SR',
  'Av_SE_2_Corrected_Optimal',
  'Av_SE_2_FE_R',
  'Av_SE_2_FE_SR',
  'Av_SE_2_FE_Optimal',
  
  'Av_SE_3_True',
  'Av_SE_3_R',
  'Av_SE_3_SR',
  'Av_SE_3_Optimal',
  'Av_SE_3_Corrected_R',
  'Av_SE_3_Corrected_SR',
  'Av_SE_3_Corrected_Optimal',
  'Av_SE_3_FE_R',
  'Av_SE_3_FE_SR',
  'Av_SE_3_FE_Optimal',
  
  'Mx_Bias1_1',
  'Mx_Bias2_1',
  'Mx_Bias3_1',
  'Mx_Bias4_1',
  'Mx_Bias5_1',
  'Mx_Bias6_1',
  'Mx_Bias7_1',
  'Mx_Bias8_1',
  'Mx_Bias9_1',
  
  'Mx_Bias1_2',
  'Mx_Bias2_2',
  'Mx_Bias3_2',
  'Mx_Bias4_2',
  'Mx_Bias5_2',
  'Mx_Bias6_2',
  'Mx_Bias7_2',
  'Mx_Bias8_2',
  'Mx_Bias9_2',
  
  'Mx_Bias1_3',
  'Mx_Bias2_3',
  'Mx_Bias3_3',
  'Mx_Bias4_3',
  'Mx_Bias5_3',
  'Mx_Bias6_3',
  'Mx_Bias7_3',
  'Mx_Bias8_3',
  'Mx_Bias9_3',
  
  'Alt1_1',
  'Alt1_2',
  'Alt1_3',
  'Alt1_4',
  'Alt1_5',
  'Alt1_6',
  
  'Alt2_1',
  'Alt2_2',
  'Alt2_3',
  'Alt2_4',
  'Alt2_5',
  'Alt2_6',
  
  'Alt3_1',
  'Alt3_2',
  'Alt3_3',
  'Alt3_4',
  'Alt3_5',
  'Alt3_6',
  
  'Alt7_1',
  'Alt8_1',
  'Alt9_1',
  
  'Alt7_2',
  'Alt8_2',
  'Alt9_2',
  
  'Alt7_3',
  'Alt8_3',
  'Alt9_3',
  
  'E_1',
  'E_2',
  'E_3',
  
  'Noise_30_100'
  
  
)

end <- Sys.time()


#Total computer time
end-start


saveRDS(Simulation.Output.Data, file= "Path:/Output/Simulation.Output.Batch.Effect.Revised.Rda")

