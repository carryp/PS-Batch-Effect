
Batch.Pheno = readRDS(file="Path:/Data/Pheno.Batch.Effect.Rda")

Expression = readRDS(file="Path:/Data/Expression.Data.Rda")


    #Check Must Sum to 0
    Check.1 = Batch.Pheno$ID
    Check.2 = colnames(Expression)
    sum(ifelse(Check.1==Check.2, 0, 1))


#Covariates
HBA1C = Batch.Pheno$HBA1C.Z
Age = Batch.Pheno$Age.Z


meQTM.lm = function(a){
  
  LM.C1 = lm (a ~ Age + HBA1C)
  
  B.Age = summary(LM.C1)$coefficients[2, 1]
  P.Age = summary(LM.C1)$coefficients[2, 4]
  
  B.A1C = summary(LM.C1)$coefficients[3, 1]
  P.A1C = summary(LM.C1)$coefficients[3, 4]
  
  Output = c(B.Age, P.Age, B.A1C, P.A1C)
  
  
  names(Output) = c("Beta.Age", "P.Age", "Beta.A1C", "P.A1C")
  
  return(Output)
}

Output.meQTM = data.frame(t(apply(Expression, 1, meQTM.lm)))

#Estimating fold change
Output.meQTM$Age.Fold = 2^abs(Output.meQTM$Beta.Age)
Output.meQTM$A1C.Fold = 2^abs(Output.meQTM$Beta.A1C)

Output.meQTM$FDR.Age = p.adjust(Output.meQTM$P.Age, method="BH")
Output.meQTM$FDR.A1C = p.adjust(Output.meQTM$P.A1C, method="BH")

#Hit defined as nominal p value <0.01 
Output.meQTM$Age.Hit = ifelse(Output.meQTM$Age.Fold>1.15 & Output.meQTM$P.Age<0.01, 1, 0)
Output.meQTM$A1C.Hit = ifelse(Output.meQTM$A1C.Fold>1.15 & Output.meQTM$P.A1C<0.01, 1, 0)

Output.meQTM$Hit.Both = ifelse(Output.meQTM$Age.Hit==1 & Output.meQTM$A1C.Hit==1, 1,0)
Output.meQTM$Hit.Either = ifelse(Output.meQTM$Age.Hit==1 | Output.meQTM$A1C.Hit==1, 1,0)


#This dataset represents the 'true' associations
saveRDS(Output.meQTM, "Path:/Output/Age.HBA1C.Output.Data.Rda")


#identifying probe associated with both age and HbA1c
View(Output.meQTM[Output.meQTM$Hit.Both==1, ])

#8051275
#This probe will be used to test performance of algorithm under the alternative hypothesis
