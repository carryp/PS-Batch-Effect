
#FROM: "GSE50398", queried using "GEOquery" R package

#Loading raw phenotype file
meQTM.Phenotype = readRDS(file = "Path:/Data/meQTM.Phenotype.Array.Rda")   

#Selecting relevant covariates
Batch.Pheno.0 = meQTM.Phenotype[, c("ID", "Age", "HBA1C")]

#Limiting to complete data only
Batch.Pheno = na.omit(Batch.Pheno.0)

#Selecting subjects for analysis (removing outliers/subjects associated with extreme values)
Batch.Pheno$Age.Z = (Batch.Pheno$Age-mean(Batch.Pheno$Age))/sd(Batch.Pheno$Age)
Batch.Pheno$HBA1C.Z = (Batch.Pheno$HBA1C-mean(Batch.Pheno$HBA1C))/sd(Batch.Pheno$HBA1C)

Batch.Pheno$P = exp(2.5 + log(1.5)*Batch.Pheno$Age.Z + log(1.5)*Batch.Pheno$HBA1C.Z) / (1 + exp(2.5 + log(1.5)*Batch.Pheno$Age.Z + log(1.5)*Batch.Pheno$HBA1C.Z))

Batch.Pheno.Order = Batch.Pheno[order(Batch.Pheno$P), ]

Batch.Pheno.Final = Batch.Pheno.Order[1:60, ]
Batch.Pheno.Final$Group = c(rep(1, 30), rep(2, 30))
Batch.Pheno.Final$Group.1 = ifelse(Batch.Pheno.Final$Group==1, 1, 0)
Batch.Pheno.Final$Group.2 = ifelse(Batch.Pheno.Final$Group==2, 1, 0)


Batch.Pheno.Final.Order = Batch.Pheno.Final[order(Batch.Pheno.Final$ID), ]


#NOTE: loading the 'raw' expression set where batch effects have already been removed (per info on GEO,  see GEO_GSM1218048_Phenotype.csv in the data folder)
Expression = readRDS(file = "Path:/Data/meQTM.Expression.Array.Rda")   

Expression.Include = data.frame(Expression[, which(colnames(Expression) %in% Batch.Pheno.Final.Order$ID)])

Expression.Include.Order = Expression.Include[, order(Batch.Pheno.Final.Order$ID)]

    #Check Must Sum to 0
    Check.1 = Batch.Pheno.Final.Order$ID
    Check.2 = colnames(Expression.Include.Order)
    sum(ifelse(Check.1==Check.2, 0, 1))


#Seleting top 10,000 most variable genes
SD = apply(Expression.Include.Order, 1, sd)

Summary = data.frame(row.names(Expression.Include.Order), SD)
colnames(Summary) = c("Gene", "SD")

Summary.Order = Summary[order(Summary$SD, decreasing = TRUE), ]

Genes.Include = Summary.Order[1:10000, ]

Expression.Final = Expression.Include.Order[which(row.names(Expression.Include.Order) %in% Genes.Include$Gene), ]


#Saving datasets for analysis
saveRDS(Batch.Pheno.Final.Order, file = "Path:/Data/Pheno.Batch.Effect.Rda")

saveRDS(Expression.Final, file = "Path:/Data/Expression.Data.Rda")

