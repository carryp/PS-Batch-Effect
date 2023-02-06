#Balancing function
Almost.Optimal  = function(Covariates, Pheno.Dataset, Randomization.Dataset){
  
  #Setting to NA
  Av.S1 = NA
  Av.S2 = NA
  Av.S3 = NA
  Av.S4 = NA
  
  Av.nS1 = NA
  Av.nS2 = NA
  Av.nS3 = NA
  Av.nS4 = NA
  
  Diff.1 = NA
  Diff.2 = NA
  Diff.3 = NA
  Diff.4 = NA
  
  Av.B = NA
  Mx.B = NA
  Min.B = NA
  
  
  Randomization.Options = Randomization.Dataset
  
  Randomization.Options$Iteration = 1:nrow(Randomization.Options)
  
  
  Batch.Propensity.Score = function(b){
    
    n = b[ ncol(Randomization.Options)   ]
    a = b[ 2:(ncol(Randomization.Options)-1)]
    
    
    N.Batch = length( unique(as.list(a)) )
    Batch.List.Order =  sort( unique(as.character(a)) )
    
    
    for(i in 1:N.Batch){
      
      B <- paste("B", i, sep = "")
      assign(B, as.vector(ifelse(a==paste( Batch.List.Order[i] ), 1, 0)))
      
      LR = paste("LR.", i, sep = "")
      assign(LR, glm( as.formula( paste(  paste(Batch.List.Order[i]), 
                                          " ~ ", 
                                          paste(  gsub(" ","", gsub(",","+",toString(Covariates)) ) ))),  
                      family = 'binomial', 
                      data = Pheno.Dataset ) ) 
      
      
      P = paste("P", i, sep = "")
      assign(P, predict( eval(as.name(paste("LR.", i, sep = ""))), type="response"))
      
      Av.S = paste("Av.S", i, sep="")
      assign(Av.S,  mean( eval(as.name(paste("P", i, sep = "")))[ eval(as.name(paste(Batch.List.Order[i])))==1] ) )
      
      Av.nS = paste("Av.nS", i, sep="")
      assign(Av.nS,  mean( eval(as.name(paste("P", i, sep = "")))[ eval(as.name(paste(Batch.List.Order[i])))==0] ) )
      
      Diff = paste("Diff.", i, sep="")
      assign(Diff,  abs( eval(as.name(paste( paste("Av.S", i, sep="") )))- eval(as.name(paste( paste("Av.nS", i, sep="") )))     ))
      
    }
    
    
    
    Av.B = mean( c(na.omit(Diff.1), na.omit(Diff.2), na.omit(Diff.3), na.omit(Diff.4)) )
    Mx.B = max( c(na.omit(Diff.1), na.omit(Diff.2), na.omit(Diff.3), na.omit(Diff.4)) )
    Min.B = min( c(na.omit(Diff.1), na.omit(Diff.2), na.omit(Diff.3), na.omit(Diff.4)) )
    
    output =  c(Av.B, Mx.B, Min.B)
    
    names(output) = c("Av.B", "Mx.B", "Min.B")
    
    cat("\rRandomization Iterations Evaluated", n, "of", nrow(Randomization.Options))
    
    return(output)
    
    
  }
  
  
  Out.Data = data.frame( t(apply(Randomization.Options, 1, Batch.Propensity.Score) )) 
  
  #Creating new phenotype file
  Pheno.Final = Pheno.Dataset
  
  #Batch assignment based on average difference
  Pheno.Final$Batch.Assignment.Av = as.character(t(Randomization.Options[row.names(Out.Data[order(Out.Data$Av.B ), ][1,]), 2:(ncol(Randomization.Options)-1)]))
  Pheno.Final$Av.Rand.N = row.names(Out.Data[order(Out.Data$Av.B ), ][1,])
  
  #unbalanced - for checking integrity of randomization
  Pheno.Final$Batch.Assignment.Imbalance.Av = as.character(t(Randomization.Options[row.names(Out.Data[order(Out.Data$Av.B ), ][nrow(Out.Data),]), 2:(ncol(Randomization.Options)-1)]))
  Pheno.Final$Imbalance.Av.N = row.names(Out.Data[order(Out.Data$Av.B ), ][nrow(Out.Data),])
  
  return(Pheno.Final)
  
}