#R.Allocation Function (Randomiation Only)
R.Allocation =  function(Pheno.Data,
                         N.Iter,
                         Initial.Seed,
                         B1, B2, B3, B4) {
  
  set.seed(Initial.Seed)
  
  Rand.Space.R = matrix(NA, nrow = N.Iter, ncol = nrow(Pheno.Data)+1)
  
  for (n in 1:N.Iter){
    new_seed <- sample(1:2147483647, 1)
    
    set.seed( new_seed )
    
    Rand.Space.R[n,] = c(new_seed, 
                         sample(c(rep("B1", B1), rep("B2", B2), rep("B3", B3), rep("B4", B4))) )
  }
  
  return(data.frame(Rand.Space.R))
  
}