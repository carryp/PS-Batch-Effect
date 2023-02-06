#SR.Allocation Function (Stratified Randomization)
SR.Allocation =  function(Pheno.Data,
                          N.Iter,
                          Initial.Seed,
                          S1.B1, S1.B2, S1.B3, S1.B4, 
                          S2.B1, S2.B2, S2.B3, S2.B4) {
  
  
  Rand.Space.SR = matrix(NA, nrow = N.Iter, ncol = nrow(Pheno.Data)+1)
  
  set.seed(Initial.Seed)
  
  for (n in 1:N.Iter){
    new_seed <- sample(1:2147483647, 1)
    
    set.seed( new_seed )
    
    Rand.Space.SR[n,] = c(new_seed, 
                          sample(c(rep("B1", S1.B1), rep("B2", S1.B2), rep("B3", S1.B3), rep("B4", S1.B4))), 
                          sample(c(rep("B1", S2.B1), rep("B2", S2.B2), rep("B3", S2.B3), rep("B4", S2.B4))))
  }
  
  return(data.frame(Rand.Space.SR))
  
}