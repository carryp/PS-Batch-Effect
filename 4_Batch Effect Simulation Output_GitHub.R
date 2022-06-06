
#Beta 1 = case/control 
#Beta 2 = Age
#Beta 3 = HBA1C

#Loading simulation output
Batch.Effect = readRDS(file = "Path:/Output/Simulation.Output.Batch.Effect.Revised.Rda")

#############################################################################################################
########################### Defining the functions for Summarizing the Output ###############################
#############################################################################################################

Function.1 = function(a){
                      want= c(
                      mean(a), 
                      sd(a),
                      min(a),
                      max(a))
  return(want)
}

row.names.1 = c("True Condition", "Randomization Condition", "Stratified Randomization Condition",
"Optimal Condition", "Randomization Condition Combat Corrected", "Stratified Randomization Combat Corrected Condition",
"Optimal Combat Corrected Condition", "Randomization Condition Regression Corrected", 
"Stratified Randomization Regression Corrected Condition", "Optimal Regression Corrected Condition")

colnames.1 = c("Mean", "Stdev", "Min", "Max")


row.names.2 = c("Randomization vs True",
                "Stratified Randomization vs True",
                "Optimal vs True",
                "Randomization Combat vs True",
                "Stratified Randomization Combat vs True",
                "Optimal Combat vs True",
                "Randomization Regression vs True",
                "Stratified Randomization Regression vs True",
                "Optimal Regression vs True")




#NOTE: Table 1 = study design table


#############################################################################################################
############################# Tables 2 & 3: Max Bias for Parameter Estimates ###############################
#############################################################################################################

T3A = Batch.Effect[, c("Mx_Bias1_1", "Mx_Bias2_1","Mx_Bias3_1","Mx_Bias4_1","Mx_Bias5_1",
                       "Mx_Bias6_1","Mx_Bias7_1", "Mx_Bias8_1", "Mx_Bias9_1")]


T3A.M = data.frame(t(apply(T3A, 2, Function.1)))  
colnames(T3A.M) = colnames.1
T3A.M$Label = row.names.2
T3A.M$Variable = "Case"
T3A.M$Condition = "Batch Effect"

T3C = Batch.Effect[, c("Mx_Bias1_2", "Mx_Bias2_2","Mx_Bias3_2","Mx_Bias4_2","Mx_Bias5_2",
                       "Mx_Bias6_2","Mx_Bias7_2", "Mx_Bias8_2", "Mx_Bias9_2")]

T3C.M = data.frame(t(apply(T3C, 2, Function.1)))  
colnames(T3C.M) = colnames.1
T3C.M$Label = row.names.2
T3C.M$Variable = "Age"
T3C.M$Condition = "Batch Effect"


T3E = Batch.Effect[, c("Mx_Bias1_3", "Mx_Bias2_3","Mx_Bias3_3","Mx_Bias4_3","Mx_Bias5_3",
                       "Mx_Bias6_3","Mx_Bias7_3", "Mx_Bias8_3", "Mx_Bias9_3")]

T3E.M = data.frame(t(apply(T3E, 2, Function.1)))  
colnames(T3E.M) = colnames.1
T3E.M$Label = row.names.2
T3E.M$Variable = "HBA1C"
T3E.M$Condition = "Batch Effect"


T3 = rbind(T3A.M, T3C.M, T3E.M)
T3$Range = T3$Max-T3$Min

write.csv(T3, 
          file = "Path:/Output/T2_T3_Maximum Bias_All Probes_GitHub.csv",
          row.names = FALSE,
          na = "")


##################################################################################################
################# Table 4: bias under the alternative hypothesis #################################
##################################################################################################

Abs.Alt1_1 = abs(Batch.Effect$Alt1_1)
Abs.Alt1_2 = abs(Batch.Effect$Alt1_2)
Abs.Alt1_3 = abs(Batch.Effect$Alt1_3)
Abs.Alt1_4 = abs(Batch.Effect$Alt1_4)
Abs.Alt1_5 = abs(Batch.Effect$Alt1_5)
Abs.Alt1_6 = abs(Batch.Effect$Alt1_6)

Abs.Alt7_1 = abs(Batch.Effect$Alt7_1)
Abs.Alt8_1 = abs(Batch.Effect$Alt8_1)
Abs.Alt9_1 = abs(Batch.Effect$Alt9_1)


Abs.T2F =  data.frame(cbind(
  Abs.Alt1_1,
  Abs.Alt1_2,
  Abs.Alt1_3,
  Abs.Alt1_4,
  Abs.Alt1_5,
  Abs.Alt1_6,
  Abs.Alt7_1,
  Abs.Alt8_1,
  Abs.Alt9_1))

Abs.T2F.M = data.frame(t(apply(Abs.T2F, 2, Function.1)))  
colnames(Abs.T2F.M) = colnames.1
Abs.T2F.M$Label = row.names.2
Abs.T2F.M$Variable = "Case"
Abs.T2F.M$Condition = "Batch Effect"


Abs.Alt2_1 = abs(Batch.Effect$Alt2_1)
Abs.Alt2_2 = abs(Batch.Effect$Alt2_2)
Abs.Alt2_3 = abs(Batch.Effect$Alt2_3)
Abs.Alt2_4 = abs(Batch.Effect$Alt2_4)
Abs.Alt2_5 = abs(Batch.Effect$Alt2_5)
Abs.Alt2_6 = abs(Batch.Effect$Alt2_6)

Abs.Alt7_2 = abs(Batch.Effect$Alt7_2)
Abs.Alt8_2 = abs(Batch.Effect$Alt8_2)
Abs.Alt9_2 = abs(Batch.Effect$Alt9_2)


Abs.T2G =  data.frame(cbind(
  Abs.Alt2_1,
  Abs.Alt2_2,
  Abs.Alt2_3,
  Abs.Alt2_4,
  Abs.Alt2_5,
  Abs.Alt2_6,
  Abs.Alt7_2,
  Abs.Alt8_2,
  Abs.Alt9_2))


Abs.T2G.M = data.frame(t(apply(Abs.T2G, 2, Function.1)))  
colnames(Abs.T2G.M) = colnames.1
Abs.T2G.M$Label = row.names.2
Abs.T2G.M$Variable = "Age"
Abs.T2G.M$Condition = "Batch Effect"


#HbA1c
Abs.Alt3_1 = abs(Batch.Effect$Alt3_1)
Abs.Alt3_2 = abs(Batch.Effect$Alt3_2)
Abs.Alt3_3 = abs(Batch.Effect$Alt3_3)
Abs.Alt3_4 = abs(Batch.Effect$Alt3_4)
Abs.Alt3_5 = abs(Batch.Effect$Alt3_5)
Abs.Alt3_6 = abs(Batch.Effect$Alt3_6)

Abs.Alt7_3 = abs(Batch.Effect$Alt7_3)
Abs.Alt8_3 = abs(Batch.Effect$Alt8_3)
Abs.Alt9_3 = abs(Batch.Effect$Alt9_3)


Abs.T2H =  data.frame(cbind(
  Abs.Alt3_1,
  Abs.Alt3_2,
  Abs.Alt3_3,
  Abs.Alt3_4,
  Abs.Alt3_5,
  Abs.Alt3_6,
  Abs.Alt7_3,
  Abs.Alt8_3,
  Abs.Alt9_3))



Abs.T2H.M = data.frame(t(apply(Abs.T2H, 2, Function.1)))  
colnames(Abs.T2H.M) = colnames.1
Abs.T2H.M$Label = row.names.2
Abs.T2H.M$Variable = "HBA1C"
Abs.T2H.M$Condition = "Batch Effect"


Abs.T2B = rbind(Abs.T2F.M, Abs.T2G.M, Abs.T2H.M)
Abs.T2B$Range = Abs.T2B$Max-Abs.T2B$Min


write.csv(Abs.T2B, 
          file = "Path:/Output/T4_Average Absolute Bias_Significant_Probe_2x_GitHub.csv",
          row.names = FALSE,
          na = "")


#############################################################################################################
############################ Table 5: Mean SE for Parameter Estimates ######################################
#############################################################################################################


T4A = Batch.Effect[, c("Av_SE_1_True",
                       "Av_SE_1_R",
                       "Av_SE_1_SR",
                       "Av_SE_1_Optimal",
                       "Av_SE_1_Corrected_R",
                       "Av_SE_1_Corrected_SR",
                       "Av_SE_1_Corrected_Optimal",
                       "Av_SE_1_FE_R",
                       "Av_SE_1_FE_SR",
                       "Av_SE_1_FE_Optimal")]


T4A.M = data.frame(t(apply(T4A, 2, Function.1)))  
colnames(T4A.M) = colnames.1
T4A.M$Label = row.names.1
T4A.M$Variable = "Case"
T4A.M$Condition = "Batch Effect"




T4C = Batch.Effect[, c("Av_SE_2_True",
                       "Av_SE_2_R",
                       "Av_SE_2_SR",
                       "Av_SE_2_Optimal",
                       "Av_SE_2_Corrected_R",
                       "Av_SE_2_Corrected_SR",
                       "Av_SE_2_Corrected_Optimal",
                       "Av_SE_2_FE_R",
                       "Av_SE_2_FE_SR",
                       "Av_SE_2_FE_Optimal")]

T4C.M = data.frame(t(apply(T4C, 2, Function.1)))  
colnames(T4C.M) = colnames.1
T4C.M$Label = row.names.1
T4C.M$Variable = "Age"
T4C.M$Condition = "Batch Effect"


T4E = Batch.Effect[, c("Av_SE_3_True",
                       "Av_SE_3_R",
                       "Av_SE_3_SR",
                       "Av_SE_3_Optimal",
                       "Av_SE_3_Corrected_R",
                       "Av_SE_3_Corrected_SR",
                       "Av_SE_3_Corrected_Optimal",
                       "Av_SE_3_FE_R",
                       "Av_SE_3_FE_SR",
                       "Av_SE_3_FE_Optimal")]

T4E.M = data.frame(t(apply(T4E, 2, Function.1)))  
colnames(T4E.M) = colnames.1
T4E.M$Label = row.names.1
T4E.M$Variable = "HBA1C"
T4E.M$Condition = "Batch Effect"


T4 = rbind(T4A.M, T4C.M, T4E.M)
T4$Range = T4$Max - T4$Min

write.csv(T4, 
          file = "Path:/Output/T5_Average SE_2x_GitHub.csv",
          row.names = FALSE,
          na = "")
