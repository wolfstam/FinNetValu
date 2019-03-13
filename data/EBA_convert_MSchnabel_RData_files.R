library(tidyverse)
library(stringr)

load("/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/MSchnabel_preprocessing/Rworkspace_before_initial_shock.RData")

write_csv(cont.Pi, "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/MSchnabel_preprocessing/Cont_Pi_MSchnabel.csv")
write_csv(cont.Theta, "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/MSchnabel_preprocessing/Cont_Theta_MSchnabel.csv")
write_csv(as.data.frame(cap), "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/MSchnabel_preprocessing/Cont_C_MSchnabel.csv")
#write_csv(as.data.frame(lambda), "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/MSchnabel_preprocessing/lambda_MSchnabel.csv")
#write_csv(as.data.frame(lambda_target), "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/MSchnabel_preprocessing/lambda_target_MSchnabel.csv")
