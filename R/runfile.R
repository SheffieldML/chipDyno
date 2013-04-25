print("Data file save @ ResultsTu_10000Ita_Lovelace.RData")
print("Results save @ test_rdemTu_10000_lovelace.txt")
sink("test_rdemTu_10000_lovelace.txt")
source("demTu.R")
sink()

