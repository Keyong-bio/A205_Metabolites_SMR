#### Double-Check #######

load(file = "CLSA_to_NAFLD_1e6_res_20240610.Rdata")
load(file = "UK_to_NAFLD_1e6_res_0610.Rdata")

load(file = "CLSA_to_NAFLD_5e8_res_0610.Rdata")
load(file = "UK_to_NAFLD_5e8_res0610.Rdata")


X1 <- CLSA_to_NAFLD_1e6$CLSA_to_MVP_mr_res_1e6_correct$exposure
X2 <- CLSA_to_NAFLD_5e8$CLSA_to_MVP_mr_res_correct$exposure
X3 <- UK_to_NAFLD_1e6$UK_to_MVP_mr_res_1e6_full_correct$exposure
X4 <- UK_to_NAFLD_5e8$UK_to_Biopsy_mr_res_5e8_full_correct$exposure
X5 <- UK_to_NAFLD_5e8$UK_to_MVP_mr_res_5e8_full_correct$exposure
X6 <- UK_to_NAFLD_5e8$UK_to_EHR_mr_res_5e8_full_correct$exposure

unique(c(X1, X2, X3, X4, X5, X6))
