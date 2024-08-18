## CETP -> Metabolites #####

## load(file = "CETP_mr.Rdata")

library(pacman)
pacman::p_load(data.table, ieugwasr, openxlsx, dplyr,
               compareGroups, ggplot2, TwoSampleMR, LDlinkR,
               # phenoscanner,
               reshape2, # melt dataset
               ggcorrplot, # create the correlation matrix
               ggforestplot, # For forest plots of metabolites
               broom, # for nice output
               glmnet, # LASSO regression
               tableone, # Descriptive 
               # plinkbinr,
               stringr,
               tidyr, gwasrapidd)

## Read the CETP GWAS data

CETP <- fread("GWAS_NAFLD/dt_CETP_GWAS.txt", header = T)
head(CETP); dim(CETP)

CETP <- CETP[grep("^rs", CETP$SNPID),]
summary(CETP$PVALUE)

CETP_sig <- subset(CETP, CETP$PVALUE < 5E-8) %>% as.data.frame() %>% 
  mutate(Pheno = "CETP") %>% 
  rename(SNP = SNPID)

## outcome data
CLSA_outcome <- fread("GWAS_Metabolites/GWAS_Met_CLSA/GCST90200060_1_1_enyl_palmitoyl_2_linoleoyl_GPC_p_16_0_18_2_.gz") %>%  as.data.frame() %>% 
  mutate(samplesize = 8260) %>% 
  relocate("variant_id")

colnames(CLSA_outcome)<- c("SNP", "chrom","position","A1", "A2", "EAF", "BETA", "SE", "pvalue", "sample_size") 
# colnames(CLSA_outcome)[1] <- "SNP"
# 
# CETP_CLSA_outcome <- merge(CETP_sig, CLSA_outcome, by = "SNP", all.x = T)
# dim(CETP_CLSA_outcome)
# 
# colnames(CETP_CLSA_outcome)[4:5] = c("CETP_effectallel","CETP_noeffectallel")
# colnames(CETP_CLSA_outcome)[7] = c("CETP_eaf")
# colnames(CETP_CLSA_outcome)[c(8,9,10)] = c("CETP_beta","CETP_se","CETP_pval")
# 
# colnames(CETP_CLSA_outcome)[16:17] = c("CLSA_effectallel","CLSA_noeffectallel")
# colnames(CETP_CLSA_outcome)[18] = c("CLSA_eaf")
# colnames(CETP_CLSA_outcome)[c(19,20,21)] = c("CLSA_beta","CLSA_se","CLSA_pval")
# 
# table(toupper(CETP_CLSA_outcome$CETP_effectallel) == CETP_CLSA_outcome$CLSA_effectallel) # 253 TRUE
# table(toupper(CETP_CLSA_outcome$CETP_effectallel) == CETP_CLSA_outcome$CLSA_noeffectallel)
# 
# inconsistent=which(is.na(CETP_CLSA_outcome$CETP_effectallel != CETP_CLSA_outcome$CLSA_effectallel) 
#                    & is.na(CETP_CLSA_outcome$CETP_effectallel != CETP_CLSA_outcome$CLSA_noeffectallel))
# 
# CETP_CLSA_outcome = CETP_CLSA_outcome[-inconsistent,]
# 
# CETP_CLSA_outcome$CLSA_beta = ifelse(CETP_CLSA_outcome$CETP_effectallel == CETP_CLSA_outcome$CLSA_effectallel, 
#                                      CETP_CLSA_outcome$CLSA_beta, -1*CETP_CLSA_outcome$CLSA_beta)
# 
# CETP_CLSA_outcome$CLSA_eaf = ifelse(CETP_CLSA_outcome$CETP_effectallel == CETP_CLSA_outcome$CLSA_effectallel, 
#                                     CETP_CLSA_outcome$CLSA_eaf, 1-CETP_CLSA_outcome$CETP_eaf)
# 
# colnames(CETP_CLSA_outcome)[1] <- "rsid"
# colnames(CETP_CLSA_outcome)[10] <- "pval"
# 
# CETP_CLSA_outcome_clp <- ld_clump(dat = CETP_CLSA_outcome, 
#                                   clump_kb = 10000, ## 10000kb
#                                   clump_r2 = 0.001, 
#                                   plink_bin = genetics.binaRies::get_plink_binary(),
#                                   bfile = "./EUR/EUR")
# 
# rs = CETP_CLSA_outcome_clp$rsid
# CETP_beta = CETP_CLSA_outcome_clp$CETP_beta
# CETP_se = CETP_CLSA_outcome_clp$CETP_se
# 
# CLSA_beta = CETP_CLSA_outcome_clp$CLSA_beta
# CLSA_se = CETP_CLSA_outcome_clp$CLSA_se
# 
# library(MendelianRandomization)
# mr.input = mr_input(bx = CETP_beta, bxse = CETP_se, by = CLSA_beta, byse = CLSA_se, 
#                     exposure = "CETP", outcome = "CLSA_meta", snps = rs)
# mr_ivw(mr.input)
# 
# png(filename = "2-sample plots/CETP-Meta.png")
# mr_plot(mr.input, interactive=FALSE)
# dev.off()
# 
# mr_allmethods(mr.input)
# 
# png(filename = "2-sample plots/CETP-Meta-all.png")
# mr_plot(mr_allmethods(mr.input, method = "main"))
# dev.off()

## CETP Exposure
IV_CETP <- TwoSampleMR::format_data(CETP_sig, # should be data.frame 
                                          type= "exposure",
                                          phenotype_col = "Pheno",
                                          snp_col = "SNP",
                                          beta_col = "BETA",
                                          se_col = "SE",
                                          eaf_col = "FRQ_A1",
                                          pval_col = "PVALUE",
                                          effect_allele_col = "A1",
                                          other_allele_col = "A0",
                                          samplesize_col = "N",
                                          chr_col = "CHR",
                                          pos_col = "POS") 
# table(IV_CETP$mr_keep.exposure)

## Clumping step
IV_CETP <- IV_CETP %>% rename(rsid = SNP, pval = pval.exposure)
  
IV_CETP_clp <- ld_clump(dat = IV_CETP, clump_kb = 10000, ## 10000kb
         clump_r2 = 0.001, 
         plink_bin = genetics.binaRies::get_plink_binary(),
         bfile = "./EUR/EUR")

IV_CETP_clp <- IV_CETP_clp %>% rename(SNP = rsid, pval.exposure = pval)


source("Functions.R")

IV_CETP_clp$SNP %in% CLSA_outcome$SNP

find_missing_SNP(exp_dat = IV_CETP_clp, out_dat = CLSA_outcome) # -> 0

outcome_var <- format_data(CLSA_outcome,
                           type = "outcome",
                           snps = IV_CETP_clp$SNP,
                           snp_col = "SNP",
                           # phenotype_col = "PHENO",
                           beta_col = "BETA", 
                           se_col = "SE", 
                           effect_allele_col = "A1", 
                           other_allele_col = "A2",
                           chr_col = "chrom",
                           pos_col = "position",
                           # eaf_col = "maf", 
                           pval_col = "pvalue",
                           # ncase_col = "cases",
                           # ncontrol_col = "controls",
                           samplesize_col = "sample_size")
# head(outcome_var)

#### filter out those SNP associated with the pval.outcome < 0.05
outcome_var <- subset(outcome_var, pval.outcome > 0.05)

dat_harm = harmonise_data(IV_CETP_clp, outcome_var)
# dat_harm = steiger_filtering(dat_harm)
# dat_harm = dat_harm[dat_harm$mr_keep == TRUE,]

# directionality_test(dat_harm)

dat_harm = fstat_fun(dat_harm)
res_1 <- mr(dat_harm, method_list = c("mr_ivw_mre", "mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
# res <- generate_odds_ratios(res)
res_1$Direction <- "CETP->CLSA_Met"
res_1

temp <- mr_egger(mr_input(bx = dat_harm$beta.exposure, bxse = dat_harm$se.exposure, by = dat_harm$beta.outcome, byse = dat_harm$se.outcome)) ## MR-Egger Intercept

mr_singlesnp(dat_harm)

png(filename = "singSNP.png")
mr_forest_plot(mr_singlesnp(dat_harm))
dev.off()

## Metabolites -> CETP level 
CLSA_exposure <- fread("GWAS_Metabolites/GWAS_Met_CLSA/GCST90200060_1_1_enyl_palmitoyl_2_linoleoyl_GPC_p_16_0_18_2_.gz") %>%  as.data.frame() %>% 
  mutate(samplesize = 8260) %>% 
  relocate("variant_id", .before = 1) %>% 
  filter(p_value < 5e-8)

CLSA_exposure <- TwoSampleMR::format_data(CLSA_exposure, # should be data.frame 
                                          type= "exposure",
                                          # phenotype_col = "PHENO",
                                          snp_col = "variant_id",
                                          beta_col = "beta",
                                          se_col = "standard_error",
                                          eaf_col = "effect_allele_frequency",
                                          pval_col = "p_value",
                                          effect_allele_col = "effect_allele",
                                          other_allele_col = "other_allele",
                                          samplesize_col = "Samplesize",
                                          chr_col = "chromosome",
                                          pos_col = "base_pair_location") 

CLSA_exposure_clp <- clump_data(CLSA_exposure, clump_kb = 1000, clump_p2 = 0.001) # 2 variants remained

CLSA_exposure_clp$SNP %in% CETP$SNPID

snps = find_missing_SNP(exp_dat = CLSA_exposure_clp, out_dat = CETP) # ""rs1260326" missing 

CETP_2 <- prep_func(snp_proxy = snps, out_full = CETP, input = CLSA_exposure_clp)

CETP_outcome <- format_data(as.data.frame(CETP_2),
                            type = "outcome",
                            snps = CLSA_exposure_clp$SNP,
                            snp_col = "SNPID",
                            # phenotype_col = "PHENO",
                            beta_col = "BETA", 
                            se_col = "SE", 
                            effect_allele_col = "A1", 
                            other_allele_col = "A0",
                            chr_col = "CHR",
                            pos_col = "POS",
                            eaf_col = "FRQ_A1", 
                            pval_col = "PVALUE",
                            # ncase_col = "cases",
                            # ncontrol_col = "controls",
                            samplesize_col = "N")

dat_harm2 = harmonise_data(CLSA_exposure_clp, CETP_outcome)
# dat_harm2 = steiger_filtering(dat_harm2)
# dat_harm2 = dat_harm2[dat_harm2$mr_keep == TRUE & dat_harm2$pval.outcome > 0.05,]

## dat_harm2 <- fstat_fun(dat_harm2)
dat_harm2$F_stat <- (dat_harm2$beta.exposure/dat_harm2$se.exposure)^2
dat_harm2$zval_steiger <- (abs(dat_harm2$beta.exposure)-abs(dat_harm2$beta.outcome))/sqrt(dat_harm2$se.exposure**2 + dat_harm2$se.outcome**2)
dat_harm2 <- dat_harm2[dat_harm2$zval_steiger > -1.96, ] # 0 obs

res_2 <- mr(dat_harm2, method_list = c("mr_ivw_mre", "mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression", "mr_wald_ratio"))
# res_2 <- generate_odds_ratios(res_2)
res_2$Direction <- "CLSA_Met->CETP"

### INTERVAL GWAS ####
INTERVAL_dat = fread("GWAS_Metabolites/GWAS_Met_UK/INTERVAL_MRC-Epi_M52682.tbl.gz")
INTERVAL_dat1 <- subset(INTERVAL_dat, `P-value` < 1e-5)

colnames(INTERVAL_dat)
chrpos <- readRDS(file = "./GWAS_Metabolites/GWAS_Met_UK/chrpos_1000.Rdata")
dt_exom  <- read.table(file = "./GWAS_Metabolites/GWAS_Met_UK/HumanCoreExome-24v1-0_A_rsids.txt", header = T)

INTERVAL_dat1 <- separate(INTERVAL_dat1, col = MarkerName, into = c("CHR", "BP", "A1", "A2"), sep = ":") %>% 
  separate(col = CHR, into = c("X", "CHR"), sep = "chr") %>% 
  dplyr::mutate(CHRPOS = paste0(CHR,":", BP)) %>% 
  dplyr::mutate(Allele1 = toupper(Allele1),
                Allele2 = toupper(Allele2),
                Metabolite = "1-(1-enyl-palmitoyl)-2-linoleoyl-GPC(P-16:0/18:2)") %>% 
  dplyr::select(CHRPOS, Allele1:Metabolite)

INTERVAL_dat1 <- merge(INTERVAL_dat1, chrpos, by.x = "CHRPOS", all.x = T) 

# INTERVAL_dat1 <- INTERVAL1_merge[!is.na(INTERVAL_dat1$SNP),]

INTERVAL_dat1 <- merge(INTERVAL_dat1, dt_exom, by.x = "SNP", by.y = "Name", all.x = T)

INTERVAL_dat1$SNP <- ifelse(is.na(INTERVAL_dat1$RsID), INTERVAL_dat1$SNP, INTERVAL_dat1$RsID)

INTERVAL_dat1 <- INTERVAL_dat1[grep("^rs", INTERVAL_dat1$SNP),]

source("Functions.R")
snps = find_missing_SNP(exp_dat = IV_CETP_clp, out_dat = INTERVAL_dat1) # "rs1436424" "rs158482" is missing 
INTERVAL_dat1 <- prep_func_2(snp_proxy = snps, out_full = INTERVAL_dat1, input = IV_CETP_clp)

INTERVAL_dat2 <- separate(INTERVAL_dat1, col = CHRPOS, into = c("CHR", "BP"), sep = ":")

INTERVAL_dat2 <- format_data(INTERVAL_dat2,
                           type = "outcome",
                           snps = IV_CETP_clp$SNP,
                           snp_col = "SNP",
                           phenotype_col = "Metabolite",
                           beta_col = "Effect", 
                           se_col = "StdErr", 
                           effect_allele_col = "Allele1", 
                           other_allele_col = "Allele2",
                           chr_col = "CHR",
                           pos_col = "BP",
                           eaf_col = "Freq1", 
                           pval_col = "P-value",
                           # ncase_col = "cases",
                           # ncontrol_col = "controls",
                           samplesize_col = "N")

dat_harm3 = harmonise_data(IV_CETP_clp, INTERVAL_dat2)
dat_harm3 = fstat_fun(dat_harm3)

res_3 <- mr(dat_harm3, method_list = c("mr_ivw_mre", "mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression", "mr_wald_ratio"))
res_3$Direction <- "CETP->INTERVAL_Met"
res_3

### INTERVAL -> CETP #######
INTERVAL_dat_exp <- subset(INTERVAL_dat1, `P-value` < 5e-8)

INTERVAL_dat_exp <- separate(INTERVAL_dat_exp, col = CHRPOS, into = c("CHR", "BP"), sep = ":")

INTERVAL_dat_exp <- TwoSampleMR::format_data(INTERVAL_dat_exp, # should be data.frame 
                                    type= "exposure",
                                    phenotype_col = "Metabolite",
                                    snp_col = "SNP",
                                    beta_col = "Effect",
                                    se_col = "StdErr",
                                    eaf_col = "Freq1",
                                    pval_col = "P-value",
                                    effect_allele_col = "Allele1",
                                    other_allele_col = "Allele2",
                                    samplesize_col = "N",
                                    chr_col = "CHR",
                                    pos_col = "BP") 

INTERVAL_dat_exp_clp <- clump_data(INTERVAL_dat_exp, clump_kb = 1000, clump_p2 = 0.001) # clump step

# INTERVAL_dat_exp_clp$SNP %in% CETP$SNPID

snps = find_missing_SNP(exp_dat = INTERVAL_dat_exp_clp, out_dat = CETP) # "rs2119690" "rs2414578" "rs597808"  are missing 
CETP_2 <- prep_func_2(snp_proxy = snps, out_full = CETP, input = INTERVAL_dat_exp_clp)

CETP_outcome <- format_data(as.data.frame(CETP_2),
                            type = "outcome",
                            snps = INTERVAL_dat_exp_clp$SNP,
                            snp_col = "SNPID",
                            # phenotype_col = "PHENO",
                            beta_col = "BETA", 
                            se_col = "SE", 
                            effect_allele_col = "A1", 
                            other_allele_col = "A0",
                            chr_col = "CHR",
                            pos_col = "POS",
                            eaf_col = "FRQ_A1", 
                            pval_col = "PVALUE",
                            # ncase_col = "cases",
                            # ncontrol_col = "controls",
                            samplesize_col = "N")

dat_harm4 = harmonise_data(INTERVAL_dat_exp_clp, CETP_outcome)
dat_harm4 = fstat_fun(dat_harm4)
res_4 <- mr(dat_harm4, method_list = c("mr_ivw_mre", "mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression", "mr_wald_ratio"))
res_4$Direction <- "INTERVAL_Met->CETP" # direction is False


#### cALT NAFLD #########
out_GWAS2 <- fread(file = "./GWAS_NAFLD/MVP_Cohort/NAFLD.EUR.MVP.2021.2.txt") %>% as.data.frame()
colnames(out_GWAS2)<- c("SNP", "chrom","position","A1", "A2", "EAF", "BETA", "SE", "pvalue",  "sample_size", "OR", "PHENO") 

snps = find_missing_SNP(exp_dat = IV_CETP_clp, out_dat = out_GWAS2) # 1 SNP is missing
out_GWAS2 <- prep_func_2(snp_proxy = snps, out_full = out_GWAS2, input = IV_CETP_clp)

out_GWAS2_1 <-  format_data(as.data.frame(out_GWAS2),
                            type = "outcome",
                            snps = IV_CETP_clp$SNP,
                            snp_col = "SNP",
                            phenotype_col = "PHENO",
                            beta_col = "BETA", 
                            se_col = "SE", 
                            effect_allele_col = "A1", 
                            other_allele_col = "A2",
                            chr_col = "chrom",
                            pos_col = "position",
                            # eaf_col = "FRQ_A1", 
                            pval_col = "pvalue",
                            # ncase_col = "cases",
                            # ncontrol_col = "controls",
                            samplesize_col = "sample_size")

CTEP_NAFLD <- harmonise_data(IV_CETP_clp, out_GWAS2_1)
dat_harm5 = fstat_fun(CTEP_NAFLD)
res_5 <- mr(dat_harm5, method_list = c("mr_ivw_mre", "mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression", "mr_wald_ratio"))
res_5 <- generate_odds_ratios(res_5)

res_5$Direction <- "CETP->NAFLD" # direction is False

### Biopsy'######
out.file = "./GWAS_NAFLD"
out_GWAS1 <- fread(file = paste0(out.file, "/Biopsy_study/01_Biopsy_NAFLD_GWAS/GCST90011885_buildGRCh37.tsv.gz")) %>% as.data.frame() %>%
  dplyr::mutate(cases = 1483,
                controls = 17781,
                sample_size = 1483+17781,
                PHENO = "NAFLD_Biopsy")

colnames(out_GWAS1)[1:11]<- c("SNP", "pvalue", "chrom", "position", "A1", "A2", "OR", "lowci", "highci", "BETA", "SE")

IV_CETP_clp$SNP %in% out_GWAS1$SNP

# snps = find_missing_SNP(exp_dat = IV_CETP_clp, out_dat = out_GWAS1) # 1 SNP is missing
out_GWAS1_1 <-  format_data(as.data.frame(out_GWAS1),
                            type = "outcome",
                            snps = IV_CETP_clp$SNP,
                            snp_col = "SNP",
                            phenotype_col = "PHENO",
                            beta_col = "BETA", 
                            se_col = "SE", 
                            effect_allele_col = "A1", 
                            other_allele_col = "A2",
                            chr_col = "chrom",
                            pos_col = "position",
                            # eaf_col = "FRQ_A1", 
                            pval_col = "pvalue",
                            ncase_col = "cases",
                            ncontrol_col = "controls",
                            samplesize_col = "sample_size")

CTEP_NAFLD_bio <- harmonise_data(IV_CETP_clp, out_GWAS1_1)
dat_harm6 = fstat_fun(CTEP_NAFLD_bio)
res_6 <- mr(dat_harm6, method_list = c("mr_ivw_mre", "mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression", "mr_wald_ratio"))
res_6 <- generate_odds_ratios(res_6)

p1 = mr_singlesnp(dat_harm6, all_method = c("mr_ivw"))
p2 <- mr_forest_plot(p1)
ggsave(p2[[1]], file = paste0("./2-sample plots/","temp.pdf"), width = 7, height = 10)

### EHR ###########

out_GWAS3 <- fread(file =  paste0(out.file, "/EHR_Cohort/02_Eletronic health record_NAFLD_GWAS/GCST90091033_buildGRCh37.tsv.gz")) %>%
  as.data.frame() %>% 
  mutate(cases = 8434,
         controls = 770180,
         sample_size = 8434+770180,
         PHENO = "Medical_Record_NAFLD")

colnames(out_GWAS3)[1:9]<- c("SNP", "chrom", "position", "A1", "A2", "EAF","BETA",  "SE" , "pvalue")

out_GWAS3 <- out_GWAS3 %>% dplyr::select(SNP, pvalue, chrom, position, A1, A2, BETA, SE, cases, controls, sample_size, PHENO)

IV_CETP_clp$SNP %in% out_GWAS3$SNP

# snps = find_missing_SNP(exp_dat = IV_CETP_clp, out_dat = out_GWAS1) # 1 SNP is missing
out_GWAS3_1 <-  format_data(as.data.frame(out_GWAS3),
                            type = "outcome",
                            snps = IV_CETP_clp$SNP,
                            snp_col = "SNP",
                            phenotype_col = "PHENO",
                            beta_col = "BETA", 
                            se_col = "SE", 
                            effect_allele_col = "A1", 
                            other_allele_col = "A2",
                            chr_col = "chrom",
                            pos_col = "position",
                            # eaf_col = "FRQ_A1", 
                            pval_col = "pvalue",
                            ncase_col = "cases",
                            ncontrol_col = "controls",
                            samplesize_col = "sample_size")

CTEP_NAFLD_EHR <- harmonise_data(IV_CETP_clp, out_GWAS3_1)
dat_harm7 = fstat_fun(CTEP_NAFLD_EHR)
res_7 <- mr(dat_harm7, method_list = c("mr_ivw_mre", "mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression", "mr_wald_ratio"))
res_7 <- generate_odds_ratios(res_7)
res_7

### meta-analysis
library(meta)
input = data.frame(b = c(res_5$b[1], res_6$b[1], res_7$b[1]),
                   se = c(res_5$se[1], res_6$se[1], res_7$se[1]),
                   stud = c("cALT_cohort", "biopsy_cohort", "Eletronic_health_record_cohort"))
meta_res <- metagen(TE = b, seTE = se, data = input, 
             sm = "OR",
             studlab = stud,
             hakn = TRUE, method.tau="DL", comb.fixed = T, comb.random = F) # fixed model

pdf(file = "2-sample plots/CETP_NAFLD.pdf", width = 10, height = 10)
forest(meta_res)
dev.off()

input2 = data.frame(b = c(res_1$b[1], res_3$b[1]),
                   se = c(res_1$se[1], res_3$se[1]),
                   stud = c("CLSA_cohort", "INTERVAL_cohort"))

meta_res_2 <- metagen(TE = b, seTE = se, data = input2,
                      studlab = stud,
                      sm = "MD",
                      hakn = FALSE, method.tau="DL", comb.fixed = T, comb.random = F) # fixed model

pdf(file = "2-sample plots/CETP_Met.pdf", width = 15, height = 15)
forest(meta_res_2)
dev.off()

############ About the outcome GWAS ######
#### 1. if the SNPs are also significantly in the outcome GWAS? ######

save.image(file = "./2-sample plots/CETP_mr.Rdata")
