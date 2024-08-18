########### Meta_analysis ##########
### 1e6 ####
### 1. CLSA #############
load(file = "CLSA_to_NAFLD_1e6_res_20240610.Rdata")
load(file = "UK_to_NAFLD_1e6_res_0610.Rdata")


## read the metabolite names for UK data
Metabo_list <- read.xlsx(xlsxFile = "Metabolite_list_Keyong_v2.xlsx", sheet = 2)
Metabo_list$MA_code <- paste0("M", Metabo_list$MB_ID)

###### 1. For all metabolites ##############
Met_1e6_CLSA <- bind_rows(CLSA_to_NAFLD_1e6[[1]],
                          CLSA_to_NAFLD_1e6[[3]],
                          CLSA_to_NAFLD_1e6[[5]]) %>% 
  mutate(exposure = tolower(exposure))

Met_1e6_UK <- bind_rows(UK_to_NAFLD_1e6[[1]],
                        UK_to_NAFLD_1e6[[3]],
                        UK_to_NAFLD_1e6[[5]]) %>% 
  merge(Metabo_list, by.x = "exposure", by.y = "MA_code", all.x = T) %>% 
  mutate(CHEMICAL_NAME = tolower(CHEMICAL_NAME)) %>% 
  relocate(CHEMICAL_NAME)

Vars1 = union(Met_1e6_CLSA$exposure, Met_1e6_UK$CHEMICAL_NAME) # 92 metabolites

Met_1e6_CLSA_select <- Met_1e6_CLSA %>% 
  filter(exposure %in% Vars1) %>% 
  dplyr::mutate(Source = "CLSA")

length(table(Met_1e6_CLSA_select$exposure)) # 81 metabolites

Met_1e6_UK_select <- Met_1e6_UK %>% 
  filter(CHEMICAL_NAME %in% Vars1)%>% 
  dplyr::mutate(Source = "UK")

Met_1e6_UK_select$exposure = Met_1e6_UK_select$CHEMICAL_NAME
Met_1e6_UK_select = Met_1e6_UK_select[,c(2:18, 20)]

# colnames(Met_1e6_CLSA_select) == colnames(Met_1e6_UK_select)

Met_1e6_res = rbind(Met_1e6_CLSA_select, Met_1e6_UK_select)

write.xlsx(Met_1e6_res, file = "Met_1e6_res.xlsx", rowNames = F)
# table(Met_1e6_res$method)

Met_1e6_res$method[Met_1e6_res$method == "Wald ratio"] <- "Inverse variance weighted"

### those Inverse variance weighted metabolites 
v1 <- unique(Met_1e6_res[Met_1e6_res$method %in%  c("Inverse variance weighted"),]$exposure)

### those only with Wald ratio metabolites 
v2 <- unique(Met_1e6_res[Met_1e6_res$method %in%  c("Wald ratio"),]$exposure)

### those only with Egger metabolites
v3 <- unique(Met_1e6_res[Met_1e6_res$method %in%  c("MR Egger"),]$exposure)

### those only with Mode based metabolites (weighted mode)
v4 <- unique(Met_1e6_res[Met_1e6_res$method %in%  c("Weighted mode"),]$exposure)

### those only with median based metabolites (Weighted median)
v5 <- unique(Met_1e6_res[Met_1e6_res$method %in%  c("Weighted median"),]$exposure)


#meta-analysis
meta_func <- function(data = Met_1e6_res, method_varname, exp_varname, 
                      out1 = "NAFLD_Biopsy", out2 = "NAFLD_MVP_cohort", out3 = "Medical_Record_NAFLD"){
  
  input <- Met_1e6_res[which(Met_1e6_res$method == method_varname),]
  input <- input[which(input$exposure == exp_varname),]
  input <- input[which(input$outcome == out1 | input$outcome == out2 | input$outcome == out3),]
  
  library(meta)
  m.gen <- metagen(TE = b,
                   seTE = se,
                   data = input,
                   sm = "OR",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "REML",
                   method.random.ci = "HK",
                   comb.fixed = F, comb.random = T)
  
  summary_m.gen  <- summary(m.gen)
  k = m.gen$k
  
  or.tibble <- as_tibble(exp(summary_m.gen$random$TE))
  ci_lower.tibble <- as_tibble(exp(summary_m.gen$random$lower))
  ci_upper.tibble <- as_tibble(exp(summary_m.gen$random$upper))
  
  or.fixed.tibble <- as_tibble(exp(summary_m.gen$fixed$TE))
  ci_lower.fixed.tibble <- as_tibble(exp(summary_m.gen$fixed$lower))
  ci_upper.fixed.tibble <- as_tibble(exp(summary_m.gen$fixed$upper))
  
  TE.tibble <- as_tibble(m.gen$TE.random)
  se.tibble <- as_tibble(m.gen$seTE.random)
  p.tibble <- as_tibble(m.gen$pval.random)
  
  TE.fixed.tibble <- as_tibble(m.gen$TE.fixed)
  se.fixed.tibble <- as_tibble(m.gen$seTE.fixed)
  p.fixed.tibble <- as_tibble(m.gen$pval.fixed)
  
  tau2.tibble <- as_tibble(m.gen$tau2)
  I2.tibble <- as_tibble(m.gen$I2)
  Q.tibble <- as_tibble(m.gen$Q)
  Q_pavlue <- as_tibble(m.gen$pval.Q)
  
  #combine tibbles and change column names
  tibble <- cbind(k, TE.tibble, se.tibble, or.tibble, ci_lower.tibble, ci_upper.tibble, p.tibble,
                  TE.fixed.tibble, se.fixed.tibble, or.fixed.tibble, ci_lower.fixed.tibble, ci_upper.fixed.tibble,
                  p.fixed.tibble, tau2.tibble, I2.tibble, Q.tibble, Q_pavlue)
  
  colnames(tibble) <- c("N", "b.random", "se.random", "OR.random", "OR_lower.random", "OR_upper.random", "pval.random", 
                        "b.fixed", "se.fixed", "OR.fixed", "OR_lower.fixed", "OR_upper.fixed", "pval.fixed", 
                        "tau2", "I2", "Q", "Q_pav")
  
  tibble$exposure <- exp_varname
  tibble$method <- method_varname
  
  return(tibble)
}

#cmeta_func(method_varname = 'Inverse variance weighted', exp_varname = v1[[1]])

## ivw
meta_ivw_1e6 <- data.frame()
for (i in 1:length(v1)) {
  x = meta_func(method_varname = 'Inverse variance weighted', exp_varname = v1[[i]])
  meta_ivw_1e6 = rbind(x, meta_ivw_1e6)
}

# meta_wald_1e6 <- data.frame()
# for (i in 1:length(v2)) {
#   x = meta_func(method_varname = 'Wald ratio', exp_varname = v2[[i]])
#   meta_wald_1e6 = rbind(x, meta_wald_1e6)
# }

meta_egger_1e6 <- data.frame()
for (i in 1:length(v3)) {
  x = meta_func(method_varname = 'MR Egger', exp_varname = v3[[i]])
  meta_egger_1e6 = rbind(x, meta_egger_1e6)
}

meta_w_mode_1e6 <- data.frame()
for (i in 1:length(v4)) {
  x = meta_func(method_varname = 'Weighted mode', exp_varname = v4[[i]])
  meta_w_mode_1e6 = rbind(x, meta_w_mode_1e6)
}

meta_w_median_1e6 <- data.frame()
for (i in 1:length(v5)) {
  x = meta_func(method_varname = 'Weighted median', exp_varname = v5[[i]])
  meta_w_median_1e6 = rbind(x, meta_w_median_1e6)
}

# meta_wald_1e6_2 = meta_wald_1e6[complete.cases(meta_wald_1e6),]

## Multiple correction for pvalue
meta_ivw_1e6$p.fixed.fdr<- p.adjust(meta_ivw_1e6$pval.fixed, method = "fdr")
meta_ivw_1e6$p.fixed.fdr.sig <- ifelse(meta_ivw_1e6$p.fixed.fdr < .05, "*", "")
table(meta_ivw_1e6$p.fixed.fdr.sig) # 21

meta_ivw_1e6$p.random.fdr<- p.adjust(meta_ivw_1e6$pval.random, method = "fdr")
meta_ivw_1e6$p.random.fdr.sig <- ifelse(meta_ivw_1e6$p.random.fdr < .05, "*", "")
table(meta_ivw_1e6$p.random.fdr.sig) # 3

write.xlsx(meta_ivw_1e6, file = "meta_ivw_1e6.xlsx", rowNames = F)

######### Plot volcano figure for fixed effect model ###########

library(ggrepel)

plot_fixed = meta_ivw_1e6 %>% 
  select(exposure, OR.fixed, p.fixed.fdr) %>% 
  dplyr::mutate(p.fixed.fdr2 = -log10(p.fixed.fdr)) 

plot_fixed$p.fixed.fdr.sig <- ifelse(plot_fixed$p.fixed.fdr < .05 & plot_fixed$OR.fixed > 1, "Yes-Increase", 
                                     ifelse(plot_fixed$p.fixed.fdr < .05 & plot_fixed$OR.fixed < 1, "Yes-Decrease", "No"))

plot_fixed$p.fixed.fdr.sig <- factor(plot_fixed$p.fixed.fdr.sig, levels = c("Yes-Decrease", "No", "Yes-Increase"))

label_metabolites <- subset(plot_fixed, p.fixed.fdr < 0.05)$exposure
plot_fixed$delabel <- ifelse(plot_fixed$exposure %in% label_metabolites, plot_fixed$exposure, NA)

theme_set(theme_bw() +
            theme(
              axis.title.x = element_text(hjust = 0.5, margin = margin(20,0,0,0), 
                                          size = rel(1.1), color = 'black'),
              axis.title.y = element_text(margin = margin(20,0,0,0), 
                                          size = rel(1.1), color = 'black'),
              legend.position = "bottom"
            ))

pdf(file = "Volcano_1e6_meta_analysis.pdf", width = 10, height = 10)
p1 <- ggplot(data = plot_fixed, aes(x = OR.fixed, y = p.fixed.fdr2, col = p.fixed.fdr.sig, label = delabel)) +
  geom_vline(xintercept =  1, col = "gray", linetype = "dashed") +
  geom_hline(yintercept =  -log10(0.05), col = "gray", linetype = "dashed") +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB","grey", "#bb0c00"),
                     labels = c("Significant-Decrease", "No-significant", "Significant-Increase")) +
  labs(color = "Group",
       x = "OR", y = expression("-log"[10]*"P_adj")) +
  ggtitle("Fixed-Effect model-relaxed pvalue (1e-6)") +
  geom_text_repel(max.overlaps = Inf)
  
p1
dev.off()

######### Plot volcano figure for random effect model ###########

plot_random = meta_ivw_1e6 %>% 
  select(exposure, OR.random, p.random.fdr) %>% 
  dplyr::mutate(p.random.fdr2 = -log10(p.random.fdr)) 

plot_random$p.random.fdr.sig <- ifelse(plot_random$p.random.fdr < .05 & plot_random$OR.random > 1, "Yes-Increase", 
                                     ifelse(plot_random$p.random.fdr < .05 & plot_random$OR.random < 1, "Yes-Decrease", "No"))

plot_random$p.random.fdr.sig <- factor(plot_random$p.random.fdr.sig, levels = c("Yes-Decrease", "No", "Yes-Increase"))

label_metabolites <- subset(plot_random, p.random.fdr < 0.05)$exposure
plot_random$delabel <- ifelse(plot_random$exposure %in% label_metabolites, plot_random$exposure, NA)

pdf(file = "Volcano_1e6_meta_analysis_random.pdf", width = 10, height = 10)
p2 <- ggplot(data = plot_random, aes(x = OR.random, y = p.random.fdr2, col = p.random.fdr.sig, label = delabel)) +
  geom_vline(xintercept =  1, col = "gray", linetype = "dashed") +
  geom_hline(yintercept =  -log10(0.05), col = "gray", linetype = "dashed") +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB","grey", "#bb0c00"),
                     labels = c("Significant-Decrease", "No-significant", "Significant-Increase")) +
  labs(color = "Group",
       x = "OR", y = expression("-log"[10]*"P_adj")) +
  ggtitle("Random-Effect model-relaxed pvalue (1e-6)") +
  geom_text_repel(max.overlaps = Inf)

p2
dev.off()

library(cowplot)
pg <- plot_grid(p1, p2 , 
                labels = c("A", "B"), nrow = 1, label_size = 12)

pdf(file = "1e6_Meta_analysis_ivw.pdf", height = 11, width = 16)
pg
dev.off()

# meta-analysis-5e-8 ------------------------------------------------------
load(file = "CLSA_to_NAFLD_5e8_res_0610.Rdata")
load(file = "UK_to_NAFLD_5e8_res0610.Rdata")

### 5e8
Met_5e8_CLSA <- bind_rows(CLSA_to_NAFLD_5e8[[1]],
                          CLSA_to_NAFLD_5e8[[3]],
                          CLSA_to_NAFLD_5e8[[5]]) %>% 
  mutate(exposure = tolower(exposure))

Met_5e8_UK <- bind_rows(UK_to_NAFLD_5e8[[1]],
                        UK_to_NAFLD_5e8[[3]],
                        UK_to_NAFLD_5e8[[5]]) %>% 
  merge(Metabo_list, by.x = "exposure", by.y = "MA_code", all.x = T) %>% 
  mutate(CHEMICAL_NAME = tolower(CHEMICAL_NAME)) %>% 
  relocate(CHEMICAL_NAME)

Vars1 = union(Met_5e8_CLSA$exposure, Met_5e8_UK$CHEMICAL_NAME) # 84 metabolites

Met_5e8_CLSA_select <- Met_5e8_CLSA %>% 
  filter(exposure %in% Vars1) %>% 
  dplyr::mutate(Source = "CLSA")

length(table(Met_5e8_CLSA_select$exposure)) # 63 metabolites

Met_5e8_UK_select <- Met_5e8_UK %>% 
  filter(CHEMICAL_NAME %in% Vars1)%>% 
  dplyr::mutate(Source = "UK")

Met_5e8_UK_select$exposure = Met_5e8_UK_select$CHEMICAL_NAME
Met_5e8_UK_select = Met_5e8_UK_select[,c(2:18, 20)]

Met_5e8_res = rbind(Met_5e8_CLSA_select, Met_5e8_UK_select)

write.xlsx(Met_5e8_res, file = "Met_5e8_res.xlsx", rowNames = F)

Met_5e8_res$method[Met_5e8_res$method == "Wald ratio"] <- "Inverse variance weighted"

### those Inverse variance weighted metabolites 
v1 <- unique(Met_5e8_res[Met_5e8_res$method %in%  c("Inverse variance weighted"),]$exposure)

### those only with Wald ratio metabolites 
# v2 <- unique(Met_1e6_res[Met_5e8_res$method %in%  c("Wald ratio"),]$exposure)

### those only with Egger metabolites
v3 <- unique(Met_5e8_res[Met_5e8_res$method %in%  c("MR Egger"),]$exposure)

### those only with Mode based metabolites (weighted mode)
v4 <- unique(Met_5e8_res[Met_5e8_res$method %in%  c("Weighted mode"),]$exposure)

### those only with median based metabolites (Weighted median)
v5 <- unique(Met_5e8_res[Met_5e8_res$method %in%  c("Weighted median"),]$exposure)

######### ivw 5e-8 ######
meta_ivw_5e8 <- data.frame()

for (i in 1:length(v1)) {
  x = meta_func(data = Met_5e8_res,method_varname = 'Inverse variance weighted', exp_varname = v1[[i]])
  meta_ivw_5e8 = rbind(x, meta_ivw_5e8)
}

meta_egger_5e8 <- data.frame()
for (i in 1:length(v3)) {
  x = meta_func(data = Met_5e8_res, method_varname = 'MR Egger', exp_varname = v3[[i]])
  meta_egger_5e8 = rbind(x, meta_egger_5e8)
}

meta_w_mode_5e8 <- data.frame()
for (i in 1:length(v4)) {
  x = meta_func(data = Met_5e8_res, method_varname = 'Weighted mode', exp_varname = v4[[i]])
  meta_w_mode_5e8 = rbind(x, meta_w_mode_5e8)
}

meta_w_median_5e8 <- data.frame()
for (i in 1:length(v5)) {
  x = meta_func(data = Met_5e8_res, method_varname = 'Weighted median', exp_varname = v5[[i]])
  meta_w_median_5e8 = rbind(x, meta_w_median_5e8)
}

######## Multiple correction #######

meta_ivw_5e8$p.fixed.fdr<- p.adjust(meta_ivw_5e8$pval.fixed, method = "fdr")
meta_ivw_5e8$p.fixed.fdr.sig <- ifelse(meta_ivw_5e8$p.fixed.fdr < .05, "*", "")
table(meta_ivw_5e8$p.fixed.fdr.sig) # 20

meta_ivw_5e8$p.random.fdr<- p.adjust(meta_ivw_5e8$pval.random, method = "fdr")
meta_ivw_5e8$p.random.fdr.sig <- ifelse(meta_ivw_5e8$p.random.fdr < .05, "*", "")
table(meta_ivw_5e8$p.random.fdr.sig) # 3

write.xlsx(meta_ivw_5e8, file = "meta_ivw_5e8.xlsx", rowNames = F)


plot_fixed = meta_ivw_5e8 %>% 
  select(exposure, OR.fixed, p.fixed.fdr) %>% 
  dplyr::mutate(p.fixed.fdr2 = -log10(p.fixed.fdr)) 

plot_fixed$p.fixed.fdr.sig <- ifelse(plot_fixed$p.fixed.fdr < .05 & plot_fixed$OR.fixed > 1, "Yes-Increase", 
                                     ifelse(plot_fixed$p.fixed.fdr < .05 & plot_fixed$OR.fixed < 1, "Yes-Decrease", "No"))

plot_fixed$p.fixed.fdr.sig <- factor(plot_fixed$p.fixed.fdr.sig, levels = c("Yes-Decrease", "No", "Yes-Increase"))

label_metabolites <- subset(plot_fixed, p.fixed.fdr < 0.05)$exposure
plot_fixed$delabel <- ifelse(plot_fixed$exposure %in% label_metabolites, plot_fixed$exposure, NA)

theme_set(theme_bw() +
            theme(
              axis.title.x = element_text(hjust = 0.5, margin = margin(20,0,0,0), 
                                          size = rel(1.1), color = 'black'),
              axis.title.y = element_text(margin = margin(20,0,0,0), 
                                          size = rel(1.1), color = 'black'),
              legend.position = "bottom"
            ))

pdf(file = "Volcano_5e8_meta_analysis.pdf", width = 10, height = 10)
p3 <- ggplot(data = plot_fixed, aes(x = OR.fixed, y = p.fixed.fdr2, col = p.fixed.fdr.sig, label = delabel)) +
  geom_vline(xintercept =  1, col = "gray", linetype = "dashed") +
  geom_hline(yintercept =  -log10(0.05), col = "gray", linetype = "dashed") +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB","grey", "#bb0c00"),
                     labels = c("Significant-Decrease", "No-significant", "Significant-Increase")) +
  labs(color = "Group",
       x = "OR", y = expression("-log"[10]*"P_adj")) +
  ggtitle("Fixed-Effect model-relaxed pvalue (5e-8)") +
  geom_text_repel(max.overlaps = Inf)

p3
dev.off()

######### Plot volcano figure for random effect model ###########

plot_random = meta_ivw_5e8 %>% 
  select(exposure, OR.random, p.random.fdr) %>% 
  dplyr::mutate(p.random.fdr2 = -log10(p.random.fdr)) 

plot_random$p.random.fdr.sig <- ifelse(plot_random$p.random.fdr < .05 & plot_random$OR.random > 1, "Yes-Increase", 
                                       ifelse(plot_random$p.random.fdr < .05 & plot_random$OR.random < 1, "Yes-Decrease", "No"))

plot_random$p.random.fdr.sig <- factor(plot_random$p.random.fdr.sig, levels = c("Yes-Decrease", "No", "Yes-Increase"))

label_metabolites <- subset(plot_random, p.random.fdr < 0.05)$exposure
plot_random$delabel <- ifelse(plot_random$exposure %in% label_metabolites, plot_random$exposure, NA)

pdf(file = "Volcano_5e8_meta_analysis_random.pdf", width = 10, height = 10)
p4 <- ggplot(data = plot_random, aes(x = OR.random, y = p.random.fdr2, col = p.random.fdr.sig, label = delabel)) +
  geom_vline(xintercept =  1, col = "gray", linetype = "dashed") +
  geom_hline(yintercept =  -log10(0.05), col = "gray", linetype = "dashed") +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB","grey", "#bb0c00"),
                     labels = c("Significant-Decrease", "No-significant", "Significant-Increase")) +
  labs(color = "Group",
       x = "OR", y = expression("-log"[10]*"P_adj")) +
  ggtitle("Random-Effect model-relaxed pvalue (5e-8)") +
  geom_text_repel(max.overlaps = Inf)

p4
dev.off()

library(cowplot)
pg <- plot_grid(p3, p4 , 
                labels = c("A", "B"), nrow = 1, label_size = 12)

pdf(file = "5e8_Meta_analysis_ivw.pdf", height = 11, width = 16)
pg
dev.off()

#########
#### We're more interested in those metabolites #####
#### 1. 1−(1−enyl−palmitoyl)−2−linoleoyl−gpc (p−16:0/18:2) ######
#### 2. n−acetyltyrosine ######









