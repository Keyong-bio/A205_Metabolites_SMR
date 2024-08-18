###### SNPs are not in Outcomes ###### 
find_missing_SNP <- function(exp_dat, out_dat) {
  snps_need_proxy <- setdiff(exp_dat$SNP, intersect(exp_dat$SNP, out_dat$SNP))
  return (snps_need_proxy)
}

##### Print the table #### 
kable_dt <- function(df) {
  require(kableExtra)
  df 
}


###### Preparation function to look for proxied SNPs in the outcome ##################
prep_func <- function(snp_proxy, out_full, input){
  
  for (i in 1:length(snp_proxy)) {
    
    # i = 4
    exp.UK.clump <- input[[i]]
    
    if(length(snp_proxy[[i]]) == 0) {
      
      print("All exposure IVs found in outcome GWAS.")
      
    } else {
      
      print("Some exposure IVs missing from outcome GWAS.")
      # out_fall = out_GWAS1
      
      missing_SNPs = snp_proxy[[i]]
      
      for (j in 1:length(missing_SNPs)){
        # j = 2
        
        tryCatch({    
          proxies <- LDproxy(missing_SNPs[[j]], pop = "EUR", r2d = "r2", token = "7c2357903077")
          proxies <- proxies[proxies$R2 > 0.8, ]
          proxy_present = FALSE
          
        } , error = function(e) {
          cat("An error occurred:", e$message, "\n")
        })
    
     
        tryCatch({
          if (length(proxies$RS_Number) == 0) {
            
            print(paste0("No proxy SNP available for ", missing_SNPs[[j]]))
            
          } else {
            
            for (k in 2:length(proxies$RS_Number)) {
              
              proxy_present <- proxies$RS_Number[k] %in% out_full$SNP
              
              if (proxy_present) {
                proxy_SNP = proxies$RS_Number[k]
                proxy_SNP_allele_1 = str_sub(proxies$Alleles[k], 2, 2)
                proxy_SNP_allele_2 = str_sub(proxies$Alleles[k], 4, 4)
                original_SNP_allele_1 = str_sub(proxies$Alleles[1], 2, 2)
                original_SNP_allele_2 = str_sub(proxies$Alleles[1], 4, 4)
                break
              }
            }
          }
        }, error = function(e) {
          cat("An error occurred:", e$message, "\n")
        }
        )
        
        if (proxy_present == T) {
          
          print(paste0("Proxy SNP found. ", missing_SNPs[[j]], " replaced with ", proxy_SNP))
          
          proxy_row <- out_full[1, ] 
          proxy_row$SNP = missing_SNPs[[j]] # change names to the missing SNPs, namely original SNP
          proxy_row$BETA = as.numeric(out_full[out_full$SNP == proxy_SNP, "BETA"])
          proxy_row$SE = as.numeric(out_full[out_full$SNP == proxy_SNP, "SE"])
          
          if (out_full[out_full$SNP == proxy_SNP, "A1"] == proxy_SNP_allele_1) proxy_row$A1 = original_SNP_allele_1
          if (out_full[out_full$SNP == proxy_SNP, "A1"] == proxy_SNP_allele_2) proxy_row$A1 = original_SNP_allele_2
          if (out_full[out_full$SNP == proxy_SNP, "A2"] == proxy_SNP_allele_1) proxy_row$A2 = original_SNP_allele_1
          if (out_full[out_full$SNP == proxy_SNP, "A2"] == proxy_SNP_allele_2) proxy_row$A2 = original_SNP_allele_2
          ## change other values
          
          proxy_row$pvalue = as.numeric(out_full[out_full$SNP == proxy_SNP, "pvalue"])
          proxy_row$sample_size = as.numeric(out_full[out_full$SNP == proxy_SNP, "sample_size"])
          
          if("cases" %in% colnames(out_full)) proxy_row$cases = as.numeric(out_full[out_full$SNP == proxy_SNP, "cases"])
          if("controls" %in% colnames(out_full))proxy_row$controls = as.numeric(out_full[out_full$SNP == proxy_SNP, "controls"])
          
          proxy_row$chrom = as.numeric(exp.UK.clump[exp.UK.clump$SNP == missing_SNPs[[j]], "chr.exposure"]) # SNPs
          proxy_row$position = as.numeric(exp.UK.clump[exp.UK.clump$SNP == missing_SNPs[[j]], "pos.exposure"]) # SNPs
          # if("AF1" %in% colnames(out_full)) proxy_row$eaf.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "AF1"])
          out_full <- rbind(out_full, proxy_row)
          
        }
        if(proxy_present == FALSE) {
          print(paste0("No proxy SNP available for ", missing_SNPs[[j]], " in outcome GWAS."))
        }
      }
    }
  }
  return(out_full)
}

######### out_func: format the outcome data ###############
out_func <- function(out_GWAS, snp_exposure) {
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  if (length(snp_exposure) > 1) {
    outcome_var <- format_data(out_GWAS,
                               type = "outcome",
                               snps = snp_exposure,
                               snp_col = "SNP",
                               phenotype_col = "PHENO",
                               beta_col = "BETA", 
                               se_col = "SE", 
                               effect_allele_col = "A1", 
                               other_allele_col = "A2",
                               chr_col = "chrom",
                               pos_col = "position",
                               # eaf_col = "maf", 
                               pval_col = "pvalue",
                               ncase_col = "cases",
                               ncontrol_col = "controls",
                               samplesize_col = "sample_size")
    return(outcome_var)
    
  } else if (length(snp_exposure) == 1 && snp_exposure %in% out_GWAS$SNP) {
    outcome_var <- format_data(out_GWAS,
                               type = "outcome",
                               snps = snp_exposure,
                               snp_col = "SNP",
                               phenotype_col = "PHENO",
                               beta_col = "BETA", 
                               se_col = "SE", 
                               effect_allele_col = "A1", 
                               other_allele_col = "A2",
                               chr_col = "chrom",
                               pos_col = "position",
                               # eaf_col = "maf", 
                               pval_col = "pvalue",
                               ncase_col = "cases",
                               ncontrol_col = "controls",
                               samplesize_col = "sample_size")
    return(outcome_var)
  } else {
    print("No available in Outcome")
  }
}

######## conditional F-statistic, filtering out ###############
fstat_fun <- function(snp_data){
  # snp_data <- res_df[res_df$mr_keep = T,]
  snp_data$F_stat <- (snp_data$beta.exposure/snp_data$se.exposure)^2
  snp_data$weak_snp <- ifelse(snp_data$F_stat < 10, "yes", "no")
  snp_data$maf = ifelse(snp_data$eaf.exposure > 0.5, 1 - snp_data$eaf.exposure, snp_data$eaf.exposure)
  snp_data = snp_data %>% filter(weak_snp == "no")
  
  # steiger filtering
  snp_data <- steiger_filtering(snp_data)
  snp_data <- snp_data[(snp_data$steiger_dir) | (!snp_data$steiger_dir & snp_data$steiger_pval > 0.05), ]
  # print(paste0("Number of SNPs after Steiger filter: ", nrow(dat)))
  
  return(snp_data)
}


####### mr_keyong function to get the result ##########
mr_keyong <- function(dat) {
  ## set up a list
  Met_to_NAFLD_mr <- list()
  for (i in 1:length(dat)) {
    Met_to_NAFLD_mr1 <- mr(dat[[i]], method_list = c("mr_ivw",
                                                    # "mr_ivw_mre",
                                                    "mr_ivw_fe",
                                                    "mr_weighted_mode",
                                                    "mr_simple_mode",
                                                    "mr_wald_ratio",
                                                    "mr_egger_regression",
                                                    "mr_weighted_median"))
    ### MR results to export
    mr_add_res_table <- Met_to_NAFLD_mr1 %>% select(c("exposure", "outcome", "method", "nsnp", 
                                                     "b", "se", "pval"))
    
    if(nrow(dat[[i]]) >= 2){ # If only one SNP can't compute heterogeneity
      q_stat <- mr_heterogeneity(dat[[i]])[,c("method","Q","Q_pval")]
    } else{
      q_stat <- data.frame("method" = c("MR Egger","Inverse variance weighted"), 
                           "Q" = c(NA, NA),
                           "Q_pval" = c(NA,NA))
    }
    
    mr_add_res_table <- merge(mr_add_res_table, q_stat, by="method", all.x = T)
    mr_add_res_table <- mr_add_res_table[,c("exposure", "outcome", "method", "nsnp", "b", "se", "pval","Q","Q_pval")]
    
    ### MR-Egger additional statistics
    mr_egger_res <- mr_egger_regression(b_exp = dat[[i]]$beta.exposure, 
                                        b_out = dat[[i]]$beta.outcome, 
                                        se_exp = dat[[i]]$se.exposure,
                                        se_out = dat[[i]]$se.outcome)
    
    mr_add_res_table$Egger_intercept <- NA
    mr_add_res_table$Egger_intercept[which(mr_add_res_table$method == "MR Egger")] <- mr_egger_res$b_i
    
    mr_add_res_table$E_se_i <- NA
    mr_add_res_table$E_se_i[which(mr_add_res_table$method == "MR Egger")] <- mr_egger_res$se_i
    
    mr_add_res_table$E_pval_i <- NA
    mr_add_res_table$E_pval_i[which(mr_add_res_table$method == "MR Egger")] <- mr_egger_res$pval_i
    
    Met_to_NAFLD_mr[[i]] <- generate_odds_ratios(mr_add_res_table)
  }
  
  return(Met_to_NAFLD_mr)
}

## Function to leave-one-out analysis
perform_loo_analysis <- function(dat, filenameloo){
  for (i in 1:length(dat)) {
    res_loo <- mr_leaveoneout(dat[[i]])
    p1 <- mr_leaveoneout_plot(res_loo)
    ggsave(p1[[1]], file = paste0("./2-sample plots/",i, filename_loo), width = 7, height = 10)
  }
}

## Function to plot scatter
#Create function for plotting results for different outcomes
plot_scatter <- function(results, data, base_size, legend_size) {
  
  p <- mr_scatter_plot(results, data)[[1]] +
    theme_minimal(base_size = base_size) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = legend_size),
          legend.position = "top",
          plot.background = element_rect(fill = "white", color = "white"))
  
  return(p)
}

#Make function for single-SNP analyses
perform_single_snp_analysis <- function(data, filename_singlesnp) {
  res_single <- mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
  p2 <- mr_forest_plot(res_single)
  ggsave(p2[[1]], file = paste0("./2-sample plots/",filename_singlesnp), width = 7, height = 10)
}


## MR-PRESSO
### MR-presso-deng function
round_up_1000 <- function(x) {
  return(ceiling(x / 1000) * 1000)
}

mr_presso_deng <- function(dat, EXPOSURE, OUTCOME, w = F, methods_list){
  
  if (nrow(dat) > 3) {

    ############################# RUN MR-PRESSO ##############################
    # Initialize variables
    nb_distribution <- round_up_1000(nrow(dat))  # starting value (at least upper 1000 of number of SNPs)
    max_nb_distribution <- 21000  # maximum allowed value
    warning_occurred <- TRUE  # to enter the while loop at least once
    nb_increase_dist <- 5000
    
    handle_warning <- function(w) {
      # Check if the warning message contains the expected string
      if (grepl("Outlier test unstable", conditionMessage(w))) {
        warning_occurred <<- TRUE  # set warning flag if a warning occurs (global environment)
        nb_distribution <<- nb_distribution + nb_increase_dist  # increase nb_distribution
        if (nb_distribution <= max_nb_distribution) {
          message("Warning occurred: ", conditionMessage(w), 
                  "\nRetrying with NbDistribution = ", nb_distribution, ".")
        } else {
          message("Maximum NbDistribution reached without resolving the issue. Proceeding with NbDistribution = ", nb_distribution, ".")
        }
      } else {
        # Handle other warnings if needed, or simply print them
        message("Warning occurred: ", conditionMessage(w))
      }
    }
    
    # Loop to adjust nb_distribution and re-run MR-PRESSO
    while (warning_occurred & nb_distribution <= max_nb_distribution) {
      warning_occurred <- FALSE  # reset warning flag
      
      # Try running MR-PRESSO and catch any warnings
      tryCatch({
        mr_results <- TwoSampleMR::run_mr_presso(
          dat, 
          NbDistribution = nb_distribution
        )
      },
      warning = handle_warning,
      error = function(e) {
        stop("Error occurred: ", conditionMessage(e))  # stop on errors
      })
    }
    
    ##########################################################################
    
    # Extracting main MR results
    main_results <- mr_results[[1]]$`Main MR results`
    main_results_df <- as.data.frame(main_results)
    
    # Extracting MR-PRESSO results
    global_test_pvalue <- as.numeric(sub("<", "", mr_results[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue))
    global_test_rssobs <- mr_results[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs
    distortion_test_pvalue <- mr_results[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue
    distortion_coeff <- mr_results[[1]]$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`
    outliers <- length(mr_results[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
    
    # Adding MR-PRESSO results to the main results dataframe
    main_results_df$global_test_pvalue <- global_test_pvalue
    main_results_df$global_test_rssobs <- global_test_rssobs
    main_results_df$distortion_test_pvalue <- distortion_test_pvalue
    main_results_df$distortion_coeff <- distortion_coeff
    main_results_df$nb_distribution_set <- nb_distribution
    
    # Rename column
    main_results_df <- main_results_df %>% 
      rename(exposure = "Exposure",
             method = "MR Analysis",
             b = "Causal Estimate",
             sd = "Sd",
             pval = "P-value",
             tstat = "T-stat")
    
    # Add number of snps (remove outliers)
    main_results_df$nsnp <- ifelse(main_results_df$method == "Raw", nrow(dat), nrow(dat) - outliers)
    
    # Add se
    main_results_df$se <- abs(main_results_df$b) / main_results_df$tstat
    
    # Add exposure and outcome
    main_results_df$exposure <- EXPOSURE
    main_results_df$outcome <- OUTCOME
    
    # Reordering columns
    main_results_df <- main_results_df %>% 
      select(exposure, outcome, method, nsnp, b, se, pval, everything())
    
    # Displaying the comprehensive dataframe
    print(main_results_df)
    
  } else {
    methods_list <- c("Raw","Outlier-corrected")
    main_results_df <- data.frame(
      exposure = rep(EXPOSURE, length(methods_list)),
      outcome = rep(OUTCOME, length(methods_list)),
      method = methods_list,
      nsnp = rep(nrow(dat), length(methods_list)),
      b = rep(NA, length(methods_list)),
      se = rep(NA, length(methods_list)),
      pval = rep(NA, length(methods_list)),
      sd = rep(NA, length(methods_list)),
      tstat = rep(NA, length(methods_list)),
      global_test_pvalue = rep(NA, length(methods_list)),
      global_test_rssobs = rep(NA, length(methods_list)),
      distortion_test_pvalue = rep(NA, length(methods_list)),
      distortion_coeff = rep(NA, length(methods_list)),
      nb_distribution_set = rep(NA, length(methods_list))
    )
  }
}

