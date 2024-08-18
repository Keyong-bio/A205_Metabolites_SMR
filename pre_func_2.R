prep_func_2 <- function(snp_proxy, out_full, input){
  
    for (i in 1:length(snp_proxy)) {
      
      exp.UK.clump <- input
      
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