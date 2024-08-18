library(grid)
library(forestploter)

#### 1. 1−(1−enyl−palmitoyl)−2−linoleoyl−gpc (p−16:0/18:2) ######
#### 2. n−acetyltyrosine ######

#### Make the plots for two metabolites ######
plot.f <- read.xlsx("Two_Metabolites_5e_8.xlsx", colNames = T)

plot.f$X1 <- ifelse(is.na(plot.f$N_SNP), 
                    plot.f$X1,
                    paste0("   ", plot.f$X1))

plot.f$X1 <- ifelse(is.na(plot.f$Padj), 
                    plot.f$X1,
                    paste0("   ", plot.f$X1))

# NA to blank or NA will be transformed to character
plot.f$N_SNP <- ifelse(is.na(plot.f$N_SNP ), "", plot.f$N_SNP)

plot.f$Padj <- ifelse(is.na(plot.f$Padj), "", 
                      ifelse(!is.na(plot.f$Padj), sprintf("%.3f", plot.f$Padj),NA))

plot.f$` ` <- paste(rep(" ", 10), collapse = " ")

# Create confidence interval column to display
plot.f$`OR (95% CI)` <- ifelse(is.na(plot.f$OR), "",
                               sprintf("%.2f (%.2f to %.2f)",
                                       plot.f$OR, plot.f$lowOR, plot.f$upOR))
colnames(plot.f)


tm <- forest_theme(base_size = 10,
                   refline_lwd = 1,
                   refline_lty = "dashed", 
                   refline_col = "grey20",
                   ci_col = "#762a83",
                   ci_fill = "black",
                   ci_alpha = 0.8,
                   ci_pch = 15,
                   ci_lty = 1,
                   ci_lwd = 1.5,
                   ci_Theight = 0.2,
                   footnote_cex = 0.5)

colnames(plot.f)[1] <- "Name"
plot.f$Name <- gsub("NAFLD", "MASLD", plot.f$Name)

p <- forest(plot.f[,c(1:2, 6, 7:8)],
            est = plot.f$OR,
            lower = plot.f$lowOR, 
            upper = plot.f$upOR,
            # sizes = dt$se,
            ci_column = 4,
            ref_line = 1,
            arrow_lab = c("Decreased risk", "Increased risk"),
            xlim = c(0.1, 2),
            ticks_at = c(0.5, 1, 1.5, 2),
            footnote = "Biospy: (Anstee QM et al, 2020)\ncALT: (Vujkovic M et al, 2022)\nEHR: Electronic health record-based (Ghodsian N et al, 2021)",
            theme = tm)

# Print plot
# plot(p)

# Change font face
g <- edit_plot(p, row = c(1, 9), 
               gp = gpar(fontface = "bold"))


# Get width and height
p_wh <- get_wh(plot = g, unit = "in")
png('Two_representative_metabolites_5e08.png', res = 300, width = p_wh[1] + 1, height = p_wh[2], units = "in")
g
dev.off()