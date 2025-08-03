setwd("E:/bioinformatics diploma/integ/Project/prooject refined 1")
if (!requireNamespace("limma")) install.packages("limma")
if (!requireNamespace("pheatmap")) install.packages("pheatmap")
library(limma)
library(pheatmap)

# Load expression data (genes x samples)
expr <- read.table("exp",
                   header = TRUE, row.names = 1,
                   check.names = FALSE, quote = "\"", sep = "",
                   stringsAsFactors = FALSE)
colnames(expr) <- toupper(gsub("\\.", "-", colnames(expr)))

# Load survival data
survival <- read.table("survival",
                       header = TRUE, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE)
survival$PatientID <- toupper(gsub("\\.", "-", survival$PatientID))
survival_nodup <- survival[!duplicated(survival$PatientID), ]

# Load clinical data (for later/optional use)
breast <- read.table("breast",
                     header = TRUE, sep = "\t", row.names = 1,
                     check.names = FALSE, stringsAsFactors = FALSE)
rownames(breast) <- toupper(rownames(breast))

# Helper function to extract patient code
get_patient_id <- function(x) sapply(strsplit(x, "-"), function(parts) paste(parts[1:3], collapse = "-"))

expr_patients   <- get_patient_id(colnames(expr))
surv_patients   <- survival_nodup$PatientID
breast_patients <- get_patient_id(rownames(breast))

# Intersect patients present in all datasets
patients_to_keep <- Reduce(intersect, list(expr_patients, surv_patients, breast_patients))

# Subset and re-order columns/rows to exactly match patients_to_keep
# One expression sample per patient!
keep_expr_idx    <- which(expr_patients %in% patients_to_keep & !duplicated(expr_patients))
ordered_expr_idx <- order(match(get_patient_id(colnames(expr)[keep_expr_idx]), patients_to_keep))
expr_matched     <- expr[, keep_expr_idx[ordered_expr_idx]]

breast_matched   <- breast[match(patients_to_keep, breast_patients), ]
survival_matched <- survival_nodup[match(patients_to_keep, surv_patients), ]

# Checks: All orders must match
stopifnot(all(get_patient_id(colnames(expr_matched)) == patients_to_keep))
stopifnot(all(get_patient_id(rownames(breast_matched)) == patients_to_keep))
stopifnot(all(survival_matched$PatientID == patients_to_keep))

# Keep samples with complete survival data
keep <- complete.cases(survival_matched)
expr_final     <- expr_matched[, keep]
breast_final   <- breast_matched[keep, ]
survival_final <- survival_matched[keep, ]

# Initial Plots
log_expr <- log2(as.matrix(expr_final) + 1)
hist(as.numeric(log_expr), breaks=50, col="orange", main="Expression Value Distribution")
plot(density(as.numeric(log_expr)), main="Density of Log2-Transformed Expression")

# Remove zero-variance genes
var_genes <- apply(log_expr, 1, function(x) sd(x, na.rm=TRUE)) > 0
log_expr_filtered <- log_expr[var_genes, ]

# PCA
pca_out <- prcomp(t(log_expr_filtered), scale.=TRUE)

group_color <- ifelse(survival_final$Death == 1, "red", "blue")
plot(pca_out$x[, 1], pca_out$x[, 2], col=group_color, pch=19,
     xlab="PC1", ylab="PC2", main="PCA: Survival Status (Red=Deceased, Blue=Alive)")
legend("topright", legend=c("Deceased","Alive"), col=c("red","blue"), pch=19)

# Differential Expression - Fix logFC Bug
group <- factor(survival_final$Death, levels=c(0,1), labels=c("Alive","Deceased"))
design <- model.matrix(~group)
fit <- lmFit(log_expr_filtered, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust.method="BH", number=Inf)
sig_genes_1 <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5, ]
write.csv(sig_genes_1, "DE_genes_deceased_vs_alive.csv")

# Volcano Plot (improved coloring)
plot(results$logFC, -log10(results$adj.P.Val),
     pch=20, cex=0.6, col=ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5, "red", "grey"),
     xlab="log2 Fold Change", ylab="-log10(adjusted p-value)",
     main="Volcano Plot: Deceased vs Alive")
abline(h=-log10(0.05), col="blue", lty=2)
abline(v=c(-0.5,0.5), col="blue", lty=2)

# Heatmap of Top Genes (Handles <50 Sig Genes)
top_n <- min(30, nrow(sig_genes_1))
top_gene_names <- rownames(sig_genes_1)[order(sig_genes_1$adj.P.Val)][1:top_n]
heat_expr <- log_expr_filtered[top_gene_names, ]

annotation_col <- data.frame(SurvivalStatus = group)
rownames(annotation_col) <- colnames(log_expr_filtered)   # Make sure rownames match columns

pheatmap(
  heat_expr,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  main = paste("Top", top_n, "DE Genes: Deceased (Red) vs Alive (Blue)")
)
########################################################################
#DE_Analysis_After_Counfounders

# First, identify numeric columns in your clinical data
is_num <- sapply(breast_final, is.numeric) | sapply(breast_final, is.integer)

# Get all numeric clinical variable names
num_vars <- names(breast_final)[is_num]

# For each numeric variable, correlate with PC1 and PC2
cat("Correlation with PC1 and PC2 for numeric variables:\n")
for(v in num_vars){
  var <- breast_final[[v]]
  # Remove missing in both var and PC
  comp <- complete.cases(pca_out$x[,1], var)
  if(sum(comp) > 4){   # needs at least 5 samples
    cor1 <- cor(pca_out$x[comp,1], var[comp])
    cor2 <- cor(pca_out$x[comp,2], var[comp])
    cat(sprintf("%-50s PC1: %6.3f PC2: %6.3f\n", v, cor1, cor2))
  }
}

####################
# Find 'factor-like' (categorical) columns: <=12 unique values, not all numeric, not NA
is_categ <- sapply(breast_final, function(x) length(unique(x[!is.na(x)])) <= 12 & !is.numeric(x))
cat_vars <- names(breast_final)[is_categ]

cat("\nANOVA p-values with PC1 and PC2 for categorical variables:\n")
for(v in cat_vars){
  f <- as.factor(breast_final[[v]])
  if (nlevels(f) > 1) {
    comp <- complete.cases(pca_out$x[,1], f)
    if(sum(comp) > 4){
      p1 <- summary(aov(pca_out$x[comp,1] ~ f[comp]))[[1]][1,"Pr(>F)"]
      p2 <- summary(aov(pca_out$x[comp,2] ~ f[comp]))[[1]][1,"Pr(>F)"]
      cat(sprintf("%-50s PC1: %8.4g  PC2: %8.4g\n", v, p1, p2))
    }
  }
}
#################
# Prepare Confounders
age   <- as.numeric(breast_final$Age_at_Initial_Pathologic_Diagnosis_nature2012)
stage <- as.factor(breast_final$AJCC_Stage_nature2012)
pam50 <- as.factor(breast_final$PAM50Call_RNAseq)
group <- factor(survival_final$Death, levels = c(0,1), labels = c("Alive", "Deceased"))

summary(age)
table(stage, useNA="ifany")
table(pam50, useNA="ifany")
# Keep only samples with all needed variables
conf_cc <- complete.cases(group, age, stage, pam50)
cat("Number of samples with all DE confounders present:", sum(conf_cc), "\n")

# Subset ALL matrices
log_expr_use <- log_expr_filtered[, conf_cc]
group_use    <- group[conf_cc]
age_use      <- age[conf_cc]
stage_use    <- stage[conf_cc]
pam50_use    <- pam50[conf_cc]

# Build design matrix for limma
design <- model.matrix(~ group_use + age_use + stage_use + pam50_use)

# Run limma DE analysis
fit <- lmFit(log_expr_use, design)
fit <- eBayes(fit)
results <- topTable(fit, coef="group_useDeceased", adjust.method="BH", number=Inf)
sig_genes_2 <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5, ]

# Output
write.csv(sig_genes_2, "DE_genes_adjusted_mainconfounders.csv")
#########################
# Visualize main confounders on PCA
cols <- as.factor(age[conf_cc])
plot(pca_out$x[conf_cc,1], pca_out$x[conf_cc,2], col=cols, pch=19,
     xlab="PC1", ylab="PC2", main="PC1/PC2 colored by Stage")
legend("topright", legend=levels(cols), col=1:length(levels(cols)), pch=19)
dim(sig_genes_2)
##############

# Read both gene lists
simple <- read.csv("DE_genes_deceased_vs_alive.csv", row.names=1)
adjusted <- read.csv("DE_genes_adjusted_mainconfounders.csv", row.names=1)

# Venn diagram!
install.packages("VennDiagram")
library(VennDiagram)
venn.plot <- venn.diagram(list(Simple=rownames(simple), Adjusted=rownames(adjusted)),
                          NULL, fill=c("lightblue", "pink"), alpha=0.5, cex=2)
grid.draw(venn.plot)
############################3
# Heatmap for adjusted results
pheatmap(
  heat_expr,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(30),
  main = paste("Top", top_n, "DE Genes (Adjusted): Deceased (Red) vs Alive (Blue)")
)

# Volcano plot for adjusted results
plot(results$logFC, -log10(results$adj.P.Val),
     pch=20, cex=0.6, col=ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5, "red", "grey"),
     xlab="log2 Fold Change", ylab="-log10(adjusted p-value)",
     main="Volcano Plot (Confounder-Adjusted): Deceased vs Alive")
abline(h=-log10(0.05), col="blue", lty=2)
abline(v=c(-0.5,0.5), col="blue", lty=2)




nrow(sig_genes_1)
nrow(sig_genes_2)
head(sig_genes_1)
head(sig_genes_2)

pca_out <- prcomp(t(log_expr_filtered), scale.=TRUE)
summary(pca_out)
var_explained <- cols$sdev^2 / sum(cols$sdev^2)
round(var_explained[1:2] * 100, 2)
summary(pca_out)$importance[2, 1:2]


library(ggplot2)

# Get percent variance
var_explained <- pca_out$sdev^2 / sum(pca_out$sdev^2)
pc1_var <- round(var_explained[1] * 100, 1)
pc2_var <- round(var_explained[2] * 100, 1)

# Create data frame for ggplot
df_pca <- data.frame(PC1 = pca_out$x[,1],
                     PC2 = pca_out$x[,2],
                     PAM50 = pam50)

ggplot(df_pca, aes(x = PC1, y = PC2, color = PAM50)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = "PCA Colored by PAM50",
       x = paste0("PC1 (", pc1_var, "%)"),
       y = paste0("PC2 (", pc2_var, "%)")) +
  theme_minimal()

