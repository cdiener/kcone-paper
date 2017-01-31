#  gene_expression.R
#
#  Copyright 2015 Felipe Muñoz, Christian Diener <resendislab.inmegen.gob.mx>
#
#  MIT license. See LICENSE for more information.


library(GEOquery, quietly = T)
library(affy, quietly = T)
library(limma, quietly = T)
library(genefilter, quietly = T)
library(hgu133plus2frmavecs, quietly = T)
library(frma, quietly = T)

# Read sample list and create download directory
sample_info <- read.csv("ge_samples.csv")
dir.create("ma")

# Download raw data for samples that don't exist yet
already_there <- dir.exists(paste0("ma/", as.character(sample_info$geoID)))
file_info <- lapply(sample_info$geoID[!already_there], getGEOSuppFiles,
                    baseDir = "ma")

# Create a list of all raw data files and assign cell lines and condition
celfiles <- list.files("ma", pattern = "*.cel*", recursive = T,
                       ignore.case = T)
names(celfiles) <- sapply(celfiles, dirname)
condition <- rep.int("disease", length(celfiles))
condition[sample_info$cell_line %in% c("HaCaT", "keratinocyte")] <- "normal"
condition <- factor(condition)

# Read the raw data files and extract expression values
raw_data <- ReadAffy(filenames = paste0("ma/", celfiles[sample_info$geoID]),
                     compress = T)
pData(raw_data)$cell_line <- sample_info$cell_line
pData(raw_data)$condition <- condition

# Normalize across arrays
eset <- frma(raw_data)

# Group sanity checks
pca <- prcomp(t(exprs(eset)))
library(ggplot2)

pca_plot <- ggplot(data.frame(pca$x), aes(x = PC1, y = PC2)) + theme_bw() +
    geom_point(aes(col = condition, shape = sample_info$cell_line))
ggsave(pca_plot, file = "pca.svg", width = 5, height = 3)
write(sprintf("Array sample pca info in first two PCs: %f%%",
    sum(pca$sdev[1:2]) / sum(pca$sdev) * 100), file = "log.txt", append = T)

cl <- kmeans(t(exprs(eset)), 2)
normal_cluster <- cl$cluster[1] # first condition is normal
cl_err <- sum(cl$cluster[condition == "normal"] !=
          normal_cluster) / length(cl$cluster)
write(sprintf("Clustering error between normal/disease: %f%%", cl_err * 100),
    file = "log.txt", append = T)

# Differential expression
mean_max <- findLargest(rownames(eset), rowMeans(exprs(eset)))
gset <- eset[mean_max, ]
rownames(gset) <- as.character(hgu133plus2ENTREZID)[mean_max]

design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)
fit <- lmFit(gset, design)
contrast.matrix <- makeContrasts(disease - normal, levels = design)
cfit <- contrasts.fit(fit, contrast.matrix)
ebfit <- eBayes(cfit)

ma_lfcs <- topTable(ebfit, number = Inf)
save(ma_lfcs, file = "gene_expression.Rd")
