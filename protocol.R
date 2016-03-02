#  hela.R
#  
#  Copyright 2015 Christian Diener <ch.diener[a]gmail.com>
#  
#  MIT license. See LICENSE for more information.


start_t <- proc.time()

VOL <- 1.54e-12  # cell volume in L
NC <- 1e6  # number of cells in each sample
C_UNIT_M <- 1e-12  # concentration unit measured (picomolar)
C_UNIT_U <- 1e-6 # used concentration unit (micromolar)
d_idx <- 5:10  # columns where data entries start in the original data file

###### 

#library(dycone)
devtools::load_all("~/code/dycone")

id_map <- read.csv("id_map.csv", stringsAsFactors=FALSE)
r <- read_reactions("reactions.csv")
metab <- read.csv("metabolome.csv")
# convert to uM concentrations
metab[, d_idx] <- C_UNIT_M / C_UNIT_U * metab[, d_idx]/(NC * VOL)

matches <- sapply(id_map$kegg, grep_id, x = metab$kegg_id)
miss <- is.na(matches)

# Assemble complete data set

full <- id_map[, 1:2]
m <- matrix(NA, nrow = nrow(full), ncol = length(d_idx))
colnames(m) <- names(metab)[d_idx]
full <- cbind(full, m)
matched_idx <- !is.na(matches)
full[matched_idx, 3:8] <- metab[matches[matched_idx], d_idx]

{
    concs <- hmdb_concentration(id_map$hmdb, add = id_map[, 1:2])
} %c% "scraped_concs.Rd"

m_concs <- as.vector(by(concs, concs$name, priority_mean))
names(m_concs) <- levels(factor(concs$name))

# Patching
scraped <- data.frame(name = names(m_concs), normal = m_concs)
rownames(scraped) <- NULL
patched <- patch(full, id = 1, normal = 3:5, treatment = 6:8, ref_data = scraped)

# Get basis
S <- stoichiometry(r)
mats <- ma_terms(S, patched[, c(1, 3:8)])

# Get missing fraction
Smod <- t((S[id_map$name, ] > 0) + 0)
n_miss <- Smod %*% miss
n_measured <- Smod %*% !miss
write(sprintf("%d reactions missing at least one metabolite", sum(n_miss > 0)),
    file="log.txt")
write(sprintf("%d reactions having at least one measured metabolite\n", sum(n_measured > 0)), 
    file="log.txt", append=T)

{
    V <- polytope_basis(S)
} %c% "basis.Rd"

library(doParallel)
registerDoParallel(cl = 6)

prolif <- c(atp = -20.7045, prpp = -0.053446, pyr = -0.50563, oaa = -0.35261, 
    glu = -0.38587, cit = -0.15446, `3pg` = -0.39253, adp = 20.6508, 
    pi = 20.6508)

# Generate hypothesis
h <- hyp(mats[, 1:3], mats[, 4:6], r, type = "fva", 
	obj = prolif, v_min = 1e-12, full = T)
pw <- rp(make_irreversible(r), "pathway")[,2]
r_ids <- rp(make_irreversible(r), "KEGG_reaction")[,2]
h$hyp <- cbind(h$hyp, pathway = pw[h$hyp$idx])
h$hyp <- cbind(h$hyp, reaction_id = r_ids[h$hyp$idx])
write.csv(h$hyp, file = "regulatory_bias.csv", quote = F, row.names = F)

# Check for reaction-order correlation
ho <- h$hyp[order(h$hyp$idx), ]
ro <- r_order(make_irreversible(r))
write(sprintf("Bias correlation\nlog-fold: %g\npval: %g\n", cor(ro, abs(ho$mean_log_fold)), 
    cor(ro, ho$pval)), file = "log.txt", append=T)

# Stability analysis
{
    stab <- foreach(i = 1:6, .combine = rbind) %dopar% {
        concs <- patched[, 2 + i]
        names(concs) <- patched$name
        s <- stability_analysis(kcone(V, mats[, i], normalize=T), S, concs)
        cbind(s, basis_idx = i)
    }
} %c% "stab.Rd"

cell_line <- rep(c("HaCaT", "HeLa"), each = 3)
stab$cell_line <- factor(cell_line[stab$basis_idx], levels = c("HaCaT", "HeLa"))

## Comparison with gene expression data
library(AnnotationDbi)
library(hgu133plus2.db)
library(dplyr)

source("gene_expression.R")

ecs <- rp(make_irreversible(r), "KEGG_enzyme")
ens <- AnnotationDbi::select(hgu133plus2.db, keys = ecs[, 2], keytype = "ENZYME", columns = c("SYMBOL", 
    "ENSEMBL", "ENTREZID"))
load("gene_expression.Rd")
info <- ecs %>% group_by(r_idx) %>% do(ens[ens$ENZYME %in% .$KEGG_enzyme,])
good <- sapply(info$ENTREZID, function(eid) !is.na(eid) & 
    (eid %in% rownames(ma_lfcs)))
write(sprintf("%f%% of enzymes found on array.", sum(good)/length(good)*100), 
    file="log.txt", append=T)
info <- info[good, ]
info <- info %>% group_by(r_idx) %>% mutate(met_lfc = 
    h$hyp$mean_log_fold[h$hyp$idx %in% r_idx], met_pval = 
    h$hyp$corr_pval[h$hyp$idx %in% r_idx], pathway = 
    h$hyp$pathway[h$hyp$idx %in% r_idx])

get_ma_lfc <- function(eid) {
    found <- which(eid[1] == rownames(ma_lfcs))

    return(data.frame(ge_lfc = ma_lfcs$logFC[found], 
        ge_pval = ma_lfcs$adj.P.Val[found]))
}
    
ge <- lapply(info$ENTREZID, get_ma_lfc)
info <- cbind(info, do.call(rbind, ge))
write.csv(info, "all_lfcs.csv")


### Plots ###
library(ggplot2)
library(tidyr)
library(pheatmap)

comb <- expand.grid(a = 1:3, b = 4:6)
mlfc <- apply(metab[, 5:10], 1, function(m) log(m[comb$b], 2) - log(m[comb$a], 2))
colnames(mlfc) <- as.character(metab$name)
mlfc <- gather(data.frame(mlfc, check.names = F), name, lfc, factor_key=T)
metab_plot <- ggplot(mlfc, aes(x = name, y = lfc)) + geom_boxplot() + theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
	xlab("") + ylab("log-fold change")
ggsave("metabs.svg", metab_plot, width = 10, height = 6)

write("Plotting log-fold changes...", file="")
x <- h$lfc_disease
m <- max(range(x))
ann <- data.frame("flux variability" = h$hyp$fva_log_fold, check.names = FALSE)

pheatmap(t(x[h$hyp$idx, ]), breaks = seq(-m, m, length.out = 102), col = dycone:::DC_DIVCOL(101), 
    cluster_rows = F, cluster_cols = F, cellwidth = 10, cellheight = 10, labels_col = h$hyp$name, 
    annotation_col = ann, show_rownames = FALSE, filename = "lfcs.pdf", width = 18, height = 3.5)

# Normality tests
write("Plotting normality tests ...", file = "")

log_all <- log(metab[, 5:10], 2)
colnames(log_all) <- c(paste0("HaCaT_", 1:3), paste0("HeLa_", 1:3))
log_hacat <- unlist(log_all[, 1:3])
log_hela <- unlist(unlist(log_all[, 4:6]))

svg("norm.svg", width=12, height=6)
par(mfrow = 1:2)
qqnorm(log_hacat, pch = 1, col = "royalblue", main = "")
qqline(log_hacat, lwd = 2, col = "darkblue")
points(qqnorm(log_hela, plot.it = F), pch = 2, col = "red3")
qqline(log_hela, lwd = 2, col = "darkred")

x <- seq(0, 16, length.out = 256)
plot(ecdf(log_hacat), pch = 1, col = "royalblue", main = "")
lines(x, pnorm(x, mean(log_hacat, na.rm = T), sd(log_hacat, na.rm = T)), lwd = 2,
    col = "darkblue")
plot(ecdf(log_hela), pch = 2, add = T, col = "red3")
lines(x, pnorm(x, mean(log_hela, na.rm = T), sd(log_hela, na.rm = T)), lwd = 2,
    col = "darkred")
dev.off()


# Plots
write("Preparing k-cones...", file="")
K <- lapply(1:6, function(i) kcone(V, mats[, i]))
combs <- expand.grid(4:6, 1:3)
lfcs <- apply(combs, 1, function(i) log(mats[, i[2]], 2) - log(mats[, i[1]], 2))
large <- abs(rowMeans(lfcs)) > 1 
K_subset <- lapply(K, function(ki) ki[large, ])

rn <- rp(make_irreversible(r), "abbreviation")[,2]
png("red_all.png", width = 6, height = 6, units = "in", res = 600)
write("All:", file = "log.txt", append = T)
par(mgp = c(1.75, 0.6, 0), cex.lab=1.25)
m <- capture.output(plot_red(K, col = rep(c("blue", "tomato2"), each=3), 
	r_names = rn))
write(m, file = "log.txt", append = T)
dev.off()

png("red_sig.png", width = 6, height = 6, units = "in", res = 600)
write("\ncutoff:", file = "log.txt", append = T)
par(mgp = c(1.75, 0.6, 0), cex.lab=1.25)
m <- capture.output(plot_red(K_subset, col = rep(c("blue", "tomato2"), each=3), 
    r_names = rn[large]))
write(m, file = "log.txt", append = T)
dev.off()

write("Plotting stability...", file="")
stab_plot <- ggplot(stab, aes(x = what, fill = cell_line, group = basis_idx)) + 
    stat_count(position = "dodge", col = "black") + theme_bw() + 
    scale_fill_manual(values = c("royalblue", "red3")) + 
    theme(legend.position = c(0.75, 0.75))
ggsave(stab_plot, file = "stab.svg", width = 4, height = 4)

write("Plotting heterogeneity and correlations...", file="")
sd_plot <- ggplot(h$hyp, aes(x = sd_normal, y = sd_disease, col = pathway)) + 
    geom_polygon(data = data.frame(x = c(0, 9, 9), y = c(0, 3, 27)), 
    aes(x = x, y = y), fill = "blue", alpha = 0.1, col = NA) + 
    geom_abline(colour="blue", alpha=0.75) + geom_point() + 
    theme_bw() + coord_cartesian(xlim = c(-0.1, 2), ylim = c(-0.1, 2)) + 
    xlab(expression(HaCaT/HaCaT ~ sigma)) + ylab(expression(HeLa/HaCaT ~ sigma))
ggsave(sd_plot, file = "sds.svg", width = 5.5, height = 3)

fc_sd <- h$hyp$sd_disease/h$hyp$sd_normal
sig <- fc_sd > 3 & is.finite(fc_sd)
type <- paste0(h$hyp$name[sig], " (", h$hyp$type[sig], ")")
d <- cor(t(h$lfc_disease[h$hyp$idx[sig], ]))
pheatmap(d, col = dycone:::DC_DIVCOL(101), breaks = seq(-1, 1, length.out = 102), 
    cellwidth = 10, cellheight = 10, labels_row = type, labels_col = type,
    filename = "cors.pdf", height = 4, width = 4)

write("Plotting pathways...", file="")
pathway_plot <- ggplot(h$hyp, aes(y = pathway, x = mean_log_fold, col = pathway)) + 
    geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_point(aes(shape = abs(mean_log_fold) > fva_log_fold), 
    position = position_jitter(height = 0.2), size = 3) + theme_bw() + 
    scale_shape_manual(values = c(1, 17)) + theme(legend.position = "none") + 
    xlab("mean log-fold change") + ylab("")
ggsave(pathway_plot, file = "pathways.svg", width = 6, height = 4)

write("Plotting expression vs. enzyme activity...", file="")
info$significant <- info$ge_pval<0.05 & info$met_pval<0.05
explot <- ggplot(info, aes(x=ge_lfc, y=met_lfc, color=pathway, shape=significant)) + 
    geom_hline(yintercept=0, linetype="dashed") + 
    geom_vline(xintercept=0, linetype="dashed") + geom_point() + 
    theme_bw() + xlab("gene expression") + scale_color_discrete(drop=FALSE) + 
    ylab("metabolome")
ggsave(explot, file = "gene_ex.svg", width=6, height=3)

write("----------\nUsed time:\n----------", file = "")
print(proc.time() - start_t) 
