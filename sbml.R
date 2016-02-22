#  sbml.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

library(dycone)

idcomb <- function(df) {
    nd <- names(df)
    comb <- apply(df[, -1], 2, paste0, collapse=", ")
    x <- data.frame(df[1,1], t(comb))
    names(x) <- nd
    x
}

prolif <- c(atp = -20.7045, prpp = -0.053446, pyr = -0.50563, oaa = -0.35261, 
    glu = -0.38587, cit = -0.15446, `3pg` = -0.39253, adp = 20.6508, 
    pi = 20.6508)

subs <- prolif < 0
prods <- prolif > 0

r_prolif <- list(S=names(prolif)[subs], P=names(prolif)[prods], N_S=-prolif[subs],
    N_P=prolif[prods], rev=FALSE, KEGG_reaction=NA, KEGG_enzyme=NA, 
    abbreviation="proliferation")

s <- read.csv("id_map.csv", header=T)
s <- do.call(rbind, tapply(1:nrow(s), s$name, function(i) idcomb(s[i, ])))
names(s) <- c("id", "kegg.compound", "hmdb")
r <- read_reactions("reactions.csv")
r[[length(r)+1]] <- r_prolif
am <- c(KEGG_reaction="kegg.reaction", KEGG_enzyme="ec-code")

write_sbml(r, s, obj=length(r), annmap=am, out="cemet.xml")
