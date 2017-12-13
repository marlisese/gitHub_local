### This script simulates the proportion of generation-biased genes using the data reported by XX. 


library(ape)
library(phytools)
library(geiger)
library(picante)

t=read.tree("all_genes.fasta_phyml_tree.txt")
t_root<- root(t, outgroup="Desmarestia")
t_4sp <- drop.tip(t_root, "Desmarestia")
t_4sp$node.label <-NULL

names_sp <- c("Ectocarpus","Scytosiphon","Macrocystis","Saccharina")
spo_tips <- c(2097,3306,4082,3594)
names(spo_tips)<- names_sp

gam_tips <- c(4105,7527,4497,3224) 
names(gam_tips)<- names_sp

total_genes <-  c(17426,27891,19112,19673)
names(total_genes)<- names_sp

# may use better the proportion of genes 

spo_level2 <- c(156, 322) / c(6438,5061)
spo_root <- 3 / 3290

gam_level2 <- c(722, 354) / c(6438,5061)
gam_root <- 35/3290

spo_tips<- spo_tips/total_genes
gam_tips<- gam_tips/total_genes

## Calculate the rates of evolution

spo_trait <- match.phylo.data (t_4sp, spo_tips)
spo_anc_states <- fastAnc(spo_trait$phy, spo_trait$data)
spo_fit_rates <- fitContinuous(spo_trait$phy, spo_trait$data)

gam_trait <- match.phylo.data (t_4sp, gam_tips)
gam_anc_states <- fastAnc(gam_trait$phy, gam_trait$data)
gam_fit_rates <- fitContinuous(gam_trait$phy, gam_trait$data)


# simulation using the empirical equal rates, and the empirical ancestral state

#simul_spo_equal <- fastBM(spo_trait$phy, a=spo_fit_rates$opt$z0, sig2=spo_fit_rates$opt$sigsq, nsim=100, internal=T)

#simul_gam_equal <- fastBM(gam_trait$phy, a=gam_fit_rates$opt$z0, sig2=gam_fit_rates$opt$sigsq, nsim=100, internal=T)

########################
# with the true root

simul_spo_equal <- fastBM(spo_trait$phy, a=spo_root, sig2=spo_fit_rates$opt$sigsq, nsim=100, internal=T)

simul_gam_equal <- fastBM(gam_trait$phy, a=gam_root, sig2=gam_fit_rates$opt$sigsq, nsim=100, internal=T)


# Simulation settings

# ES is the clade with Ectocarpus and Scytosiphon, while MS is Saccharina and Macrocystis, where sporophyte is bigger than gametophyte and maybe some genes evolve faster. 

# edgelabels 4,5,6 are ES
# edgelabels 1,2,3 are MS
 
 MS_edges <- t_4sp$edge[c(1,2,3),2]
 tree_painted <- paintBranches (t_4sp, edge=MS_edges, "MS", anc.state="ES")
 
 pdf ("Plot_simulation_MS_ES.pdf")
 plotSimmap (tree_painted)
 dev.off()
 
 # simulation bounded by 0-1 which is the percentage of genes
sim.rates_bound <- function (tree, sig2, anc = 0, nsim = 1, internal = F, plot = F) {
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (is.null(tree$mapped.edge)) {
    message("tree does not contain a mapped discrete character history, using fastBM")
    X <- fastBM(tree, a = anc, sig2 = sig2[1], nsim = nsim, 
                internal = internal)
  }
  else {
    if (is.null(names(sig2))) {
      message("names absent from sig2: assuming same order as $mapped.edge")
      if (length(sig2) == ncol(tree$mapped.edge)) 
        names(sig2) <- colnames(tree$mapped.edge)
      else stop("the number of elements in sig2 should match the number of rows in mapped.edge")
    }
    sig2 <- sig2[colnames(tree$mapped.edge)]
    edge.length <- rep(0, nrow(tree$edge))
    for (i in 1:ncol(tree$mapped.edge)) edge.length <- edge.length + 
      sig2[i] * tree$mapped.edge[, i]
    names(edge.length) <- NULL
    tree_new <- list(Nnode = tree$Nnode, edge = tree$edge, tip.label = tree$tip.label, 
                 edge.length = edge.length)
    class(tree_new) <- "phylo"
    if (plot) 
      plot(tree_new)
    X <- fastBM(tree_new, a = anc, nsim = nsim, internal = internal, bounds=c(0,1))
  }
  X
}

## 3x faster in MS

sig2 = spo_fit_rates$opt$sigsq
fast=3
rate_vector <- c(sig2, sig2*fast) 
names(rate_vector) <- c("ES", "MS")
trait_spo_x3<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=spo_root, nsim=100, internal=T)

sig2 = gam_fit_rates$opt$sigsq
fast= 3
rate_vector <- c(sig2, sig2*fast) 
names(rate_vector) <- c("ES", "MS")
trait_gam_x3<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=gam_root, nsim=100, internal=T)

## 5x faster in MS

sig2 = spo_fit_rates$opt$sigsq
fast= 5
rate_vector <- c(sig2, sig2*fast) 
names(rate_vector) <- c("ES", "MS")
trait_spo_x5<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=spo_root, nsim=100, internal=T)

sig2 = gam_fit_rates$opt$sigsq
fast= 5
rate_vector <- c(sig2, sig2*fast) 
names(rate_vector) <- c("ES", "MS")
trait_gam_x5<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=gam_root, nsim=100, internal=T)

## 10x faster in MS

sig2 = spo_fit_rates$opt$sigsq
fast=10
rate_vector <- c(sig2, sig2*fast) 
names(rate_vector) <- c("ES", "MS")
trait_spo_x10<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=spo_root, nsim=100, internal=T)

sig2 = gam_fit_rates$opt$sigsq
fast= 10
rate_vector <- c(sig2, sig2*fast) 
names(rate_vector) <- c("ES", "MS")
trait_gam_x10<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=gam_root, nsim=100, internal=T)

## 20x faster in MS

sig2 = spo_fit_rates$opt$sigsq
fast=20
rate_vector <- c(sig2, sig2*fast) 
names(rate_vector) <- c("ES", "MS")
trait_spo_x20<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=spo_root, nsim=100, internal=T)

sig2 = gam_fit_rates$opt$sigsq
fast= 20
rate_vector <- c(sig2, sig2*fast) 
names(rate_vector) <- c("ES", "MS")
trait_gam_x20<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=gam_root, nsim=100, internal=T)


## Simulation where tips are equal, but t2 is slower.

# t1 is where the tips are, and t2 is where the division of the two groups ocurred, in time  

# edgelabels 2,3,5,6 are t1
# edgelabels 1,4 are t2
 
 MS_edges <- t_4sp$edge[c(1,4),2]
 tree_painted <- paintBranches (t_4sp, edge=MS_edges, "t2", anc.state="t1")
 
 pdf ("Plot_simulation_t1_t2.pdf")
 plotSimmap (tree_painted)
 dev.off()
 
## 5x faster in t1
 
 sig2 = spo_fit_rates$opt$sigsq
 fast= 5
 rate_vector <- c(sig2*fast, sig2) 
 names(rate_vector) <- c("t1", "t2")
 trait_spo_time_x5<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=spo_root, nsim=100, internal=T)

 sig2 = gam_fit_rates$opt$sigsq
 fast= 5
 rate_vector <- c(sig2*fast, sig2) 
 names(rate_vector) <- c("t1", "t2")
 trait_gam_time_x5<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=gam_root, nsim=100, internal=T)
 
## 10x faster in t1
 
 sig2 = spo_fit_rates$opt$sigsq
 fast= 10
 rate_vector <- c(sig2*fast, sig2) 
 names(rate_vector) <- c("t1", "t2")
 trait_spo_time_x10<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=spo_root, nsim=100, internal=T)

 sig2 = gam_fit_rates$opt$sigsq
 fast= 10
 rate_vector <- c(sig2*fast, sig2) 
 names(rate_vector) <- c("t1", "t2")
 trait_gam_time_x10<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=gam_root, nsim=100, internal=T)
 
## 20x faster in t1
 
 sig2 = spo_fit_rates$opt$sigsq
 fast= 20
 rate_vector <- c(sig2*fast, sig2) 
 names(rate_vector) <- c("t1", "t2")
 trait_spo_time_x20<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=spo_root, nsim=100, internal=T)

 sig2 = gam_fit_rates$opt$sigsq
 fast= 20
 rate_vector <- c(sig2*fast, sig2) 
 names(rate_vector) <- c("t1", "t2")
 trait_gam_time_x20<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=gam_root, nsim=100, internal=T) 
 
 ## 50x faster in t1
 
 sig2 = spo_fit_rates$opt$sigsq
 fast= 50
 rate_vector <- c(sig2*fast, sig2) 
 names(rate_vector) <- c("t1", "t2")
 trait_spo_time_x50<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=spo_root, nsim=100, internal=T)

 sig2 = gam_fit_rates$opt$sigsq
 fast= 50
 rate_vector <- c(sig2*fast, sig2) 
 names(rate_vector) <- c("t1", "t2")
 trait_gam_time_x50<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=gam_root, nsim=100, internal=T) 
 
  ## 10x faster in t1 and 10x slower in t2, from the estimation
 
 sig2 = spo_fit_rates$opt$sigsq
 fast= 10
 rate_vector <- c(sig2*fast, sig2/10) 
 names(rate_vector) <- c("t1", "t2")
 trait_spo_time_x10_x10<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=spo_root, nsim=100, internal=T)

 sig2 = gam_fit_rates$opt$sigsq
 fast= 10
 rate_vector <- c(sig2*fast, sig2/10) 
 names(rate_vector) <- c("t1", "t2")
 trait_gam_time_x10_x10<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=gam_root, nsim=100, internal=T) 
 
  ## 10x faster in t1 and 10x slower in t2, from the estimation
 
 sig2 = spo_fit_rates$opt$sigsq
 fast= 20
 rate_vector <- c(sig2*fast, sig2*2) 
 names(rate_vector) <- c("t1", "t2")
 trait_spo_time_x20_x2<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=spo_root, nsim=100, internal=T)

 sig2 = gam_fit_rates$opt$sigsq
 fast= 20
 rate_vector <- c(sig2*fast, sig2*2) 
 names(rate_vector) <- c("t1", "t2")
 trait_gam_time_x20_x2<-sim.rates_bound(tree_painted, sig2=rate_vector,anc=gam_root, nsim=100, internal=T) 
 
 
 
## plotting all results

pdf("Prop_genes_SPO_GAM_biased.pdf", width=15, height=8) 

par(mfrow=c(1,2))

plot(t_4sp, main="Sporophyte specific % genes")
add.scale.bar()
text (x=0.115, y=c(1,2,3,4), labels= round(spo_tips, digits=2))
text (x=c(0.008,0.042, 0.05), y=c(2.5,3.5, 1.5) ,labels=round(c(spo_root, spo_level2), digits=3) )


plot(t_4sp, main="Gametophyte specific % genes")
add.scale.bar()
text (x=0.115, y=c(1,2,3,4), labels= round(gam_tips, digits=2))
text (x=c(0.008,0.042, 0.05), y=c(2.5,3.5, 1.5) ,labels=round(c(gam_root, gam_level2), digits=3) )

dev.off()

# Expected ancestral states using the tips information vs observed values
# Expected rates 

pdf("Observed_Expected_SPO_GAM_biased.pdf", width=9, height=4) 

par(mfrow=c(1,3))

plot(t_4sp)
nodelabels()

barplot (rbind(as.numeric(spo_anc_states), c(spo_root, spo_level2)), beside=T, main="Exp vs Obs % genes Spo-biased", names=c(5,6,7), ylim=c(0,0.25))
text(x=4, y=0.22, labels=paste("Est Rate ", round(spo_fit_rates$opt$sigsq, digits=3), sep=""))

barplot (rbind(as.numeric(gam_anc_states), c(gam_root, gam_level2)), beside=T, main="Exp vs Obs % genes Gam-biased", names=c(5,6,7),ylim=c(0,0.25))
text(x=4, y=0.22, labels=paste("Est Rate ", round(gam_fit_rates$opt$sigsq, digits=3), sep=""))

dev.off()

## simulations with equal rates between clades, and taken from  empirical with root empirical value 


pdf("Simulations_rates_violin.pdf", width=6, height=12) 

## Equal rates (from estimated rates for SPO and GAM)


library(vioplot)

par(mfrow=c(5,1))

## lineage variation sporophyte

vioplot ( simul_spo_equal [1,], simul_spo_equal [2,], simul_spo_equal [3,], simul_spo_equal [4,], simul_spo_equal [6,],simul_spo_equal [7,], names=rownames(simul_spo_equal)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.12, 0.3))
title(main="Equal rates: sporophyte", ylab="Simulated prop. genes")
points(c(1:6), c(spo_tips[c(3,4,1,2)], spo_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_spo_x3 [1,], trait_spo_x3 [2,], trait_spo_x3 [3,], trait_spo_x3 [4,], trait_spo_x3 [6,],trait_spo_x3 [7,], names=rownames(trait_spo_x3)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.3))
title(main="3X MS rates: sporophyte", ylab="Simulated prop. genes")
points(c(1:6), c(spo_tips[c(3,4,1,2)], spo_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_spo_x5 [1,], trait_spo_x5 [2,], trait_spo_x5 [3,], trait_spo_x5 [4,], trait_spo_x5 [6,],trait_spo_x5 [7,], names=rownames(trait_spo_x5)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.3))
title(main="5X MS rates: sporophyte", ylab="Simulated prop. genes")
points(c(1:6), c(spo_tips[c(3,4,1,2)], spo_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_spo_x10 [1,], trait_spo_x10 [2,], trait_spo_x10 [3,], trait_spo_x10 [4,], trait_spo_x10 [6,],trait_spo_x10 [7,], names=rownames(trait_spo_x10)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.4))
title(main="10X MS rates: sporophyte", ylab="Simulated prop. genes")
points(c(1:6), c(spo_tips[c(3,4,1,2)], spo_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_spo_x20 [1,], trait_spo_x20 [2,], trait_spo_x20 [3,], trait_spo_x20 [4,], trait_spo_x20 [6,],trait_spo_x20 [7,], names=rownames(trait_spo_x20)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.5))
title(main="20X MS rates: sporophyte", ylab="Simulated prop. genes")
points(c(1:6), c(spo_tips[c(3,4,1,2)], spo_level2[c(2,1)]), col="red", pch=19)


## time variation 

par(mfrow=c(4,1))

vioplot ( trait_spo_time_x5 [1,], trait_spo_time_x5 [2,], trait_spo_time_x5 [3,], trait_spo_time_x5 [4,], trait_spo_time_x5 [6,],trait_spo_time_x5 [7,], names=rownames(trait_spo_time_x5)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.25))
title(main="5X faster rates in t1: sporophyte", ylab="Simulated prop. genes")
points(c(1:6), c(spo_tips[c(3,4,1,2)], spo_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_spo_time_x10 [1,], trait_spo_time_x10 [2,], trait_spo_time_x10 [3,], trait_spo_time_x10 [4,], trait_spo_time_x10 [6,],trait_spo_time_x10 [7,], names=rownames(trait_spo_time_x10)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.30))
title(main="10X faster rates in t1: sporophyte", ylab="Simulated prop. genes")
points(c(1:6), c(spo_tips[c(3,4,1,2)], spo_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_spo_time_x20 [1,], trait_spo_time_x20 [2,], trait_spo_time_x20 [3,], trait_spo_time_x20 [4,], trait_spo_time_x20 [6,],trait_spo_time_x20 [7,], names=rownames(trait_spo_time_x20)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.5))
title(main="20X faster rates in t1: sporophyte", ylab="Simulated prop. genes")
points(c(1:6), c(spo_tips[c(3,4,1,2)], spo_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_spo_time_x50 [1,], trait_spo_time_x50 [2,], trait_spo_time_x50 [3,], trait_spo_time_x50 [4,], trait_spo_time_x50 [6,],trait_spo_time_x50 [7,], names=rownames(trait_spo_time_x50)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.65))
title(main="50X faster rates in t1: sporophyte", ylab="Simulated prop. genes")
points(c(1:6), c(spo_tips[c(3,4,1,2)], spo_level2[c(2,1)]), col="red", pch=19)

## variation gametophyte 

par(mfrow=c(5,1))

vioplot ( simul_gam_equal [1,], simul_gam_equal [2,], simul_gam_equal [3,], simul_gam_equal [4,], simul_gam_equal [6,],simul_gam_equal [7,], names=rownames(simul_gam_equal)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.12, 0.3))
title(main="Equal rates: gametophyte", ylab="Simulated prop. genes")
points(c(1:6), c(gam_tips[c(3,4,1,2)], gam_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_gam_x3 [1,], trait_gam_x3 [2,], trait_gam_x3 [3,], trait_gam_x3 [4,], trait_gam_x3 [6,],trait_gam_x3 [7,], names=rownames(trait_gam_x3)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.35))
title(main="3X MS rates: gametophyte", ylab="Simulated prop. genes")
points(c(1:6), c(gam_tips[c(3,4,1,2)], gam_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_gam_x5 [1,], trait_gam_x5 [2,], trait_gam_x5 [3,], trait_gam_x5 [4,], trait_gam_x5 [6,],trait_gam_x5 [7,], names=rownames(trait_gam_x5)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.4))
title(main="5X MS rates: gametophyte", ylab="Simulated prop. genes")
points(c(1:6), c(gam_tips[c(3,4,1,2)], gam_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_gam_x10 [1,], trait_gam_x10 [2,], trait_gam_x10 [3,], trait_gam_x10 [4,], trait_gam_x10 [6,],trait_gam_x10 [7,], names=rownames(trait_gam_x10)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.4))
title(main="10X MS rates: gametophyte", ylab="Simulated prop. genes")
points(c(1:6), c(gam_tips[c(3,4,1,2)], gam_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_gam_x20 [1,], trait_gam_x20 [2,], trait_gam_x20 [3,], trait_gam_x20 [4,], trait_gam_x20 [6,],trait_gam_x20 [7,], names=rownames(trait_gam_x20)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.7))
title(main="20X MS rates: gametophyte", ylab="Simulated prop. genes")
points(c(1:6), c(gam_tips[c(3,4,1,2)], gam_level2[c(2,1)]), col="red", pch=19)

par(mfrow=c(4,1))


vioplot ( trait_gam_time_x5 [1,], trait_gam_time_x5 [2,], trait_gam_time_x5 [3,], trait_gam_time_x5 [4,], trait_gam_time_x5 [6,],trait_gam_time_x5 [7,], names=rownames(trait_gam_time_x5)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.35))
title(main="5X faster rates in t1: gametophyte", ylab="Simulated prop. genes")
points(c(1:6), c(gam_tips[c(3,4,1,2)], gam_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_gam_time_x10 [1,], trait_gam_time_x10 [2,], trait_gam_time_x10 [3,], trait_gam_time_x10 [4,], trait_gam_time_x10 [6,],trait_gam_time_x10 [7,], names=rownames(trait_gam_time_x10)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.35))
title(main="10X faster rates in t1: gametophyte", ylab="Simulated prop. genes")
points(c(1:6), c(gam_tips[c(3,4,1,2)], gam_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_gam_time_x20 [1,], trait_gam_time_x20 [2,], trait_gam_time_x20 [3,], trait_gam_time_x20 [4,], trait_gam_time_x20 [6,],trait_gam_time_x20 [7,], names=rownames(trait_gam_time_x20)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.45))
title(main="20X faster rates in t1: gametophyte", ylab="Simulated prop. genes")
points(c(1:6), c(gam_tips[c(3,4,1,2)], gam_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_gam_time_x50 [1,], trait_gam_time_x50 [2,], trait_gam_time_x50 [3,], trait_gam_time_x50 [4,], trait_gam_time_x50 [6,],trait_gam_time_x50 [7,], names=rownames(trait_gam_time_x50)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.7))
title(main="50X faster rates in t1: gametophyte", ylab="Simulated prop. genes")
points(c(1:6), c(gam_tips[c(3,4,1,2)], gam_level2[c(2,1)]), col="red", pch=19)


## time twice variation 

vioplot ( trait_spo_time_x10_x10 [1,], trait_spo_time_x10_x10 [2,], trait_spo_time_x10_x10 [3,], trait_spo_time_x10_x10 [4,], trait_spo_time_x10_x10 [6,],trait_spo_time_x10_x10 [7,], names=rownames(trait_spo_time_x10_x10)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.3))
title(main="10X faster t1-t2, 10X slower t2-t3: sporophyte", ylab="Simulated prop. genes")
points(c(1:6), c(spo_tips[c(3,4,1,2)], spo_level2[c(2,1)]), col="red", pch=19)


vioplot ( trait_spo_time_x20_x2 [1,], trait_spo_time_x20_x2 [2,], trait_spo_time_x20_x2 [3,], trait_spo_time_x20_x2 [4,], trait_spo_time_x20_x2 [6,],trait_spo_time_x20_x2 [7,], names=rownames(trait_spo_time_x20_x2)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.5))
title(main="20X faster t1-t2, 2X faster t2-t3: sporophyte", ylab="Simulated prop. genes")
points(c(1:6), c(spo_tips[c(3,4,1,2)], spo_level2[c(2,1)]), col="red", pch=19)


vioplot ( trait_gam_time_x10_x10 [1,], trait_gam_time_x10_x10 [2,], trait_gam_time_x10_x10 [3,], trait_gam_time_x10_x10 [4,], trait_gam_time_x10_x10 [6,],trait_gam_time_x10_x10 [7,], names=rownames(trait_gam_time_x10_x10)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.4))
title(main="10X faster t1-t2, 10X slower t2-t3: gametophyte", ylab="Simulated prop. genes")
points(c(1:6), c(gam_tips[c(3,4,1,2)], gam_level2[c(2,1)]), col="red", pch=19)

vioplot ( trait_gam_time_x20_x2 [1,], trait_gam_time_x20_x2 [2,], trait_gam_time_x20_x2 [3,], trait_gam_time_x20_x2 [4,], trait_gam_time_x20_x2 [6,],trait_gam_time_x20_x2 [7,], names=rownames(trait_gam_time_x20_x2)[c(1,2,3,4,6,7)], col="gray", ylim=c(-0.01, 0.55))
title(main="20X faster t1-t2, 2X faster t2-t3: gametophyte", ylab="Simulated prop. genes")
points(c(1:6), c(gam_tips[c(3,4,1,2)], gam_level2[c(2,1)]), col="red", pch=19)


dev.off()
