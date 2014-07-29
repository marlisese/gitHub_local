# This function simulates two binary trait for testing significance in correlations 
# the inputs are the phylogeny, the number of simulations, the rates of trait evolution for trait 1 and 2
# the ancestral states, and the names of the output simulated traits
# Function returns a table with species names and traits 1 and 2

require (phytools)
setwd   ("path_to_place_simulations")

tree = rtree(10)
nsim = 10
rates1 = c(0.03, 0.025)   # rates 0>1 and 1>0
rates2 = c(0.0284, 0.0284) # equal rates
anc1 = 1
anc2 = 1
name = "-GPsim.txt"


simul_trait <- function (tree, nsim, rates1, rates2, anc1, anc2, name) { 

    for (i in 1:nsim) {

     print(i)
    # simulating the confidence of dependent , independent_evolution

    Q1 <- matrix(c(-rates1[1],rates1[1],rates1[2],-rates1[2]),2,2)
  
    tree_simG <- sim.history(tree,Q1, anc=anc1)

    Q2 <- matrix(c(-rates2[1],rates2[1],rates2[2],-rates2[2]),2,2)
    
    tree_simP<-sim.history(tree, Q2, anc=anc2)

    
    outfile <- paste(i, name, sep="")

    write.table(cbind(names(tree_simG$states), (as.numeric(tree_simG$states)-1), (as.numeric(tree_simP$states)-1)), file=outfile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

    }
    }
    
    
simul_trait (tree, nsim, rates1, rates2, anc1, anc2, name)
