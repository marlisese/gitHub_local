# This function creates a richness matrix for Medusa when you have missing species values for genera, 
# but they have to be assign to species tips. 
# Function takes input of a tree and table of species richness where columns are genus, species names 
# (as in the tree) and richness. Richness values are amount of missing species in the genus, repeated over
# all species belonging to that genus. 
# Function produces a even distribution of the missing species: If # missing species is smaller 
# than sampled species it assigns randomly to the tips, one single missing species. If # missing species 
# is bigger than sampled, it distributes # species evenly (as integers) and randomly assign the residuals. 
# Martha Serrano 
# July 24, 2014 


richness_sp <- function(tree, data) { 

# creating the richness table 

resamp <- function(x,...){if(length(x)==1) x else sample(x,...)} 
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

rich <- matrix(nrow=length(tree$tip.label), ncol=2)
rich[,1]=tree$tip.label
rich[,2]=0

for (i in 1:length(table(data$genus))) {

p = which(data$genus==names(table(data$genus)[i]))

# number species missing is bigger than # sp present 

    if (length(p)<data[p[1],3]) {
      
	  q = which(sapply(rich[,1],function(m){any(m==data$species_tip[p])}),  arr.ind = TRUE)
    
	  if (is.wholenumber(data[p[1],3] / length(p)) ==TRUE) { 
	  
		  w = data[p[1],3] / length(p)
		  rich[q,2] = w
	      
			  } else { 
		  
		  v = round(data[p[1],3] / length(p))
	  	  d = data[p[1],3] - (v*length(p))
		  rich[q,2] = v
	  	 
		    if (d <0) { 
		  	rich[resamp(q, size=abs(d)),2]= (v-1)
			      }  else { 
	  		rich[resamp(q, size=d),2] = (v+1) 
				      }
				  }
	  } else { 
      
	  q = which(sapply(rich[,1],function(m){any(m==data$species_tip[p])}),  arr.ind = TRUE)
    
	  rich[resamp(q, size=data[p[1],3])  , 2] = "1" 
      		  }
					 }
# 
# for (i in 1:length(table(data$genus))) {
# 	  y = which(data$genus==names(table(data$genus)[i]))
# 	  x = resamp(y, size=1)
# 	  z = which(rich[,1]==data$species_tip[x])
# 	  rich[z,2]= data[x,3]
# 	  }

rich=as.data.frame(rich)
colnames(rich)=c("taxon", "n.taxa")

return (rich)

}

