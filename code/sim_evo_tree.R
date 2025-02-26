# Simulate the evolution of editing-type modifications, and obtain a distribution of modification level across genes/modification events (different modification events assumed to evolve independently)

library(expm) # Required to do matrix exponential (in A %^% T form)
library(ggplot2)

# Calculate normalized cis-genotypic value
# Input: genotype (a vector of 0 and 1), effect sizes of effector allele at all loci (a vector)
cis.gv <- function(gt,effect){
	v=sum(gt*effect)
	vmax=sum(effect)
	return(v/vmax)
}

# Calculate modification rate beta
# beta as a linear function of cis-genotypic value
# Input: normalized cis-genotypic value (a number, calculated using cis.gv()), trans- genotypic value, shape parameter characterizing relationship between cis- loci and enzyme-substrate affinity (a number), background/non-specific modification activity (a number)
beta.calc <- function(cis,trans,C,epsilon=0){
	beta=trans*(C*cis+epsilon)
	return(beta)
}

# Genotype-phenotype mapping (2 isoforms)
# Input: genotypic value of overall expression level, rate of modification (calculated using beta.calc()), decay rates (numbers)
g2p <- function(alpha,beta,gamma_0,gamma_1){
	P_0=alpha/(beta+gamma_0)
	P_1=alpha*beta/(gamma_1*(beta+gamma_0))
	return(c(P_0,P_1))
}

# Calculate fitness w/ respect to an isoform under stabilizing selection
# Input: isoform abundance (a number), optimum (a number), SD of fitness function
# Input values are not log-transformed and should be log-transformed by the function
fitness.stab <- function(z,opt,sig){
	d=log(z/opt) # Convert to log scale
	if(sig>0){
		w=dnorm(d,mean=0,sd=sig)/dnorm(0,mean=0,sd=sig)
	}else{
		w=1
	}
	return(w)
}

# Calculate fitness w/ respect to a deleterious isoform
# Input: isoform abundance, shape parameter (a positive number)
fitness.del <- function(z,par){
	w=exp(-par*z)
	return(w)
}

# Overall fitness as a product of fitness w/ respect to each isoform
# Input: isoform abundances (a vector), selection parameter (a list; vector-->stabilizing, number-->deleterious, NULL=neutral/equivalent to I_0), whether each modified isoform is functionally equivalent to I_0 (a vector, 1 means equivalent to I_0)
fitness.prod <- function(z,selection,eq=NULL){
	w=rep(0,length(z))
	type=lengths(selection) # Type of selection, which corresponds to number of selection parameters
	if(length(eq)>0){
		z[1]=z[1]+sum(z[which(eq==1)])
	}
	for(i in 1:length(z)){
		if(type[i]==2){
			w[i]=fitness.stab(z[i],selection[[i]][1],selection[[i]][2])
		}else{
			if(type[i]==1){
				w[i]=fitness.del(z[i],selection[[i]])
			}else{
				w[i]=1
			}
		}
	}
	return(prod(w))
}

# Fixation probability given ancestral and mutant fitness
fix.prob <- function(wa,wm,Ne){
	if(wa>0){
		s=(wm/wa)-1 # Coefficient of selection
		if(s==0){
			p=1/(2*Ne)
		}else{
			p=(1-exp(-2*s))/(1-exp(-4*Ne*s))
		}
	}else{ # The ancestral fitness is close to 0 (close enough to be recognized as zero by R)
		if(wm==0){ # Both ancestral and mutant fitness are close to 0
			p=1/(2*Ne) # The mutation is considered neutral
		}else{ # Ancestral fitness is close to 0 while mutant fitness isn't
			p=1 # The mutation is considered strongly beneficial
		}
	}
	return(p)
}

# Obtain transition matrix given other parameters
# Simplest version: 2 isoforms, constant transition matrix (no trans- substitutions, no environmental change), all cis- loci have equal effect size (i.e., no. of states = no. of loci+1) and mutational spectrum
# Input: number of cis- loci, mutation rates (a vector, elements being u01 and u10, respectively; 2*Ne*u should not exceed 1), selection parameter (a list), effective population size, other parameters (alpha,trans,C,gamma_0,gamma_1,epsilon,eq)
tm <- function(nloci,u,selection,Ne,alpha,trans,C,epsilon=0,gamma_0,gamma_1,eq=NULL){
	mat=matrix(0,nrow=nloci+1,ncol=nloci+1)
	for(i in 1:nrow(mat)){
		v1=(i-1)/nloci
		beta1=beta.calc(v1,trans,C,epsilon)
		z1=g2p(alpha,beta1,gamma_0,gamma_1)
		w1=fitness.prod(z1,selection,eq)
		for(j in 1:ncol(mat)){
			if(abs(i-j)==1){ # Consider adjescent states only
				v2=(j-1)/nloci
				beta2=beta.calc(v2,trans,C,epsilon)
				z2=g2p(alpha,beta2,gamma_0,gamma_1)
				w2=fitness.prod(z2,selection,eq)
				fp=fix.prob(w1,w2,Ne)
				if(i<j){
					uij=2*Ne*(nloci-(i-1))*u[1]
				}else{
					uij=2*Ne*(i-1)*u[2]
				}
				mat[i,j]=uij*fp
			}
		}
	}
	for(i in 1:nrow(mat)){
		mat[i,i]=1-sum(mat[i,])
	}
	return(mat)
}

# Read or simulate phylogenetic tree
tr=read.tree(text="((oct:5,bim:5):270,(squ:170,sep:170):105);") # Coleoidea tree used in the study

# Rescale tree such that tree height is of the desired value
scale=1e6
tr$edge.length=tr$edge.length*scale

# Label external branches among all branches for convenience of later analysis
sp=c()
for(i in 1:nrow(tr$edge)){
	des=which(tr$edge[,1]==tr$edge[i,2])
	if(length(des)==0){
		sp=c(sp,i)
	}
}

pgv=list()
for(l in 1:10){
	pgv[[l]]=rep(0,l+1)
	for(i in 0:l){
		pgv[[l]][i+1]=choose(l,i)
	}
	pgv[[l]]=pgv[[l]]/(sum(pgv[[l]]))
}

# Parameter values (All can be custom set)
# Gene-specific parameters sampled from pre-chosen distributions
Ne=1e3
u=c(1e-9,1e-9) # Mutation rates at cis-loci (to and from effector allele, respectively)
trans=1 # trans-genotypic value characterizing modification activity

# Number of different types of modification events
n1=0 # Number of modification events where the I_0 is functional
n2=0 # Number of modification events where the I_1 is functional
n3=1e4 # Number of neutral modification events, where both isoforms are functional
ngene=n1+n2+n3 # Total number of modification events
type.all=c(rep(1,n1),rep(2,n2),rep(3,n3)) # Assign 'types' to modification events
alpha.all=exp(rnorm(ngene,mean=0,sd=1)) # Expression rate, sampled from log-normal distribution
nloci.all=sample(1:10,ngene,replace=TRUE) # Number of cis- loci, sampled uniformly from a pre-chosen set
C.all=rexp(ngene,rate=10) # Parameter characterizing relationship between cis- loci and enzyme-substrate affinity, sampled from log-normal distribution
gamma_0.all=rep(1,ngene) # Decay rate of I_0
gamma_1.all=rep(1,ngene) # Decay rate of I_1
sig.all=rep(10,ngene) # Width of Gaussian fitness function characterizing strength of stabilizing selection on abundance of the functional isoform
lambda.all=rep(1e-3,ngene) # Strength of selection on abundance of deleterious isoform (for neutral modifications, the value is assigned by default but unused for fitness calculation)
epsilon.all=rep(1e-4,ngene) # Rate of non-specific, cis-genotype-independent modification

comb.all=data.frame(alpha.all,nloci.all,C.all,gamma_0.all,gamma_1.all,type.all,sig.all,lambda.all,epsilon.all) # Data matrix containing each gene's parameter

# Matrix for writing output
# Each row for a modification event, each coloumn for a terminal node (a species), values are modification levels
out=matrix(0,nrow=ngene,ncol=length(sp))
for(i in 1:ngene){
	# Read paramters for the gene under consideration
	alpha=comb.all[i,1]
	nloci=comb.all[i,2]
	C=comb.all[i,3]
	gamma_0=comb.all[i,4]
	gamma_1=comb.all[i,5]
	type=comb.all[i,6]
	sig=comb.all[i,7]
	epsilon=comb.all[i,9]
	
	# Read modification type, and arrange selection parameters into the right format
	if(type==3){ # Modification is neutral
		eq=c(0,1)
		par=list(c(alpha/gamma_0,sig),NULL)
	}else{ # Modification is non-neutral
		eq=c(0,0)
		lambda=comb.all[i,8]
		if(type==1){ # I_0 is functional
			par=list(c(alpha/gamma_0,sig),lambda)
		}else{ # I_1 is functional
			par=list(lambda,c(alpha/gamma_1,sig))
		}
	}
	
	mat=tm(nloci,u,par,Ne,alpha,trans,C,epsilon,gamma_0,gamma_1,eq) # Get transition matrix
	start=rep(0,nloci+1);v0=sample(0:(nloci),1,prob=pgv[[nloci]]);start[v0+1]=1 # Initial cis-genotype (no effector allele)
	
	v.end=rep(0,nrow(tr$edge)) # Vector to store end-point cis-genotypic value for each branch
	lv.end=rep(0,nrow(tr$edge)) # Vector to store end-point modification level for each branch
	
	# Go through each branch of the tree and simulate evolution along each branch
	for(j in 1:nrow(tr$edge)){
		ances=which(tr$edge[,2]==tr$edge[j,1]) # Find ancestral node
		if(length(ances)==0){ # If the ancestral node is the root, assign the initial cis-genotype as the branch's starting point
			start.edge=start
		}else{ # If the ancestral node is not the root, assign the end-point cis-genotype of the ancestral branch as the focal branch's starting point
			start.edge=rep(0,nloci+1)
			start.edge[v.end[ances]+1]=1
		}
		T=as.integer(tr$edge.length[j]) # Convert branch length to an integer
		distr=start.edge%*%(mat %^% T) # Get probability distribution of cis-genotypic value at the end
		v.end[j]=sample(0:nloci,1,prob=distr) # Sample a cis-genotypic value from the distribution
		beta=beta.calc(v.end[j]/nloci,trans,C,epsilon) # Calculate the corresponding modification rate
		z=g2p(alpha,beta,gamma_0,gamma_1) # Calculate isoform abundances
		lv.end[j]=z[2]/sum(z) # Get modification level
	}
	out[i,]=lv.end[sp] # Write end-point modification levels of terminal branches
}

# Write data file
colnames(out)=c("lv_bim","lv_oct","lv_squ","lv_sep")
d=data.frame(comb.all,out)
fn=paste("out_sim_tree_",n1,"_",n2,"_",n3,".txt",sep="")
write.table(d,file=fn,sep="\t")

# Re-read for analysis, if data needed are not in the R environment
d1<-read.table("out_sim_tree_0_0_10000.txt",sep="\t") # All modifications are neutral
d2<-read.table("out_sim_tree_10000_0_0.txt",sep="\t") # All modifications are deleterious
d=rbind(d1,d2) # Merge into one data frame

# A histogram for modification level distribution
g<-ggplot(d,aes(x=lv_oct))
g=g+geom_histogram()
g=g+theme_classic()+xlab("Editing level (octopus)")+ylab("Number of editing sites")
ggsave(g,file="distr_lv_sim.pdf",width=6,height=5)

# Plot the phylogenetic tree and the editome tree to compare their structures
ld=log(d[,10:13]) # Log-transform modification levels
ed=dist(t(ld),method="euclidean");edm=as.matrix(ed) # Euclidean distances between species, arranged into a distance matrix
tr.lv=nj(ed) # Neighbor-joining tree
plot(tr.lv,type="unrooted",no.margin=TRUE,show.tip.label=FALSE,edge.width=4) # Plot modification level tree
plot(tr,type="unrooted",no.margin=TRUE,show.tip.label=FALSE,edge.width=4) # Plot phylogenetic tree






