# Mapping between genotype and phenotype

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

# Calculate modification rate beta
# cis- genotypic value's effect on beta diminishes
# Input: normalized cis-genotypic value (a number, calculated using cis.gv()), trans- genotypic value, shape parameter characterizing relationship between cis- loci and enzyme-substrate affinity (a number), background/non-specific modification activity (a number)
beta.calc <- function(cis,trans,C,epsilon=0){
	beta=trans*(cis+epsilon)/(cis+C)
	return(beta)
}

# Calculate modification rate beta
# Two effector alleles, one increasing beta_1, the other increasing beta_2 (e.g., alternative splicing)
# Input: normalized cis-genotypic value for beta_1 (that for beta_2 is 1-cis), trans- genotypic value, shape parameter (a vector), background/non-specific modification activity (a vector; epsilon[1] should be high enough to ensure low mis-splicing rate)
beta.calc.type.2 <- function(cis,trans,C,epsilon=c(0,0)){
	beta1=trans*(C[1]*cis+epsilon[1])
	beta2=trans*(C[2]*(1-cis)+epsilon[2])
	return(c(beta1,beta2))
}

# Calculate modification rate beta
# Quantity of interest is the number of modifications in an RNA/protein molecule (e.g., m6A, protein phosphorylation)
# cis- loci shared by all modifications, effect on declines as modification number increases
# Input: number of modification sites (a number), normalized cis-genotypic value (a number, calculated using cis.gv()), trans- genotypic value, shape parameter characterizing relationship between cis- loci and enzyme-substrate affinity (a number), diminishing resturn parameter characterizing the effect of existing modifications on beta (a number), background/non-specific modification activity (a number)
# Output is a vector of length n
beta.calc.type.3 <- function(n,cis,trans,C,epsilon=0){
	beta=rep(0,n)
	for(i in 1:n){
		beta[i]=trans*((C*cis+epsilon)^((n-i)/n))
	}
	return(beta)
}

# Genotype-phenotype mapping (2 isoforms)
# Input: genotypic value of overall expression rate, rate of modification (calculated using beta.calc()), decay rates (numbers)
g2p <- function(alpha,beta,gamma_0,gamma_1){
	P_0=alpha/(beta+gamma_0)
	P_1=alpha*beta/(gamma_1*(beta+gamma_0))
	return(c(P_0,P_1))
}

# Genotype-phenotype mapping (>2 isoforms)
# Input: genotypic value of overall expression rate, rate of modification (a vector, each element calculated using beta.calc()), decay rate of unmodified isoform (a number), decay rates of modified isoforms (a vector)
g2p <- function(alpha,beta,gamma_0,gamma){
	ni=length(beta) # Number of modified isoforms
	A=matrix(0,nrow=ni+1,ncol=ni+1)
	A[1,1]=sum(beta)+gamma_0
	A[2:nrow(A),1]=beta
	diag(A)[2:nrow(A)]=-gamma
	b=c(alpha,rep(0,ni))
	return(solve(A,b))
}

# Genotype-phenotype mapping (type-3 modification)
# Input: Input: genotypic value of overall expression level, rate of modification (a vector, each element calculated using beta.calc()), decay rate of unmodified isoform (a number), decay rates of modified isoforms (a vector)
g2p.type.3 <- function(alpha,beta,gamma_0,gamma){
	ni=length(beta) # Number of modified isoforms
	A=matrix(0,nrow=ni+1,ncol=ni+1)
	A[1,1]=beta[1]+gamma_0
	A[ni+1,ni]=beta[ni];A[ni+1,ni+1]=-gamma[ni]
	for(i in 2:ni){
		A[i,i-1]=beta[i-1]
		A[i,i]=-beta[i]-gamma[i-1]
	}
	b=c(alpha,rep(0,ni))
	return(solve(A,b))
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

# Obtain transition matrix given other parameters (type-2 modification)
# Input: number of cis- loci, mutation rates (a vector, elements being u01 and u10, respectively; 2*Ne*u should not exceed 1), selection parameter (a list), effective population size, other parameters (alpha,trans,C,gamma_0,gamma_1,epsilon,eq)
tm.type.2 <- function(nloci,u,selection,Ne,alpha,trans,C,epsilon=c(0,0),gamma_0,gamma,eq=NULL){
	mat=matrix(0,nrow=nloci+1,ncol=nloci+1)
	for(i in 1:nrow(mat)){
		v1=(i-1)/nloci # cis- genotypic value for beta_1
		beta1=beta.calc.type.2(v1,trans,C,epsilon)
		z1=g2p(alpha,beta1,gamma_0,gamma) # Remember to use the generic version of g2p()
		w1=fitness.prod(z1,selection,eq)
		for(j in 1:ncol(mat)){
			if(abs(i-j)==1){ # Consider adjescent states only
				v2=(j-1)/nloci
				beta2=beta.calc.type.2(v2,trans,C,epsilon)
				z2=g2p(alpha,beta2,gamma_0,gamma)
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

# Simulate cis-trans coevolution
# Two isoforms per gene, one functional and one deleterious
sim.evo<-function(nloci.all,alpha.all,gamma_0.all,gamma_1.all,selection.all,T,Ne,u.cis,u.trans,C,epsilon,start,opt.trans,sig.trans){
	nmut=rpois(1,lambda=u.trans[1]*2*Ne*T)
	if(nmut>=T){
		return("error")
	}else{
		ngene=length(nloci.all)
		tm.all=list()
		distr.all=list()
		v.all=start[[1]]
		trans=start[[2]]
		beta.all=rep(0,ngene)
		z.all=matrix(0,nrow=ngene,ncol=2)
		for(i in 1:ngene){
			tm.all[[i]]=tm(nloci.all[i],u.cis,selection.all[[i]],Ne,alpha.all[i],trans,C,epsilon,gamma_0.all[i],gamma_1.all[i])
			distr.all[[i]]=rep(0,nloci.all[i]+1);distr.all[[i]][v.all[i]+1]=1 # Starting distribution of cis- gv; probability of the assigned starting value is 1
			beta.all[i]=beta.calc(v.all[i]/nloci.all[i],trans,C,epsilon)
			z.all[i,]=g2p(alpha.all[i],beta.all[i],gamma_0.all[i],gamma_1.all[i])
		}
		if(nmut>0){
			tmut=sort(sample(1:T,nmut))
			t.last=0
			for(j in 1:nmut){
				t=tmut[j]
				w0=1
				for(i in 1:ngene){
					distr.all[[i]]=distr.all[[i]]%*%(tm.all[[i]] %^% (t-t.last)) # Calculate distribution at the time point under consideration
					v.all[i]=sample(0:nloci.all[i],1,prob=distr.all[[i]])
					w0=w0*fitness.prod(z.all[i,],selection.all[[i]])
				}
				w0.trans=fitness.stab(trans,opt.trans,sig.trans)
				w0=w0*w0.trans
				trans.mut=exp(log(trans)+rnorm(1,mean=0,sd=u.trans[2]))
				w1=1
				beta.mut=rep(0,ngene)
				z.mut=matrix(0,nrow=ngene,ncol=2)
				for(i in 1:ngene){
					beta.mut[i]=beta.calc(v.all[i]/nloci.all[i],trans.mut,C,epsilon)
					z.mut[i,]=g2p(alpha.all[i],beta.mut[i],gamma_0.all[i],gamma_1.all[i])
					w1=w1*fitness.prod(z.mut[i,],selection.all[[i]])
				}
				w1.trans=fitness.stab(trans.mut,opt.trans,sig.trans)
				w1=w1*w1.trans
				fp=fix.prob(w0,w1,Ne)
				if.fix=rbinom(n=1,size=1,prob=fp)
				if(if.fix==1){
					trans=trans.mut
					beta.all=beta.mut
					z.all=z.mut
					for(i in 1:ngene){
						tm.all[[i]]=tm(nloci.all[i],u.cis,selection.all[[i]],Ne,alpha.all[i],trans.mut,C,epsilon,gamma_0.all[i],gamma_1.all[i])
					}
				}
				t.last=t
			}
			for(i in 1:ngene){
				distr.all[[i]]=distr.all[[i]]%*%(tm.all[[i]] %^% (T-t.last)) # Calculate distribution at the time point under consideration
				v.all[i]=sample(0:nloci.all[i],1,prob=distr.all[[i]])
				beta.all[i]=beta.calc(v.all[i]/nloci.all[i],trans,C,epsilon)
				z.all[i,]=g2p(alpha.all[i],beta.all[i],gamma_0.all[i],gamma_1.all[i])
			}
		}
		return(list(z.all,trans,v.all))
	}
}

# Simulate cis-trans coevolution
# Two modified isoforms per gene, one functional and one deleterious
sim.evo.type.2<-function(nloci.all,alpha.all,gamma_0.all,gamma_1.all,selection.all,T,Ne,u.cis,u.gamma,trans,C,epsilon,start){
	nmut=rpois(1,lambda=u.gamma[1]*2*Ne*T)
	if(nmut>=T){
		return("error")
	}else{
		ngene=length(nloci.all)
		tm.all=list()
		distr.all=list()
		v.all=start[[1]]
		gamma_2=start[[2]]
		beta.all=matrix(0,nrow=ngene,ncol=2)
		z.all=matrix(0,nrow=ngene,ncol=3)
		for(i in 1:ngene){
			tm.all[[i]]=tm.type.2(nloci.all[i],u.cis,selection.all[[i]],Ne,alpha.all[i],trans,C,epsilon,gamma_0.all[i],c(gamma_1.all[i],gamma_2))
			distr.all[[i]]=rep(0,nloci.all[i]+1);distr.all[[i]][v.all[i]+1]=1 # Starting distribution of cis- gv; probability of the assigned starting value is 1
			beta.all[i,]=beta.calc.type.2(v.all[i]/nloci.all[i],trans,C,epsilon)
			z.all[i,]=g2p(alpha.all[i],beta.all[i,],gamma_0.all[i],c(gamma_1.all[i],gamma_2))
		}
		if(nmut>0){
			tmut=sort(sample(1:T,nmut))
			t.last=0
			for(j in 1:nmut){
				t=tmut[j]
				w0=1
				for(i in 1:ngene){
					distr.all[[i]]=distr.all[[i]]%*%(tm.all[[i]] %^% (t-t.last)) # Calculate distribution at the time point under consideration
					v.all[i]=sample(0:nloci.all[i],1,prob=distr.all[[i]])
					beta.all[i,]=beta.calc.type.2(v.all[i]/nloci.all[i],trans,C,epsilon)
					z.all[i,]=g2p(alpha.all[i],beta.all[i,],gamma_0.all[i],c(gamma_1.all[i],gamma_2))
					w0=w0*fitness.prod(z.all[i,],selection.all[[i]])
				}
				gamma_2.mut=exp(log(gamma_2)+rnorm(1,mean=0,sd=u.gamma[2]))
				w1=1
				z.mut=matrix(0,nrow=ngene,ncol=3)
				for(i in 1:ngene){
					z.mut[i,]=g2p(alpha.all[i],beta.all[i,],gamma_0.all[i],c(gamma_1.all[i],gamma_2.mut))
					w1=w1*fitness.prod(z.mut[i,],selection.all[[i]])
				}
				fp=fix.prob(w0,w1,Ne)
				if.fix=rbinom(n=1,size=1,prob=fp)
				if(if.fix==1){
					z.all=z.mut
					gamma_2=gamma_2.mut
					for(i in 1:ngene){
						tm.all[[i]]=tm.type.2(nloci.all[i],u.cis,selection.all[[i]],Ne,alpha.all[i],trans,C,epsilon,gamma_0.all[i],c(gamma_1.all[i],gamma_2))
					}
				}
				t.last=t
			}
			for(i in 1:ngene){
				distr.all[[i]]=distr.all[[i]]%*%(tm.all[[i]] %^% (T-t.last)) # Calculate distribution at the time point under consideration
				v.all[i]=sample(0:nloci.all[i],1,prob=distr.all[[i]])
				beta.all[i,]=beta.calc.type.2(v.all[i]/nloci.all[i],trans,C,epsilon)
				z.all[i,]=g2p(alpha.all[i],beta.all[i,],gamma_0.all[i],c(gamma_1.all[i],gamma_2))
			}
		}
		return(list(z.all,gamma_2,v.all))
	}
}

# Simulate evolution of the cis- genotypic value as a Markov process
# Simplest version: 2 isoforms, constant transition matrix (no trans- substitutions, no environmental change), all cis- loci have equal effect size (i.e., no. of states = no. of loci+1) and mutational spectrum
# Input: starting genotypic value, transition matrix, time
mc <- function(v0,mat,T){
	v=rep(0,nrow(mat));v[v0+1]=1;v=t(v) # Convert starting genotypic value to a distribution (corresponding element=1. others=0)
	for(t in 1:T){
		v=v%*%mat
	}
	phe=sample(0:(nrow(mat)-1),1,prob=v) # Sample phenotype at the end
	return(list(v,phe))
}

# Solve for stationary distribution
ss <- function(mat){
	e=eigen(mat)
	lvec=ginv(e$vectors)
	s=lvec[1,]/sum(lvec[1,])
	return(s)
}

