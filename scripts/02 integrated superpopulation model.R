# Modification of Jim Lyons's code (Lyons et al 2018)

library(R2WinBUGS)
library(jagsUI)
library(MCMCvis)

options(scipen=999)

# set location of WinBUGS
bugs.dir <- "C:/WinBUGS14"


# ----- Load marked ration scan samples -----
Scans <- read.csv('raw_data/SC scans.csv')

# ----- Load encounter histories -----
Eh <- read.csv('raw_data/SC eh cleaned.csv')
colnames(Eh)[1] <- "X1"

m <- as.numeric(Scans$marked)
K <- as.numeric(Scans$checked)
n.scans = length(m)

n.per <- length(unique(Scans$week))  # missing scan data for weeks 5 & 6
per <- as.numeric(factor(Scans$week))


# Augment
nz <- 600
CH <- as.matrix(Eh)
n.ch <- dim(CH)[1]
CH.aug <- rbind(CH, matrix(0, ncol=dim(CH)[2], nrow=nz))


# Known latent state
# Function to create a matrix with information about known latent state z for Jolly-Seber model with data augmentation
known.state.js <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    if (sum (ch[i,]) == 0) next       # Modified for jolly-seber w/ data aug
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
  }
  state[state==0] <- NA
  return(state)
}  

(str(bugs.data <- list(y=CH.aug, 
                       n.occasions=dim(CH.aug)[2], 
                       M=dim(CH.aug)[1], 
                       n.ch=dim(CH)[1],
                       m=m,
                       K=K, 
                       n=n.scans, 
                       per=per, n.per=n.per, 
                       z=known.state.js(CH.aug))))



# Occupancy parameterization (adapted from Kery and Schaub)
# Four different models are listed below -- 
# Corresponding to the four models in Lyons et al 2018 + one with added random effect

# ----- Model 1 t-t-t -----

## BUGS
cat(file='scripts/Modified Lyons model1.txt', "
# Restricted Dynamic Occupancy Model Parameterization
# Adapted from Kery and Schaub 2012

model {

# Priors and constraints
for (i in 1:M) {
  for (t in 1:(n.occasions-1)) {
    phi[i,t] <- phi_t[t]
    } # t
  for (t in 1:n.occasions) {
    p[i,t] <- psight[t]
    } 
  }

for (t in 1:(n.occasions-1)) {
  phi_t[t] ~ dunif(0,1)
  }

p.first ~ dunif(0,1)
p.last ~ dunif(0,1)
for (i in 1:2) {
  psight[i] <- p.first
  }
for (i in 3:11) {
  psight[i] ~ dunif(0,1)
  }
psight[12] <- p.last
psight[13] <- p.last
  
for (t in 1:n.occasions) {
  gamma[t] ~ dunif(0,1)
  } 

# Likelihood
for (i in 1:M) {
  # First occasion
  # State process
  z[i,1] ~ dbern(gamma[1])
  mu1[i] <- z[i,1]*p[i,1]
  # Observation process
  y[i,1] ~ dbern(mu1[i])

  # Subsequent occasions
  for (t in 2:n.occasions) {
    # State process
    q[i,t-1] <- 1 - z[i,t-1] # Availability for recruitment
    mu2[i, t] <- phi[i,t-1] * z[i,t-1] + gamma[t]*prod(q[i,1:(t-1)])
    z[i,t] ~ dbern(mu2[i,t])
    # Observation process
    mu3[i,t] <- z[i,t] * p[i,t]
    y[i,t] ~ dbern(mu3[i,t])
    } # t
  } # i
  
# GLMM for scan samples
alpha ~ dnorm(0,0.001)
for (i in 1:n.per) {
  delta[i] ~ dnorm(0, tau.per)
  }
tau.per <- pow(sigma.per, -2)
sigma.per ~ dunif(0,10)

for (i in 1:n) {
  m[i] ~ dbin(lp[i],K[i])
  logit(lp[i]) <- alpha + delta[ per[i] ]
  }
	
# ------------- scan samples ---------------------------------
for (i in 1:n.per) {
  pflag[i] <- exp(alpha + delta[i])/(1 + exp(alpha + delta[i]))
  }
# ------------- end : scan samples ---------------------------------
meanpflag <- exp(alpha)/(1 + exp(alpha))

# Calculate derived population parameters
for (t in 1:n.occasions) {
  qgamma [t] <- 1 - gamma[t]
  }
cprob[1] <- gamma[1]
for (t in 2:n.occasions) {
  cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
  } # t
psi <- sum(cprob[])                                 # Inclusion probability
for(t in 1:n.occasions) {
  b[t] <- cprob[t]/psi                              # Entry probability
  } # t
for(i in 1:M) {
  recruit[i,1] <- z[i,1]
  for (t in 2:n.occasions) {
    recruit[i,t] <- (1-z[i,t-1]) * z[i,t]
    }
  }

for(t in 1:n.occasions) {
  Nflag[t]   <- sum( z[1:M,t] )
  Nstop_t[t] <- sum( z[1:M,t] )/pflag[t]
  Bflag[t]   <- sum(recruit[1:M,t])                     
  Bstop_t[t] <- (sum(recruit[1:M,t]))/pflag[t]             
  } 
  
for(i in 1:M) {
  Nind[i] <- sum(z[i,1:n.occasions])
  Nalive[i] <- 1 - equals(Nind[i],0)
  }
   
Nsuperflag <- sum(Nalive[])  
Nsuperstop_Bt <- sum(Bstop_t[])


# ------------------- Stopover duration -----------------------------------
for (i in 1:n.ch) {
  zstop[i] <- sum(z[i,1:n.occasions])*7
  }
S.z <- mean(zstop[])
# ------------------- end : Stopover duration ---------------------------------


}
")


# ----- Model 2 t-t-c -----

## BUGS
cat(file='scripts/Modified Lyons model2.txt', "
# Restricted Dynamic Occupancy Model Parameterization
# Adapted from Kery and Schaub 2012

model {

# Priors and constraints
for (i in 1:M) {
  for (t in 1:(n.occasions-1)) {
    phi[i,t] <- phi_t[t]
    } # t
  for (t in 1:n.occasions) {
    p[i,t] <- mean.p
    } 
  }

for (t in 1:(n.occasions-1)) {
  phi_t[t] ~ dunif(0,1)
  }

mean.p ~ dunif(0, 1)

for (t in 1:n.occasions) {
  gamma[t] ~ dunif(0,1)
  } 

# Likelihood
for (i in 1:M) {
  # First occasion
  # State process
  z[i,1] ~ dbern(gamma[1])
  mu1[i] <- z[i,1]*p[i,1]
  # Observation process
  y[i,1] ~ dbern(mu1[i])

  # Subsequent occasions
  for (t in 2:n.occasions) {
    # State process
    q[i,t-1] <- 1 - z[i,t-1] # Availability for recruitment
    mu2[i, t] <- phi[i,t-1] * z[i,t-1] + gamma[t]*prod(q[i,1:(t-1)])
    z[i,t] ~ dbern(mu2[i,t])
    # Observation process
    mu3[i,t] <- z[i,t] * p[i,t]
    y[i,t] ~ dbern(mu3[i,t])
    } # t
  } # i
  
# GLMM for scan samples
alpha ~ dnorm(0,0.001)
for (i in 1:n.per) {
  delta[i] ~ dnorm(0, tau.per)
  }
tau.per <- pow(sigma.per, -2)
sigma.per ~ dunif(0,10)

for (i in 1:n) {
  m[i] ~ dbin(lp[i],K[i])
  logit(lp[i]) <- alpha + delta[ per[i] ]
  }
	
# ------------- scan samples ---------------------------------
for (i in 1:n.per) {
  pflag[i] <- exp(alpha + delta[i])/(1 + exp(alpha + delta[i]))
  }
# ------------- end : scan samples ---------------------------------
meanpflag <- exp(alpha)/(1 + exp(alpha))

# Calculate derived population parameters
for (t in 1:n.occasions) {
  qgamma [t] <- 1 - gamma[t]
  }
cprob[1] <- gamma[1]
for (t in 2:n.occasions) {
  cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
  } # t
psi <- sum(cprob[])                                 # Inclusion probability
for(t in 1:n.occasions) {
  b[t] <- cprob[t]/psi                              # Entry probability
  } # t
for(i in 1:M) {
  recruit[i,1] <- z[i,1]
  for (t in 2:n.occasions) {
    recruit[i,t] <- (1-z[i,t-1]) * z[i,t]
    }
  }

for(t in 1:n.occasions) {
  Nflag[t]   <- sum( z[1:M,t] )
  Nstop_t[t] <- sum( z[1:M,t] )/pflag[t]
  Bflag[t]   <- sum(recruit[1:M,t])                     
  Bstop_t[t] <- (sum(recruit[1:M,t]))/pflag[t]             
  } 
  
for(i in 1:M) {
  Nind[i] <- sum(z[i,1:n.occasions])
  Nalive[i] <- 1 - equals(Nind[i],0)
  }
   
Nsuperflag <- sum(Nalive[])  
Nsuperstop_Bt <- sum(Bstop_t[])


# ------------------- Stopover duration -----------------------------------
for (i in 1:n.ch) {
  zstop[i] <- sum(z[i,1:n.occasions])*7
  }
S.z <- mean(zstop[])
# ------------------- end : Stopover duration ---------------------------------


}
")
# end : BUGS




# ----- Model 3 t-c-c -----

## BUGS
cat(file='scripts/Modified Lyons model3.txt', "
# Restricted Dynamic Occupancy Model Parameterization
# Adapted from Kery and Schaub 2012

model {

# Priors and constraints
for (i in 1:M) {
  for (t in 1:(n.occasions-1)) {
    phi[i,t] <- mean.phi
    } # t
  for (t in 1:n.occasions) {
    p[i,t] <- mean.p
    } 
  }

mean.phi ~ dunif(0, 1)
mean.p ~ dunif(0, 1)

for (t in 1:n.occasions) {
  gamma[t] ~ dunif(0,1)
  } 

# Likelihood
for (i in 1:M) {
  # First occasion
  # State process
  z[i,1] ~ dbern(gamma[1])
  mu1[i] <- z[i,1]*p[i,1]
  # Observation process
  y[i,1] ~ dbern(mu1[i])

  # Subsequent occasions
  for (t in 2:n.occasions) {
    # State process
    q[i,t-1] <- 1 - z[i,t-1] # Availability for recruitment
    mu2[i, t] <- phi[i,t-1] * z[i,t-1] + gamma[t]*prod(q[i,1:(t-1)])
    z[i,t] ~ dbern(mu2[i,t])
    # Observation process
    mu3[i,t] <- z[i,t] * p[i,t]
    y[i,t] ~ dbern(mu3[i,t])
    } # t
  } # i
  
# GLMM for scan samples
alpha ~ dnorm(0,0.001)
for (i in 1:n.per) {
  delta[i] ~ dnorm(0, tau.per)
  }
tau.per <- pow(sigma.per, -2)
sigma.per ~ dunif(0,10)

for (i in 1:n) {
  m[i] ~ dbin(lp[i],K[i])
  logit(lp[i]) <- alpha + delta[ per[i] ]
  }
	
# ------------- scan samples ---------------------------------
for (i in 1:n.per) {
  pflag[i] <- exp(alpha + delta[i])/(1 + exp(alpha + delta[i]))
  }
# ------------- end : scan samples ---------------------------------
meanpflag <- exp(alpha)/(1 + exp(alpha))

# Calculate derived population parameters
for (t in 1:n.occasions) {
  qgamma [t] <- 1 - gamma[t]
  }
cprob[1] <- gamma[1]
for (t in 2:n.occasions) {
  cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
  } # t
psi <- sum(cprob[])                                 # Inclusion probability
for(t in 1:n.occasions) {
  b[t] <- cprob[t]/psi                              # Entry probability
  } # t
for(i in 1:M) {
  recruit[i,1] <- z[i,1]
  for (t in 2:n.occasions) {
    recruit[i,t] <- (1-z[i,t-1]) * z[i,t]
    }
  }

for(t in 1:n.occasions) {
  Nflag[t]   <- sum( z[1:M,t] )
  Nstop_t[t] <- sum( z[1:M,t] )/pflag[t]
  Bflag[t]   <- sum(recruit[1:M,t])                     
  Bstop_t[t] <- (sum(recruit[1:M,t]))/pflag[t]             
  } 
  
for(i in 1:M) {
  Nind[i] <- sum(z[i,1:n.occasions])
  Nalive[i] <- 1 - equals(Nind[i],0)
  }
   
Nsuperflag <- sum(Nalive[])  
Nsuperstop_Bt <- sum(Bstop_t[])


# ------------------- Stopover duration -----------------------------------
for (i in 1:n.ch) {
  zstop[i] <- sum(z[i,1:n.occasions])*7
  }
S.z <- mean(zstop[])
# ------------------- end : Stopover duration ---------------------------------


}
")
# end : BUGS



# ----- Model 4 t-c-t -----

## BUGS
cat(file='scripts/Modified Lyons model4.txt', "
# Restricted Dynamic Occupancy Model Parameterization
# Adapted from Kery and Schaub 2012

model {

# Priors and constraints
for (i in 1:M) {
  for (t in 1:(n.occasions-1)) {
    phi[i,t] <- mean.phi
    } # t
  for (t in 1:n.occasions) {
    p[i,t] <- psight[t]
    } 
  }

mean.phi ~ dunif(0, 1)

p.first ~ dunif(0,1)
p.last ~ dunif(0,1)
for (i in 1:2) {
  psight[i] <- p.first
  }
for (i in 3:11) {
  psight[i] ~ dunif(0,1)
  }
psight[12] <- p.last
psight[13] <- p.last

for (t in 1:n.occasions) {
  gamma[t] ~ dunif(0,1)
  } 

# Likelihood
for (i in 1:M) {
  # First occasion
  # State process
  z[i,1] ~ dbern(gamma[1])
  mu1[i] <- z[i,1]*p[i,1]
  # Observation process
  y[i,1] ~ dbern(mu1[i])

  # Subsequent occasions
  for (t in 2:n.occasions) {
    # State process
    q[i,t-1] <- 1 - z[i,t-1] # Availability for recruitment
    mu2[i, t] <- phi[i,t-1] * z[i,t-1] + gamma[t]*prod(q[i,1:(t-1)])
    z[i,t] ~ dbern(mu2[i,t])
    # Observation process
    mu3[i,t] <- z[i,t] * p[i,t]
    y[i,t] ~ dbern(mu3[i,t])
    } # t
  } # i
  
# GLMM for scan samples
alpha ~ dnorm(0,0.001)
for (i in 1:n.per) {
  delta[i] ~ dnorm(0, tau.per)
  }
tau.per <- pow(sigma.per, -2)
sigma.per ~ dunif(0,10)

for (i in 1:n) {
  m[i] ~ dbin(lp[i],K[i])
  logit(lp[i]) <- alpha + delta[ per[i] ]
  }
	
# ------------- scan samples ---------------------------------
for (i in 1:n.per) {
  pflag[i] <- exp(alpha + delta[i])/(1 + exp(alpha + delta[i]))
  }
# ------------- end : scan samples ---------------------------------
meanpflag <- exp(alpha)/(1 + exp(alpha))

# Calculate derived population parameters
for (t in 1:n.occasions) {
  qgamma [t] <- 1 - gamma[t]
  }
cprob[1] <- gamma[1]
for (t in 2:n.occasions) {
  cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
  } # t
psi <- sum(cprob[])                                 # Inclusion probability
for(t in 1:n.occasions) {
  b[t] <- cprob[t]/psi                              # Entry probability
  } # t
for(i in 1:M) {
  recruit[i,1] <- z[i,1]
  for (t in 2:n.occasions) {
    recruit[i,t] <- (1-z[i,t-1]) * z[i,t]
    }
  }

for(t in 1:n.occasions) {
  Nflag[t]   <- sum( z[1:M,t] )
  Nstop_t[t] <- sum( z[1:M,t] )/pflag[t]
  Bflag[t]   <- sum(recruit[1:M,t])                     
  Bstop_t[t] <- (sum(recruit[1:M,t]))/pflag[t]             
  } 
  
for(i in 1:M) {
  Nind[i] <- sum(z[i,1:n.occasions])
  Nalive[i] <- 1 - equals(Nind[i],0)
  }
   
Nsuperflag <- sum(Nalive[])  
Nsuperstop_Bt <- sum(Bstop_t[])


# ------------------- Stopover duration -----------------------------------
for (i in 1:n.ch) {
  zstop[i] <- sum(z[i,1:n.occasions])*7
  }
S.z <- mean(zstop[])
# ------------------- end : Stopover duration ---------------------------------


}
")
# end : BUGS


# ------ Model 5 t-t-t w/ random effects -------------


## BUGS
cat(file='scripts/Modified Lyons model5.txt', "
# Restricted Dynamic Occupancy Model Parameterization
# Adapted from Kery and Schaub 2012

model {

# Priors and constraints
for (i in 1:M) {
  for (t in 1:(n.occasions-1)) {
    phi[i,t] <- phi_t[t]
    } # t
  for (t in 1:n.occasions) {
    p[i,t] <- psight[t]
    } 
  }

# random effect for phi

for (t in 1:(n.occasions-1)) {
  logit(phi_t[t]) <- mu_phi + epsilon_phi[t]
  epsilon_phi[t] ~ dnorm(0, tau_phi)
  } 

mean.phi ~ dunif(0, 1) # prior for mean survival
mu_phi <- log(mean.phi / (1-mean.phi)) # logit transformation
sigma_phi ~ dunif(0, 5) # prior for standard deviation
tau_phi <- pow(sigma_phi, -2)
sigma_phi2 <- pow(sigma_phi, 2) # temporal variance
sigma_phi2.real <- sigma_phi2 * pow(mean.phi, 2) * pow((1-mean.phi), 2) # Temporal variance on real scale

# random effect for p

for (t in 1:n.occasions) {
  logit(psight[t]) <- mu_p + epsilon_p[t]
  epsilon_p[t] ~ dnorm(0, tau_p)
  } 

# Hyperpriors
mean.p ~ dunif(0, 1) # prior for mean resight probability
mu_p <- log(mean.p / (1-mean.p)) # logit transformation
sigma_p ~ dunif(0, 5) # prior for standard deviation
tau_p <- pow(sigma_p, -2)
sigma_p2 <- pow(sigma_p, 2) # temporal variance
sigma_p2.real <- sigma_p2 * pow(mean.p, 2) * pow((1-mean.p), 2) # Temporal variance on real scale


for (t in 1:n.occasions) {
  gamma[t] ~ dunif(0,1)
  } 

# Likelihood
for (i in 1:M) {
  # First occasion
  # State process
  z[i,1] ~ dbern(gamma[1])
  mu1[i] <- z[i,1]*p[i,1]
  # Observation process
  y[i,1] ~ dbern(mu1[i])

  # Subsequent occasions
  for (t in 2:n.occasions) {
    # State process
    q[i,t-1] <- 1 - z[i,t-1] # Availability for recruitment
    mu2[i, t] <- phi[i,t-1] * z[i,t-1] + gamma[t]*prod(q[i,1:(t-1)])
    z[i,t] ~ dbern(mu2[i,t])
    # Observation process
    mu3[i,t] <- z[i,t] * p[i,t]
    y[i,t] ~ dbern(mu3[i,t])
    } # t
  } # i
  
# GLMM for scan samples
alpha ~ dnorm(0,0.001)
for (i in 1:n.per) {
  delta[i] ~ dnorm(0, tau.per)
  }
tau.per <- pow(sigma.per, -2)
sigma.per ~ dunif(0,10)

for (i in 1:n) {
  m[i] ~ dbin(lp[i],K[i])
  logit(lp[i]) <- alpha + delta[ per[i] ]
  }
	
# ------------- scan samples ---------------------------------
for (i in 1:n.per) {
  pflag[i] <- exp(alpha + delta[i])/(1 + exp(alpha + delta[i]))
  }
# ------------- end : scan samples ---------------------------------
meanpflag <- exp(alpha)/(1 + exp(alpha))

# Calculate derived population parameters
for (t in 1:n.occasions) {
  qgamma [t] <- 1 - gamma[t]
  }
cprob[1] <- gamma[1]
for (t in 2:n.occasions) {
  cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
  } # t
psi <- sum(cprob[])                                 # Inclusion probability
for(t in 1:n.occasions) {
  b[t] <- cprob[t]/psi                              # Entry probability
  } # t
for(i in 1:M) {
  recruit[i,1] <- z[i,1]
  for (t in 2:n.occasions) {
    recruit[i,t] <- (1-z[i,t-1]) * z[i,t]
    }
  }

for(t in 1:n.occasions) {
  Nflag[t]   <- sum( z[1:M,t] )
  Nstop_t[t] <- sum( z[1:M,t] )/pflag[t]
  Bflag[t]   <- sum(recruit[1:M,t])                     
  Bstop_t[t] <- (sum(recruit[1:M,t]))/pflag[t]             
  } 
  
for(i in 1:M) {
  Nind[i] <- sum(z[i,1:n.occasions])
  Nalive[i] <- 1 - equals(Nind[i],0)
  }
   
Nsuperflag <- sum(Nalive[])  
Nsuperstop_Bt <- sum(Bstop_t[])


# ------------------- Stopover duration -----------------------------------
for (i in 1:n.ch) {
  zstop[i] <- sum(z[i,1:n.occasions])*7
  }
S.z <- mean(zstop[])
# ------------------- end : Stopover duration ---------------------------------

}
")


# ----------- Prepare & Run Model --------------------------



# Initial values

# Create vector with occasion of marking (Adapted for Jolly-Seber--jel)
js.get.first <- function(x) ifelse(sum(x)==0, 0, min(which(x != 0)))
f <- apply(CH.aug, 1, js.get.first)
## Do not give initial values for elements of z specified as data; they get NA.
js.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (f[i] == 0) next
    if (sum(ch[i,]) == 1) next
    n2 <- max(which(ch[i,] == 1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    if (f[i] == 0) next
    ch[i,f[i]] <- NA
  }
  return(ch)
}

# Initial values Model4
inits <- function() {list(mean.phi=runif(1, 0, 1),
                          p.first=runif(1,0,1),
                          p.last=runif(1,0,1),
                          gamma=c(runif(dim(CH.aug)[2],0,0.2)),
                          z=js.init.z(CH.aug,f), 
                          alpha=rnorm(1,0,1), 
                          sigma.per=runif(1))
}




# MCMC settings
# ni <- 10 ; nb <- 2 ; nt <- 1 ; nc <- 3
# ni <- 10000 ; nb <- 5000 ; nt <- 1 ; nc <- 3  # about 25 min.
# ni <- 20000 ; nb <- 10000 ; nt <- 1 ; nc <- 3 # about 45 min.
# ni <- 40000 ; nb <- 20000 ; nt <- 1 ; nc <- 2
ni <- 90000 ; nb <- 30000 ; nt <- 3 ; nc <- 3

# Parameters
# mean.phi or phi_t
# mean.p or psight

parameters <- c("psi",
                "Nsuperflag","Nsuperstop_Bt",
                "b","mean.phi","psight", "alpha","sigma.per",
                "meanpflag","pflag","S.z",
                "Nflag","Nstop_t",
                "Bflag","Bstop_t")


# Run in JAGS


mod4 <- jagsUI::jags(bugs.data, inits, parameters, 'scripts/Modified Lyons model4.txt', n.chains=nc,
             n.thin=nt, n.iter=ni, n.burnin=nb, n.cores = 7, parallel=TRUE)
mod4


# Inspect model results

jagsUI::traceplot(mod1)
jagsUI::traceplot(mod4, parameters = c("Nsuperflag", "S.z"))
jagsUI::traceplot(mod4, parameters = c("zes"))
jagsUI::densityplot(mod4, parameters = c("Nsuperflag", "zes", "Nsuperstop", "phi.mean"))

jagsUI::summary(mod1)


# Save results

mod4_df <- MCMCsummary(mod4)
options(scipen = 999)
mod4_df

write.csv(mod4_df, "processed_data/Final Results.csv", row.names = TRUE)

