
#--- CITATION ------------------------------------------------------------------

#*******************************************************************************
# Hawkins, BL. 2019. Exploring the role of temperature and possible alternative 
#   stable states in brook trout (Salvelinus fontinalis) population structure.
#   M.S. Thesis. Univeristy of California San Diego.
#*******************************************************************************

#--- MODELS --------------------------------------------------------------------

# model for one normal component
one.peak.model <- "model {
sdf[1]~dunif(1,550)
m1~dunif(1,550)

m[1]<-m1

omega<-1 # prob a fish is in group 1

for(i in 1:nf) {
z[i]~dcat(omega) # stochastic switch for fish group assignment
d[i,1]~dnorm(m[z[i]],(1/ (sdf[z[i]]*sdf[z[i]]) )) # likelihood
}

}"

# model for two normal components
two.peak.model <- "model {
sdf[1]~dunif(1,550)
sdf[2]~dunif(1,550)
m1~dunif(1,550)
aj~dunif(50,200)

m[1]<-m1
m[2]<-m1+aj

omega~ddirich(c(1,1)) # prob a fish is in group 1 vs 2

for(i in 1:nf) {
z[i]~dcat(omega) # stochastic switch for fish group assignment
d[i,1]~dnorm(m[z[i]],(1/ (sdf[z[i]]*sdf[z[i]]) )) # likelihood
}

}"

# model for three normal components
three.peak.model <- "model {
sdf[1]~dunif(1,550)
sdf[2]~dunif(1,550)
sdf[3]~dunif(1,550)
m1~dunif(1,550)
aj1~dunif(50,200)
aj2~dunif(50,200)

m[1]<-m1
m[2]<-m1+aj1
m[3]<-m1+aj2

omega~ddirich(c(1,1,1)) # prob a fish is in group 1 vs 2 vs 3

for(i in 1:nf) {
z[i]~dcat(omega) # stochastic switch for fish group assignment
d[i,1]~dnorm(m[z[i]],(1/ (sdf[z[i]]*sdf[z[i]]) )) # likelihood
}

}"

# model for four normal components
four.peak.model <- "model {
sdf[1]~dunif(1,550)
sdf[2]~dunif(1,550)
sdf[3]~dunif(1,550)
sdf[4]~dunif(1,550)
m1~dunif(1,550)
aj1~dunif(50,200)
aj2~dunif(50,200)
aj3~dunif(50,200)

m[1]<-m1
m[2]<-m1+aj1
m[3]<-m1+aj2
m[4]<-m1+aj3

omega~ddirich(c(1,1,1,1)) # prob a fish is in group 1 vs 2 vs 3 vs 4

for(i in 1:nf) {
z[i]~dcat(omega) # stochastic switch for fish group assignment
d[i,1]~dnorm(m[z[i]],(1/ (sdf[z[i]]*sdf[z[i]]) )) # likelihood
}

}"

# model for five normal components
five.peak.model <- "model {
sdf[1]~dunif(1,550)
sdf[2]~dunif(1,550)
sdf[3]~dunif(1,550)
sdf[4]~dunif(1,550)
sdf[5]~dunif(1,550)
m1~dunif(1,550)
aj1~dunif(50,200)
aj2~dunif(50,200)
aj3~dunif(50,200)
aj4~dunif(50,200)

m[1]<-m1
m[2]<-m1+aj1
m[3]<-m1+aj2
m[4]<-m1+aj3
m[5]<-m1+aj4

omega~ddirich(c(1,1,1,1,1)) # prob a fish is in group 1 vs 2 vs 3 vs 4 vs 5

for(i in 1:nf) {
z[i]~dcat(omega) # stochastic switch for fish group assignment
d[i,1]~dnorm(m[z[i]],(1/ (sdf[z[i]]*sdf[z[i]]) )) # likelihood
}

}"

# model for six normal components
six.peak.model <- "model {
sdf[1]~dunif(1,550)
sdf[2]~dunif(1,550)
sdf[3]~dunif(1,550)
sdf[4]~dunif(1,550)
sdf[5]~dunif(1,550)
sdf[6]~dunif(1,550)
m1~dunif(1,550)
aj1~dunif(50,200)
aj2~dunif(50,200)
aj3~dunif(50,200)
aj4~dunif(50,200)
aj5~dunif(50,200)

m[1]<-m1
m[2]<-m1+aj1
m[3]<-m1+aj2
m[4]<-m1+aj3
m[5]<-m1+aj4
m[6]<-m1+aj5

omega~ddirich(c(1,1,1,1,1,1)) # prob a fish is in group 1 vs 2 vs 3 vs 4 vs 5 vs 6

for(i in 1:nf) {
z[i]~dcat(omega) # stochastic switch for fish group assignment
d[i,1]~dnorm(m[z[i]],(1/ (sdf[z[i]]*sdf[z[i]]) )) # likelihood
}

}"
