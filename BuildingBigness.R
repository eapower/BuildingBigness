# This includes code to:
# (1) generate the networks, 
# (2) run the ERGMs, 
# (3) visualize the results (plots of fitted probabilities), 
# (4) run the binomial regressions predicting external alters of "high position", 
# (5) produce scatter plots of reputation and weighted PageRank centrality, 
# (6) run PCAs across the reputational qualities.

require(igraph)

Ten <- read.csv("TenMetadata.csv",header=TRUE)
Ala <- read.csv("AlaMetadata.csv",header=TRUE)


##############
## Tenpatti ##
##############

elTen <- read.csv("TenEdgeList.csv")

elSocTen <- elTen[,c(1:14,23)]
elSocTen <- subset(elSocTen, (ImpIss > 0 | Work > 0 | Errand > 0 | Borrow > 0 | Cash > 0 | Loan > 0 | Babysit > 0 | Talk > 0 | Defend > 0 | Position > 0 | Close > 0 | Advice > 0))

snaTen <- graph.data.frame(elSocTen)
V(snaTen)$label <- V(snaTen)$name
IndivIDTen <- V(snaTen)$label
IndivIDTen <- data.frame(IndivIDTen)

attTen <- merge(IndivIDTen,Ten,by.x="IndivIDTen",by.y="IndivID",sort=FALSE,all.x=TRUE)
colnames(attTen)[1] <- "IndivID"

snaTenSupFull <- graph.data.frame(d = elSocTen, vertices = attTen, directed=TRUE)

## reduce the network down to include only those who completed the survey in the village
snaTenSup <- delete.vertices(snaTenSupFull,V(snaTenSupFull)[degree(snaTenSupFull,mode="out")==0])

## To get networks of reputational nominations
snaTenFull <- graph.data.frame(d = elTen, directed=TRUE)
snaTenFull <- delete.vertices(snaTenFull,V(snaTenFull)[degree(snaTenFull,mode="out")==0])

## This represents a tie from A to B, if A nominated B
snaTenGenerous <- delete.edges(snaTenFull, E(snaTenFull)[get.edge.attribute(snaTenFull, name = "Generous") == 0])
snaTenAllAdvice <- delete.edges(snaTenFull, E(snaTenFull)[get.edge.attribute(snaTenFull, name = "AllAdvice") == 0])
snaTenInfluence <- delete.edges(snaTenFull, E(snaTenFull)[get.edge.attribute(snaTenFull, name = "Influence") == 0])
snaTenCharacter <- delete.edges(snaTenFull, E(snaTenFull)[get.edge.attribute(snaTenFull, name = "Character") == 0])
snaTenStrength <- delete.edges(snaTenFull, E(snaTenFull)[get.edge.attribute(snaTenFull, name = "Strength") == 0])

## This represents a tie from A to B, if B nominated A as 
snaTenGenerousI <- graph.adjacency(t(get.adjacency(snaTenGenerous,sparse=FALSE)))
snaTenAllAdviceI <- graph.adjacency(t(get.adjacency(snaTenAllAdvice,sparse=FALSE)))
snaTenInfluenceI <- graph.adjacency(t(get.adjacency(snaTenInfluence,sparse=FALSE)))
snaTenCharacterI <- graph.adjacency(t(get.adjacency(snaTenCharacter,sparse=FALSE)))
snaTenStrengthI <- graph.adjacency(t(get.adjacency(snaTenStrength,sparse=FALSE)))

#################
## Azhagapuram ##
#################

elAla <- read.csv("AlaEdgeList.csv")

elSocAla <- elAla[,c(1:14,23)]
elSocAla <- subset(elSocAla, (ImpIss > 0 | Work > 0 | Errand > 0 | Borrow > 0 | Cash > 0 | Loan > 0 | Babysit > 0 | Talk > 0 | Defend > 0 | Position > 0 | Close > 0 | Advice > 0))


snaAla <- graph.data.frame(elSocAla)
V(snaAla)$label <- V(snaAla)$name
IndivIDAla <- V(snaAla)$label
IndivIDAla <- data.frame(IndivIDAla)

attAla <- merge(IndivIDAla,Ala,by.x="IndivIDAla",by.y="IndivID",sort=FALSE,all.x=TRUE)
colnames(attAla)[1] <- "IndivID"

snaAlaSupFull <- graph.data.frame(d = elSocAla, vertices = attAla, directed=TRUE)

## reduce the network down to include only those who completed the survey in the village
snaAlaSup <- delete.vertices(snaAlaSupFull,V(snaAlaSupFull)[degree(snaAlaSupFull,mode="out")==0])

## To get networks of reputational nominations
snaAlaFull <- graph.data.frame(d = elAla, directed=TRUE)
snaAlaFull <- delete.vertices(snaAlaFull,V(snaAlaFull)[degree(snaAlaFull,mode="out")==0])

snaAlaGenerous <- delete.edges(snaAlaFull, E(snaAlaFull)[get.edge.attribute(snaAlaFull, name = "Generous") == 0])
snaAlaAllAdvice <- delete.edges(snaAlaFull, E(snaAlaFull)[get.edge.attribute(snaAlaFull, name = "AllAdvice") == 0])
snaAlaInfluence <- delete.edges(snaAlaFull, E(snaAlaFull)[get.edge.attribute(snaAlaFull, name = "Influence") == 0])
snaAlaCharacter <- delete.edges(snaAlaFull, E(snaAlaFull)[get.edge.attribute(snaAlaFull, name = "Character") == 0])
snaAlaStrength <- delete.edges(snaAlaFull, E(snaAlaFull)[get.edge.attribute(snaAlaFull, name = "Strength") == 0])

snaAlaGenerousI <- graph.adjacency(t(get.adjacency(snaAlaGenerous,sparse=FALSE)))
snaAlaAllAdviceI <- graph.adjacency(t(get.adjacency(snaAlaAllAdvice,sparse=FALSE)))
snaAlaInfluenceI <- graph.adjacency(t(get.adjacency(snaAlaInfluence,sparse=FALSE)))
snaAlaCharacterI <- graph.adjacency(t(get.adjacency(snaAlaCharacter,sparse=FALSE)))
snaAlaStrengthI <- graph.adjacency(t(get.adjacency(snaAlaStrength,sparse=FALSE)))


######################
## Kinship Networks ##
######################

closekinEL <- read.csv("CloseKinEdgeList.csv",header=TRUE)

require(car)

kinship<-function(g){
  kin <- V(g)$name
  kinel <- expand.grid(kin,kin)
  colnames(kinel) <- c("Ego","Alter")
  kinel$Edge <- paste(kinel[,1],kinel[,2],sep="_")
  kin1 <- merge(kinel,closekinEL[,c(3,4)],by="Edge",all.x=TRUE)
  kin1$Kin <- recode(kin1$Kin,"NA=0")
  kin1$Edge<- NULL
  knet <- graph.data.frame(kin1)
  return(knet)
}

knetTen<-kinship(snaTenSup)
knetAla<-kinship(snaAlaSup)


########################
## Household Distance ##
########################

TenIndivDist<-read.csv("TenIndivDist.csv",as.is=TRUE)
AlaIndivDist<-read.csv("AlaIndivDist.csv",as.is=TRUE)

## first making empty matrix with names as the House ID
Tendistancemat <-matrix(nrow=length(V(snaTenSup)$name),ncol=length(V(snaTenSup)$name),dimnames=list(as.vector(V(snaTenSup)$name),as.vector(V(snaTenSup)$name)))

## fill in the values from the Distance file. This will take a moment.
for(i in 1:nrow(TenIndivDist)){Tendistancemat[rownames(Tendistancemat)==TenIndivDist[i,1],colnames(Tendistancemat)==TenIndivDist[i,2]] <- as.numeric(TenIndivDist[i,3])}

## Same for Alagapuram
Aladistancemat <-matrix(nrow=length(V(snaAlaSup)$name),ncol=length(V(snaAlaSup)$name),dimnames=list(as.vector(V(snaAlaSup)$name),as.vector(V(snaAlaSup)$name)))
for(i in 1:nrow(AlaIndivDist)){Aladistancemat[rownames(Aladistancemat)==AlaIndivDist[i,1],colnames(Aladistancemat)==AlaIndivDist[i,2]] <- as.numeric(AlaIndivDist[i,3])}

########################
## PORTING TO STATNET ##
########################

detach(package:igraph)
require(intergraph)
require(network)
require(ergm)

Net_snaTenSup <- asNetwork(snaTenSup)
Net_kinTen <- asNetwork(knetTen)
Net_snaTenGenerous <- asNetwork(snaTenGenerous)
Net_snaTenAllAdvice <- asNetwork(snaTenAllAdvice)
Net_snaTenInfluence <- asNetwork(snaTenInfluence)
Net_snaTenCharacter <- asNetwork(snaTenCharacter)
Net_snaTenStrength <- asNetwork(snaTenStrength)

Net_snaTenGenerousI <- asNetwork(snaTenGenerousI)
Net_snaTenAllAdviceI <- asNetwork(snaTenAllAdviceI)
Net_snaTenInfluenceI <- asNetwork(snaTenInfluenceI)
Net_snaTenCharacterI <- asNetwork(snaTenCharacterI)
Net_snaTenStrengthI <- asNetwork(snaTenStrengthI)


Net_snaAlaSup <- asNetwork(snaAlaSup)
Net_kinAla <- asNetwork(knetAla)
Net_snaAlaGenerous <- asNetwork(snaAlaGenerous)
Net_snaAlaAllAdvice <- asNetwork(snaAlaAllAdvice)
Net_snaAlaInfluence <- asNetwork(snaAlaInfluence)
Net_snaAlaCharacter <- asNetwork(snaAlaCharacter)
Net_snaAlaStrength <- asNetwork(snaAlaStrength)

Net_snaAlaGenerousI <- asNetwork(snaAlaGenerousI)
Net_snaAlaAllAdviceI <- asNetwork(snaAlaAllAdviceI)
Net_snaAlaInfluenceI <- asNetwork(snaAlaInfluenceI)
Net_snaAlaCharacterI <- asNetwork(snaAlaCharacterI)
Net_snaAlaStrengthI <- asNetwork(snaAlaStrengthI)


#########################################################################
#########################################################################
########## FROM HERE RUNNING FULL MODELS THAT THAT A LONG TIME ##########
#########################################################################
#########################################################################

#modelTenSup <- ergm(Net_snaTenSup ~ edges + nodecov("Age") + nodematch("Gender") + nodeifactor("Gender") + edgecov(Net_kinTen,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + nodeicov("HouseholdWealthCalc") + absdiff("EduYears") + edgecov(Tendistancemat/10)+ mutual + gwesp(0.5,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
#summary(modelTenSup)

#modelTenSup.WithGen <- ergm(Net_snaTenSup ~ edges + nodecov("Age") + nodematch("Gender") + nodeifactor("Gender") + edgecov(Net_kinTen,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + nodeicov("HouseholdWealthCalc") + absdiff("EduYears") + edgecov(Tendistancemat/10) + nodeicov("GenerousScaled") + nodeocov("GenerousScaled") + edgecov(Net_snaTenGenerous,attrname="Generous") + edgecov(Net_snaTenGenerousI) + mutual + gwesp(0.5,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
#summary(modelTenSup.WithGen)
#saveRDS(modelTenSup.WithGen,"modelTenSupWithGen.rds")

#modelTenSup.WithChar <- ergm(Net_snaTenSup ~ edges + nodecov("Age") + nodematch("Gender") + nodeifactor("Gender") + edgecov(Net_kinTen,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + nodeicov("HouseholdWealthCalc") + absdiff("EduYears") + edgecov(Tendistancemat/10) + nodeicov("CharacterScaled") + nodeocov("CharacterScaled") + edgecov(Net_snaTenCharacter,attrname="Character") + edgecov(Net_snaTenCharacterI) + mutual + gwesp(0.5,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
#summary(modelTenSup.WithChar)
#saveRDS(modelTenSup.WithChar,"modelTenSupWithChar.rds")

#modelTenSup.WithAdv <- ergm(Net_snaTenSup ~ edges + nodecov("Age") + nodematch("Gender") + nodeifactor("Gender") + edgecov(Net_kinTen,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + nodeicov("HouseholdWealthCalc") + absdiff("EduYears") + edgecov(Tendistancemat/10) + nodeicov("AllAdviceScaled") + nodeocov("AllAdviceScaled") + edgecov(Net_snaTenAllAdvice,attrname="AllAdvice") + edgecov(Net_snaTenAllAdviceI) + mutual + gwesp(0.5,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
#summary(modelTenSup.WithAdv)
#saveRDS(modelTenSup.WithAdv,"modelTenSupWithAdv.rds")

#modelTenSup.WithInfl <- ergm(Net_snaTenSup ~ edges + nodecov("Age") + nodematch("Gender") + nodeifactor("Gender") + edgecov(Net_kinTen,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + nodeicov("HouseholdWealthCalc") + absdiff("EduYears") + edgecov(Tendistancemat/10) + nodeicov("InfluenceScaled") + nodeocov("InfluenceScaled") + edgecov(Net_snaTenInfluence,attrname="Influence") + edgecov(Net_snaTenInfluenceI) + mutual + gwesp(0.5,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
#summary(modelTenSup.WithInfl)
#saveRDS(modelTenSup.WithInfl,"modelTenSupWithInfl.rds")





#modelAlaSup <- ergm(Net_snaAlaSup ~ edges + nodecov("Age") + nodematch("Gender") + nodeifactor("Gender") + edgecov(Net_kinAla,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + nodeicov("HouseholdWealthCalc") + absdiff("EduYears") + edgecov(Aladistancemat/10)+ mutual + gwesp(0.6,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
#summary(modelAlaSup)

#modelAlaSup.WithGen <- ergm(Net_snaAlaSup ~ edges + nodecov("Age") + nodematch("Gender") + nodeifactor("Gender") + edgecov(Net_kinAla,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + nodeicov("HouseholdWealthCalc") + absdiff("EduYears") + edgecov(Aladistancemat/10) + nodeicov("GenerousScaled") + nodeocov("GenerousScaled") + edgecov(Net_snaAlaGenerous,attrname="Generous") + edgecov(Net_snaAlaGenerousI) + mutual + gwesp(0.6,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
#summary(modelAlaSup.WithGen)
#saveRDS(modelAlaSup.WithGen,"modelAlaSupWithGen.rds")

#modelAlaSup.WithChar <- ergm(Net_snaAlaSup ~ edges + nodecov("Age") + nodematch("Gender") + nodeifactor("Gender") + edgecov(Net_kinAla,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + nodeicov("HouseholdWealthCalc") + absdiff("EduYears") + edgecov(Aladistancemat/10) + nodeicov("CharacterScaled") + nodeocov("CharacterScaled") + edgecov(Net_snaAlaCharacter,attrname="Character") + edgecov(Net_snaAlaCharacterI) + mutual + gwesp(0.6,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
#summary(modelAlaSup.WithChar)
#saveRDS(modelAlaSup.WithChar,"modelAlaSupWithChar.rds")

#modelAlaSup.WithAdv <- ergm(Net_snaAlaSup ~ edges + nodecov("Age") + nodematch("Gender") + nodeifactor("Gender") + edgecov(Net_kinAla,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + nodeicov("HouseholdWealthCalc") + absdiff("EduYears") + edgecov(Aladistancemat/10) + nodeicov("AllAdviceScaled") + nodeocov("AllAdviceScaled") + edgecov(Net_snaAlaAllAdvice,attrname="AllAdvice") + edgecov(Net_snaAlaAllAdviceI) + mutual + gwesp(0.6,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
#summary(modelAlaSup.WithAdv)
#saveRDS(modelAlaSup.WithAdv,"modelAlaSupWithAdv.rds")

#modelAlaSup.WithInfl <- ergm(Net_snaAlaSup ~ edges + nodecov("Age") + nodematch("Gender") + nodeifactor("Gender") + edgecov(Net_kinAla,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + nodeicov("HouseholdWealthCalc") + absdiff("EduYears") + edgecov(Aladistancemat/10) + nodeicov("InfluenceScaled") + nodeocov("InfluenceScaled") + edgecov(Net_snaAlaInfluence,attrname="Influence") + edgecov(Net_snaAlaInfluenceI) + mutual + gwesp(0.6,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
#summary(modelAlaSup.WithInfl)
#saveRDS(modelAlaSup.WithInfl,"modelAlaSupWithInfl.rds")


#modelTenSup.WithGenInfl <- ergm(Net_snaTenSup ~ edges + nodecov("Age") + nodematch("Gender") + nodeifactor("Gender") + edgecov(Net_kinTen,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + nodeicov("HouseholdWealthCalc") + absdiff("EduYears") + edgecov(Tendistancemat/10) + nodeicov("GenerousScaled") + nodeocov("GenerousScaled") + edgecov(Net_snaTenGenerous,attrname="Generous") + edgecov(Net_snaTenGenerousI) + nodeicov("InfluenceScaled") + nodeocov("InfluenceScaled") + edgecov(Net_snaTenInfluence,attrname="Influence") + edgecov(Net_snaTenInfluenceI) + mutual + gwesp(0.5,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
#summary(modelTenSup.WithGenInfl)
#saveRDS(modelTenSup.WithGenInfl,"modelTenSupWithGenInfl.rds")

#modelAlaSup.WithGenInfl <- ergm(Net_snaAlaSup ~ edges + nodecov("Age") + nodematch("Gender") + nodeifactor("Gender") + edgecov(Net_kinAla,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + nodeicov("HouseholdWealthCalc") + absdiff("EduYears") + edgecov(Aladistancemat/10) + nodeicov("GenerousScaled") + nodeocov("GenerousScaled") + edgecov(Net_snaAlaGenerous,attrname="Generous") + edgecov(Net_snaAlaGenerousI) + nodeicov("InfluenceScaled") + nodeocov("InfluenceScaled") + edgecov(Net_snaAlaInfluence,attrname="Influence") + edgecov(Net_snaAlaInfluenceI) + mutual + gwesp(0.6,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
#summary(modelAlaSup.WithGenInfl)
#saveRDS(modelAlaSup.WithGenInfl,"modelAlaSupWithGenInfl.rds")


#### odds ratio, reporting
require(xtable)
require(texreg)

texreg(list(T_gen, T_char, T_alladv, T_infl),digits=3)
texreg(list(A_gen, A_char, A_alladv, A_infl),digits=3)

# specify model
model = T_infl

or <- exp( model$coef ) 
ste <- sqrt( diag( model$covar ) ) 
lci <- exp( model$coef-1.96*ste ) 
uci <- exp( model$coef+1.96*ste ) 
oddsratios <- rbind( round( lci,digits = 4 ),round( or,digits = 4 ),round( uci,digits = 4 ) ) 
oddsratios <- t( oddsratios ) 
colnames( oddsratios ) <- c( "Lower","OR","Upper" )
teststat <- model$coef/ste 
teststats <- rbind( round( teststat,digits = 4 )) 
teststats <- t( teststats ) 
colnames( teststats ) <- c("Wald")

results<-cbind(summary(model)$coefs[1],summary(model)$coefs[2],or,summary(model)$coefs[4])
xtable(results,digits=c(3,3,3,3,4))





#####PLOTTING

## plots of model predictions that show the shift in probability with different parameters
# TERGMpredict_building_giving.csv has two women, both aged 42, same edu, both of the Pallar caste, etc.
# nothing but indirect reputation, plus direct reputation, plus mutual, plus mutual reputation nomination
# TERGMpredict_building_requesting.csv changes the nodeicov and nodeocov reputation terms
# For adding sharing partners; a tie between two nodes with no shared partners has a GWESP change statistic of 0; a tie that closes one triangle has a GWESP change statistic of 1 (one triangle contains three pairs of nodes with one shared partner)

Tpredergm<-read.csv("TERGMpredict_building_requesting.csv",header=TRUE,as.is=TRUE)
#Tpredergm<-read.csv("TERGMpredict_building_giving.csv",header=TRUE,as.is=TRUE)


#plots of model predictions

Tgen_coefs <- T_gen$coef

Tinfl_coefs <- T_infl$coef

Tchar_coefs <- T_char$coef

Talladv_coefs <- T_alladv$coef

estoprob <- function(b) {
  exp(b)/(1+exp(b))
}


Tgen_pred_vect = Tchar_pred_vect = Tinfl_pred_vect = Talladv_pred_vect = rep(0,20)
for (i in 2:52) {
  Tgen_pred_vect[i-1]=estoprob(sum(as.numeric(Tpredergm[,i])*Tgen_coefs))
  Tinfl_pred_vect[i-1]=estoprob(sum(as.numeric(Tpredergm[,i])*Tinfl_coefs))
  Tchar_pred_vect[i-1]=estoprob(sum(as.numeric(Tpredergm[,i])*Tchar_coefs))
  Talladv_pred_vect[i-1]=estoprob(sum(as.numeric(Tpredergm[,i])*Talladv_coefs))
}

require(RColorBrewer)
colors <- brewer.pal(8,"Dark2")
ygen <- c(Tgen_pred_vect[1:11], NA, Tgen_pred_vect[12:21], NA, Tgen_pred_vect[22:31], NA, Tgen_pred_vect[32:41], NA, Tgen_pred_vect[42:51])
yinfl <- c(Tinfl_pred_vect[1:11], NA, Tinfl_pred_vect[12:21], NA, Tinfl_pred_vect[22:31], NA, Tinfl_pred_vect[32:41], NA, Tinfl_pred_vect[42:51])
ychar <- c(Tchar_pred_vect[1:11], NA, Tchar_pred_vect[12:21], NA, Tchar_pred_vect[22:31], NA, Tchar_pred_vect[32:41], NA, Tchar_pred_vect[42:51])
yalladv <- c(Talladv_pred_vect[1:11], NA, Talladv_pred_vect[12:21], NA, Talladv_pred_vect[22:31], NA, Talladv_pred_vect[32:41], NA, Talladv_pred_vect[42:51])
x <- c(0:54)


pdf("Tenpatti_requesting.pdf", height=5, width=10)
#pdf("Tenpatti_giving.pdf", height=5, width=10)
plot(x, ygen, type="l",col=colors[1],lwd=3,ylab="Probability of Tie",xlab="Model Specifications",
     main="Predicted Probability of a Support Tie for an Older Pallar Woman in Tenpatti",
     ylim=c(0,1),xaxt="n")
abline(0,0,lty=3,col="gray")
abline(0.2,0,lty=3,col="gray")
abline(0.4,0,lty=3,col="gray")
abline(0.6,0,lty=3,col="gray")
abline(0.8,0,lty=3,col="gray")
abline(1,0,lty=3,col="gray")
lines(x,yinfl,type="l",col=colors[4],lwd=3,ylim=c(0,1))
lines(x,ychar,type="l",col=colors[2],lwd=3,ylim=c(0,1))
lines(x,yalladv,type="l",col=colors[3],lwd=3,ylim=c(0,1))
legend("topleft",c("Generous","Good Character","Gives Good Advice","Influential"),cex=0.6,col=colors[1:4],lwd=3)
labloc=c(0:54)
lab=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,NA,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,NA,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,NA,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,NA,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
text(x=labloc,y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[4]/3),labels=lab,cex=0.6,pos=3,adj=2,xpd=TRUE)
dev.off()


## same for Alagapuram

Apredergm<-read.csv("AERGMpredict_building_requesting.csv",header=TRUE,as.is=TRUE)
#Apredergm<-read.csv("AERGMpredict_building_giving.csv",header=TRUE,as.is=TRUE)

#plots of model predictions

Agen_coefs <- A_gen$coef

Ainfl_coefs <- A_infl$coef

Achar_coefs <- A_char$coef

Aalladv_coefs <- A_alladv$coef

Agen_pred_vect = Achar_pred_vect = Ainfl_pred_vect = Aalladv_pred_vect = rep(0,20)
for (i in 2:52) {
  Agen_pred_vect[i-1]=estoprob(sum(as.numeric(Apredergm[,i])*Agen_coefs))
  Ainfl_pred_vect[i-1]=estoprob(sum(as.numeric(Apredergm[,i])*Ainfl_coefs))
  Achar_pred_vect[i-1]=estoprob(sum(as.numeric(Apredergm[,i])*Achar_coefs))
  Aalladv_pred_vect[i-1]=estoprob(sum(as.numeric(Apredergm[,i])*Aalladv_coefs))
}

ygen <- c(Agen_pred_vect[1:11], NA, Agen_pred_vect[12:21], NA, Agen_pred_vect[22:31], NA, Agen_pred_vect[32:41], NA, Agen_pred_vect[42:51])
yinfl <- c(Ainfl_pred_vect[1:11], NA, Ainfl_pred_vect[12:21], NA, Ainfl_pred_vect[22:31], NA, Ainfl_pred_vect[32:41], NA, Ainfl_pred_vect[42:51])
ychar <- c(Achar_pred_vect[1:11], NA, Achar_pred_vect[12:21], NA, Achar_pred_vect[22:31], NA, Achar_pred_vect[32:41], NA, Achar_pred_vect[42:51])
yalladv <- c(Aalladv_pred_vect[1:11], NA, Aalladv_pred_vect[12:21], NA, Aalladv_pred_vect[22:31], NA, Aalladv_pred_vect[32:41], NA, Aalladv_pred_vect[42:51])
x <- c(0:54)

#pdf("Alaga_giving.pdf", width=10, height=5)
pdf("Alaga_requesting.pdf", width=10, height=5)
plot(x, ygen, type="l",col=colors[1],lwd=3,ylab="Probability of Tie",xlab="Model Specifications",
     main="Predicted Probability of a Support Tie for an Older Pallar Woman in Alagapuram",
     ylim=c(0,1),xaxt="n")
abline(0,0,lty=3,col="gray")
abline(0.2,0,lty=3,col="gray")
abline(0.4,0,lty=3,col="gray")
abline(0.6,0,lty=3,col="gray")
abline(0.8,0,lty=3,col="gray")
abline(1,0,lty=3,col="gray")
lines(x,yinfl,type="l",col=colors[4],lwd=3,ylim=c(0,1))
lines(x,ychar,type="l",col=colors[2],lwd=3,ylim=c(0,1))
lines(x,yalladv,type="l",col=colors[3],lwd=3,ylim=c(0,1))
legend("topleft",c("Generous","Good Character","Gives Good Advice","Influential"),cex=0.6,col=colors[1:4],lwd=3)
labloc=c(0:54)
lab=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,NA,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,NA,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,NA,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,NA,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
text(x=labloc,y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[4]/3),labels=lab,cex=0.6,pos=3,adj=2,xpd=TRUE)
dev.off()


################
## HIGH COUNT ##
################

detach(package:ergm)
detach(package:network)
detach(package:intergraph)
library(igraph)

library(scales)
library(rethinking)

#get total individual degree
TendFullSupout <- degree(snaTenSupFull, mode="out")
AladFullSupout <- degree(snaAlaSupFull, mode="out")
TenFullSup<-data.frame(V(snaTenSupFull)$name,TendFullSupout)
colnames(TenFullSup)<-c("IndivID","FullOutDegree")
AlaFullSup<-data.frame(V(snaAlaSupFull)$name,AladFullSupout)
colnames(AlaFullSup)<-c("IndivID","FullOutDegree")

Ten<-merge(Ten,TenFullSup,by="IndivID",all.x=TRUE,sort=FALSE)
Ala<-merge(Ala,AlaFullSup,by="IndivID",all.x=TRUE,sort=FALSE)

Ten$Age2 <- Ten$Age^2
Ala$Age2 <- Ala$Age^2

#create dataframes for rethinking
participdataTen <- data.frame(Ten$IndivID, as.integer(Ten$HighCountOut), as.integer(Ten$HighCountOutNoAbr), as.integer(Ten$FullOutDegree))
colnames(participdataTen) <- c("indiv","n_highout","n_highoutNA","n_alters")
participdataTen$is_SC <- Ten$Caste=="Pallar" | Ten$Caste == "Arundhathiyar"
participdataTen$is_male <- Ten$Gender=="Male"
participdataTen$eduyears <- (Ten$EduYears-mean(Ten$EduYears))/sd(Ten$EduYears)
participdataTen$age <- (Ten$Age-mean(Ten$Age))/sd(Ten$Age)
participdataTen$age2 <- (Ten$Age2-mean(Ten$Age2))/sd(Ten$Age2)
participdataTen$eduyears <- (Ten$EduYears-mean(Ten$EduYears))/sd(Ten$EduYears)
participdataTen$wealth <- (Ten$HouseholdWealthCalc-mean(Ten$HouseholdWealthCalc))/sd(Ten$HouseholdWealthCalc)
participdataTen$infl <- Ten$InfluenceScaled 

participdataAla <- data.frame(Ala$IndivID, as.integer(Ala$HighCountOut), as.integer(Ala$HighCountOutNoAbr), as.integer(Ala$FullOutDegree))
colnames(participdataAla) <- c("indiv","n_highout","n_highoutNA","n_alters")
participdataAla$is_SC <- Ala$Caste=="Pallar" | Ala$Caste == "Arundhathiyar" | Ala$Caste == "Paraiyar"
participdataAla$is_male<-Ala$Gender=="Male"
participdataAla$eduyears <- (Ala$EduYears-mean(Ala$EduYears))/sd(Ala$EduYears)
participdataAla$age <- (Ala$Age-mean(Ala$Age))/sd(Ala$Age)
participdataAla$age2 <- (Ala$Age2-mean(Ala$Age2))/sd(Ala$Age2)
participdataAla$wealth <- (Ala$HouseholdWealthCalc-mean(Ala$HouseholdWealthCalc))/sd(Ala$HouseholdWealthCalc)
participdataAla$infl <- Ala$InfluenceScaled

#store in list so can run villages with same code
dat <- list(Ten=participdataTen, Ala=participdataAla)
for (dfs in 1:2){
  village <- names(dat)[dfs]
  villagedat <- dat[[dfs]]
  m_infl<-map2stan(
    alist(
      n_highoutNA ~ dbinom(n_alters, p),
      logit(p) <- alpha + tau*z[indiv] + beta_infl*infl + beta_edu*eduyears + beta_age*age + beta_age2*age2 + beta_male*is_male + beta_SC*is_SC+ beta_wealth*wealth,
      z[indiv] ~ dnorm(0, 1),
      alpha ~ dnorm(0, 1),
      tau ~ dcauchy(0, 1),
      beta_infl ~ dnorm(0, 1),
      beta_edu ~ dnorm(0, 1),
      beta_age ~ dnorm(0, 1),
      beta_age2 ~ dnorm(0, 1),
      beta_male ~ dnorm(0, 1),
      beta_SC ~ dnorm(0, 1),
      beta_wealth ~ dnorm(0, 1)
    ),
    data=villagedat, chains=4, iter=2000, warmup=1000, control=list(max_treedepth=20),constraints=list(tau="lower=0")
  )
  pdf(paste(village, "_m_infl_pairs.pdf", sep=""))
  pairs(m_infl, pars=c("alpha", "tau", "beta_infl",
                       "beta_edu", "beta_age", "beta_age2", 
                       "beta_male", "beta_SC", "beta_wealth"))
  dev.off()
  models <- list(m_infl=m_infl)
  for (i in 1:length(models)) {
    model <- models[[i]]
    name <- names(models)[i]
    pdf(paste(village, "_", name, "_trace.pdf", sep=""))
    plot(model, ask=FALSE)
    dev.off()
    modelsum <- precis(model, digits=3, prob=0.95)
    write.csv(modelsum, paste(village, "_", name, ".csv", sep=""))
    line=c("WAIC", WAIC(model))
    write(line,file=paste(village, "_", name, ".csv", sep=""), append=TRUE)
    pdf(paste(village, "_", name, ".pdf", sep=""))
    plot(modelsum)
    mtext(paste(village, name), 3, 1)
    dev.off()
    if (village == "Ten") {m_infl_Ten <- m_infl}
    if (village == "Ala") {m_infl_Ala <- m_infl}
  }  
}

pdf("influence_models.pdf", width=6, height=3, pointsize=9)
b <- precis(m_infl_Ala, prob=0.65)
b2 <- precis(m_infl_Ala, prob=0.95)
points <- cbind(b[,1:4], b2[,3:4])
par(mar=c(5,7,4,1), mfrow=c(1,2), lend=3)
plot(points[,1], c(9:1), xlim=c(-2.5, 1.25), ylim=c(0.5, 9.5), xlab=("Estimate"),ylab="", yaxt="n", main="Alakapuram", cex.axis=0.8, type="n")
abline(v=0, lty=1, lwd=0.75)
for (i in 9:1) {
  abline(h=10-i, lty=3, lwd=0.5)
  lines(c(points[10-i,3], points[10-i,4]), c(i, i), lwd=5)
  lines(c(points[10-i,5], points[10-i,6]), c(i, i), lwd=2)
}
points(points[,1], c(9:1), pch=3, cex=1.5)
axis(2, at = c(9:1), labels = c("Alpha", "Tau", "Influential", "Education", "Age", "Age squared", "Gender (F=0)", "Caste (SC=1)", "Wealth"), cex.axis=0.9, tick = FALSE, line=FALSE, las=TRUE)
b <- precis(m_infl_Ten, prob=0.65)
b2 <- precis(m_infl_Ten, prob=0.95)
points <- cbind(b[,1:4], b2[,3:4])
par(mar=c(5,1,4,7), lend=3)
plot(points[,1], c(9:1), xlim=c(-2.5, 1.25), ylim=c(0.5, 9.5), xlab=("Estimate"),ylab="", yaxt="n", main="Tenpatti", cex.axis=0.8, type="n")
abline(v=0, lty=1, lwd=0.75)
for (i in 9:1) {
  abline(h=10-i, lty=3, lwd=0.5)
  lines(c(points[10-i,3], points[10-i,4]), c(i, i), lwd=5)
  lines(c(points[10-i,5], points[10-i,6]), c(i, i), lwd=2)
}
points(points[,1], c(9:1), pch=3, cex=1.5)
dev.off()



#########################################################################
########### NETWORK MEASURES OF SOCIAL CAPITAL vs. REPUTATION ###########
#########################################################################


TenSupname <- V(snaTenSup)$name

Tsupsum<-E(snaTenSup)$ImpIss+E(snaTenSup)$Work+E(snaTenSup)$Errand+E(snaTenSup)$Borrow+E(snaTenSup)$Cash+E(snaTenSup)$Loan+E(snaTenSup)$Babysit+E(snaTenSup)$Talk+E(snaTenSup)$Defend+E(snaTenSup)$Position+E(snaTenSup)$Close+E(snaTenSup)$Advice

TenprwSup <- page.rank(snaTenSup,directed=TRUE,weights=Tsupsum)

Ten$WeightedPageRank <- TenprwSup$vector



Asupsum<-E(snaAlaSup)$ImpIss+E(snaAlaSup)$Work+E(snaAlaSup)$Errand+E(snaAlaSup)$Borrow+E(snaAlaSup)$Cash+E(snaAlaSup)$Loan+E(snaAlaSup)$Babysit+E(snaAlaSup)$Talk+E(snaAlaSup)$Defend+E(snaAlaSup)$Position+E(snaAlaSup)$Close+E(snaAlaSup)$Advice

AlaprwSup <- page.rank(snaAlaSup,directed=TRUE,weights=Asupsum)

Ala$WeightedPageRank <- AlaprwSup$vector


Ala$WPRscaled<-Ala$WeightedPageRank/max(Ala$WeightedPageRank)
Ten$WPRscaled<-Ten$WeightedPageRank/max(Ten$WeightedPageRank)


palette(c("orchid", "forestgreen"))
palette(adjustcolor(palette(), alpha.f = 0.5))


## Influence, Generous
# WITHOUT GENDER, LOESS
require(msir)
pdf("centrality_results_wprinflgen.pdf", height=3, width=3, pointsize=6)
par(mfrow=c(2,2), mar=c(0.5,4,4,0.5))
plot(Ten$WPRscaled, Ten$InfluenceScaled, pch=19, cex=0.5, main="Tenpatti", axes=F, ylab="GivesGoodAdvice", xlab="")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
axis(2, lwd=0.5)
grid(NULL, col="grey90", lty=3, lwd=0.5)
points(Ten$WPRscaled, Ten$InfluenceScaled, pch=19, cex=0.5)
#lines(loess.smooth(Ten$WPRscaled, Ten$InfluenceScaled, span = 0.8, degree = 1, family = c("gaussian")), lwd=1)
lines(l<-loess.sd(Ten$InfluenceScaled~Ten$WPRscaled))
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("black",alpha.f=0.25),border=NA)

par(mar=c(0.5,0.5,4,4))
plot(Ala$WPRscaled, Ala$InfluenceScaled,  main="Azhagapuram", axes=F, ylab="", xlab="")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
grid(NULL, col="grey90", lty=3, lwd=0.5)
points(Ala$WPRscaled, Ala$InfluenceScaled, pch=19, cex=0.5)
#lines(loess.smooth(Ala$WPRscaled, Ala$InfluenceScaled, span = 0.8, degree = 1, family = c("gaussian")), lwd=1)
lines(l<-loess.sd(Ala$InfluenceScaled~Ala$WPRscaled))
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("black",alpha.f=0.25),border=NA)

par(mar=c(4,4,0.5,0.5)) 
plot(Ten$WPRscaled, Ten$GenerousScaled, axes=F, ylab="GoodGenerous", xlab="Weighted PageRank centrality")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
axis(2, lwd=0.5)
axis(1, lwd=0.5)
grid(NULL, col="grey90", lty=3, lwd=0.5)
points(Ten$WPRscaled, Ten$GenerousScaled, pch=19, cex=0.5)
lines(l<-loess.sd(Ten$GenerousScaled~Ten$WPRscaled))
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("black",alpha.f=0.25),border=NA)

par(mar=c(4,0.5,0.5,4))
plot(Ala$WPRscaled, Ala$GenerousScaled, axes=F, ylab="", xlab="Weighted PageRank centrality")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
grid(NULL, col="grey90", lty=3, lwd=0.5)
axis(1, lwd=0.5)
points(Ala$WPRscaled, Ala$GenerousScaled, pch=19, cex=0.5)
lines(l<-loess.sd(Ala$GenerousScaled~Ala$WPRscaled))
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("black",alpha.f=0.25),border=NA)
dev.off()


# WITH GENDER, LOESS
pdf("centrality_results_wprinflgenGEN.pdf", height=3, width=3, pointsize=6)
par(mfrow=c(2,2), mar=c(0.5,4,4,0.5))
plot(Ten$WPRscaled, Ten$InfluenceScaled, pch=19, cex=0.5, main="Tenpatti", axes=F, ylab="GivesGoodAdvice", xlab="")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
axis(2, lwd=0.5)
grid(NULL, col="grey90", lty=3, lwd=0.5)
points(Ten$WPRscaled, Ten$InfluenceScaled, col=Ten$Gender, pch=19, cex=0.5)
lines(l<-loess.sd(Ten$InfluenceScaled[which(Ten$Gender=="Female")]~Ten$WPRscaled[which(Ten$Gender=="Female")]),col="orchid")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("orchid",alpha.f=0.25),border=NA)
lines(l<-loess.sd(Ten$InfluenceScaled[which(Ten$Gender=="Male")]~Ten$WPRscaled[which(Ten$Gender=="Male")]),col="forestgreen")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("forestgreen",alpha.f=0.25),border=NA)


par(mar=c(0.5,0.5,4,4))
plot(Ala$WPRscaled, Ala$InfluenceScaled,  main="Azhagapuram", axes=F, ylab="", xlab="")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
grid(NULL, col="grey90", lty=3, lwd=0.5)
points(Ala$WPRscaled, Ala$InfluenceScaled, col=Ala$Gender, pch=19, cex=0.5)
lines(l<-loess.sd(Ala$InfluenceScaled[which(Ala$Gender=="Female")]~Ala$WPRscaled[which(Ala$Gender=="Female")]),col="orchid")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("orchid",alpha.f=0.25),border=NA)
lines(l<-loess.sd(Ala$InfluenceScaled[which(Ala$Gender=="Male")]~Ala$WPRscaled[which(Ala$Gender=="Male")]),col="forestgreen")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("forestgreen",alpha.f=0.25),border=NA)

par(mar=c(4,4,0.5,0.5)) 
plot(Ten$WPRscaled, Ten$GenerousScaled, axes=F, ylab="GoodGenerous", xlab="Weighted PageRank centrality")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
axis(2, lwd=0.5)
axis(1, lwd=0.5)
grid(NULL, col="grey90", lty=3, lwd=0.5)
points(Ten$WPRscaled, Ten$GenerousScaled, col=Ten$Gender, pch=19, cex=0.5)
lines(l<-loess.sd(Ten$GenerousScaled[which(Ten$Gender=="Female")]~Ten$WPRscaled[which(Ten$Gender=="Female")]),col="orchid")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("orchid",alpha.f=0.25),border=NA)
lines(l<-loess.sd(Ten$GenerousScaled[which(Ten$Gender=="Male")]~Ten$WPRscaled[which(Ten$Gender=="Male")]),col="forestgreen")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("forestgreen",alpha.f=0.25),border=NA)

par(mar=c(4,0.5,0.5,4))
plot(Ala$WPRscaled, Ala$GenerousScaled, axes=F, ylab="", xlab="Weighted PageRank centrality")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
grid(NULL, col="grey90", lty=3, lwd=0.5)
axis(1, lwd=0.5)
points(Ala$WPRscaled, Ala$GenerousScale, col=Ala$Gender, pch=19, cex=0.5)
lines(l<-loess.sd(Ala$GenerousScaled[which(Ala$Gender=="Female")]~Ala$WPRscaled[which(Ala$Gender=="Female")]),col="orchid")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("orchid",alpha.f=0.25),border=NA)
lines(l<-loess.sd(Ala$GenerousScaled[which(Ala$Gender=="Male")]~Ala$WPRscaled[which(Ala$Gender=="Male")]),col="forestgreen")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("forestgreen",alpha.f=0.25),border=NA)
legend("bottomright", c("Men", "Women"), pch=19, cex=0.75, pt.lwd=0, col=c( "forestgreen","orchid"), bty="n")
dev.off()

## Good Advice, Character
# WITHOUT GENDER, LOESS
require(msir)
pdf("centrality_results_wpradvchar.pdf", height=3, width=3, pointsize=6)
par(mfrow=c(2,2), mar=c(0.5,4,4,0.5))
plot(Ten$WPRscaled, Ten$AllAdviceScaled, pch=19, cex=0.5, main="Tenpatti", axes=F, ylab="GivesGoodAdvice", xlab="")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
axis(2, lwd=0.5)
grid(NULL, col="grey90", lty=3, lwd=0.5)
points(Ten$WPRscaled, Ten$AllAdviceScaled, pch=19, cex=0.5)
#lines(loess.smooth(Ten$WPRscaled, Ten$AllAdviceScaled, span = 0.8, degree = 1, family = c("gaussian")), lwd=1)
lines(l<-loess.sd(Ten$AllAdviceScaled~Ten$WPRscaled))
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("black",alpha.f=0.25),border=NA)

par(mar=c(0.5,0.5,4,4))
plot(Ala$WPRscaled, Ala$AllAdviceScaled,  main="Azhagapuram", axes=F, ylab="", xlab="")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
grid(NULL, col="grey90", lty=3, lwd=0.5)
points(Ala$WPRscaled, Ala$AllAdviceScaled, pch=19, cex=0.5)
#lines(loess.smooth(Ala$WPRscaled, Ala$AllAdviceScaled, span = 0.8, degree = 1, family = c("gaussian")), lwd=1)
lines(l<-loess.sd(Ala$AllAdviceScaled~Ala$WPRscaled))
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("black",alpha.f=0.25),border=NA)

par(mar=c(4,4,0.5,0.5)) 
plot(Ten$WPRscaled, Ten$CharacterScaled, axes=F, ylab="GoodCharacter", xlab="Weighted PageRank centrality")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
axis(2, lwd=0.5)
axis(1, lwd=0.5)
grid(NULL, col="grey90", lty=3, lwd=0.5)
points(Ten$WPRscaled, Ten$CharacterScaled, pch=19, cex=0.5)
lines(l<-loess.sd(Ten$CharacterScaled~Ten$WPRscaled))
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("black",alpha.f=0.25),border=NA)

par(mar=c(4,0.5,0.5,4))
plot(Ala$WPRscaled, Ala$CharacterScaled, axes=F, ylab="", xlab="Weighted PageRank centrality")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
grid(NULL, col="grey90", lty=3, lwd=0.5)
axis(1, lwd=0.5)
points(Ala$WPRscaled, Ala$CharacterScaled, pch=19, cex=0.5)
lines(l<-loess.sd(Ala$CharacterScaled~Ala$WPRscaled))
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("black",alpha.f=0.25),border=NA)
dev.off()


# WITH GENDER, LOESS
pdf("centrality_results_wpradvcharGEN.pdf", height=3, width=3, pointsize=6)
par(mfrow=c(2,2), mar=c(0.5,4,4,0.5))
plot(Ten$WPRscaled, Ten$AllAdviceScaled, pch=19, cex=0.5, main="Tenpatti", axes=F, ylab="GivesGoodAdvice", xlab="")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
axis(2, lwd=0.5)
grid(NULL, col="grey90", lty=3, lwd=0.5)
points(Ten$WPRscaled, Ten$AllAdviceScaled, col=Ten$Gender, pch=19, cex=0.5)
lines(l<-loess.sd(Ten$AllAdviceScaled[which(Ten$Gender=="Female")]~Ten$WPRscaled[which(Ten$Gender=="Female")]),col="orchid")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("orchid",alpha.f=0.25),border=NA)
lines(l<-loess.sd(Ten$AllAdviceScaled[which(Ten$Gender=="Male")]~Ten$WPRscaled[which(Ten$Gender=="Male")]),col="forestgreen")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("forestgreen",alpha.f=0.25),border=NA)


par(mar=c(0.5,0.5,4,4))
plot(Ala$WPRscaled, Ala$AllAdviceScaled,  main="Azhagapuram", axes=F, ylab="", xlab="")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
grid(NULL, col="grey90", lty=3, lwd=0.5)
points(Ala$WPRscaled, Ala$AllAdviceScaled, col=Ala$Gender, pch=19, cex=0.5)
lines(l<-loess.sd(Ala$AllAdviceScaled[which(Ala$Gender=="Female")]~Ala$WPRscaled[which(Ala$Gender=="Female")]),col="orchid")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("orchid",alpha.f=0.25),border=NA)
lines(l<-loess.sd(Ala$AllAdviceScaled[which(Ala$Gender=="Male")]~Ala$WPRscaled[which(Ala$Gender=="Male")]),col="forestgreen")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("forestgreen",alpha.f=0.25),border=NA)

par(mar=c(4,4,0.5,0.5)) 
plot(Ten$WPRscaled, Ten$CharacterScaled, axes=F, ylab="GoodCharacter", xlab="Weighted PageRank centrality")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
axis(2, lwd=0.5)
axis(1, lwd=0.5)
grid(NULL, col="grey90", lty=3, lwd=0.5)
points(Ten$WPRscaled, Ten$CharacterScaled, col=Ten$Gender, pch=19, cex=0.5)
lines(l<-loess.sd(Ten$CharacterScaled[which(Ten$Gender=="Female")]~Ten$WPRscaled[which(Ten$Gender=="Female")]),col="orchid")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("orchid",alpha.f=0.25),border=NA)
lines(l<-loess.sd(Ten$CharacterScaled[which(Ten$Gender=="Male")]~Ten$WPRscaled[which(Ten$Gender=="Male")]),col="forestgreen")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("forestgreen",alpha.f=0.25),border=NA)

par(mar=c(4,0.5,0.5,4))
plot(Ala$WPRscaled, Ala$CharacterScaled, axes=F, ylab="", xlab="Weighted PageRank centrality")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white", lty=0)
grid(NULL, col="grey90", lty=3, lwd=0.5)
axis(1, lwd=0.5)
points(Ala$WPRscaled, Ala$GenerousScale, col=Ala$Gender, pch=19, cex=0.5)
lines(l<-loess.sd(Ala$CharacterScaled[which(Ala$Gender=="Female")]~Ala$WPRscaled[which(Ala$Gender=="Female")]),col="orchid")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("orchid",alpha.f=0.25),border=NA)
lines(l<-loess.sd(Ala$CharacterScaled[which(Ala$Gender=="Male")]~Ala$WPRscaled[which(Ala$Gender=="Male")]),col="forestgreen")
polygon(c(l$x,rev(l$x)),c(l$upper,rev(l$lower)),col=adjustcolor("forestgreen",alpha.f=0.25),border=NA)
legend("bottomright", c("Men", "Women"), pch=19, cex=0.75, pt.lwd=0, col=c( "forestgreen","orchid"), bty="n")
dev.off()


#########
## PCA ##
#########

## Factor analysis of the 4 reputations under study

ten.pca <- prcomp(Ten[,20:23], center = TRUE,scale. = FALSE)
summary(ten.pca)
print(ten.pca)

ala.pca <- prcomp(Ala[,20:23], center = TRUE,scale. = FALSE)
summary(ala.pca)
print(ala.pca)

library(factoextra)

t_pca<-fviz_pca_biplot(ten.pca, repel = TRUE,
                        col.var = "#2E9FDF",
                        col.ind = "#696969",
                        geom="point",
                        title="Tenpatti PCA"
)

t_pca+scale_y_reverse() ## since the signs are arbitrary, reversing this so it aligns with that for the other village

a_pca<-fviz_pca_biplot(ala.pca, repel = TRUE,
                        col.var = "#2E9FDF",
                        col.ind = "#696969",
                        geom="point",
                        title="Alakapuram PCA"
)

a_pca
