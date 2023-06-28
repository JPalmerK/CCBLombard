## Re-do Lombard 
rm(list=ls())

library(sp) # for point.in.polygon
library(mgcv)  
library(boot)
library(ggplot2)
library(formula.tools) # for converting formulas
if (dir.exists(
  'C:\\Users\\Kaitlin\\OneDrive - SOI Group Limited\\KJP\\CCB_Comparison\\SECR')){
    setwd('C:\\Users\\Kaitlin\\OneDrive - SOI Group Limited\\KJP\\CCB_Comparison\\SECR')

}

if (dir.exists(
  'C:\\Users\\Kaitlin Palmer\\OneDrive - SOI Group Limited\\KJP\\CCB_Comparison\\SECR')){
  setwd('C:\\Users\\Kaitlin Palmer\\OneDrive - SOI Group Limited\\KJP\\CCB_Comparison\\SECR')
  
}

# Run Mick's predict VCV
source("predictvcv.R")



# MARU locations
MARUloc=read.table("Mass Bay MARU_locs.TXT",sep="",header=T)


# Read in all raven logs
SelTables = list.files(pattern =  "\\FreqLim_RLAdded.txt$")
sel_table_out= data.frame()

# Load selection tables add pertinant dates and aggregate NL values
for (ii in 1:length(SelTables)){
  
  idx_start = nrow(sel_table_out)
  
  sel_table =  read.csv( SelTables[ii], header = T, sep = '\t')
  sel_table$ID = as.numeric(sel_table$ID + nrow(sel_table_out)) 
  
  # Separate the RL and NL
  sel_table_NL = subset(sel_table, RLNLtag != 'RL')
  sel_table_RL = subset(sel_table, RLNLtag == 'RL')
  sel_table$RLNLtag[sel_table$RLNLtag == 'NL2'] = as.factor('NL1')
  
  
  # Get minimum noise level and add it
  mm= aggregate(data = sel_table, Inband.Power..dB.~ ID + Channel, FUN = min)
  colnames(mm)[3] = 'NL'
  sel_table_RL = merge(sel_table_RL, mm, all.x = T)
  
  # Add date and change RLNL tag to RL
  colnames(sel_table_RL)[colnames(sel_table_RL)=="Inband.Power..dB."] = 'RL'
  sel_table_RL$Date = as.factor(substr(SelTables[ii], 12,14))
  sel_table_RL$SNR = sel_table_RL$RL-sel_table_RL$NL
  
  if (ii==1){
    sel_table_out = sel_table_RL
  }else{
  
  # Make the column names match
    for (jj in 1:ncol(sel_table_out)){
      if (colnames(sel_table_out)[jj] %in% colnames(sel_table_RL)){1+1}
      else{
        sel_table_RL[colnames(sel_table_out)[jj]] = NaN}
    }
    
    for (jj in 1:ncol(sel_table_RL)){
      if (colnames(sel_table_RL)[jj] %in% colnames(sel_table_out)){1+1}
      else{
        sel_table_out[colnames(sel_table_RL)[jj]] = NaN}
    }
  

  sel_table_out = rbind(sel_table_out, sel_table_RL)
  }
}

rm(sel_table_RL, sel_table, sel_table_NL, mm)

sel_table_out$date = as.factor(sel_table_out$Date)
sel_table_out$Capture[sel_table_out$Capture>1]=1
sel_table_out$Channel = as.factor(sel_table_out$Channel)
colnames(sel_table_out)[12:13]=c('Lat', 'Lon')
levels(sel_table_out$date) <- c("Feb-18", "Feb-19", "Feb-20",'Feb-21','Feb-22','Apr-17')



################################################################################
# Data Cleaning  #
################################################################################

# Remove arrivals with noise levels at the tails of the distribution
quants = quantile(sel_table_out$NL, c(.025,.975))
sel_table_out = sel_table_out[sel_table_out$NL>quants[1],]
sel_table_out = sel_table_out[sel_table_out$NL<quants[2],]
sel_table_out = sel_table_out[!is.na(sel_table_out$Capture),]



################################################################################
# Restric datapoints to only those in the within the bounds of the array #
################################################################################


ip=which(point.in.polygon(sel_table_out$Lat, sel_table_out$Lon,
                          MARUloc$x[c(1,2,3,4,5,9,10,8,1)],
                          MARUloc$y[c(1,2,3,4,5,9,10,8,1)])==1)
sel_table_out$Included=0
sel_table_out$Included[ip]=1
sel_table_IP=sel_table_out[ip,]
rownames(sel_table_IP)=seq(1, nrow(sel_table_IP))
# 
# ## Check the plot ##
# plot(sel_table_out$Lon, sel_table_out$Lat, pch=19)
# points(MARUloc$x, MARUloc$y, pch=19, col='red')
# points(dat$Lat[131], dat$Lon[131], pch=19, col='pink')
# text(MARUloc$x+900, MARUloc$y+500, labels = as.character(MARUloc$phone), col='green')

#########################################################################
# TL coefficient
########################################################################
library(lme4)

# Remove calls where the start time was estimated
TL_data = sel_table_IP#[sel_table_out$StartTimeEstimated==0,]

# Set the ID as a factor and put distance on the log scale for estimating recieved level
TL_data$ID = as.factor(TL_data$ID)
TL_data$Distance = log10(TL_data$DistanceKm*1000)
TL_data$Distancem = TL_data$DistanceKm*1000
TL_data=TL_data[TL_data$DistanceKm>.03,]
TL_data = TL_data[TL_data$Channel !=3,] 
TL_data = TL_data[TL_data$SNR>0,]
TL_data= TL_data[TL_data$NCaptures>4,]

# Step through the calls and calculate the distance and delta rl values
unique_calls = unique(TL_data$ID)


TL_data$Lp = 10* log10(10^(TL_data$RL/10)-10^(TL_data$NL/10))  
TL_data$Lp[TL_data$SNR>10] =TL_data$RL[TL_data$SNR>10] 


fits <- lmList(Lp ~ Distance | ID, data=TL_data)
slope = coef(fits)[,2]
hist(slope[slope<0])

# Get the median value of the slope, only consider values below -1 and make positive for ease later
TL_exp = abs(median(slope[slope<0], na.rm = TRUE))

#########################################################################
# Estimate Source Levels and when the call was produced
########################################################################
# use TL to estimate SL for all calls and the time at which the call was produced


# For testing channel 5 effect
#sel_table_IP=sel_table_IP[sel_table_IP$Channel !=5,]
sel_table_IP=droplevels(sel_table_IP)
# Add autocorrelation values for each channel and time period


sel_table_IP$SL = NA
unique_calls = unique(sel_table_IP$ID)


for (ii in 1:length(unique_calls)){
  
  # Pull out the call
  call_idx = which(sel_table_IP$ID == unique_calls[ii])
  
  # Pull out the call
  call = sel_table_IP[call_idx,]
  
  # Identify calls where the SNR is greater than 10 or less than 2
  SNR = call$RL-call$NL
  
  
  Lp = 10* log10(10^(call$RL/10)-10^(call$NL/10))  
  Lp[SNR>10]= call$RL
  Lp[SNR<3]= NaN
  
  
  
  # Use tl exponent to estimate RL
  SL = Lp + TL_exp*log10(sel_table_IP$DistanceKm[call_idx]*1000)
  
  SL =10*log10(mean(10^(SL/10), na.rm=TRUE))
  
  sel_table_IP$SL[call_idx] = SL
  sel_table_IP$CallTimeProduced[call_idx] = 
   mean(call$Begin.Time..s.-
    (call$DistanceKm*1000/1432))
  
  
}


#
sel_table_IP = sel_table_IP[!is.nan(sel_table_IP$SL),]
sel_table_IP = sel_table_IP[!is.nan(sel_table_IP$NL),]
sel_table_IP = sel_table_IP[!is.na(sel_table_IP$Capture),]

ggplot(sel_table_IP)+geom_point(aes(NL,SL, 
                                    color = as.character(Capture)))+
theme_bw()

# Add correlation analysis

# #Investigate source level autocorrelation
  dataCh5 = sel_table_IP[sel_table_IP$Channel==5,]
 bb = lm(SL ~ Begin.Time..s.,data =dataCh5 )
 plot(acf(rstandard(bb), lag = 900), main='Source Level Autocorrelation')



#plot(acf(rstandard(sel_table_IP$SL), lag = 400))
timeBins = seq(0,max(sel_table_IP$End.Time..s.)+200, 200)



# Sort by start time so R doesn't pitch a fit
sel_table_IP[order(sel_table_IP$End.Time..s.),]

sel_table_IP$BinnedTime <- as.factor(paste(sel_table_IP$Date, sel_table_IP$Channel, 
                                 cut(sel_table_IP$Begin.Time..s.,
                               breaks = timeBins,
                               labels = timeBins[2:length(timeBins)])))

# ggplot(sel_table_IP[sel_table_IP$Channel==1,])+
#   geom_point(aes(Lon, Lat, color = BinnedTime))+
#   theme(legend.position = "none") 
dataSubset = sel_table_IP[sel_table_IP$Channel==1,]

ggplot(dataSubset)+
  geom_point(aes(Lon, Lat, color = BinnedTime))+
  theme_bw()+
  theme(legend.position = "none") 
  


######################################################################
## NL FIG
######################################################################


ggplot(sel_table_IP, aes(x=Channel, y=NL)) + 
  geom_boxplot(notch = TRUE) +
  facet_wrap(~date, nrow = 2)
  theme(text=element_text(family="serif"))

######################################################################
## Build Likely GEEGLM models and select from canidates using QIC   ##
######################################################################

require(sp)
require(rgeos)
library(ggplot2)
library(geepack)
library(splines)
library(MuMIn)
library(MRSea)






MOD1.geeglm <- geeglm(Capture ~ poly(DistanceKm, 2)+
                        bs(NL, knots = mean(NL))+
                        Channel+date, corstr = 'exchangeable',
                      family = binomial, # leave out constrains
                      scale.fix=T, id=BinnedTime, data = sel_table_IP)

MOD2.geeglm <- geeglm(Capture ~poly(DistanceKm, 2)+
                        bs(NL, knots = mean(NL))+
                        date, corstr = 'exchangeable',
                      family = binomial, # leave out constrains
                      scale.fix=T, id=BinnedTime, data = sel_table_IP )

MOD3.geeglm <- geeglm(Capture ~poly(DistanceKm, 2)+
                        Channel+date, corstr = 'exchangeable',
                      family = binomial, # leave out constrains
                      scale.fix=T, id=BinnedTime, data = sel_table_IP )


MOD4.geeglm <- geeglm(Capture ~poly(DistanceKm, 2)+
                        bs(NL, knots = mean(NL))+
                        Channel+date, 
                      family = binomial, # leave out constrains
                      scale.fix=T, id=BinnedTime, data = sel_table_IP)


MOD5.geeglm <- geeglm(Capture ~poly(DistanceKm, 2)+
                        bs(NL, knots = mean(NL))+
                        Channel, corstr = 'exchangeable',
                      family = binomial, # leave out constrains
                      scale.fix=T, id=BinnedTime, data = sel_table_IP)


MOD6.geeglm <- geeglm(Capture ~poly(DistanceKm, 2)+
                        bs(NL, knots = mean(NL)), 
                      family = binomial,corstr = 'exchangeable', # leave out constrains
                      scale.fix=T, id=BinnedTime, data = sel_table_IP)

MOD7.geeglm <- geeglm(Capture ~poly(DistanceKm, 2)+
                        Channel, 
                      family = binomial, corstr = 'exchangeable',# leave out constrains
                      scale.fix=T, id=BinnedTime, data = sel_table_IP)

MOD8.geeglm <- geeglm(Capture ~ bs(NL, knots = mean(NL))+
                        Channel, 
                      family = binomial, corstr = 'exchangeable',# leave out constrains
                      scale.fix=T, id=BinnedTime, data = sel_table_IP)


MOD9.geeglm <- geeglm(Capture ~poly(DistanceKm, 2), 
                      family = binomial, corstr = 'exchangeable',# leave out constrains
                      scale.fix=T, id=BinnedTime, data = sel_table_IP)

MOD10.geeglm <- geeglm(Capture ~ bs(NL, knots = mean(NL)), 
                       family = binomial,corstr = 'exchangeable', # leave out constrains
                       scale.fix=T, id=BinnedTime, data = sel_table_IP)

# knock out the contstraint
MOD11.geeglm <- geeglm(Capture ~poly(DistanceKm, 2)+
                         bs(NL, knots = mean(NL))+
                         Channel, 
                       family = binomial, # leave out constrains
                       scale.fix=T, id=BinnedTime, data = sel_table_IP )


MOD12.geeglm<- geeglm(Capture ~ poly(DistanceKm, 2)+
                        NL+
                        Channel+date, corstr = 'exchangeable',
                      family = binomial, # leave out constrains
                      scale.fix=T, id=BinnedTime, data = sel_table_IP)

MOD13.geeglm<- geeglm(Capture ~ poly(DistanceKm, 2)+
                        NL+
                        Channel+date, 
                      family = binomial, # leave out constrains
                      scale.fix=T, id=BinnedTime, data = sel_table_IP)



#acf(residuals(MOD1.geeglm, type='pearson')) # to tell if there is an auto
#getPvalues(MOD1.geeglm)

Qicvals=data.frame(Mod=as.character(seq(1,13)), QIC=rep(0,13), Constr=rep(0,13),
                   deltaqic=rep(0,13))
Qicvals$QIC[1]=QIC(MOD1.geeglm)
Qicvals$QIC[2]=QIC(MOD2.geeglm)
Qicvals$QIC[3]=QIC(MOD3.geeglm)
Qicvals$QIC[4]=QIC(MOD4.geeglm)
Qicvals$QIC[5]=QIC(MOD5.geeglm)
Qicvals$QIC[6]=QIC(MOD6.geeglm)
Qicvals$QIC[7]=QIC(MOD7.geeglm)
Qicvals$QIC[8]=QIC(MOD8.geeglm)
Qicvals$QIC[9]=QIC(MOD9.geeglm)
Qicvals$QIC[10]=QIC(MOD10.geeglm)
Qicvals$QIC[11]=QIC(MOD11.geeglm)
Qicvals$QIC[12]=QIC(MOD12.geeglm)
Qicvals$QIC[13]=QIC(MOD13.geeglm)
#Qicvals$QIC[14]=QIC(MOD14.geeglm)


Qicvals$Constr=unlist(c(MOD1.geeglm$modelInfo$corstr, MOD2.geeglm$modelInfo$corstr, MOD3.geeglm$modelInfo$corstr,
                        MOD4.geeglm$modelInfo$corstr,MOD5.geeglm$modelInfo$corstr, MOD6.geeglm$modelInfo$corstr,
                        MOD7.geeglm$modelInfo$corstr,MOD8.geeglm$modelInfo$corstr,MOD9.geeglm$modelInfo$corstr,
                        MOD10.geeglm$modelInfo$corstr,MOD11.geeglm$modelInfo$corstr,MOD12.geeglm$modelInfo$corstr,
                        MOD13.geeglm$modelInfo$corstr))#,MOD14.geeglm$modelInfo$corstr))

Qicvals$Formula=(c(as.character(MOD1.geeglm$formula), as.character(MOD2.geeglm$formula), as.character(MOD3.geeglm$formula), 
                         as.character(MOD4.geeglm$formula), as.character(MOD5.geeglm$formula), as.character(MOD6.geeglm$formula),
                         as.character(MOD7.geeglm$formula), as.character(MOD8.geeglm$formula), as.character(MOD9.geeglm$formula), 
                         as.character(MOD10.geeglm$formula), as.character(MOD11.geeglm$formula),  as.character(MOD12.geeglm$formula),
                         as.character(MOD13.geeglm$formula)))#, as.character(MOD14.geeglm$formula)))

Qicvals$deltaqic=Qicvals$QIC-min(Qicvals$QIC)

# Linear model
MOD1Linear = update(MOD1.geeglm, ~poly(DistanceKm, 2) + NL+ Channel + date )


# Model summary
library(jtools)
summ(MOD1.geeglm)

#############################################################
# Plot model outcomes
#############################################################



runPartialPlots(MOD1Linear, includeB0=TRUE,savedata = TRUE,
                              varlist = c("NL", "DistanceKm"), 
                             factorlist = c("Channel", 'date'),
                              data=sel_table_IP )


runPartialPlots(MOD1.geeglm, savedata = TRUE,includeB0=TRUE,
                varlist = c("NL", "DistanceKm"), 
                factorlist = c("Channel", 'date'),
                data=sel_table_IP )


# Load the data
 load('PartialData_DistanceKm.RData')
ggplot(partialdata, aes(x = newX))+
  geom_line(aes(y= partialfit))+
  geom_ribbon(aes(ymin=cis...1.,ymax=cis...2.),
              alpha=0.2)+
  xlab('Range (km)')+
  ylab('Probability of Detection')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  # geom_rug(data=subset(sel_table_IP, Capture==0), 
  #          aes(x=DistanceKm, y=Capture), 
  #          sides = 'b', alpha=.4) +
  # geom_rug(data=subset(sel_table_IP, Capture==1), 
  #          aes(x=DistanceKm, y=Capture),
  #          sides = 't', alpha=.4) +
  xlim(c(0,35))
  

# Load the data
load('PartialData_NL.RData')
ggplot(partialdata, aes(x = newX))+
  geom_line(aes(y= partialfit))+
  geom_ribbon(aes(ymin=cis...1.,ymax=cis...2.),
              alpha=0.2)+
  xlab('Noise Level')+
  ylab('Probability of Detection')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  # geom_rug(data=subset(sel_table_IP, Capture==0), 
  #          aes(x=NL, y=Capture), 
  #         sides = 'b', alpha=.4) +
  # geom_rug(data=subset(sel_table_IP, Capture==1), 
  #          aes(x=NL, y=Capture),
  #          sides = 't', alpha=.4) 
  xlim(c(100,130))


ggplot(sel_table_IP[sel_table_IP$Channel==1,])+
  geom_point(aes(Lon, Lat, color = BinnedTime))+
  theme(legend.position = "none")

################################################################
# Model performance and Effective Detection Area/Radius
###################################################################
library(boot) # just for inv.logit

# Model Prediction matrix #

nlVals = round(as.numeric(quantile(sel_table_IP$NL, seq(.05,.95,.05))),1)
distRes = .1
maxRange = 25
rangeVals = seq(0,maxRange, by=distRes)

# Prediction data
newdat <- expand.grid( 
  DistanceKm = rangeVals,
  NL = nlVals,
  Channel=as.character(seq(1,10)),
  date=as.character(221),
  Capture=0)


levels(newdat$date)=as.factor(unique(sel_table_IP$date))
levels(newdat$Channel)=levels(sel_table_IP$Channel)

# Store linear data
newdatLinear<-newdat

# Model selected using QIC
CIvals=predictvcv(MOD1.geeglm, newdata = newdat) 
newdat$pdet=inv.logit(CIvals$fit)
newdat$lcl=inv.logit(CIvals$lwr)
newdat$ucl=inv.logit(CIvals$upr)
newdat$pdf=newdat$pdet*(pi*newdat$DistanceKm^2-pi*(newdat$DistanceKm-distRes)^2)/
  (pi*maxRange^2)
newdat$pdf.low=newdat$lcl*(pi*newdat$DistanceKm^2-pi*(newdat$DistanceKm-distRes)^2)/
  (pi*maxRange^2)
newdat$pdf.hi=newdat$ucl*(pi*newdat$DistanceKm^2-pi*(newdat$DistanceKm-distRes)^2)/
  (pi*maxRange^2)
newdat$NL=as.factor(newdat$NL)

# Calculate effective detection radius - areas under curve/area monitored
newdat$AnulusArea = 2*pi*(newdat$DistanceKm^2- 
                           (newdat$DistanceKm-distRes)^2)
newdat$grR = newdat$AnulusArea*newdat$pdet
  
effArea = aggregate(data=newdat, grR~NL+Channel+date, FUN =sum)
effArea$Radius = sqrt(effArea$grR/pi)
effArea$NL =  as.numeric(as.character(effArea$NL))
ggplot(effArea)+
  geom_point(aes(NL,Radius, color=Channel))+
  geom_line(aes(NL,Radius, color=Channel))+
  scale_color_brewer(type = "div", 
                     name= expression('Ambient Sound Level (dB[rms] 20-200Hz)'))+
  theme(legend.position = "none") +ylim(c(20,37))+
  ggtitle('smooth')
  


# Linear Model

CIvals=predictvcv(MOD1Linear, newdata = newdatLinear) # Model selected using QIC
newdatLinear$pdet=inv.logit(CIvals$fit)
newdatLinear$lcl=inv.logit(CIvals$lwr)
newdatLinear$ucl=inv.logit(CIvals$upr)
newdatLinear$pdf=newdatLinear$pdet*pi*(newdatLinear$DistanceKm^2-(newdatLinear$DistanceKm-distRes)^2)/
  (pi*max(maxRange)^2)
newdatLinear$pdf.low=newdatLinear$lcl*(pi*newdatLinear$DistanceKm^2-pi*(newdatLinear$DistanceKm-distRes)^2)/
  (pi*max(maxRange)^2)
newdatLinear$pdf.hi=newdatLinear$ucl*(pi*newdatLinear$DistanceKm^2-pi*(newdatLinear$DistanceKm-distRes)^2)/
  (pi*max(maxRange)^2)
newdatLinear$NL=as.factor(newdatLinear$NL)


# Calculate effective detection radius - areas under curve/area monitored
newdatLinear$AnulusArea = 2*pi*(newdatLinear$DistanceKm^2- 
                            (newdatLinear$DistanceKm-distRes)^2)
newdatLinear$grR = newdat$AnulusArea*newdatLinear$pdet
effAreaLinear = aggregate(data=newdatLinear, grR~NL+Channel+date, FUN =sum)
effAreaLinear$Radius = sqrt(effAreaLinear$grR/pi)
effAreaLinear$NL =  as.numeric(as.character(effAreaLinear$NL))
ggplot(effAreaLinear)+
  geom_point(aes(NL,Radius, color=Channel))+
  geom_line(aes(NL,Radius, color=Channel))+
  scale_color_brewer(type = "div",
                     name= 'Linear RMS Noise(dB 20-200Hz)')+
  theme(legend.position = "none") +ylim(c(20,37))+
  ggtitle('Linear')


ggplot(effAreaLinear)+
  geom_point(aes(NL,diffR, color=Channel))+
  geom_line(aes(NL,diffR, color=Channel))+
  scale_color_brewer(type = "div",
                     name= 'Linear RMS Noise(dB 20-200Hz)')+
  theme(legend.position = "none") 



##########################################################
# Make Prediction Plots
##############################################################


# Model Prediction matrix #
newdat <- expand.grid(
  DistanceKm = seq(min(sel_table_IP$DistanceKm), max(sel_table_IP$DistanceKm), by=.2),
  NL = round(as.numeric(quantile(sel_table_IP$NL, seq(0.25,1, .25))),1),
  Channel=as.character(seq(1,10)),
  date=as.character(221),
  Capture=0)
levels(newdat$date)=as.factor(unique(sel_table_IP$date))
levels(newdat$Channel)=levels(sel_table_IP$Channel)
CIvals=predictvcv(MOD1.geeglm, newdata = newdat) # Model selected using QIC


newdat$pdet=inv.logit(CIvals$fit)
newdat$lcl=inv.logit(CIvals$lwr)
newdat$ucl=inv.logit(CIvals$upr)

# Probability Density Function 
newdat$pdf=newdat$pdet*(pi*newdat$DistanceKm^2-pi*(newdat$DistanceKm-.1)^2)/(pi*max(sel_table_out$DistanceKm)^2)
newdat$pdf.low=newdat$lcl*(pi*newdat$DistanceKm^2-pi*(newdat$DistanceKm-.1)^2)/(pi*max(sel_table_out$DistanceKm)^2)
newdat$pdf.hi=newdat$ucl*(pi*newdat$DistanceKm^2-pi*(newdat$DistanceKm-.1)^2)/(pi*max(sel_table_out$DistanceKm)^2)
newdat$NL=as.factor(newdat$NL)
newdat$Model = 'Smooth'


# Prediction model matrix for the linear version
newdatLM <- expand.grid(
  DistanceKm = seq(min(sel_table_IP$DistanceKm), max(sel_table_IP$DistanceKm), by=.2),
  NL = round(as.numeric(quantile(sel_table_IP$NL, seq(0.25,1, .25))),1),
  Channel=as.character(seq(1,10)),
  date=as.character(221),
  Capture=0)
levels(newdatLM$date)=as.factor(unique(sel_table_IP$date))
levels(newdatLM$Channel)=levels(sel_table_IP$Channel)
CIvals=predictvcv(MOD1Linear, newdata = newdatLM) # Model selected using QIC
newdatLM$pdet=inv.logit(CIvals$fit)
newdatLM$lcl=inv.logit(CIvals$lwr)
newdatLM$ucl=inv.logit(CIvals$upr)
newdatLM$pdf=newdatLM$pdet*(pi*newdatLM$DistanceKm^2-pi*(newdatLM$DistanceKm-.1)^2)/(pi*max(sel_table_out$DistanceKm)^2)
newdatLM$pdf.low=newdatLM$lcl*(pi*newdatLM$DistanceKm^2-pi*(newdatLM$DistanceKm-.1)^2)/(pi*max(sel_table_out$DistanceKm)^2)
newdatLM$pdf.hi=newdatLM$ucl*(pi*newdatLM$DistanceKm^2-pi*(newdatLM$DistanceKm-.1)^2)/(pi*max(sel_table_out$DistanceKm)^2)
newdatLM$NL=as.factor(newdatLM$NL)
newdatLM$Model = 'Linear'



newdat=rbind(newdat,newdatLM)

# Lists to save figures #
p=list()
p1=list()
AreaAll=list()

# Make Pdet and PDF plots for all channels on on one day, only use 2 for publication 
# because 10 is simply too many 
legendTitle = expression(paste("Ambient Sound Level\n(dB[rms] re: 1\u03BC Pa 20-200Hz)", 
                           sep = "")) 
for (ii in 1:10){
  
  datsub=subset(newdat, Channel==as.character(ii))
  ccbsub=subset(sel_table_out, Channel==ii, select=c('Capture', 'DistanceKm'))
  
  datsub$AnulusArea = pi*(datsub$DistanceKm^2-(datsub$DistanceKm-.2)^2)
  datsub$aa = datsub$pdet* datsub$AnulusArea
  Area =aggregate(data=datsub, aa~NL+Model, FUN = sum)
  Area$Radius = sqrt(Area$aa/pi)
  Area$Channel = ii 
  AreaAll[[ii]] =Area
  
  p[[ii]]=ggplot()+
    #geom_rug(data=subset(ccbsub, Capture==0), aes(x=DistanceKm, y=Capture), 
    #        sides = 'b', alpha=.4) +
    #geom_rug(data=subset(ccbsub, Capture==1), aes(x=DistanceKm, y=Capture),
    #         sides = 't', alpha=.4) +
    
    geom_line(data=datsub, size=1, aes(DistanceKm, 
                                       pdet, colour=NL)) +
    geom_ribbon(data=datsub, aes(x=DistanceKm, ymin=lcl, ymax=ucl, color=NL), 
                alpha=.2,linetype='blank') +
    facet_grid(~Model)+
    scale_color_brewer(type = "div", name= legendTitle) +         
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()  +        
    labs(x='Range (km)', y='Detection Probability') +
    theme(axis.title=element_text(size=12, face="bold")) +
    ggtitle(paste('MARU', as.character(ii))) +
    scale_fill_manual("",values="grey12") +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  
  p1[[ii]]=ggplot(datsub, aes(DistanceKm, pdf, colour=NL)) + 
    scale_colour_brewer(type = "div", 
                        name= 'RMs Noise Level\n (dB 20-200Hz)') +
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw() +
    geom_line(size=1) +
    geom_ribbon(aes(ymin=pdf.low, ymax=pdf.hi), 
                alpha=.2,linetype='blank') +
    labs(x='Range (km)', y='Probability Density') +
    theme(axis.title=element_text(size=12,face="bold"))+
    scale_fill_manual("",values="grey12") +
    ggtitle(paste('MARU', as.character(ii))) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  #print(p)
  #print(p1)
  
}




######################################################
# Noise Level Plots
######################################################
legendTitle = expression(paste("Noise Level (dB re: 1\u03BC Pa)", 
                               sep = "")) 

ggplot(sel_table_IP, aes(x=Channel, y=NL)) + 
   geom_boxplot(notch = TRUE) +
  facet_wrap(~ date, nrow = 2)+
  theme_bw()+
  theme(text=element_text(family="serif"))+
  ylab(legendTitle)+
  xlab('MARU')



#############################################################
# Evidence of Lombard Long
#############################################################

unique_calls = unique(sel_table_IP$ID)

LomData = data.frame(SL = numeric(length(unique_calls)))
LomData$NL = 0
LomData$Time = 0
LomData$date = sel_table_IP$date[1]
LomData$DistanceKm = 0
LomData$ID = 0


# For each call pull out the NL associated with the shortest range
for(ii in 1:length(unique_calls)){
  call= sel_table_IP[sel_table_IP$ID==unique_calls[ii],]
  
  LomData$SL[ii] = call$SL[1]
  LomData$date[ii]= call$date[1]
  LomData$Time[ii] = call$CallTimeProduced[1]
  LomData$ID = call$ID[1]
  
  nlIdx = which.min(call$DistanceKm)
  LomData$NL[ii] = call$NL[nlIdx]
  LomData$DistanceKm[ii] = call$DistanceKm[nlIdx]
  LomData$RL[ii] = call$RL[nlIdx]
  
  
  
}

# Restric this analysis to less than 5km
LomData= LomData[LomData$DistanceKm<5,]

# Trim the upper 2.5%
LomData=LomData[LomData$SL<quantile(LomData$SL,.975),]
LomData=LomData[LomData$SL>quantile(LomData$SL,.025),]
LomData[with(LomData, order(Date, Time)), ]      

LomData$BinnedTime <- 
  as.factor(paste(LomData$Date, cut(LomData$Time,
                                               breaks = timeBins,
                                               labels = timeBins[2:length(timeBins)])))



LomMod1 <- geeglm(SL ~ poly(NL,2), corstr = 'exchangeable',
                      family = gaussian, 
                      scale.fix=T, id=BinnedTime, data = LomData)

LomMod2 <- geeglm(SL ~ NL, corstr = 'exchangeable',
                  family = gaussian, 
                  scale.fix=T, id=BinnedTime, data = LomData)

LomMod3 <- geeglm(SL ~ NL + date, corstr = 'exchangeable',
                  family = gaussian,
                  scale.fix=T, id=BinnedTime, data = LomData)

LomMod4<- geeglm(SL ~ NL, 
                 family = gaussian, 
                 scale.fix=T, id=BinnedTime, data = LomData)




QIC(LomMod1,LomMod2,LomMod3, LomMod4)
QicvalsLom=data.frame(Mod=as.character(seq(1,4)), QIC=rep(0,1), 
                      Constr=rep(0,4),
                   deltaqic=rep(0,4))
QicvalsLom$QIC[1]=QIC(LomMod1)
QicvalsLom$QIC[2]=QIC(LomMod2)
QicvalsLom$QIC[3]=QIC(LomMod3)
QicvalsLom$QIC[4]=QIC(LomMod4)


QicvalsLom$Constr=unlist(c(LomMod1$modelInfo$corstr, 
                        LomMod2$modelInfo$corstr, 
                        LomMod3$modelInfo$corstr,
                        LomMod4$modelInfo$corstr))

QicvalsLom$Formula=c(as.character(LomMod1$formula),
                   as.character(LomMod2$formula),
                   as.character(LomMod3$formula),
                   as.character(LomMod4$formula))
#as.character(MOD13.geeglm$formula))), as.character(MOD14.geeglm$formula)))

QicvalsLom$deltaqic=QicvalsLom$QIC-min(QicvalsLom$QIC)


####################################################################
# Lombard Plotting and prediction
######################################################################

CIvals=predictvcv(LomMod4, newdata = LomData) # Model selected using QIC
CIvals2=predictvcv(LomMod2, newdata = LomData) # Model selected using QIC




ggplot(data =LomData)+
  geom_point(aes(NL, SL, color = date))+
  ylim(c(145,180))+
  theme_bw()+
  scale_color_viridis_d()+
  theme(legend.title = element_blank())+
  xlab(expression(paste("Noise Level (dB re: 1\u03BC Pa)", 
                      sep = "")))+
  ylab(expression(paste("Source Level (dB re: 1\u03BC Pa)", 
                        sep = "")))+
  geom_line(aes(x = NL, y=CIvals$fit), size=1) +
  geom_ribbon( aes(x = NL, y=SL,ymin=CIvals$lwr, ymax=CIvals$upr), 
              alpha=.2,linetype='blank')


ggplot(data =LomData)+
  geom_point(aes(DistanceKm, SL-NL))

####################################################################
# SLNR PLOT`
######################################################################

sel_table_IP$SLNR = sel_table_IP$SL- sel_table_IP$NL

DistBreak = seq(0, 30, by=2) 

SLNLbreaks = round(seq(min(sel_table_IP$SLNR), max(sel_table_IP$SLNR), 
                     by = 2))

sel_table_IP$SLNLbin = cut(sel_table_IP$SLNR, 
                           breaks=SLNLbreaks, 
                           include.lowest=TRUE, 
                           right=FALSE)
sel_table_IP$DistBin = cut(sel_table_IP$DistanceKm, 
                           breaks=DistBreak, 
                           include.lowest=TRUE, 
                           right=FALSE)

SLNRden = aggregate(data = sel_table_IP, 
                    Capture~SLNLbin+DistBin, FUN =mean)
SLNRden$Nobs = aggregate(data = sel_table_IP, 
                         Capture~SLNLbin+DistBin, FUN =length)[3]

SLNRden$Dist =DistBreak[as.numeric(SLNRden$DistBin)]
SLNRden$SLNLBinNumeric =SLNLbreaks[as.numeric(SLNRden$SLNLbin)]

SLNRden =SLNRden[SLNRden$Nobs>5,]
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))


df = data.frame(x = seq(0,max(DistBreak)), 
                y = TL_exp*log10(seq(0,max(DistBreak))*1000))
ggplot() +
  geom_tile(data = SLNRden, mapping = aes(x = Dist,
                                          y = SLNLBinNumeric,
                                          fill = Capture))+
  scale_fill_gradientn(name  ="g(SLNR,r)", colours =  jet.colors(7))+
  theme_bw()+xlab('Range (km)') + ylab('Signal Excess (SLNL)')+
  geom_line(data = df, aes(x = x,  y = y), size=2)



# Aggregate over SNLN 
ggplot(SLNRden)+
  geom_point(aes(x=Dist, y=Capture, color= SLNLBinNumeric))+
  scale_color_viridis_c()


ggplot(SLNRden)+geom_point(aes(x = SLNLbin, y = Capture))
SLNRden = aggregate(data = sel_table_IP, 
                    Capture~SLNLbin+DistBin, FUN =mean)

######################################################
# Correlation in noise levels across the
######################################################
tempdf=subset(sel_table_out, 
              Channel==unique(sel_table_out$Channel)[1],
              select=c( ID, NL))
colnames(tempdf)[2]=unique(sel_table_out$Channel)[1]

# Reshape the dataframe... undoubtedy a better way to do this
for (ii in 2:length(unique(sel_table_out$Channel))){
  x1=subset(sel_table_out, 
            Channel==ii,
            select=c(ID, NL))
  colnames(x1)[2]= as.character(ii)
  tempdf=merge(tempdf, x1, by='ID', all = TRUE)
  rm(list=c('x1'))  
}


Pairwise_Correlation=cor(tempdf[,2:11],
                         use="pairwise.complete.obs")
library(corrgram)


corrgram(Pairwise_Correlation,
         main="Noise Level pairwise  correlations",
         lower.panel=panel.pts, upper.panel=panel.conf,
         diag.panel=panel.density)






###########################################################
# Model cross validation
##############################################################

library(caret)
# Define training control
train.control <- trainControl(method = "cv", number = 10)

# Train the model


cost <- function(r, pi = 0) mean(abs(r-pi) > 0.5)


mycost <- function(r, pi){
  weight1 = 1 #cost for getting 1 wrong
  weight0 = 1 #cost for getting 0 wrong
  c1 = (r==1)&(pi<0.5) #logical vector - true if actual 1 but predict 0
  c0 = (r==0)&(pi>=0.5) #logical vector - true if actual 0 but predict 1
  return(mean(weight1*c1+weight0*c0))
}

cvErr = cv.gamMRSea(data=sel_table_IP, cost = mycost, 
                    modelobject = MOD1.geeglm, K=5)$delta
cvErrLinear = cv.gamMRSea(data=sel_table_IP,cost = mycost,
                          modelobject = MOD1Linear, K=5)






