#################
# COVID-19 time-series infection date double adjustment
# Citation: Goldstein ND, Burstyn I. Further Improving Analysis of Date-Based COVID-19 Surveillance Data. Manuscript in preparation.
# 12/14/21 -- Neal Goldstein
#################


### FUNCTIONS ###

library("R2jags") #JAGS interface
library("MCMCvis") #visualizations

#function to add CIs to plot
add.ci <- function(x, est, se, lvl=.95, trnc.lo=0) {
  for (i in 1:length(x)) {
    points(rep(x[i],2),
           pmax(est[i]+qnorm(1-(1-lvl)/2)*c(-1,1)*se[i],rep(trnc.lo,2)),
           type="l")
  }
}

#Bayesian model specification for true counts based on reported positives; see: https://pubmed.ncbi.nlm.nih.gov/32505172/
genmod.JAGS = function()
{
  
  ### prior distribution
  
  ### how much weight on the linear component for sens
  ### (earlier version corresponds to no weight)
  
  ### sn.wt <- 0 no smoothness (as per earlier versions)
  ### sn.wt <- 1 perfectly linear
  sn.wt ~ dunif(0.5,0.9)  ### some unspecified amount of smoothing
  
  ### linear component endpoints
  sn.LL ~ dunif(sn.wt*sn.lo[1], sn.wt*sn.hi[1])
  sn.LR ~ dunif(sn.wt*sn.lo[num.kn], sn.wt*sn.hi[num.kn])
    
  ### prev, sn piecewise linear
  ### parameterized by value at knots
  for (i in 1:num.kn) {
    r.kn[i] ~ dunif(0, r.hi[i])
    
    ### linear component
    sn.L[i] <- ((knts[i]-knts[1])*sn.LR+(knts[num.kn]-knts[i])*sn.LL)/
               (knts[num.kn]-knts[1])
    
    ### jumpy component
    sn.J[i] ~ dunif((1-sn.wt)*sn.lo[i],(1-sn.wt)*sn.hi[i]) 
    
    ### two together
    sn.kn[i] ~ dsum(sn.L[i], sn.J[i])
  }
  
  sp ~ dunif(sp.lo,1)
  
  ### these imply the daily values
  for (i in 1:(num.kn-1)) {
    for (j in 0:(spc.kn[i]-1)) {
      r[knts[i]+j] <- ((spc.kn[i]-j)*r.kn[i]+j*r.kn[i+1])/(spc.kn[i])
    }  
  }    
  r[knts[num.kn]] <- r.kn[num.kn]
  
  ### these imply the daily values
  for (i in 1:(num.kn-1)) {
    for (j in 0:(spc.kn[i]-1)) {
      sn[knts[i]+j] <- ((spc.kn[i]-j)*sn.kn[i]+j*sn.kn[i+1])/(spc.kn[i])
    }  
  }    
  sn[knts[num.kn]] <- sn.kn[num.kn]
  
  for (i in 1:(knts[num.kn])) {
    y[i] ~ dbinom(r[i], n[i])            ### true positives
    ystr1[i] ~ dbinom(sn[i], y[i])       ### correct positives
    ystr2[i] ~ dbinom(1-sp, n[i]-y[i])   ### false positives
    ystr[i] ~ sum(ystr1[i], ystr2[i])
  }

}


### READ DATA ###

#Data may be downloaded from https://dph.georgia.gov/covid-19-daily-status-report
#See link "Download the data (CSV)", and then select the epicurve_symptom_date.csv 
georgia_data = read.csv("epicurve_symptom_date.csv", as.is=T, stringsAsFactors=F)

#retain only relevant data and variables
georgia_data = georgia_data[georgia_data$measure=="state_total", c("symptom.date","cases")]

#set date type
georgia_data$symptom.date = as.Date(georgia_data$symptom.date, format="%Y-%m-%d")

#trim to date ranges that are reported in Hennessee et al.: December 1, 2020–April 11, 2021
georgia_data = georgia_data[305:436, ]


### BAYESIAN MODELING for TRUE CASE COUNTS ###

#active data set
dta = georgia_data

#population size of GA from 2020 census: https://www.census.gov/library/stories/state-by-state/georgia-population-change-between-census-decade.html
pop_size = 10711908

#number of time points
T.end = nrow(dta)

#total observed positive results
ystr_onset = dta$cases[1:T.end]

#susceptible population at each time point
n_onset = c(pop_size - (cumsum(dta$cases) - dta$cases[1]))

#observed prevalence per time point
q.hat_onset = ystr_onset/n_onset

#error for observed test positivity
se_onset = sqrt(q.hat_onset*(1-q.hat_onset)/n_onset)

#select interior knots based on linear trends
#knts = c(1, 15, 26, 55, 89, 140, T.end)
knts = c(1, T.end)
num.kn = length(knts)
spc.kn = knts[-1] - knts[-num.kn]

#plot observed data and knots (report)
plot(1:T.end, q.hat_onset, 
     xlab=paste("Dates (",dta$symptom.date[1]," - ",dta$symptom.date[T.end],")",sep=""),
     ylab="Proportion Positive", ylim=c(0, max(q.hat_onset+2.1*se_onset)))
add.ci((1:T.end), q.hat_onset, se_onset)
points(knts, rep(0,num.kn),pch=17,col="red")

#specify hyperparameters
sn.lo.onset = rep(0.6, num.kn)  # same bound at each knot
sn.hi.onset = rep(0.9, num.kn)  # same upper-bound at each knot
r.hi.onset = rep(0.05, num.kn)   # same upper-bound at each knot
sp.lo.onset = 0.95
line_data_onset = list(knts=knts, num.kn=num.kn, spc.kn=spc.kn, sp.lo=sp.lo.onset, sn.lo=sn.lo.onset, sn.hi=sn.hi.onset, r.hi=r.hi.onset, ystr=ystr_onset, n=n_onset)

#initial values (assigned within JAGS)
line_inits = function() {list(y=round(1.2*ystr), ystr1=ystr, ystr2=rep(0,length(ystr)))}
#line_inits_onset = function() {list(y=round(1.2*ystr_onset), ystr1=ystr_onset, ystr2=rep(0,length(ystr_onset)))}

#run JAGS with benchmarking
time1 = Sys.time()
model_onset_posterior = jags.parallel(data=line_data_onset, inits=line_inits, parameters.to.save=c("sp","sn.wt","sn.kn","r.kn","y"), model.file=genmod.JAGS, n.chains=4, n.thin=200, n.iter=400000, n.burnin=1000)
time2 = Sys.time()
time2-time1
rm(time1, time2)

#diagnostics
options(max.print=9999)
print(model_onset_posterior)
MCMCtrace(model_onset_posterior, params=c("sp","sn.wt","sn.kn","r.kn","y"), wd="~/Downloads/", filename="Appendix.pdf")

#extract true counts and add to original (dated) data
y_post_onset = data.frame(model_onset_posterior$BUGSoutput$summary[grep("y",substr(row.names(model_onset_posterior$BUGSoutput$summary),1,1)),], stringsAsFactors=F)

georgia_data$true_onset_median = round(y_post_onset$X50.)
georgia_data$true_onset_lo = round(y_post_onset$X2.5.)
georgia_data$true_onset_hi = round(y_post_onset$X97.5.)

#extract posterior counts from all chains for reported date
count_posterior = function(covid_posterior) {
  #extract each chain as a dataframe
  #nrow = (n.iter - n.burnin) / n.thin
  chain1 = as.data.frame(as.mcmc(covid_posterior)[[1]])
  chain2 = as.data.frame(as.mcmc(covid_posterior)[[2]])
  chain3 = as.data.frame(as.mcmc(covid_posterior)[[3]])
  chain4 = as.data.frame(as.mcmc(covid_posterior)[[4]])

  #add zip code specific posteriors to dataframe
  count_posterior = data.frame(matrix(NA,ncol=nrow(georgia_data),nrow=nrow(chain1)*4))
  col_mapping = colnames(count_posterior)[order(names(count_posterior))]
  for (i in 1:ncol(count_posterior)) {
    
    #determine correct column (offset is based on first y.t column in mcmc output)
    col_index = which(col_mapping==colnames(count_posterior)[i]) + (which(colnames(chain1)=="y[1]") - 1)
    
    #add data
    count_posterior[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index])
  }
  rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4)
  
  return(count_posterior)
}

model_onset_posterior_counts = count_posterior(model_onset_posterior)

#clean up
rm(dta, line_data_onset, y_post_onset, knts, n, n_onset, num.kn, q.hat_onset, r.hi, r.hi.onset, se_onset, sn.hi, sn.hi.onset, sn.lo, sn.lo.onset, sp.lo, sp.lo.onset, spc.kn, T.end, ystr, ystr_onset, add.ci, genmod.JAGS, line_inits, count_posterior)

#save
save.image("bayes_posterior.RData")


### READ DATA ###

load("bayes_posterior.RData")


### FUNCTIONS ###

library("EpiNow2") #deconvolution; see: https://epiforecasts.io/EpiNow2/; tested on v1.1.0


### EPINOW2 APPROACH TO RECOVERING INFECTION DATE ###

#https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.17.2000257
generation_time = list(mean = covid_generation_times[1, ]$mean,
                       mean_sd = covid_generation_times[1, ]$mean_sd,
                       sd = covid_generation_times[1, ]$sd,
                       sd_sd = covid_generation_times[1, ]$sd_sd,
                       max = 30)

#https://annals.org/aim/fullarticle/2762808/incubation-period-coronavirus-disease-2019-covid-19-from-publicly-reported
incubation_period = list(mean = covid_incubation_period[1, ]$mean,
                         mean_sd = covid_incubation_period[1, ]$mean_sd,
                         sd = covid_incubation_period[1, ]$sd,
                         sd_sd = covid_incubation_period[1, ]$sd_sd,
                         max = 30)

#random sample from posteriors
set.seed(777)
rand_samples = sample(1:nrow(model_onset_posterior_counts), 100, replace=F)

#final posterior samples
epinow_posterior = data.frame("sample"=NA, "date"=as.Date("2020-01-01"), "value"=NA, "model"=NA, stringsAsFactors=F)

for (i in 1:length(rand_samples)) {

  cat("\n\n************** ","Observation: ",i," **************\n",sep="")
  
  #data frame for projections
  reported_cases = data.frame("date"=georgia_data$symptom.date, "confirm"=as.integer(model_onset_posterior_counts[rand_samples[i], ]))
  
  #projections
  #time1 = Sys.time()
  #options(mc.cores = ifelse(interactive(), 4, 1))
  model_report_posterior_epinow = epinow(reported_cases = reported_cases, generation_time = generation_time,
                                         delays = list(incubation_period),
                                         rt_prior = list(mean = 1, sd = 1),
                                         samples = 2000, warmup = 200, cores = ifelse(interactive(), 4, 1), chains = 4,
                                         verbose = TRUE, return_fit = TRUE)
  
  #time2 = Sys.time()
  #time2-time1
  #rm(time1, time2)
  
  #extract posteriors
  model_report_posterior_epinow_counts = model_report_posterior_epinow$estimates$samples[model_report_posterior_epinow$estimates$samples$variable=="infections"]
  model_report_posterior_epinow_counts = model_report_posterior_epinow_counts[model_report_posterior_epinow_counts$sample %in% sample(1:2000, length(rand_samples), replace=F), c("sample","date","value")]
  model_report_posterior_epinow_counts$model = rand_samples[i]
  
  #add to dataframe
  epinow_posterior = rbind(epinow_posterior, model_report_posterior_epinow_counts)
  
  #clean up
  rm(model_report_posterior_epinow, reported_cases, model_report_posterior_epinow_counts)
  gc()
  
}
rm(i, incubation_period, generation_time, model_onset_posterior_counts)
epinow_posterior = epinow_posterior[-1, ]

#unique time-series indicator
epinow_posterior$sample_model = paste(epinow_posterior$sample, epinow_posterior$model, sep="-")

#save
save.image("bayes_posterior_EpiNow.RData")


### ADJUSTED CURVE ###

load("bayes_posterior_EpiNow.RData")

#impute infection date via median incubation (5.1 days) subtraction as reported in Hennessee et al.
georgia_data$infection.date = georgia_data$symptom.date - 5.1

#compute 2.5, 50, 97.5 quantiles
plot_med = aggregate(epinow_posterior$value, list(epinow_posterior$date), median)
plot_lo = aggregate(epinow_posterior$value, list(epinow_posterior$date), quantile, probs=0.025)
plot_hi = aggregate(epinow_posterior$value, list(epinow_posterior$date), quantile, probs=0.975)
plot_dates = unique(epinow_posterior$date)

palette = c("#000000", "#000000", "#666666")
plot(plot_dates, plot_med$x, type="l", lwd=3, xaxt="n", xlab="Date", ylab="Cases", ylim=c(0,10000), col=palette[3], lwed=3, lty=2)
points(plot_dates[seq(1, length(plot_dates), 10)], plot_med$x[seq(1, length(plot_med$x), 10)], col=palette[3], pch=19)
axis(side=1, at=seq(min(plot_dates), max(plot_dates), by=7), labels=format(seq(min(plot_dates), max(plot_dates), by=7), "%m-%d"), las=2, cex.axis=1)
polygon(x=c(plot_dates,rev(plot_dates)), y=c(plot_lo$x,rev(plot_hi$x)), col=rgb(t(col2rgb(palette[3])), alpha=100, maxColorValue=255), border=NA)
lines(georgia_data$infection.date, georgia_data$cases, col=palette[2], lwd=2, lty=1)
points(georgia_data$infection.date[seq(1, length(georgia_data$infection.date), 10)], georgia_data$cases[seq(1, length(georgia_data$cases), 10)], col=palette[2], pch=17)
legend("topright", legend=c("Hennessee et al.", "Goldstein et al."), lty=c(1,2), col=palette[c(2,3)], pch=c(17,19), lwd=3, cex=1)

rm(plot_med,plot_lo,plot_hi,plot_dates)

