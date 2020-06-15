require(deSolve)
require(ggplot2)
require(benford.analysis)
require(dplyr)
require(reshape2)
require(gridExtra)
library(grid)
library(circular)
library(BenfordTests)
library(kuiper.2samp)
library(parallel)
library(foreach)
library(doSNOW)

options(digits = 9)




################################################################################################################
simulate_benford <- function(howoften) {
  SIR_fn <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      N  <- S + I + R
      dS <- -beta * S * I / N
      dI <- beta * S * I / N - gamma * I - mu * I
      dR <- gamma * I
      dD <- mu * I
      return(list(c(dS, dI, dR, dD)))
    })
  }
  

  
# start with 100000 susceptibles and 1 infected. 
  initial_state_values <- c(S = 100000,
                            I = 1,
                            R = 0,
                            D = 0)
  
  # choose values to start your optimisation
  # Data-based analysis, modelling and forecasting of the COVID-19 outbreak
  # https://doi.org/10.1371/journal.pone.0230405
  # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0230405
  beta_start  <- 0.191
  gamma_start <- 0.064  # 13 to 20
  mu_start    <- 0.01
  # set the time long so we catch the disease outbreak for certain combinations as quite long time frame 
  times <- seq(from = 0, to = 1000, by = 1)
  #how often is the simulation to be run
  howoften=howoften
  # start or parallel frame work
  raus <- foreach(zerun=1:howoften,.combine=rbind,.packages=c('circular', 'BenfordTests','kuiper.2samp','deSolve','benford.analysis')) %dopar% { 
 
    # choose a random beta/gamma/mu with the mean around the expected beta of Covid19 and a not too extreme sd to get valid disease curves
    rbeta <- rnorm(1, mean = beta_start, sd = 0.05)
    #catch unlikely negative beta/gamma/mu events
    if (rbeta < 0) {rbeta = beta_start}
    rgamma <- rnorm(1, mean = gamma_start, sd = 0.005)
    if (rgamma < 0) {rgamma = gamma_start}
    rmu <- rnorm(1, mean = mu_start, sd = 0.005)
    if (rmu < 0) {rmu = mu_start}
    #set the parameters to the randomised disease parameters
    parameters <- c(beta = rbeta,gamma = rgamma,mu = rmu)
    # run the ode solver
    output <- as.data.frame(ode(y = initial_state_values,times = times,func = SIR_fn,parms = parameters))
    # get the dimensions of the output vector (Infected) 
    n = dim(output["I"])[1]
    #declare an empty vector to catch the increasing curve for the Benford consition > 10% rise and magnitude change > 3
    tobenbool = vector()
    # check the values that are bigger than 10% than the one before
    for (i in 1:(n - 1)) {
      dazu <- (output["I"][[i + 1, 1]] - output["I"][[i, 1]]) > 0.1
      tobenbool <- append(tobenbool, dazu)
    }
    # add an extra FALSE value to make the correct vector length 
    tobenbool <- append(tobenbool, FALSE, after = 0)
    # subset the Infected output vector 
    toben <- subset(output["I"], tobenbool)
    #get the length of the vector
    bn = length(toben[[1]])
    #if the length of the infected incidence is bigger that 49 proceed
    if (bn > 49) {
      #check if there is a magnitude change > 3 in the rising infected cases
      ismag = round(log10(toben[[bn, 1]]), 0) - round(log10(toben[[2, 1]]), 0) >= 3
      # if there is proceed with Benford tests
      if (ismag & bn > 1) {
        #check Benford adherence with benford.analysis package
        isitben <- benford(toben[[1]],number.of.digits = 1,sign = "positive",discrete = FALSE,round = 3)
        #prepare a data frame with number of incidences and the time step
        dfforlm <-  data.frame(n = seq(1, bn), seq = toben[[1]])
        #make a random Benford sequence for BenfordTests package tests
        benny <- rbenf(n = length(toben[[1]]))
        #get the significant digits from the randomly created series and the SIR created series
        X <- signifd(x = toben[[1]], digits = 1)
        Y <- signifd(x = benny, digits = 1)
        #execute a Kuipers 2 sample test on both as per Okamura paper
        kuip <- kuiper.2samp(X, Y)$p.value
        #execute a chisquared test from the BenfordTests package
        csop <- chisq.benftest(toben[[1]])[[2]]
        #Euclidean Distance for Benford Distribution
        ed <-edist.benftest(x = toben[[1]],digits = 1,pvalmethod = "simulate",pvalsims = 10000)[[2]]
        #A Hotelling T-square Type Test for Benford's Law
        htstt <-jointdigit.benftest(x = toben[[1]],digits = 1,eigenvalues = "all",tol = 1e-15,pvalmethod = "asymptotic",pvalsims = 10000)[[2]][[1]]
        #Joenssen's JP-square Test for Benford's Law
        jjpt <-jpsq.benftest(x = toben[[1]],digits = 1,pvalmethod = "simulate",pvalsims = 10000)[[2]]
        #Kolmogorov-Smirnov Test for Benford's Law
        kst <-ks.benftest(x = toben[[1]],digits = 1,pvalmethod = "simulate",pvalsims = 10000)[[2]]
        #Chebyshev Distance Test (maximum norm) for Benford's Law
        cdtmn <-mdist.benftest(x = toben[[1]],digits = 1,pvalmethod = "simulate",pvalsims = 10000)[[2]]
        #Judge-Schechter Mean Deviation Test for Benford's Law
        jsmdt <-meandigit.benftest(x = toben[[1]],digits = 1,pvalmethod = "asymptotic",pvalsims = 10000)[[2]]
        #Freedman-Watson U-square Test for Benford's Law
        fwust <-usq.benftest(x = toben[[1]],digits = 1,pvalmethod = "simulate",pvalsims = 10000)[[2]]
        #take all Benford's law test values and generate a mean 
        meanassess <-mean(c(chisq(isitben)[[3]][[1]],kuip,csop,ed,htstt,jjpt,kst,cdtmn,jsmdt,fwust))
        #take all Benford's law test values and generate a median
        medianassess <-median(c(chisq(isitben)[[3]][[1]],kuip,csop,ed,htstt,jjpt,kst,cdtmn,jsmdt,fwust))
        #take all Benford's law test values and generate a Standard Deviation 
        sdassess <-sd(c(chisq(isitben)[[3]][[1]],kuip,csop,ed,htstt,jjpt,kst,cdtmn,jsmdt,fwust))
        # throw that out for the parallel algo to collect
        c(zerun,isitben[[7]],chisq(isitben)[[3]][[1]],parameters["beta"][[1]],parameters["gamma"][[1]],parameters["mu"][[1]],bn,lm(seq ~ n, data = dfforlm)[[1]][[2]],kuip,csop,ed,htstt,jjpt,kst,cdtmn,jsmdt,fwust,meanassess,medianassess,sdassess)
      } else {
        #print("no benford")
      }
     }
  }
  #return all simulation runs 
  return(raus)
}

#how often the calculations have to be done
wieoft = 100000

# get maximun cores - one for OS
cl<-makeCluster(detectCores() - 1) 
#register the Cluster for Windows
registerDoSNOW(cl)
# get start time
old <- Sys.time() 
#run the 
wasisses <- simulate_benford(wieoft)
# calculate difference of times after the runs finished
new <- Sys.time() - old 
# print in nice format
print(new) 
stopCluster(cl)  
#Time difference of 5.15469722 mins 1000 runs one core
#Time difference of 33.5846059 secs 1000 runs 15 cores
#[1] 9.20903565 overall speed up through paralisation

#create data frame
dfb <- data.frame(run = as.integer(integer()),Nigrini = character(),Xsquared = as.numeric(double()),beta = as.numeric(double()),gamma = as.numeric(double()),mu = as.numeric(double()),laenge = as.integer(integer()),slope = as.numeric(double()),kuiper = as.numeric(double()),csop = as.numeric(double()),ed = as.numeric(double()),htstt = as.numeric(double()),jjpt = as.numeric(double()),kst = as.numeric(double()),cdtmn = as.numeric(double()),jsmdt = as.numeric(double()),fwust = as.numeric(double()),meanassess = as.numeric(double()),medianassess = as.numeric(double()),sdassess = as.numeric(double()), stringsAsFactors = F)
#fill it up from parallel runs simulations
for (i in 1:dim(wasisses)[[1]]) {
  dfb[i,] = c(wasisses[i,])
}
#through out NA cases 
dfb <- na.omit(dfb)

#delete the old file
unlink("E:/datafakts/covid/yohan/montecarlo-parallel.csv")
#write csv
write.csv(dfb, "E:/datafakts/covid/yohan/montecarlo-parallel.csv", row.names = FALSE)

#some formatting in case data frame made a mess of the types

dfb["slope"][[1]] = as.numeric(as.character(dfb["slope"][[1]]))
dfb["laenge"][[1]] = as.numeric(as.character(dfb["laenge"][[1]]))
dfb["beta"][[1]] = as.numeric(as.character(dfb["beta"][[1]]))
dfb["gamma"][[1]] = as.numeric(as.character(dfb["gamma"][[1]]))
dfb["mu"][[1]] = as.numeric(as.character(dfb["mu"][[1]]))
dfb["Xsquared"][[1]] = as.numeric(as.character(dfb["Xsquared"][[1]]))
dfb["kuiper"][[1]] = as.numeric(as.character(dfb["kuiper"][[1]]))
dfb["csop"][[1]] = as.numeric(as.character(dfb["csop"][[1]]))
dfb["ed"][[1]] = as.numeric(as.character(dfb["ed"][[1]]))
dfb["htstt"][[1]] = as.numeric(as.character(dfb["htstt"][[1]]))
dfb["jjpt"][[1]] = as.numeric(as.character(dfb["jjpt"][[1]]))
dfb["kst"][[1]] = as.numeric(as.character(dfb["kst"][[1]]))
dfb["cdtmn"][[1]] = as.numeric(as.character(dfb["cdtmn"][[1]]))
dfb["jsmdt"][[1]] = as.numeric(as.character(dfb["jsmdt"][[1]]))
dfb["fwust"][[1]] = as.numeric(as.character(dfb["fwust"][[1]]))
dfb["meanassess"][[1]] = as.numeric(as.character(dfb["meanassess"][[1]]))
dfb["medianassess"][[1]] = as.numeric(as.character(dfb["medianassess"][[1]]))
dfb["sdassess"][[1]] = as.numeric(as.character(dfb["sdassess"][[1]]))

# get the longest and shortest of the runs for info
min(dfb["laenge"][[1]])
max(dfb["laenge"][[1]])

# code for text analysis if wanted
dfb %>%
  group_by(Nigrini) %>%
  summarise(slope = mean(slope))

# make the graphs
slopeplot <- ggplot(dfb,   aes(x=Nigrini, y=slope))+
  geom_violin(aes(fill=Nigrini))  +
  theme(legend.position="bottom",axis.text.x = element_text(angle = -6))+
  #geom_violin(trim=TRUE)+
  stat_summary(fun=median, geom="point", size=2, color="black")+
  geom_boxplot(width=0.1)+ 
  scale_fill_brewer(palette="Blues")

laengeplot <- ggplot(dfb,   aes(x=Nigrini, y=laenge))+
  geom_violin(aes(fill=Nigrini))  +
  theme(legend.position="bottom",axis.text.x = element_text(angle = -6))+
  #geom_violin(trim=TRUE)+
  stat_summary(fun=median, geom="point", size=2, color="black")+
  geom_boxplot(width=0.1)+ 
  scale_fill_brewer(palette="Blues")


betaplot <- ggplot(dfb,   aes(x=Nigrini, y=beta))+
  geom_violin(aes(fill=Nigrini))  +
  theme(legend.position="bottom",axis.text.x = element_text(angle = -6))+
  #geom_violin(trim=TRUE)+
  stat_summary(fun=median, geom="point", size=2, color="black")+
  geom_boxplot(width=0.1)+ 
  scale_fill_brewer(palette="Blues")

gammaplot <- ggplot(dfb,   aes(x=Nigrini, y=gamma))+
  geom_violin(aes(fill=Nigrini))  +
  theme(legend.position="bottom",axis.text.x = element_text(angle = -6))+
  #geom_violin(trim=TRUE)+
  stat_summary(fun=median, geom="point", size=2, color="black")+
  geom_boxplot(width=0.1)+ 
  scale_fill_brewer(palette="Blues")

muplot <- ggplot(dfb,   aes(x=Nigrini, y=mu))+
  geom_violin(aes(fill=Nigrini))  +
  theme(legend.position="bottom",axis.text.x = element_text(angle = -6))+
  #geom_violin(trim=TRUE)+
  stat_summary(fun=median, geom="point", size=2, color="black")+
  geom_boxplot(width=0.1)+ 
  scale_fill_brewer(palette="Blues")

meanplot <- ggplot(dfb,   aes(x=Nigrini, y=meanassess))+
  geom_violin(aes(fill=Nigrini))  +
  theme(legend.position="bottom",axis.text.x = element_text(angle = -6))+
  #geom_violin(trim=TRUE)+
  stat_summary(fun=median, geom="point", size=2, color="black")+
  geom_boxplot(width=0.1)+ 
  scale_fill_brewer(palette="Blues")

medianplot <- ggplot(dfb,   aes(x=Nigrini, y=medianassess))+
  geom_violin(aes(fill=Nigrini))  +
  theme(legend.position="bottom",axis.text.x = element_text(angle = -6))+
  #geom_violin(trim=TRUE)+
  stat_summary(fun=median, geom="point", size=2, color="black")+
  geom_boxplot(width=0.1)+ 
  scale_fill_brewer(palette="Blues")

sdplot <- ggplot(dfb,   aes(x=Nigrini, y=sdassess))+
  geom_violin(aes(fill=Nigrini))  +
  theme(legend.position="bottom",axis.text.x = element_text(angle = -6))+
  #geom_violin(trim=TRUE)+
  stat_summary(fun=median, geom="point", size=2, color="black")+
  geom_boxplot(width=0.1)+ 
  scale_fill_brewer(palette="Blues")

chi2plot <- ggplot(dfb,   aes(x=Nigrini, y=csop))+
  geom_violin(aes(fill=Nigrini))  +
  theme(legend.position="bottom",axis.text.x = element_text(angle = -6))+
  #geom_violin(trim=TRUE)+
  stat_summary(fun=median, geom="point", size=2, color="black")+
  geom_boxplot(width=0.1)+ 
  scale_fill_brewer(palette="Blues")


kuipplot <- ggplot(dfb,   aes(x=Nigrini, y=kuiper))+
  geom_violin(aes(fill=Nigrini))  +
  theme(legend.position="bottom",axis.text.x = element_text(angle = -6))+
  #geom_violin(trim=TRUE)+
  stat_summary(fun=median, geom="point", size=2, color="black")+
  geom_boxplot(width=0.1)+ 
  scale_fill_brewer(palette="Blues")


Xsquaredplot <- ggplot(dfb,   aes(x=Nigrini, y=Xsquared))+
  geom_violin(aes(fill=Nigrini))  +
  theme(legend.position="bottom",axis.text.x = element_text(angle = -6))+
  #geom_violin(trim=TRUE)+
  stat_summary(fun=median, geom="point", size=2, color="black")+
  geom_boxplot(width=0.1)+ 
  scale_fill_brewer(palette="Blues")

edplot <- ggplot(dfb,   aes(x=Nigrini, y=ed))+
  geom_violin(aes(fill=Nigrini))  +
  theme(legend.position="bottom",axis.text.x = element_text(angle = -6))+
  #geom_violin(trim=TRUE)+
  stat_summary(fun=median, geom="point", size=2, color="black")+
  geom_boxplot(width=0.1)+ 
  scale_fill_brewer(palette="Blues")

# put it all in one graph
grid.arrange(
  slopeplot,
  laengeplot,
  betaplot,
  gammaplot,
  muplot,
  Xsquaredplot,
  chi2plot,
  meanplot,
  medianplot,
  sdplot,
  kuipplot,
  edplot,
  ncol = 3,
  top = "Benford's Law Covid-19 Monte Carlo SIRD Simulations n=100000"
  ,
  bottom = textGrob("",
                    gp = gpar(
                      fontface = 3, fontsize = 9
                    ),
                    just = "bottom"))

