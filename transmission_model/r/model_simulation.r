library(rstan)
library(matrixStats)
library(data.table)
library(lubridate)
library(gdata)
library(dplyr)
library(tidyr)
library(EnvStats)
library(scales)
library(tidyverse)
library(dplyr)
library(abind)
library(xtable)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(bayesplot)
library(cowplot)
library(openxlsx)
library(loo)
library(posterior)
library(here)
library(rstan)
library(matrixStats)
library(data.table)
library(lubridate)
library(gdata)
library(dplyr)
library(tidyr)
library(EnvStats)
library(scales)
library(tidyverse)
library(dplyr)
library(abind)
library(xtable)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(bayesplot)
library(cowplot)
library(openxlsx)
library(loo)
library(cmdstanr)
library(ks)
library(ggpubr)
library(here)
library(optparse)
library(ggExtra)



source(here("transmission_model/utils/geom-stepribbon.r"))
source(here("transmission_model/utils/gammaAlt.r"))

# check if sivep data exists if not download and read else just read
RELEASE_DATE_i = "08-02-2021"
SIVEP_ORIGINAL20 <- here(paste0("transmission_model/data/INFLUD-",RELEASE_DATE_i,".csv"))
SIVEP_ORIGINAL21 <- here(paste0("transmission_model/data/INFLUD21-",RELEASE_DATE_i,".csv"))

if(!file.exists(SIVEP_ORIGINAL20)){
  download.file(paste0("https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2020/INFLUD-",RELEASE_DATE_i,".csv"),SIVEP_ORIGINAL20)
}
if(!file.exists(SIVEP_ORIGINAL21)){
  download.file(paste0("https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2021/INFLUD21-",RELEASE_DATE_i,".csv"),SIVEP_ORIGINAL21)
}

df_SIVEP_ORIGINAL20 <- suppressMessages(read_csv2(SIVEP_ORIGINAL20))
df_SIVEP_ORIGINAL21 <- suppressMessages(read_csv2(SIVEP_ORIGINAL21))

df_SIVEP_original <- rbind(df_SIVEP_ORIGINAL20,df_SIVEP_ORIGINAL21)
deaths_nowcast_original = read.csv(here("transmission_model/data/nowcasted_daily_Manaus_class4and5_2021-02-11.csv"),
                                   stringsAsFactors = FALSE)
df_pop_original = read.csv(here("transmission_model/data/brazil_population.csv"),sep=";", stringsAsFactors = FALSE)
serial.interval = read.csv(here("transmission_model/data/serial_interval.csv"))
onset_paras_original = read.csv(here("transmission_model/data/statelevel_OnsetDeathParams.csv"))
pcr_and_sero <-read.csv(here("transmission_model/data/pcr_and_sero_pos.csv"))
pcr_genome_fraction <- read.csv(here("transmission_model/data/pcr_genome_fraction.csv"))

back_date = "2021-02-07" 
p1back_date = "2021-02-07" 
end_date = "2021-02-07" 
selected_state <- "MANAUS"
deaths_source <- "Nowcast4and5"
job = "base"
name="non_P1_P1_model"

rayleigh_par = 310
rel_IFR1 = 1/100
T2_date = ymd("2020-11-09")
ITER = 400
WARM = 10
CORES = 1
THIN = 4
SEED = 4444
DELTA = 0.95
TREE = 12
JOBID=1

pcr_genome_fraction$date <- dmy(pcr_genome_fraction$date)
pcr_genome_fraction$negative <- pcr_genome_fraction$number - pcr_genome_fraction$positive
pcr_genome_fraction <- pcr_genome_fraction[pcr_genome_fraction$date > T2_date,]
pcr_genome_fraction <- pcr_genome_fraction[pcr_genome_fraction$date <= p1back_date,]

onset_paras = onset_paras_original[c("state","mean","cv")]  %>% drop_na()

deaths_nowcast <- deaths_nowcast_original
deaths_nowcast$"DateRep" = ymd(as.Date(deaths_nowcast$"Date"))
deaths_nowcast = deaths_nowcast[c("DateRep","Deaths","Var")]
deaths_nowcast$region = as.factor(selected_state)
deaths_nowcast_back = deaths_nowcast[deaths_nowcast$Date <= ymd(back_date),]
deaths_nowcast = deaths_nowcast[deaths_nowcast$Date <= ymd(end_date),]

if (deaths_source == "Nowcast4and5")
{
  df_SIVEP_original -> df_SIVEP
  df_SIVEP = df_SIVEP[df_SIVEP$ID_MUNICIP == selected_state,]
  df_SIVEP = df_SIVEP[which(df_SIVEP$CLASSI_FIN==5) || which(df_SIVEP$CLASSI_FIN==4),]
  df_SIVEP = df_SIVEP[which(df_SIVEP$EVOLUCAO==2),]
  df_SIVEP$DT_EVOLUCA = dmy(df_SIVEP$DT_EVOLUCA)
  filter_date = head(sort(df_SIVEP[df_SIVEP$CLASSI_FIN==5,]$DT_EVOLUCA),1)
  df_SIVEP = df_SIVEP %>% filter(DT_EVOLUCA >= ymd(filter_date))
  df_SIVEP = df_SIVEP[,c("DT_EVOLUCA", "ID_MUNICIP")]
  dim(df_SIVEP)
  df_SIVEP$Deaths = 1
  df_SIVEP = aggregate(. ~ID_MUNICIP+DT_EVOLUCA, data=df_SIVEP, sum, na.rm=TRUE)
  colnames(df_SIVEP) = c("region","DateRep","Deaths")
  df_SIVEP = rbind(df_SIVEP[df_SIVEP$DateRep < ymd(head(deaths_nowcast[c("DateRep","Deaths","region")],1)$DateRep),], deaths_nowcast[c("DateRep","Deaths","region")])
  regions <- unique(df_SIVEP$region)
  aux_df_zeros_brazil = NULL
  for( jj in 1:length(regions)){
    aux_sub = subset(df_SIVEP,region==regions[jj])
    aux_sub = aux_sub[order(aux_sub$DateRep),]
    last_date = aux_sub$DateRep[1]
    zeroes<-data.frame(DateRep=seq(as.Date('2020/01/01'),as.Date(last_date-1),"days"),
                       region=regions[jj],
                       Deaths=0)
    zeroes$DateRep <- ymd(zeroes$DateRep)
    aux_df_zeros_brazil = rbind.data.frame(zeroes,aux_df_zeros_brazil)
  }
  df_SIVEP_deaths=rbind.data.frame(aux_df_zeros_brazil,df_SIVEP)
} else
{
  df_SIVEP_original -> df_SIVEP
  df_SIVEP = df_SIVEP[df_SIVEP$ID_MUNICIP == selected_state,]
  df_SIVEP = df_SIVEP[which(df_SIVEP$CLASSI_FIN==5) || which(df_SIVEP$CLASSI_FIN==4),]
  df_SIVEP = df_SIVEP[which(df_SIVEP$EVOLUCAO==2),]
  df_SIVEP$DT_EVOLUCA = dmy(df_SIVEP$DT_EVOLUCA)
  df_SIVEP = df_SIVEP[,c("DT_EVOLUCA", "ID_MUNICIP")] #Date, and Federative unit patient residence.
  df_SIVEP$Deaths = 1
  df_SIVEP = aggregate(. ~ID_MUNICIP+DT_EVOLUCA, data=df_SIVEP, sum, na.rm=TRUE)
  colnames(df_SIVEP) = c("region","DateRep","Deaths")
  regions <- unique(df_SIVEP$region)
  aux_df_zeros_brazil = NULL
  for( jj in 1:length(regions)){
    aux_sub = subset(df_SIVEP,region==regions[jj])
    aux_sub = aux_sub[order(aux_sub$DateRep),]
    last_date = aux_sub$DateRep[1]
    zeroes<-data.frame(DateRep=seq(as.Date('2020/01/01'),as.Date(last_date-1),"days"),
                       region=regions[jj],
                       Deaths=0)
    zeroes$DateRep <- ymd(zeroes$DateRep)
    aux_df_zeros_brazil = rbind.data.frame(zeroes,aux_df_zeros_brazil)
  }
  df_SIVEP_deaths=rbind.data.frame(aux_df_zeros_brazil,df_SIVEP)
}
df_SIVEP_original -> df_SIVEP
df_SIVEP = df_SIVEP[df_SIVEP$ID_MUNICIP == selected_state,]
df_SIVEP = df_SIVEP[which(df_SIVEP$CLASSI_FIN==5),] 
df_SIVEP = df_SIVEP[,c("DT_SIN_PRI", "ID_MUNICIP")]
df_SIVEP$DT_SIN_PRI = dmy(df_SIVEP$DT_SIN_PRI)
df_SIVEP = df_SIVEP[order(df_SIVEP$DT_SIN_PRI),] %>% filter(DT_SIN_PRI >= ymd('2020-02-26')) #painel index case
df_SIVEP$Cases = 1
df_SIVEP = aggregate(. ~ID_MUNICIP+DT_SIN_PRI, data=df_SIVEP, sum, na.rm=TRUE)
colnames(df_SIVEP) = c("region","DateRep","Cases")
regions <- unique(df_SIVEP$region)
aux_df_zeros_brazil = NULL
for( jj in 1:length(regions)){
  aux_sub = subset(df_SIVEP,region==regions[jj])
  aux_sub = aux_sub[order(aux_sub$DateRep),]
  last_date = aux_sub$DateRep[1]
  zeroes<-data.frame(DateRep=seq(as.Date('2020/01/01'),as.Date(last_date-1),"days"),
                     region=regions[jj],
                     Cases=0)
  zeroes$DateRep <- ymd(zeroes$DateRep)
  aux_df_zeros_brazil = rbind.data.frame(zeroes,aux_df_zeros_brazil)
}
df_SIVEP_cases=rbind.data.frame(aux_df_zeros_brazil,df_SIVEP)
regions <- unique(df_SIVEP_deaths$region)
aux_df_zeros_brazil = NULL
for( jj in 1:length(regions)){
  zeroes<-data.frame(DateRep=seq(as.Date('2020/01/01'),as.Date(tail(sort(unique(df_SIVEP_deaths$DateRep)),1)),"days"),
                     region=regions[jj],
                     AUX=0)
  zeroes$DateRep <- ymd(zeroes$DateRep)
  aux_df_zeros_brazil = rbind.data.frame(zeroes,aux_df_zeros_brazil)
}
aux_df_zeros_brazil$key = paste0(aux_df_zeros_brazil$DateRep,aux_df_zeros_brazil$region)
df_SIVEP_deaths$key = paste0(df_SIVEP_deaths$DateRep,df_SIVEP_deaths$region)
df_SIVEP_cases$key = paste0(df_SIVEP_cases$DateRep,df_SIVEP_cases$region)
df_SIVEP_deaths = merge(aux_df_zeros_brazil, df_SIVEP_deaths, by = "key", all=TRUE)
df_SIVEP_cases = merge(aux_df_zeros_brazil, df_SIVEP_cases, by = "key", all=TRUE)
df_SIVEP = merge(df_SIVEP_deaths,df_SIVEP_cases, by = "key")
df_SIVEP = df_SIVEP[,c("DateRep.y.x","region.y.x","Deaths","Cases")]
colnames(df_SIVEP) = c("DateRep","region","Deaths","Cases")
df_SIVEP[c("Cases")] = df_SIVEP[c("Cases")] %>% replace(is.na(.), 0)
df_SIVEP = df_SIVEP %>% group_by(region) %>% mutate(cumulativeCases = cumsum(Cases))
df_SIVEP = df_SIVEP %>% group_by(region) %>% mutate(cumulativeDeaths = cumsum(Deaths))

df_pop <- df_pop_original[c("region","population")]
df_pop <- df_pop[1:(nrow(df_pop)-1),]
df_SIVEP <- merge(x = df_SIVEP, y = df_pop, by = "region", all = TRUE)
df_SIVEP <- df_SIVEP[order(df_SIVEP$DateRep),]

df=df_SIVEP %>% filter(DateRep <= ymd(back_date))
START_TIME=ymd(as.Date("2020-03-03"))
END_TIME=ymd(as.Date(end_date))
RANGE_TIME=seq(START_TIME,END_TIME,by = '1 day')
StanModel = as.character(job)
models = job
writeDir = job
d<-df
countries = selected_state
N2 = length(RANGE_TIME)
dates = list()
reported_cases = list()
stan_data = list(M=length(countries),
                 N=NULL,
                 deaths=NULL,
                 deaths_combined=NULL,
                 f=NULL,
                 N0=6,
                 cases=NULL,
                 SI=NULL,
                 EpidemicStart = NULL,
                 pop = NULL,
                 par = NULL,
                 T2 = NULL, 
                 phylo_N = NULL, 
                 phylo_PSamples = NULL, 
                 phylo_NSamples = NULL,
                 sero_N = NULL, 
                 sero_sigma = NULL, 
                 sero_prev = NULL,
                 len_excess = NULL, 
                 excess_N = NULL, 
                 excess = NULL,
                 nowcast_sd_N = NULL, 
                 nowcast_sd = NULL, 
                 nowcast_sd_len = NULL
) 
deaths_by_country = list()
deaths_by_country_combined = list()
mean1 = 5.1; cv1 = 0.86;
x1 = rgammaAlt(1e6,mean1,cv1)
aux.epidemicStart = NULL
Country = selected_state
mean2 = onset_paras[onset_paras$state == Country,]$mean
cv2 = onset_paras[onset_paras$state == Country,]$cv
x2 = rgammaAlt(1e6,mean2,cv2)
ecdf.saved = ecdf(x1+x2)
d1=d[d$region==Country,]
d1$DateRep = as.Date(d1$DateRep)
d1_pop = df_pop[df_pop$region==Country,]
d1 = d1[order(d1$DateRep),]
index = which(d1$Cases>0)[1]
index1 = which(cumsum(d1$Deaths)>=5)[1] # 5, 10
index2 = index1-26 #26
d1=d1[index2:nrow(d1),]

aux.epidemicStart = c(aux.epidemicStart,d1$DateRep[index1+1-index2])
stan_data$EpidemicStart = c(stan_data$EpidemicStart,index1+1-index2)
stan_data$EpidemicStart = as.array(stan_data$EpidemicStart)
stan_data$pop = c(stan_data$pop, d1_pop$population)
stan_data$pop = as.array(stan_data$pop)
stan_data$par = rayleigh_par
dates[[Country]] = d1$DateRep
N = length(d1$Cases)
forecast = N2 - N
if(forecast < 0) {
  N2 = N
  forecast = N2 - N
}
convolution1 = function(u) (rel_IFR1 * ecdf.saved(u))
f = rep(0,N2) 
f[1] = (convolution1(1.5) - convolution1(0))
for(i in 2:N2)
{
  f[i] = (convolution1(i+.5) - convolution1(i-.5))
}
cases = as.vector(as.numeric(d1$Cases))
deaths = as.vector(as.numeric(d1$Deaths))
stan_data$N = c(stan_data$N,N)
stan_data$f = cbind(stan_data$f,f)
stan_data$deaths = cbind(stan_data$deaths,deaths)
stan_data$cases = cbind(stan_data$cases,cases)
stan_data$N2 = N2
stan_data$x=1:N2
if(length(stan_data$N) == 1) {
  stan_data$N = as.array(stan_data$N)
}
stan_data$W <- ceiling(stan_data$N2/7)
stan_data$week_index <- matrix(1,stan_data$M,stan_data$N2)
for(state.i in 1:stan_data$M) {
  stan_data$week_index[state.i,] <- rep(2:(stan_data$W+1),each=7)[1:stan_data$N2]
  last_ar_week = which(dates[[state.i]]==max(df$Date) -  7)
  stan_data$week_index[state.i,last_ar_week:ncol(stan_data$week_index)] <-
    stan_data$week_index[state.i,last_ar_week]
}
stan_data$AR_SD_MEAN = 0.2
stan_data$P <- 0
stan_data$SI = serial.interval$fit[1:N2]
stan_data$T2 = seq(1, length(dates[[Country]]))[dates[[Country]] %in%  T2_date]

if(length(stan_data$deaths < N2))
{
  stan_data$deaths = rbind(stan_data$deaths,as.matrix(rep(-1,N2-N)))
}

stan_data$CasesStart <- 20 #40
stan_data$compute_likelihood <- 0

length_pcr_sero_data <- length(pcr_and_sero$p_PCR_positive)
padding <- stan_data$N2 - length_pcr_sero_data
stan_data$PCR_pos_prob <- as.matrix(c(pcr_and_sero$p_PCR_positive, rep(0, padding)))
stan_data$seroconv_cdf <- as.matrix(c(pcr_and_sero$cum_seropositive, rep(1, padding)))
stan_data$serorev_surv <- as.matrix(1 - pweibull(seq(1, N2), shape = 2.933, scale = 208))

stan_data$phylo_N = seq(1, length(dates[[Country]]))[dates[[Country]] %in%  pcr_genome_fraction$date]
stan_data$phylo_N_len = length(stan_data$phylo_N)
stan_data$phylo_PSamples = pcr_genome_fraction$positive
stan_data$phylo_NSamples = pcr_genome_fraction$negative

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
######simulation
######you need to run the model and get the mean posterior parameters
a<-rep(0,342)
b<-rep(0,342)
weekly_sd =0.63;
weekly_rho =0.75;
weekly_rho1 = 0.06;
weekly_effect <- matrix(ncol = 1,nrow=50)
#simulation_first
weekly_effect[,1]<-c(-0.000045,-0.110813,-0.244019,-0.235055,-0.115224,0.272611,0.925742,1.710220,1.860155,1.521232,1.133688,0.886009,0.780019,0.774908,0.779310,0.791570,0.722043,0.570189,0.403496,0.295753,0.276571,0.294532,0.258221,0.185230,0.099074,0.104312,0.117610,0.065354,-0.032622,-0.117900,-0.182262,-0.156510,-0.081098,-0.030074,-0.010139,-0.055340,-0.189456,-0.438358,-0.631959,-0.743087,-0.716995,-0.501209,-0.130450,0.276777,0.632271,0.844858,0.681098,0.474862,0.383688,0.327117)
#simulation_second
#weekly_effect[,1]<-rep(0,50) 
prediction_v1 <- matrix(0,ncol = 1,nrow=342)
prediction_v2 <- matrix(0,ncol = 1,nrow=342)
Rt_adj_immune_v1<-matrix(0,ncol = 1,nrow=342)
Rt_adj_immune_v2<-matrix(0,ncol = 1,nrow=342)
alpha_sus1<-matrix(0,ncol = 1,nrow=342)
immune_v1<-matrix(0,ncol = 1,nrow=342)
immune_v2<-matrix(0,ncol = 1,nrow=342)
Rt_v1 <- matrix(0,ncol = 1,nrow=342)
Rt_v2 <- matrix(0,ncol = 1,nrow=342)
n1<- matrix(0,ncol = 1,nrow=342)
n2<- matrix(0,ncol = 1,nrow=342)
E_deaths_v1<- matrix(0,ncol = 1,nrow=342)
E_deaths_v2<- matrix(0,ncol = 1,nrow=342)
E_deaths<- matrix(0,ncol = 1,nrow=342)
alpha_sus2<-matrix(0,ncol = 1,nrow=342)
#weekly_effect[1,1 ] = rnorm(1,0, 0.01);
#weekly_effect[2,1 ] = rnorm(1,0,weekly_sd *sqrt(1-((weekly_rho)^2)-((weekly_rho1)^2) - 2 *((weekly_rho)^2) * weekly_rho1/(1-weekly_rho1)));
#for (w in 3:50)
#{
#  weekly_effect[w, 1] = rnorm(1, weekly_effect[w-1,1]* weekly_rho + weekly_effect[w-2,1]* weekly_rho1,weekly_sd *sqrt(1-((weekly_rho)^2)-((weekly_rho1)^2) - 2 *((weekly_rho)^2) * weekly_rho1/(1-weekly_rho1)));
#}
y_v1<-148.14
y_v2<-2.38
N0<-stan_data$N0
N2<-stan_data$N2
T2<-stan_data$T2
SI<-stan_data$SI
par<-stan_data$par
pop<-stan_data$pop
f<-stan_data$f
week_index<-stan_data$week_index
R_difference<-2.07
cross<-0.65
ifr1<-0.26
ifr2<-0.42
m<-1
prediction_v1[1:N0,m] = rep(y_v1,N0); 
prediction_v2[T2:(T2+N0-1),m] = rep(y_v2,N0); 
Rt_v1[1:342,m] = 3.28 * 2 * exp(- weekly_effect[week_index[1:342],m])/(1+exp(- weekly_effect[week_index[1:342],m]));
#Rt_v1[1:342,m] = 3.28 * 2 * exp(- weekly_effect[week_index[1:342]])/(1+exp(- weekly_effect[week_index[1:342]]));
Rt_v2[T2:N2,m] = Rt_v1[T2:N2,m] * (R_difference); 
Rt_adj_immune_v1[1:N0,m] = Rt_v1[1:N0,m]; 
for (i in (N0+1):N2) 
{ 
  convolution_v1 = 1e-15;
  for (j in 1:(i-1)) 
  {
    convolution_v1 = convolution_v1+prediction_v1[j, m] * SI[i-j];
    immune_v1[i,m] = immune_v1[i,m]+prediction_v1[j, m] * exp( - 0.5 * (i-j) * (i-j) / (par * par) );
  }
  if ( i > (T2) ) 
  {
    convolution_v2 = 1e-15;
    for(j in T2:(i-1)) 
    {
      convolution_v2 = convolution_v2+prediction_v2[j, m] * SI[i-j];
      immune_v2[i,m] = immune_v2[i,m]+prediction_v2[j, m] * exp( - 0.5 * (i-j) * (i-j) / ( par * par) );
    }
    alpha_sus2[i,m] = (1 - cross) * immune_v2[i,m] / pop;
    n2[i,m] = immune_v2[i,m] + cross * (immune_v1[i,m] * ( 1 - alpha_sus2[i,m] ));
    Rt_adj_immune_v2[i,m] = ( 1 - n2[i,m] / pop) * Rt_v2[i,m]; 
    prediction_v2[i, m] = Rt_adj_immune_v2[i,m] * convolution_v2;
  }
  alpha_sus1[i,m] = (1 - cross) * immune_v1[i,m] / pop;
  n1[i,m] = immune_v1[i,m] + cross * (immune_v2[i,m] * ( 1 - alpha_sus1[i,m] ));
  Rt_adj_immune_v1[i,m] = ( 1 - n1[i,m] / pop) * Rt_v1[i,m];
  prediction_v1[i, m] = Rt_adj_immune_v1[i,m] * convolution_v1;
}
E_deaths_v1[1, m]= 1e-15 * prediction_v1[1,m];
E_deaths_v2[1, m]= 1e-15 * prediction_v2[1,m];
E_deaths[1, m]= 1e-15 * (prediction_v1[1,m] + prediction_v2[1,m]);
for (i in 2:N2)
{
  for(j in 1:(i-1))
  {	
    E_deaths_v1[i,m] = E_deaths_v1[i,m]+ prediction_v1[j,m] * f[i-j,m] * ifr1;
  }
  if (i > T2)  
  {
    for(j in T2:(i-1)) 
    {
      E_deaths_v2[i,m] = E_deaths_v2[i,m]+prediction_v2[j,m] * f[i-j,m] * ifr2;
    }
  }
  E_deaths[i,m] = E_deaths_v1[i,m] + E_deaths_v2[i,m];
}

E_deaths_v1_1<-E_deaths_v1[,1]
E_deaths_v1_2<-E_deaths_v2[,1]

excess= read.csv(here("transmission_model/data/excess_burials.csv"))
colnames(excess) = c("Date","Deaths_raw","Mean7_Deaths","X")
excess = excess[c("Date","Deaths_raw")]
excess$Date = dmy(excess$Date)
excess$Deaths_raw = as.numeric(excess$Deaths_raw)
df_SIVEP_original -> df_SIVEP
df_SIVEP = df_SIVEP[df_SIVEP$ID_MUNICIP == selected_state,]
df_SIVEP = df_SIVEP[which(df_SIVEP$CLASSI_FIN==5) || which(df_SIVEP$CLASSI_FIN==4),]
df_SIVEP = df_SIVEP[which(df_SIVEP$EVOLUCAO==2),]
df_SIVEP$DT_EVOLUCA = dmy(df_SIVEP$DT_EVOLUCA)
filter_date = head(sort(df_SIVEP[df_SIVEP$CLASSI_FIN==5,]$DT_EVOLUCA),1)#only after first confirmed
df_SIVEP = df_SIVEP %>% filter(DT_EVOLUCA >= ymd(filter_date))
df_SIVEP = df_SIVEP[,c("DT_EVOLUCA", "ID_MUNICIP")]
dim(df_SIVEP)
df_SIVEP$Deaths = 1
df_SIVEP = aggregate(. ~ID_MUNICIP+DT_EVOLUCA, data=df_SIVEP, sum, na.rm=TRUE)
sum(df_SIVEP$Deaths)
df_SIVEP$type=rep("SARI mortality data (SIVEP-Gripe)", length(df_SIVEP$Deaths))
colnames(df_SIVEP) = c("region","DateRep","Deaths","type")

Sys.setlocale("LC_TIME", "English")

E_deaths_non_P.1_111 = data.frame(
  "time" = dates[[Country]],
  "type" = rep("Non-P.1", length(dates[[Country]])),
  "rt" = E_deaths_v1_1)

E_deaths_P.1_111 = data.frame(
  "time" = dates[[Country]],
  "type" = rep("P.1", length(dates[[Country]])),
  "rt" = E_deaths_v1_2)

rt_df = rbind(E_deaths_non_P.1_111,E_deaths_P.1_111)
plt_deaths <- function(){
  ggplot() +
    geom_line(data = rt_df,aes(x = time, y = rt, color = type)) +
    scale_color_manual(name = "", labels = c("Non-P.1 (Expected deaths)","P.1 (Expected deaths)"),
                       values = c("Deepskyblue4","tan2"))+
    theme_pubr(base_size = 26) +
    theme(legend.position = c(0.4, 0.9)) +
    xlab("Date") +
    ylab("Daily no. of deaths")+
    ylim(-5,160)+
    geom_point(data = df_SIVEP, aes(x=DateRep, y=Deaths),col="firebrick2",alpha=0.5) +
    geom_point(aes(x=df_SIVEP$DateRep[142], y=144),col="firebrick2",alpha=0.5,size=3) +
    geom_text(aes(label = 'SARI mortality data (SIVEP-Gripe)'), x=df_SIVEP$DateRep[169], y=144,size=4) +
    scale_x_date(date_breaks = "1 months", labels = date_format("%b/%y")) +
    theme(legend.text = element_text(size = 11))+
    #theme(legend.title=element_text(face="italic", family="Times", colour="red",
    #                                    size=14))
    theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))
}

p_F<-plt_deaths()
p_F
ggsave(here(paste0("transmission_model/figures/fig_PG_contour.pdf")),plot=p_F)


