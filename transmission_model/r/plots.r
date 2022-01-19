###############model_normal_dis
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
Sys.setlocale("LC_TIME", "English")
option_list <- list(
    make_option(c("--fileName"),action="store", default=here("transmission_model/results/base/non_P1_P1_model_normal.Rdata"),help="File to be loaded [default \"%default\"]"),
    make_option(c("--beta"),action="store", default="0.2",help="Cross immunity factor [default \"%default\"]")
)

opt <- parse_args(OptionParser(option_list=option_list))

load(opt$fileName)
out <- rstan::extract(fit)
posterior <- as.matrix(fit)

# posteriors of interest
######################################################### 
posteriors_of_interest <- 
    bind_rows(
        tibble(parameter = "Transmissibility Increase(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(out$R_difference - 1,probs = c(0.25)),0), round(100*quantile(out$R_difference - 1,probs = c(0.75)),0))
        ),
        tibble(parameter = "Relative Risk Increase(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(out$RR - 1,probs = c(0.25)),0), round(100*quantile(out$RR - 1,probs = c(0.75)),0))
        ),
        tibble(parameter = "Cross Immunity(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(1 - out$cross,probs = c(0.25)),0), round(100*quantile( 1 - out$cross,probs = c(0.75)),0))
        )
    )
print(posteriors_of_interest)
grid.table(posteriors_of_interest)
############################################################
# cross immunity vs transmissibility
############################################################
x <- data.frame(R_difference = out$R_difference, cross =out$cross)
kd <- ks::kde(x, compute.cont=TRUE,H=matrix(c(0.01,0,0,0.0033),2))
contour_95 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["5%"])[[1]]))
contour_75 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["25%"])[[1]]))
contour_50 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["50%"])[[1]]))
contour_25 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["75%"])[[1]]))
contour_5 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                              z=estimate, levels=cont["95%"])[[1]]))
pA <- mcmc_scatter(posterior, pars = c("R_difference","cross"), size = 3.5, alpha = 0.1) +
    ylab("Cross-immunity") +
    xlab("Transmissibility increase") +
    scale_y_continuous(expand=c(0,0)) +
    theme_pubr(base_size = 26) +
    theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
    ggplot2::geom_path(aes(x, y), data=contour_95) +
    ggplot2::geom_path(aes(x, y), data=contour_75) +
    ggplot2::geom_path(aes(x, y), data=contour_50) +
    ggplot2::geom_path(aes(x, y), data=contour_25) +
    ggplot2::geom_path(aes(x, y), data=contour_5) 
pA <- ggMarginal(pA, R_difference, cross, type = c("density"),
                 margins = 'both',
                 size = 4,
                 colour = "#589e73",
                 fill = '#589e73',
                 alpha = 0.25,
                 xparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)),
                 yparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)))
pA
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PA_contour.pdf")),plot=pA)
############################################################
# cross immunity vs relative risk increase
############################################################
x <- data.frame(R_difference = out$RR, cross = out$cross)
kd <- ks::kde(x, compute.cont=TRUE,H=matrix(c(0.015,0,0,0.0075),2))
contour_95 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["5%"])[[1]]))
contour_75 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["25%"])[[1]]))
contour_50 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["50%"])[[1]]))
contour_25 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["75%"])[[1]]))
contour_5 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                              z=estimate, levels=cont["95%"])[[1]]))

pB <-  mcmc_scatter(posterior, pars = c("RR[1]","cross"), size = 3.5, alpha = 0.1) +
    ylab("Cross-immunity") +
    xlab("Relative risk of mortality") +
    scale_y_continuous(expand=c(0,0)) +
    theme_pubr(base_size = 26) +
    theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
    xlim(c(0,4)) +
    ggplot2::geom_path(aes(x, y), data=contour_95) +
    ggplot2::geom_path(aes(x, y), data=contour_75) +
    ggplot2::geom_path(aes(x, y), data=contour_50) +
    ggplot2::geom_path(aes(x, y), data=contour_25) +
    ggplot2::geom_path(aes(x, y), data=contour_5) 
pB <- ggMarginal(pB, R_difference, cross, type = c("density"),
                 margins = 'both',
                 size = 4,
                 colour = "#589e73",
                 fill = '#589e73',
                 alpha = 0.25,
                 xparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)),
                 yparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)))
theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))

pB
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PB_contour.pdf")),plot=pB)
############################################################
# Deaths Curve
############################################################
plt_deaths <- function(outcheck)
{
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
    
    deaths_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v1", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.025),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.975)
    )
    
    deaths_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v2", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.025),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.975)
    )
    
    deaths_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v1", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.25),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.75)
    )
    
    deaths_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v2", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.25),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.75)
    )
    
    deaths_df = rbind(deaths_df95v1,deaths_df95v2,deaths_df50v1,deaths_df50v2)
    
    
    ggplot() + 
        
        geom_ribbon(data = deaths_df, 
                    aes(x=time, ymax=deaths_l, ymin=deaths_u, fill = key)) +
        scale_fill_manual(name = "", labels = c("Non-P.1 (50% CI)","P.1 (50% CI)","Non-P.1 (95% CI)","P.1 (50% CI)"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("tan2",0.85),
                                     alpha("Deepskyblue4",0.45),
                                     alpha("tan2",0.45))) +
        guides(fill=guide_legend(ncol=2)) +
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.54, 0.95)) +
        xlab("Date") +
        ylab("Daily no. of deaths")+
        ylim(-5,160)+
        geom_point(data = df_SIVEP, aes(x=DateRep, y=Deaths),col="firebrick2",alpha=0.5) +
        geom_point(aes(x=df_SIVEP$DateRep[87], y=115),col="firebrick2",alpha=0.5,size=3) +
        geom_text(aes(label = 'SARI mortality data (SIVEP-Gripe)'), x=df_SIVEP$DateRep[161], y=115,size=4) +
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(legend.text = element_text(size = 11))+
        #theme(legend.title=element_text(face="italic", family="Times", colour="red",
        #                                    size=14))
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))

}

p_C <- plt_deaths(out)
p_C
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PC_contour.pdf")),plot=p_C)
############################################################
# P1 fraction curve
############################################################
plt_P1fraction <- function(outcheck)
{
    datapoints = readRDS(here("transmission_model/data/datapoints.rds"))
    colnames(datapoints)[1] <- "date"
    datapoints = datapoints[datapoints$date > ymd("2020-11-01"),]
    #     dates[[Country]] = c(dates[[Country]],seq(tail(dates[[Country]],1)+1,tail(dates[[Country]],1)+342-311,by = '1 day'))
    E_fraction_df_95 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("nintyfive", length(dates[[Country]])),
        "E_fraction" = colMeans(outcheck$E_fraction),
        "E_fraction_ui" = colQuantiles(outcheck$E_fraction[,,], probs=.975),
        "E_fraction_li" = colQuantiles(outcheck$E_fraction[,,], probs=.025)
    )
    E_fraction_df_50 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50", length(dates[[Country]])),
        "E_fraction" = colMeans(outcheck$E_fraction),
        "E_fraction_ui" = colQuantiles(outcheck$E_fraction[,,], probs=.75),
        "E_fraction_li" = colQuantiles(outcheck$E_fraction[,,], probs=.25)
    )
    E_fraction_df = rbind(E_fraction_df_95,E_fraction_df_50)
    late_E_fraction_df = E_fraction_df[E_fraction_df$time >= ymd("2020-11-06"),]
    
    ggplot() + 
        geom_point(data = datapoints, aes(x = date, y = prop)) +
        geom_point(aes(x = datapoints$date[1], y = 1.15), size=3)+
        geom_text(aes(label = 'Independent estimates from sequence data from Manaus'), x=ymd("2020-12-07"), y=1.15,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2021-02-01"), y=0.9,size=4)+
        geom_errorbar(data = datapoints, aes(x = date, ymin=lower_CI, ymax=upper_CI,width=2))+
        geom_ribbon(data = late_E_fraction_df, aes(x=time, ymax=E_fraction_ui, ymin=E_fraction_li, fill = key))+
        
        scale_fill_manual(name = "", labels = c("50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("Deepskyblue4",0.55))) +
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.2, 0.9)) +
        xlab("Date") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab("P1 fraction") +
        ylim(-0.05,1.2)+
        scale_x_date(date_breaks = "1 month", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
    
}

p_D <- plt_P1fraction(out)
p_D
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PD_contour.pdf")),plot=p_D)
############################################################
# Seropervalence Curve
############################################################
plt_sero_conv <- function(outcheck)
{
    manausPopulation = 2219580
    seroconv_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v1", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.025)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.975)/stan_data$pop[[1]]
    )
    seroconv_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v2", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.025)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.975)/stan_data$pop[[1]]
    )
    seroconv_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v1", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.25)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.75)/stan_data$pop[[1]]
    )
    seroconv_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v2", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v2)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.25)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.75)/stan_data$pop[[1]]
    )
    seroconv_df = rbind(seroconv_df95v1,seroconv_df95v2,seroconv_df50v1,seroconv_df50v2)
    
    
    manaus_seroprev = read.csv(here("transmission_model/data/manaus_seroprev.csv"))
    manaus_seroprev$date = ymd(manaus_seroprev$date)
    manaus_seroprev = manaus_seroprev[manaus_seroprev$upper < 150 & manaus_seroprev$upper > 0 ,] 
    manaus_seroprev$sero_sigma = (manaus_seroprev$prevalence - manaus_seroprev$lower)/2 
    tail(manaus_seroprev,1)
    
    ggplot() +
        geom_point(data=manaus_seroprev,aes(x=date,y=prevalence/100)) +
        geom_point(aes(x = manaus_seroprev$date[1], y = 1.1), size=3)+
        geom_text(aes(label = 'Independent seroprevalence estimates from Manaus'), x=ymd("2020-06-29"), y=1.1,size=4)+
        geom_text(aes(label = 'non-P.1'), x=ymd("2021-01-07"), y=1.07,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2021-01-07"), y=0.45,size=4)+
        geom_ribbon(data=seroconv_df,aes(x = time, ymin = seroconv_l, ymax = seroconv_u, fill =key)) +
        scale_fill_manual(name = "", labels = c("50%", "95%","50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("tan2",0.85),
                                     alpha("Deepskyblue4",0.55),
                                     alpha("tan2",0.55))) +
        geom_errorbar(data=manaus_seroprev,aes(x=date,ymin=lower/100, ymax=upper/100,width=10)) +
        
        xlab("Date") +
        ylab("Cumulative incidence per capita\n") +
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.2, 0.9)) +
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
    
}
p_E <-  plt_sero_conv(out)
p_E
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PE_contour.pdf")),plot=p_E)
############################################################
# Rt Curve
############################################################
############################################################
plt_rt <- function(outcheck)
{
    rt_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v1_95)", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v1)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.025),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.975)
    )
    rt_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v2_95", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v2)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.025),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.975)
    )
    rt_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v1_50", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v1),
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.25),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.75)
    )
    rt_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v2_50", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v2)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.25),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.75)
    )
    
    rt_df = rbind(rt_df50v1,rt_df50v2,rt_df95v1,rt_df95v2)
    tail(rt_df,1)
    ggplot() +
        geom_ribbon(data=rt_df,aes(x = time, ymin = rt_l, ymax = rt_u, fill =key)) +
        scale_fill_manual(name = "", labels = c("50%", "95%","50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("Deepskyblue4",0.55),
                                     alpha("tan2",0.85),
                                     alpha("tan2",0.55))) +
        geom_text(aes(label = 'non-P.1'), x=ymd("2020-03-12"), y=6,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2020-11-25"), y=6.1,size=4)+
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.25, 0.9)) +
        xlab("Date") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab(expression(Reproduction ~ number * "," ~ R[t]))+
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
}


p_F <- plt_rt(out)
p_F
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PF_contour.pdf")),plot=p_F)

library(cowplot)
plot_grid(plot_grid(pA, p_C, p_E, labels = c("A", "C", "E"),nrow = 1, align = "h",axis = "b",
                    rel_widths = c(0.6,1,1)), 
          plot_grid(pB, p_D, p_F, labels = c("B", "D", "F"),nrow = 1, align = "h",axis = "b",
                    rel_widths = c(0.6,1,1)),
          nrow = 2 )
ggsave(here(paste0("transmission_model/figures/",job,"/fig_4.pdf")))




###############model_original
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
Sys.setlocale("LC_TIME", "English")
option_list <- list(
    make_option(c("--fileName"),action="store", default=here("transmission_model/results/base/non_P1_P1_model_original.Rdata"),help="File to be loaded [default \"%default\"]"),
    make_option(c("--beta"),action="store", default="0.2",help="Cross immunity factor [default \"%default\"]")
)

opt <- parse_args(OptionParser(option_list=option_list))

load(opt$fileName)
out <- rstan::extract(fit)
posterior <- as.matrix(fit)

# posteriors of interest
######################################################### 
posteriors_of_interest <- 
    bind_rows(
        tibble(parameter = "Transmissibility Increase(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(out$R_difference - 1,probs = c(0.25)),0), round(100*quantile(out$R_difference - 1,probs = c(0.75)),0))
        ),
        tibble(parameter = "Relative Risk Increase(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(out$RR - 1,probs = c(0.25)),0), round(100*quantile(out$RR - 1,probs = c(0.75)),0))
        ),
        tibble(parameter = "Cross Immunity(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(1 - out$cross,probs = c(0.25)),0), round(100*quantile( 1 - out$cross,probs = c(0.75)),0))
        )
    )
print(posteriors_of_interest)
grid.table(posteriors_of_interest)
############################################################
# cross immunity vs transmissibility
############################################################
x <- data.frame(R_difference = out$R_difference, cross =out$cross)
kd <- ks::kde(x, compute.cont=TRUE,H=matrix(c(0.01,0,0,0.0033),2))
contour_95 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["5%"])[[1]]))
contour_75 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["25%"])[[1]]))
contour_50 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["50%"])[[1]]))
contour_25 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["75%"])[[1]]))
contour_5 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                              z=estimate, levels=cont["95%"])[[1]]))
pA <- mcmc_scatter(posterior[, c(1,5)], pars = c("R_difference","cross"), size = 3.5, alpha = 0.1) +
    ylab("Cross-immunity") +
    xlab("Transmissibility increase") +
    scale_y_continuous(expand=c(0,0)) +
    theme_pubr(base_size = 26) +
    theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
    ggplot2::geom_path(aes(x, y), data=contour_95) +
    ggplot2::geom_path(aes(x, y), data=contour_75) +
    ggplot2::geom_path(aes(x, y), data=contour_50) +
    ggplot2::geom_path(aes(x, y), data=contour_25) +
    ggplot2::geom_path(aes(x, y), data=contour_5) 
pA <- ggMarginal(pA, R_difference, cross, type = c("density"),
                 margins = 'both',
                 size = 4,
                 colour = "#589e73",
                 fill = '#589e73',
                 alpha = 0.25,
                 xparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)),
                 yparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)))
pA
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PA_contour.pdf")),plot=pA)
############################################################
# cross immunity vs relative risk increase
############################################################
x <- data.frame(R_difference = out$RR, cross = out$cross)
kd <- ks::kde(x, compute.cont=TRUE,H=matrix(c(0.015,0,0,0.0075),2))
contour_95 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["5%"])[[1]]))
contour_75 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["25%"])[[1]]))
contour_50 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["50%"])[[1]]))
contour_25 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["75%"])[[1]]))
contour_5 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                              z=estimate, levels=cont["95%"])[[1]]))
pB <-  mcmc_scatter(posterior[, c(5,8)], pars = c("RR[1]","cross"), size = 3.5, alpha = 0.1) +
    ylab("Cross-immunity") +
    xlab("Relative risk of mortality") +
    scale_y_continuous(expand=c(0,0)) +
    theme_pubr(base_size = 26) +
    theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
    xlim(c(0,4)) +
    ggplot2::geom_path(aes(x, y), data=contour_95) +
    ggplot2::geom_path(aes(x, y), data=contour_75) +
    ggplot2::geom_path(aes(x, y), data=contour_50) +
    ggplot2::geom_path(aes(x, y), data=contour_25) +
    ggplot2::geom_path(aes(x, y), data=contour_5) 
pB <- ggMarginal(pB, R_difference, cross, type = c("density"),
                 margins = 'both',
                 size = 4,
                 colour = "#589e73",
                 fill = '#589e73',
                 alpha = 0.25,
                 xparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)),
                 yparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)))
theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))

pB
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PB_contour.pdf")),plot=pB)
############################################################
# Deaths Curve
############################################################
plt_deaths <- function(outcheck)
{
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
    
    deaths_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v1", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.025),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.975)
    )
    
    deaths_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v2", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.025),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.975)
    )
    
    deaths_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v1", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.25),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.75)
    )
    
    deaths_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v2", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.25),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.75)
    )
    
    deaths_df = rbind(deaths_df95v1,deaths_df95v2,deaths_df50v1,deaths_df50v2)
    
    
    ggplot() + 
        
        geom_ribbon(data = deaths_df, 
                    aes(x=time, ymax=deaths_l, ymin=deaths_u, fill = key)) +
        scale_fill_manual(name = "", labels = c("Non-P.1 (50% CI)","P.1 (50% CI)","Non-P.1 (95% CI)","P.1 (50% CI)"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("tan2",0.85),
                                     alpha("Deepskyblue4",0.45),
                                     alpha("tan2",0.45))) +
        guides(fill=guide_legend(ncol=2)) +
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.54, 0.95)) +
        xlab("Date") +
        ylab("Daily no. of deaths")+
        ylim(-5,160)+
        geom_point(data = df_SIVEP, aes(x=DateRep, y=Deaths),col="firebrick2",alpha=0.5) +
        geom_point(aes(x=df_SIVEP$DateRep[87], y=115),col="firebrick2",alpha=0.5,size=3) +
        geom_text(aes(label = 'SARI mortality data (SIVEP-Gripe)'), x=df_SIVEP$DateRep[161], y=115,size=4) +
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(legend.text = element_text(size = 11))+
        #theme(legend.title=element_text(face="italic", family="Times", colour="red",
        #                                    size=14))
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))
    
}

p_C <- plt_deaths(out)
p_C
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PC_contour.pdf")),plot=p_C)
############################################################
# P1 fraction curve
############################################################
plt_P1fraction <- function(outcheck)
{
    datapoints = readRDS(here("transmission_model/data/datapoints.rds"))
    colnames(datapoints)[1] <- "date"
    datapoints = datapoints[datapoints$date > ymd("2020-11-01"),]
    #     dates[[Country]] = c(dates[[Country]],seq(tail(dates[[Country]],1)+1,tail(dates[[Country]],1)+342-311,by = '1 day'))
    E_fraction_df_95 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("nintyfive", length(dates[[Country]])),
        "E_fraction" = colMeans(outcheck$E_fraction),
        "E_fraction_ui" = colQuantiles(outcheck$E_fraction[,,], probs=.975),
        "E_fraction_li" = colQuantiles(outcheck$E_fraction[,,], probs=.025)
    )
    E_fraction_df_50 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50", length(dates[[Country]])),
        "E_fraction" = colMeans(outcheck$E_fraction),
        "E_fraction_ui" = colQuantiles(outcheck$E_fraction[,,], probs=.75),
        "E_fraction_li" = colQuantiles(outcheck$E_fraction[,,], probs=.25)
    )
    E_fraction_df = rbind(E_fraction_df_95,E_fraction_df_50)
    late_E_fraction_df = E_fraction_df[E_fraction_df$time >= ymd("2020-11-06"),]
    
    ggplot() + 
        geom_point(data = datapoints, aes(x = date, y = prop)) +
        geom_point(aes(x = datapoints$date[1], y = 1.15), size=3)+
        geom_text(aes(label = 'Independent estimates from sequence data from Manaus'), x=ymd("2020-12-07"), y=1.15,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2021-02-01"), y=0.9,size=4)+
        geom_errorbar(data = datapoints, aes(x = date, ymin=lower_CI, ymax=upper_CI,width=2))+
        geom_ribbon(data = late_E_fraction_df, aes(x=time, ymax=E_fraction_ui, ymin=E_fraction_li, fill = key))+
        
        scale_fill_manual(name = "", labels = c("50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("Deepskyblue4",0.55))) +
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.2, 0.9)) +
        xlab("Date") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab("P1 fraction") +
        ylim(-0.05,1.2)+
        scale_x_date(date_breaks = "1 month", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
    
}

p_D <- plt_P1fraction(out)
p_D
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PD_contour.pdf")),plot=p_D)
############################################################
# Seropervalence Curve
############################################################
plt_sero_conv <- function(outcheck)
{
    manausPopulation = 2219580
    seroconv_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v1", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.025)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.975)/stan_data$pop[[1]]
    )
    seroconv_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v2", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.025)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.975)/stan_data$pop[[1]]
    )
    seroconv_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v1", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.25)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.75)/stan_data$pop[[1]]
    )
    seroconv_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v2", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v2)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.25)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.75)/stan_data$pop[[1]]
    )
    seroconv_df = rbind(seroconv_df95v1,seroconv_df95v2,seroconv_df50v1,seroconv_df50v2)
    
    
    manaus_seroprev = read.csv(here("transmission_model/data/manaus_seroprev.csv"))
    manaus_seroprev$date = ymd(manaus_seroprev$date)
    manaus_seroprev = manaus_seroprev[manaus_seroprev$upper < 150 & manaus_seroprev$upper > 0 ,] 
    manaus_seroprev$sero_sigma = (manaus_seroprev$prevalence - manaus_seroprev$lower)/2 
    tail(manaus_seroprev,1)
    
    ggplot() +
        geom_point(data=manaus_seroprev,aes(x=date,y=prevalence/100)) +
        geom_point(aes(x = manaus_seroprev$date[1], y = 1.1), size=3)+
        geom_text(aes(label = 'Independent seroprevalence estimates from Manaus'), x=ymd("2020-06-29"), y=1.1,size=4)+
        geom_text(aes(label = 'non-P.1'), x=ymd("2021-01-07"), y=1.07,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2021-01-07"), y=0.45,size=4)+
        geom_ribbon(data=seroconv_df,aes(x = time, ymin = seroconv_l, ymax = seroconv_u, fill =key)) +
        scale_fill_manual(name = "", labels = c("50%", "95%","50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("tan2",0.85),
                                     alpha("Deepskyblue4",0.55),
                                     alpha("tan2",0.55))) +
        geom_errorbar(data=manaus_seroprev,aes(x=date,ymin=lower/100, ymax=upper/100,width=10)) +
        
        xlab("Date") +
        ylab("Cumulative incidence per capita\n") +
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.2, 0.9)) +
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
    
}
p_E <-  plt_sero_conv(out)
p_E
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PE_contour.pdf")),plot=p_E)
############################################################
# Rt Curve
############################################################
############################################################
plt_rt <- function(outcheck)
{
    rt_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v1_95)", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v1)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.025),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.975)
    )
    rt_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v2_95", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v2)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.025),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.975)
    )
    rt_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v1_50", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v1),
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.25),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.75)
    )
    rt_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v2_50", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v2)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.25),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.75)
    )
    
    rt_df = rbind(rt_df50v1,rt_df50v2,rt_df95v1,rt_df95v2)
    tail(rt_df,1)
    ggplot() +
        geom_ribbon(data=rt_df,aes(x = time, ymin = rt_l, ymax = rt_u, fill =key)) +
        scale_fill_manual(name = "", labels = c("50%", "95%","50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("Deepskyblue4",0.55),
                                     alpha("tan2",0.85),
                                     alpha("tan2",0.55))) +
        geom_text(aes(label = 'non-P.1'), x=ymd("2020-03-12"), y=5.1,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2020-11-25"), y=5.4,size=4)+
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.25, 0.9)) +
        xlab("Date") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab(expression(Reproduction ~ number * "," ~ R[t]))+
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
}


p_F <- plt_rt(out)
p_F
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PF_contour.pdf")),plot=p_F)

library(cowplot)
plot_grid(plot_grid(pA, p_C, p_E, labels = c("A", "C", "E"),nrow = 1, align = "h",axis = "b",
                    rel_widths = c(0.6,1,1)), 
          plot_grid(pB, p_D, p_F, labels = c("B", "D", "F"),nrow = 1, align = "h",axis = "b",
                    rel_widths = c(0.6,1,1)),
          nrow = 2 )
ggsave(here(paste0("transmission_model/figures/",job,"/fig_4.pdf")))








###############model_uniform
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
Sys.setlocale("LC_TIME", "English")
option_list <- list(
    make_option(c("--fileName"),action="store", default=here("transmission_model/results/base/non_P1_P1_model_uniform.Rdata"),help="File to be loaded [default \"%default\"]"),
    make_option(c("--beta"),action="store", default="0.2",help="Cross immunity factor [default \"%default\"]")
)

opt <- parse_args(OptionParser(option_list=option_list))

load(opt$fileName)
out <- rstan::extract(fit)
posterior <- as.matrix(fit)

# posteriors of interest
######################################################### 
posteriors_of_interest <- 
    bind_rows(
        tibble(parameter = "Transmissibility Increase(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(out$R_difference - 1,probs = c(0.25)),0), round(100*quantile(out$R_difference - 1,probs = c(0.75)),0))
        ),
        tibble(parameter = "Relative Risk Increase(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(out$RR - 1,probs = c(0.25)),0), round(100*quantile(out$RR - 1,probs = c(0.75)),0))
        ),
        tibble(parameter = "Cross Immunity(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(1 - out$cross,probs = c(0.25)),0), round(100*quantile( 1 - out$cross,probs = c(0.75)),0))
        )
    )
print(posteriors_of_interest)
grid.table(posteriors_of_interest)
############################################################
# cross immunity vs transmissibility
############################################################
x <- data.frame(R_difference = out$R_difference, cross =out$cross)
kd <- ks::kde(x, compute.cont=TRUE,H=matrix(c(0.01,0,0,0.0033),2))
contour_95 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["5%"])[[1]]))
contour_75 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["25%"])[[1]]))
contour_50 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["50%"])[[1]]))
contour_25 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["75%"])[[1]]))
contour_5 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                              z=estimate, levels=cont["95%"])[[1]]))
pA <- mcmc_scatter(posterior[, c(1,5)], pars = c("R_difference","cross"), size = 3.5, alpha = 0.1) +
    ylab("Cross-immunity") +
    xlab("Transmissibility increase") +
    scale_y_continuous(expand=c(0,0)) +
    theme_pubr(base_size = 26) +
    theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
    ggplot2::geom_path(aes(x, y), data=contour_95) +
    ggplot2::geom_path(aes(x, y), data=contour_75) +
    ggplot2::geom_path(aes(x, y), data=contour_50) +
    ggplot2::geom_path(aes(x, y), data=contour_25) +
    ggplot2::geom_path(aes(x, y), data=contour_5) 
pA <- ggMarginal(pA, R_difference, cross, type = c("density"),
                 margins = 'both',
                 size = 4,
                 colour = "#589e73",
                 fill = '#589e73',
                 alpha = 0.25,
                 xparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)),
                 yparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)))
pA
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PA_contour.pdf")),plot=pA)
############################################################
# cross immunity vs relative risk increase
############################################################
x <- data.frame(R_difference = out$RR, cross = out$cross)
kd <- ks::kde(x, compute.cont=TRUE,H=matrix(c(0.015,0,0,0.0075),2))
contour_95 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["5%"])[[1]]))
contour_75 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["25%"])[[1]]))
contour_50 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["50%"])[[1]]))
contour_25 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["75%"])[[1]]))
contour_5 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                              z=estimate, levels=cont["95%"])[[1]]))
pB <-  mcmc_scatter(posterior[, c(5,8)], pars = c("RR[1]","cross"), size = 3.5, alpha = 0.1) +
    ylab("Cross-immunity") +
    xlab("Relative risk of mortality") +
    scale_y_continuous(expand=c(0,0)) +
    theme_pubr(base_size = 26) +
    theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
    xlim(c(0,4)) +
    ggplot2::geom_path(aes(x, y), data=contour_95) +
    ggplot2::geom_path(aes(x, y), data=contour_75) +
    ggplot2::geom_path(aes(x, y), data=contour_50) +
    ggplot2::geom_path(aes(x, y), data=contour_25) +
    ggplot2::geom_path(aes(x, y), data=contour_5) 
pB <- ggMarginal(pB, R_difference, cross, type = c("density"),
                 margins = 'both',
                 size = 4,
                 colour = "#589e73",
                 fill = '#589e73',
                 alpha = 0.25,
                 xparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)),
                 yparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)))
theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))

pB
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PB_contour.pdf")),plot=pB)
############################################################
# Deaths Curve
############################################################
plt_deaths <- function(outcheck)
{
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
    
    deaths_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v1", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.025),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.975)
    )
    
    deaths_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v2", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.025),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.975)
    )
    
    deaths_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v1", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.25),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.75)
    )
    
    deaths_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v2", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.25),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.75)
    )
    
    deaths_df = rbind(deaths_df95v1,deaths_df95v2,deaths_df50v1,deaths_df50v2)
    
    
    ggplot() + 
        
        geom_ribbon(data = deaths_df, 
                    aes(x=time, ymax=deaths_l, ymin=deaths_u, fill = key)) +
        scale_fill_manual(name = "", labels = c("Non-P.1 (50% CI)","P.1 (50% CI)","Non-P.1 (95% CI)","P.1 (50% CI)"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("tan2",0.85),
                                     alpha("Deepskyblue4",0.45),
                                     alpha("tan2",0.45))) +
        guides(fill=guide_legend(ncol=2)) +
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.54, 0.95)) +
        xlab("Date") +
        ylab("Daily no. of deaths")+
        ylim(-5,160)+
        geom_point(data = df_SIVEP, aes(x=DateRep, y=Deaths),col="firebrick2",alpha=0.5) +
        geom_point(aes(x=df_SIVEP$DateRep[87], y=115),col="firebrick2",alpha=0.5,size=3) +
        geom_text(aes(label = 'SARI mortality data (SIVEP-Gripe)'), x=df_SIVEP$DateRep[161], y=115,size=4) +
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(legend.text = element_text(size = 11))+
        #theme(legend.title=element_text(face="italic", family="Times", colour="red",
        #                                    size=14))
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))
    
}

p_C <- plt_deaths(out)
p_C
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PC_contour.pdf")),plot=p_C)
############################################################
# P1 fraction curve
############################################################
plt_P1fraction <- function(outcheck)
{
    datapoints = readRDS(here("transmission_model/data/datapoints.rds"))
    colnames(datapoints)[1] <- "date"
    datapoints = datapoints[datapoints$date > ymd("2020-11-01"),]
    #     dates[[Country]] = c(dates[[Country]],seq(tail(dates[[Country]],1)+1,tail(dates[[Country]],1)+342-311,by = '1 day'))
    E_fraction_df_95 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("nintyfive", length(dates[[Country]])),
        "E_fraction" = colMeans(outcheck$E_fraction),
        "E_fraction_ui" = colQuantiles(outcheck$E_fraction[,,], probs=.975),
        "E_fraction_li" = colQuantiles(outcheck$E_fraction[,,], probs=.025)
    )
    E_fraction_df_50 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50", length(dates[[Country]])),
        "E_fraction" = colMeans(outcheck$E_fraction),
        "E_fraction_ui" = colQuantiles(outcheck$E_fraction[,,], probs=.75),
        "E_fraction_li" = colQuantiles(outcheck$E_fraction[,,], probs=.25)
    )
    E_fraction_df = rbind(E_fraction_df_95,E_fraction_df_50)
    late_E_fraction_df = E_fraction_df[E_fraction_df$time >= ymd("2020-11-06"),]
    
    ggplot() + 
        geom_point(data = datapoints, aes(x = date, y = prop)) +
        geom_point(aes(x = datapoints$date[1], y = 1.15), size=3)+
        geom_text(aes(label = 'Independent estimates from sequence data from Manaus'), x=ymd("2020-12-07"), y=1.15,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2021-02-01"), y=0.9,size=4)+
        geom_errorbar(data = datapoints, aes(x = date, ymin=lower_CI, ymax=upper_CI,width=2))+
        geom_ribbon(data = late_E_fraction_df, aes(x=time, ymax=E_fraction_ui, ymin=E_fraction_li, fill = key))+
        
        scale_fill_manual(name = "", labels = c("50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("Deepskyblue4",0.55))) +
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.2, 0.9)) +
        xlab("Date") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab("P1 fraction") +
        ylim(-0.05,1.2)+
        scale_x_date(date_breaks = "1 month", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
    
}

p_D <- plt_P1fraction(out)
p_D
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PD_contour.pdf")),plot=p_D)
############################################################
# Seropervalence Curve
############################################################
plt_sero_conv <- function(outcheck)
{
    manausPopulation = 2219580
    seroconv_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v1", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.025)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.975)/stan_data$pop[[1]]
    )
    seroconv_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v2", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.025)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.975)/stan_data$pop[[1]]
    )
    seroconv_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v1", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.25)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.75)/stan_data$pop[[1]]
    )
    seroconv_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v2", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v2)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.25)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.75)/stan_data$pop[[1]]
    )
    seroconv_df = rbind(seroconv_df95v1,seroconv_df95v2,seroconv_df50v1,seroconv_df50v2)
    
    
    manaus_seroprev = read.csv(here("transmission_model/data/manaus_seroprev.csv"))
    manaus_seroprev$date = ymd(manaus_seroprev$date)
    manaus_seroprev = manaus_seroprev[manaus_seroprev$upper < 150 & manaus_seroprev$upper > 0 ,] 
    manaus_seroprev$sero_sigma = (manaus_seroprev$prevalence - manaus_seroprev$lower)/2 
    tail(manaus_seroprev,1)
    
    ggplot() +
        geom_point(data=manaus_seroprev,aes(x=date,y=prevalence/100)) +
        geom_point(aes(x = manaus_seroprev$date[1], y = 1.1), size=3)+
        geom_text(aes(label = 'Independent seroprevalence estimates from Manaus'), x=ymd("2020-06-29"), y=1.1,size=4)+
        geom_text(aes(label = 'non-P.1'), x=ymd("2021-01-07"), y=1.07,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2021-01-07"), y=0.45,size=4)+
        geom_ribbon(data=seroconv_df,aes(x = time, ymin = seroconv_l, ymax = seroconv_u, fill =key)) +
        scale_fill_manual(name = "", labels = c("50%", "95%","50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("tan2",0.85),
                                     alpha("Deepskyblue4",0.55),
                                     alpha("tan2",0.55))) +
        geom_errorbar(data=manaus_seroprev,aes(x=date,ymin=lower/100, ymax=upper/100,width=10)) +
        
        xlab("Date") +
        ylab("Cumulative incidence per capita\n") +
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.2, 0.9)) +
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
    
}
p_E <-  plt_sero_conv(out)
p_E
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PE_contour.pdf")),plot=p_E)
############################################################
# Rt Curve
############################################################
############################################################
plt_rt <- function(outcheck)
{
    rt_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v1_95)", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v1)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.025),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.975)
    )
    rt_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v2_95", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v2)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.025),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.975)
    )
    rt_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v1_50", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v1),
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.25),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.75)
    )
    rt_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v2_50", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v2)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.25),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.75)
    )
    
    rt_df = rbind(rt_df50v1,rt_df50v2,rt_df95v1,rt_df95v2)
    tail(rt_df,1)
    ggplot() +
        geom_ribbon(data=rt_df,aes(x = time, ymin = rt_l, ymax = rt_u, fill =key)) +
        scale_fill_manual(name = "", labels = c("50%", "95%","50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("Deepskyblue4",0.55),
                                     alpha("tan2",0.85),
                                     alpha("tan2",0.55))) +
        geom_text(aes(label = 'non-P.1'), x=ymd("2020-03-12"), y=5.1,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2020-11-25"), y=5.4,size=4)+
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.25, 0.9)) +
        xlab("Date") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab(expression(Reproduction ~ number * "," ~ R[t]))+
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
}


p_F <- plt_rt(out)
p_F
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PF_contour.pdf")),plot=p_F)

library(cowplot)
plot_grid(plot_grid(pA, p_C, p_E, labels = c("A", "C", "E"),nrow = 1, align = "h",axis = "b",
                    rel_widths = c(0.6,1,1)), 
          plot_grid(pB, p_D, p_F, labels = c("B", "D", "F"),nrow = 1, align = "h",axis = "b",
                    rel_widths = c(0.6,1,1)),
          nrow = 2 )
ggsave(here(paste0("transmission_model/figures/",job,"/fig_4.pdf")))








###############model_0.8%_rr=1
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
Sys.setlocale("LC_TIME", "English")
option_list <- list(
    make_option(c("--fileName"),action="store", default=here("transmission_model/results/base/non_P1_P1_model_0.8%_rr=1.Rdata"),help="File to be loaded [default \"%default\"]"),
    make_option(c("--beta"),action="store", default="0.2",help="Cross immunity factor [default \"%default\"]")
)

opt <- parse_args(OptionParser(option_list=option_list))

load(opt$fileName)
out <- rstan::extract(fit)
posterior <- as.matrix(fit)

# posteriors of interest
######################################################### 
posteriors_of_interest <- 
    bind_rows(
        tibble(parameter = "Transmissibility Increase(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(out$R_difference - 1,probs = c(0.25)),0), round(100*quantile(out$R_difference - 1,probs = c(0.75)),0))
        ),
        tibble(parameter = "Cross Immunity(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(1 - out$cross,probs = c(0.25)),0), round(100*quantile( 1 - out$cross,probs = c(0.75)),0))
        )
    )
print(posteriors_of_interest)
grid.table(posteriors_of_interest)
############################################################
# cross immunity vs transmissibility
############################################################
x <- data.frame(R_difference = out$R_difference, cross =out$cross)
kd <- ks::kde(x, compute.cont=TRUE,H=matrix(c(0.01,0,0,0.0033),2))
contour_95 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["5%"])[[1]]))
contour_75 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["25%"])[[1]]))
contour_50 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["50%"])[[1]]))
contour_25 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["75%"])[[1]]))
contour_5 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                              z=estimate, levels=cont["95%"])[[1]]))
# Deaths Curve
############################################################
plt_deaths <- function(outcheck)
{
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
    
    deaths_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v1", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.025),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.975)
    )
    
    deaths_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v2", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.025),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.975)
    )
    
    deaths_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v1", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.25),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.75)
    )
    
    deaths_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v2", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.25),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.75)
    )
    
    deaths_df = rbind(deaths_df95v1,deaths_df95v2,deaths_df50v1,deaths_df50v2)
    
    
    ggplot() + 
        
        geom_ribbon(data = deaths_df, 
                    aes(x=time, ymax=deaths_l, ymin=deaths_u, fill = key)) +
        scale_fill_manual(name = "", labels = c("Non-P.1 (50% CI)","P.1 (50% CI)","Non-P.1 (95% CI)","P.1 (50% CI)"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("tan2",0.85),
                                     alpha("Deepskyblue4",0.45),
                                     alpha("tan2",0.45))) +
        guides(fill=guide_legend(ncol=2)) +
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.54, 0.95)) +
        xlab("Date") +
        ylab("Daily no. of deaths")+
        ylim(-5,160)+
        geom_point(data = df_SIVEP, aes(x=DateRep, y=Deaths),col="firebrick2",alpha=0.5) +
        geom_point(aes(x=df_SIVEP$DateRep[90], y=115),col="firebrick2",alpha=0.5,size=3) +
        geom_text(aes(label = 'SARI mortality data (SIVEP-Gripe)'), x=df_SIVEP$DateRep[146], y=115,size=4) +
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(legend.text = element_text(size = 11))+
        #theme(legend.title=element_text(face="italic", family="Times", colour="red",
        #                                    size=14))
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))
    
}

p_C <- plt_deaths(out)
p_C
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PC_contour.pdf")),plot=p_C)
############################################################
# P1 fraction curve
############################################################
plt_P1fraction <- function(outcheck)
{
    datapoints = readRDS(here("transmission_model/data/datapoints.rds"))
    colnames(datapoints)[1] <- "date"
    datapoints = datapoints[datapoints$date > ymd("2020-11-01"),]
    #     dates[[Country]] = c(dates[[Country]],seq(tail(dates[[Country]],1)+1,tail(dates[[Country]],1)+342-311,by = '1 day'))
    E_fraction_df_95 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("nintyfive", length(dates[[Country]])),
        "E_fraction" = colMeans(outcheck$E_fraction),
        "E_fraction_ui" = colQuantiles(outcheck$E_fraction[,,], probs=.975),
        "E_fraction_li" = colQuantiles(outcheck$E_fraction[,,], probs=.025)
    )
    E_fraction_df_50 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50", length(dates[[Country]])),
        "E_fraction" = colMeans(outcheck$E_fraction),
        "E_fraction_ui" = colQuantiles(outcheck$E_fraction[,,], probs=.75),
        "E_fraction_li" = colQuantiles(outcheck$E_fraction[,,], probs=.25)
    )
    E_fraction_df = rbind(E_fraction_df_95,E_fraction_df_50)
    late_E_fraction_df = E_fraction_df[E_fraction_df$time >= ymd("2020-11-06"),]
    
    ggplot() + 
        geom_point(data = datapoints, aes(x = date, y = prop)) +
        geom_point(aes(x = datapoints$date[1], y = 1.15), size=3)+
        geom_text(aes(label = 'Independent estimates from sequence data from Manaus'), x=ymd("2020-11-29"), y=1.15,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2021-02-01"), y=0.9,size=4)+
        geom_errorbar(data = datapoints, aes(x = date, ymin=lower_CI, ymax=upper_CI,width=2))+
        geom_ribbon(data = late_E_fraction_df, aes(x=time, ymax=E_fraction_ui, ymin=E_fraction_li, fill = key))+
        
        scale_fill_manual(name = "", labels = c("50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("Deepskyblue4",0.55))) +
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.2, 0.9)) +
        xlab("Date") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab("P1 fraction") +
        ylim(-0.05,1.2)+
        scale_x_date(date_breaks = "1 month", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
    
}

p_D <- plt_P1fraction(out)
p_D
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PD_contour.pdf")),plot=p_D)
############################################################
# Seropervalence Curve
############################################################
plt_sero_conv <- function(outcheck)
{
    manausPopulation = 2219580
    seroconv_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v1", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.025)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.975)/stan_data$pop[[1]]
    )
    seroconv_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v2", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.025)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.975)/stan_data$pop[[1]]
    )
    seroconv_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v1", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.25)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.75)/stan_data$pop[[1]]
    )
    seroconv_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v2", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v2)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.25)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.75)/stan_data$pop[[1]]
    )
    seroconv_df = rbind(seroconv_df95v1,seroconv_df95v2,seroconv_df50v1,seroconv_df50v2)
    
    
    manaus_seroprev = read.csv(here("transmission_model/data/manaus_seroprev.csv"))
    manaus_seroprev$date = ymd(manaus_seroprev$date)
    manaus_seroprev = manaus_seroprev[manaus_seroprev$upper < 150 & manaus_seroprev$upper > 0 ,] 
    manaus_seroprev$sero_sigma = (manaus_seroprev$prevalence - manaus_seroprev$lower)/2 
    tail(manaus_seroprev,1)
    
    ggplot() +
        geom_point(data=manaus_seroprev,aes(x=date,y=prevalence/100)) +
        geom_point(aes(x = manaus_seroprev$date[1], y = 1.1), size=3)+
        geom_text(aes(label = 'Independent seroprevalence estimates from Manaus'), x=ymd("2020-06-1"), y=1.1,size=4)+
        geom_text(aes(label = 'non-P.1'), x=ymd("2021-01-07"), y=0.95,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2021-01-01"), y=0.15,size=4)+
        geom_ribbon(data=seroconv_df,aes(x = time, ymin = seroconv_l, ymax = seroconv_u, fill =key)) +
        scale_fill_manual(name = "", labels = c("50%", "95%","50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("tan2",0.85),
                                     alpha("Deepskyblue4",0.55),
                                     alpha("tan2",0.55))) +
        geom_errorbar(data=manaus_seroprev,aes(x=date,ymin=lower/100, ymax=upper/100,width=10)) +
        
        xlab("Date") +
        ylab("Cumulative incidence per capita\n") +
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.2, 0.9)) +
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
    
}
p_E <-  plt_sero_conv(out)
p_E
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PE_contour.pdf")),plot=p_E)
############################################################
# Rt Curve
############################################################
############################################################
plt_rt <- function(outcheck)
{
    rt_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v1_95)", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v1)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.025),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.975)
    )
    rt_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v2_95", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v2)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.025),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.975)
    )
    rt_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v1_50", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v1),
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.25),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.75)
    )
    rt_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v2_50", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v2)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.25),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.75)
    )
    
    rt_df = rbind(rt_df50v1,rt_df50v2,rt_df95v1,rt_df95v2)
    tail(rt_df,1)
    ggplot() +
        geom_ribbon(data=rt_df,aes(x = time, ymin = rt_l, ymax = rt_u, fill =key)) +
        scale_fill_manual(name = "", labels = c("50%", "95%","50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("Deepskyblue4",0.55),
                                     alpha("tan2",0.85),
                                     alpha("tan2",0.55))) +
        geom_text(aes(label = 'non-P.1'), x=ymd("2020-03-12"), y=5.0,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2020-11-25"), y=5.7,size=4)+
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.25, 0.9)) +
        xlab("Date") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab(expression(Reproduction ~ number * "," ~ R[t]))+
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
}


p_F <- plt_rt(out)
p_F
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PF_contour.pdf")),plot=p_F)

library(cowplot)
plot_grid(plot_grid(p_C, p_E, labels = c("A", "C"),nrow = 1, align = "h",axis = "b",
                    rel_widths = c(1,1)), 
          plot_grid(p_D, p_F, labels = c("B", "D"),nrow = 1, align = "h",axis = "b",
                    rel_widths = c(1,1)),
          nrow = 2 )
ggsave(here(paste0("transmission_model/figures/",job,"/fig_4.pdf")))






###############model_ifr1=1%_rr=1_remove_reff
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
Sys.setlocale("LC_TIME", "English")
option_list <- list(
    make_option(c("--fileName"),action="store", default=here("transmission_model/results/base/non_P1_P1_model_ifr1=1%_rr=1_remove_reff.Rdata"),help="File to be loaded [default \"%default\"]"),
    make_option(c("--beta"),action="store", default="0.2",help="Cross immunity factor [default \"%default\"]")
)

opt <- parse_args(OptionParser(option_list=option_list))

load(opt$fileName)
out <- rstan::extract(fit)
posterior <- as.matrix(fit)

# posteriors of interest
######################################################### 
posteriors_of_interest <- 
    bind_rows(
        tibble(parameter = "Transmissibility Increase(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(out$R_difference - 1,probs = c(0.25)),0), round(100*quantile(out$R_difference - 1,probs = c(0.75)),0))
        )
    )
print(posteriors_of_interest)
grid.table(posteriors_of_interest)
############################################################
# Deaths Curve
############################################################
plt_deaths <- function(outcheck)
{
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
    
    deaths_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v1", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.025),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.975)
    )
    
    deaths_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v2", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.025),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.975)
    )
    
    deaths_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v1", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.25),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.75)
    )
    
    deaths_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v2", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.25),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.75)
    )
    
    deaths_df = rbind(deaths_df95v1,deaths_df95v2,deaths_df50v1,deaths_df50v2)
    
    
    ggplot() + 
        
        geom_ribbon(data = deaths_df, 
                    aes(x=time, ymax=deaths_l, ymin=deaths_u, fill = key)) +
        scale_fill_manual(name = "", labels = c("Non-P.1 (50% CI)","P.1 (50% CI)","Non-P.1 (95% CI)","P.1 (50% CI)"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("tan2",0.85),
                                     alpha("Deepskyblue4",0.45),
                                     alpha("tan2",0.45))) +
        guides(fill=guide_legend(ncol=2)) +
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.54, 0.95)) +
        xlab("Date") +
        ylab("Daily no. of deaths")+
        ylim(-5,160)+
        geom_point(data = df_SIVEP, aes(x=DateRep, y=Deaths),col="firebrick2",alpha=0.5) +
        geom_point(aes(x=df_SIVEP$DateRep[90], y=115),col="firebrick2",alpha=0.5,size=3) +
        geom_text(aes(label = 'SARI mortality data (SIVEP-Gripe)'), x=df_SIVEP$DateRep[146], y=115,size=4) +
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(legend.text = element_text(size = 11))+
        #theme(legend.title=element_text(face="italic", family="Times", colour="red",
        #                                    size=14))
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))
    
}

p_C <- plt_deaths(out)
p_C
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PC_contour.pdf")),plot=p_C)
############################################################
# P1 fraction curve
############################################################
plt_P1fraction <- function(outcheck)
{
    datapoints = readRDS(here("transmission_model/data/datapoints.rds"))
    colnames(datapoints)[1] <- "date"
    datapoints = datapoints[datapoints$date > ymd("2020-11-01"),]
    #     dates[[Country]] = c(dates[[Country]],seq(tail(dates[[Country]],1)+1,tail(dates[[Country]],1)+342-311,by = '1 day'))
    E_fraction_df_95 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("nintyfive", length(dates[[Country]])),
        "E_fraction" = colMeans(outcheck$E_fraction),
        "E_fraction_ui" = colQuantiles(outcheck$E_fraction[,,], probs=.975),
        "E_fraction_li" = colQuantiles(outcheck$E_fraction[,,], probs=.025)
    )
    E_fraction_df_50 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50", length(dates[[Country]])),
        "E_fraction" = colMeans(outcheck$E_fraction),
        "E_fraction_ui" = colQuantiles(outcheck$E_fraction[,,], probs=.75),
        "E_fraction_li" = colQuantiles(outcheck$E_fraction[,,], probs=.25)
    )
    E_fraction_df = rbind(E_fraction_df_95,E_fraction_df_50)
    late_E_fraction_df = E_fraction_df[E_fraction_df$time >= ymd("2020-11-06"),]
    
    ggplot() + 
        geom_point(data = datapoints, aes(x = date, y = prop)) +
        geom_point(aes(x = datapoints$date[1], y = 1.15), size=3)+
        geom_text(aes(label = 'Independent estimates from sequence data from Manaus'), x=ymd("2020-11-29"), y=1.15,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2021-02-01"), y=0.9,size=4)+
        geom_errorbar(data = datapoints, aes(x = date, ymin=lower_CI, ymax=upper_CI,width=2))+
        geom_ribbon(data = late_E_fraction_df, aes(x=time, ymax=E_fraction_ui, ymin=E_fraction_li, fill = key))+
        
        scale_fill_manual(name = "", labels = c("50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("Deepskyblue4",0.55))) +
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.2, 0.9)) +
        xlab("Date") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab("P1 fraction") +
        ylim(-0.05,1.2)+
        scale_x_date(date_breaks = "1 month", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
    
}

p_D <- plt_P1fraction(out)
p_D
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PD_contour.pdf")),plot=p_D)
############################################################
# Seropervalence Curve
############################################################
plt_sero_conv <- function(outcheck)
{
    manausPopulation = 2219580
    seroconv_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v1", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.025)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.975)/stan_data$pop[[1]]
    )
    seroconv_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v2", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.025)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.975)/stan_data$pop[[1]]
    )
    seroconv_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v1", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.25)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.75)/stan_data$pop[[1]]
    )
    seroconv_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v2", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v2)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.25)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.75)/stan_data$pop[[1]]
    )
    seroconv_df = rbind(seroconv_df95v1,seroconv_df95v2,seroconv_df50v1,seroconv_df50v2)
    
    
    manaus_seroprev = read.csv(here("transmission_model/data/manaus_seroprev.csv"))
    manaus_seroprev$date = ymd(manaus_seroprev$date)
    manaus_seroprev = manaus_seroprev[manaus_seroprev$upper < 150 & manaus_seroprev$upper > 0 ,] 
    manaus_seroprev$sero_sigma = (manaus_seroprev$prevalence - manaus_seroprev$lower)/2 
    tail(manaus_seroprev,1)
    
    ggplot() +
        geom_point(data=manaus_seroprev,aes(x=date,y=prevalence/100)) +
        geom_point(aes(x = manaus_seroprev$date[1], y = 1.1), size=3)+
        geom_text(aes(label = 'Independent seroprevalence estimates from Manaus'), x=ymd("2020-06-1"), y=1.1,size=4)+
        geom_text(aes(label = 'non-P.1'), x=ymd("2021-01-09"), y=0.39,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2021-01-09"), y=0.12,size=4)+
        geom_ribbon(data=seroconv_df,aes(x = time, ymin = seroconv_l, ymax = seroconv_u, fill =key)) +
        scale_fill_manual(name = "", labels = c("50%", "95%","50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("tan2",0.85),
                                     alpha("Deepskyblue4",0.55),
                                     alpha("tan2",0.55))) +
        geom_errorbar(data=manaus_seroprev,aes(x=date,ymin=lower/100, ymax=upper/100,width=10)) +
        
        xlab("Date") +
        ylab("Cumulative incidence per capita\n") +
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.2, 0.9)) +
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
    
}
p_E <-  plt_sero_conv(out)
p_E
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PE_contour.pdf")),plot=p_E)
############################################################
# Rt Curve
############################################################
############################################################
plt_rt <- function(outcheck)
{
    rt_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v1_95)", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v1)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.025),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.975)
    )
    rt_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v2_95", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v2)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.025),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.975)
    )
    rt_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v1_50", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v1),
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.25),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.75)
    )
    rt_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v2_50", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v2)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.25),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.75)
    )
    
    rt_df = rbind(rt_df50v1,rt_df50v2,rt_df95v1,rt_df95v2)
    tail(rt_df,1)
    ggplot() +
        geom_ribbon(data=rt_df,aes(x = time, ymin = rt_l, ymax = rt_u, fill =key)) +
        scale_fill_manual(name = "", labels = c("50%", "95%","50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("Deepskyblue4",0.55),
                                     alpha("tan2",0.85),
                                     alpha("tan2",0.55))) +
        geom_text(aes(label = 'non-P.1'), x=ymd("2020-03-12"), y=4.9,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2020-11-25"), y=5.6,size=4)+
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.25, 0.9)) +
        xlab("Date") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab(expression(Reproduction ~ number * "," ~ R[t]))+
        scale_x_date(date_breaks = "2 months", labels = date_format("%b/%y")) +
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
        theme(legend.position="none")
}


p_F <- plt_rt(out)
p_F
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PF_contour.pdf")),plot=p_F)

library(cowplot)
plot_grid(plot_grid(p_C, p_E, labels = c("A", "C"),nrow = 1, align = "h",axis = "b",
                    rel_widths = c(1,1)), 
          plot_grid(p_D, p_F, labels = c("B", "D"),nrow = 1, align = "h",axis = "b",
                    rel_widths = c(1,1)),
          nrow = 2 )
ggsave(here(paste0("transmission_model/figures/",job,"/fig_4.pdf")))











#############################Plot of R0
#MMM<-read_excel("C://Users/林立新/Desktop/Rt_v1_v2.xlsx",sheet=1,na="NA")
plt_rt <- function(outcheck)
{
    rt_v1_1 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("rt_v1_method1", length(dates[[Country]])),
        "rt"  = c(3.4562913,3.4562913,3.4562913,3.4562913,3.4562913,3.4562913,3.4562913,3.6649339,3.6649339,3.6649339,3.6649339,3.6649339,3.6649339,3.6649339,3.6511880,3.6511880,3.6511880,3.6511880,3.6511880,3.6511880,3.6511880,3.4642038,3.4642038,3.4642038
                  ,3.4642038,3.4642038,3.4642038,3.4642038,2.8396066,2.8396066,2.8396066,2.8396066,2.8396066,2.8396066,2.8396066,1.8708350,1.8708350,1.8708350,1.8708350,1.8708350,1.8708350,1.8708350,1.0162616,1.0162616,1.0162616,1.0162616,1.0162616,1.0162616
                  ,1.0162616,0.8964013,0.8964013,0.8964013,0.8964013,0.8964013,0.8964013,0.8964013,1.1917068,1.1917068,1.1917068,1.1917068,1.1917068,1.1917068,1.1917068,1.6153167,1.6153167,1.6153167,1.6153167,1.6153167,1.6153167,1.6153167,1.9331437,1.9331437
                  ,1.9331437,1.9331437,1.9331437,1.9331437,1.9331437,2.0790953,2.0790953,2.0790953,2.0790953,2.0790953,2.0790953,2.0790953,2.0882066,2.0882066,2.0882066,2.0882066,2.0882066,2.0882066,2.0882066,2.0815359,2.0815359,2.0815359,2.0815359,2.0815359
                  ,2.0815359,2.0815359,2.0660855,2.0660855,2.0660855,2.0660855,2.0660855,2.0660855,2.0660855,2.1634384,2.1634384,2.1634384,2.1634384,2.1634384,2.1634384,2.1634384,2.3873429,2.3873429,2.3873429,2.3873429,2.3873429,2.3873429,2.3873429,2.6398181
                  ,2.6398181,2.6398181,2.6398181,2.6398181,2.6398181,2.6398181,2.8073797,2.8073797,2.8073797,2.8073797,2.8073797,2.8073797,2.8073797,2.8392953,2.8392953,2.8392953,2.8392953,2.8392953,2.8392953,2.8392953,2.8102365,2.8102365,2.8102365,2.8102365
                  ,2.8102365,2.8102365,2.8102365,2.8679786,2.8679786,2.8679786,2.8679786,2.8679786,2.8679786,2.8679786,2.9837702,2.9837702,2.9837702,2.9837702,2.9837702,2.9837702,2.9837702,3.1214750,3.1214750,3.1214750,3.1214750,3.1214750,3.1214750,3.1214750
                  ,3.1131101,3.1131101,3.1131101,3.1131101,3.1131101,3.1131101,3.1131101,3.0920479,3.0920479,3.0920479,3.0920479,3.0920479,3.0920479,3.0920479,3.1752114,3.1752114,3.1752114,3.1752114,3.1752114,3.1752114,3.1752114,3.3318582,3.3318582,3.3318582
                  ,3.3318582,3.3318582,3.3318582,3.3318582,3.4682689,3.4682689,3.4682689,3.4682689,3.4682689,3.4682689,3.4682689,3.5694458,3.5694458,3.5694458,3.5694458,3.5694458,3.5694458,3.5694458,3.5284956,3.5284956,3.5284956,3.5284956,3.5284956,3.5284956
                  ,3.5284956,3.4085541,3.4085541,3.4085541,3.4085541,3.4085541,3.4085541,3.4085541,3.3280102,3.3280102,3.3280102,3.3280102,3.3280102,3.3280102,3.3280102,3.2953224,3.2953224,3.2953224,3.2953224,3.2953224,3.2953224,3.2953224,3.3660457,3.3660457
                  ,3.3660457,3.3660457,3.3660457,3.3660457,3.3660457,3.5755229,3.5755229,3.5755229,3.5755229,3.5755229,3.5755229,3.5755229,3.9527219,3.9527219,3.9527219,3.9527219,3.9527219,3.9527219,3.9527219,4.2323513,4.2323513,4.2323513,4.2323513,4.2323513
                  ,4.2323513,4.2323513,4.3887860,4.3887860,4.3887860,4.3887860,4.3887860,4.3887860,4.3887860,4.3573167,4.3573167,4.3573167,4.3573167,4.3573167,4.3573167,4.3573167,4.0484068,4.0484068,4.0484068,4.0484068,4.0484068,4.0484068,4.0484068,3.4850131
                  ,3.4850131,3.4850131,3.4850131,3.4850131,3.4850131,3.4850131,2.8405428,2.8405428,2.8405428,2.8405428,2.8405428,2.8405428,2.8405428,2.2972192,2.2972192,2.2972192,2.2972192,2.2972192,2.2972192,2.2972192,2.0167519,2.0167519,2.0167519,2.0167519
                  ,2.0167519,2.0167519,2.0167519,2.2608960,2.2608960,2.2608960,2.2608960,2.2608960,2.2608960,2.2608960,2.5680549,2.5680549,2.5680549,2.5680549,2.5680549,2.5680549,2.5680549,2.7032260,2.7032260,2.7032260,2.7032260,2.7032260,2.7032260,2.7032260
                  ,2.7032260,2.7032260,2.7032260,2.7032260,2.7032260,2.7032260)
    )
    rt_v2_1 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("rt_v2_method1", length(dates[[Country]])),
        "rt" = c(0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,7.366274,8.137542,8.137542,8.137542,8.137542,8.137542,8.137542,8.137542,8.701203,8.701203,8.701203,8.701203,8.701203,8.701203,8.701203,9.010355,9.010355,9.010355,9.010355
                 ,9.010355,9.010355,9.010355,8.930846,8.930846,8.930846,8.930846,8.930846,8.930846,8.930846,8.259462,8.259462,8.259462,8.259462,8.259462,8.259462,8.259462,7.079676,7.079676,7.079676,7.079676,7.079676,7.079676,7.079676,5.795990,5.795990,5.795990
                 ,5.795990,5.795990,5.795990,5.795990,4.750104,4.750104,4.750104,4.750104,4.750104,4.750104,4.750104,4.218131,4.218131,4.218131,4.218131,4.218131,4.218131,4.218131,4.739760,4.739760,4.739760,4.739760,4.739760,4.739760,4.739760,5.380523,5.380523
                 ,5.380523,5.380523,5.380523,5.380523,5.380523,5.642108,5.642108,5.642108,5.642108,5.642108,5.642108,5.642108,5.642108,5.642108,5.642108,5.642108,5.642108,5.642108)
    )
    rt_v1_2 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("rt_v1_method2", length(dates[[Country]])),
        "rt" = c(3.3727426,3.3727426,3.3727426,3.3727426,3.3727426,3.3727426,3.3727426,3.9408053,3.9408053,3.9408053,3.9408053,3.9408053,3.9408053,3.9408053,3.9195590,3.9195590,3.9195590,3.9195590,3.9195590,3.9195590,3.9195590,3.8375452,3.8375452,3.8375452
                 ,3.8375452,3.8375452,3.8375452,3.8375452,2.6930792,2.6930792,2.6930792,2.6930792,2.6930792,2.6930792,2.6930792,2.3057366,2.3057366,2.3057366,2.3057366,2.3057366,2.3057366,2.3057366,0.7951525,0.7951525,0.7951525,0.7951525,0.7951525,0.7951525
                 ,0.7951525,0.8017285,0.8017285,0.8017285,0.8017285,0.8017285,0.8017285,0.8017285,1.3588357,1.3588357,1.3588357,1.3588357,1.3588357,1.3588357,1.3588357,2.0127513,2.0127513,2.0127513,2.0127513,2.0127513,2.0127513,2.0127513,1.9799190,1.9799190
                 ,1.9799190,1.9799190,1.9799190,1.9799190,1.9799190,2.0695608,2.0695608,2.0695608,2.0695608,2.0695608,2.0695608,2.0695608,2.2336339,2.2336339,2.2336339,2.2336339,2.2336339,2.2336339,2.2336339,2.1826874,2.1826874,2.1826874,2.1826874,2.1826874
                 ,2.1826874,2.1826874,2.0633756,2.0633756,2.0633756,2.0633756,2.0633756,2.0633756,2.0633756,2.2415898,2.2415898,2.2415898,2.2415898,2.2415898,2.2415898,2.2415898,2.5364645,2.5364645,2.5364645,2.5364645,2.5364645,2.5364645,2.5364645,2.8118984
                 ,2.8118984,2.8118984,2.8118984,2.8118984,2.8118984,2.8118984,3.0849154,3.0849154,3.0849154,3.0849154,3.0849154,3.0849154,3.0849154,3.0609812,3.0609812,3.0609812,3.0609812,3.0609812,3.0609812,3.0609812,2.8181927,2.8181927,2.8181927,2.8181927
                 ,2.8181927,2.8181927,2.8181927,2.7752565,2.7752565,2.7752565,2.7752565,2.7752565,2.7752565,2.7752565,3.1080927,3.1080927,3.1080927,3.1080927,3.1080927,3.1080927,3.1080927,3.5130723,3.5130723,3.5130723,3.5130723,3.5130723,3.5130723,3.5130723
                 ,3.3843250,3.3843250,3.3843250,3.3843250,3.3843250,3.3843250,3.3843250,3.0762421,3.0762421,3.0762421,3.0762421,3.0762421,3.0762421,3.0762421,3.1089245,3.1089245,3.1089245,3.1089245,3.1089245,3.1089245,3.1089245,3.5162905,3.5162905,3.5162905
                 ,3.5162905,3.5162905,3.5162905,3.5162905,3.6596881,3.6596881,3.6596881,3.6596881,3.6596881,3.6596881,3.6596881,3.8160522,3.8160522,3.8160522,3.8160522,3.8160522,3.8160522,3.8160522,3.6057639,3.6057639,3.6057639,3.6057639,3.6057639,3.6057639
                 ,3.6057639,3.4208496,3.4208496,3.4208496,3.4208496,3.4208496,3.4208496,3.4208496,3.4823080,3.4823080,3.4823080,3.4823080,3.4823080,3.4823080,3.4823080,3.5087825,3.5087825,3.5087825,3.5087825,3.5087825,3.5087825,3.5087825,3.4856190,3.4856190
                 ,3.4856190,3.4856190,3.4856190,3.4856190,3.4856190,3.4559032,3.4559032,3.4559032,3.4559032,3.4559032,3.4559032,3.4559032,4.2555239,4.2555239,4.2555239,4.2555239,4.2555239,4.2555239,4.2555239,4.1888270,4.1888270,4.1888270,4.1888270,4.1888270
                 ,4.1888270,4.1888270,4.5022800,4.5022800,4.5022800,4.5022800,4.5022800,4.5022800,4.5022800,4.5496121,4.5496121,4.5496121,4.5496121,4.5496121,4.5496121,4.5496121,4.3403041,4.3403041,4.3403041,4.3403041,4.3403041,4.3403041,4.3403041,3.5641051
                 ,3.5641051,3.5641051,3.5641051,3.5641051,3.5641051,3.5641051,2.9089800,2.9089800,2.9089800,2.9089800,2.9089800,2.9089800,2.9089800,2.7528399,2.7528399,2.7528399,2.7528399,2.7528399,2.7528399,2.7528399,1.7014939,1.7014939,1.7014939,1.7014939
                 ,1.7014939,1.7014939,1.7014939,2.7317767,2.7317767,2.7317767,2.7317767,2.7317767,2.7317767,2.7317767,3.5739580,3.5739580,3.5739580,3.5739580,3.5739580,3.5739580,3.5739580,3.3770423,3.3770423,3.3770423,3.3770423,3.3770423,3.3770423,3.3770423
                 ,3.3770423,3.3770423,3.3770423,3.3770423,3.3770423,3.3770423)
        
    )
    rt_v2_2 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("rt_v2_method2", length(dates[[Country]])),
        "rt" = c(0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
                 ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,7.087556,8.662749,8.662749,8.662749,8.662749,8.662749,8.662749,8.662749,8.530359,8.530359,8.530359,8.530359,8.530359,8.530359,8.530359,9.169386,9.169386,9.169386,9.169386
                 ,9.169386,9.169386,9.169386,9.228455,9.228455,9.228455,9.228455,9.228455,9.228455,9.228455,8.759453,8.759453,8.759453,8.759453,8.759453,8.759453,8.759453,7.119942,7.119942,7.119942,7.119942,7.119942,7.119942,7.119942,5.862922,5.862922,5.862922
                 ,5.862922,5.862922,5.862922,5.862922,5.622718,5.622718,5.622718,5.622718,5.622718,5.622718,5.622718,3.515001,3.515001,3.515001,3.515001,3.515001,3.515001,3.515001,5.624220,5.624220,5.624220,5.624220,5.624220,5.624220,5.624220,7.334476,7.334476
                 ,7.334476,7.334476,7.334476,7.334476,7.334476,6.884059,6.884059,6.884059,6.884059,6.884059,6.884059,6.884059,6.884059,6.884059,6.884059,6.884059,6.884059,6.884059)
    )
    rt_v1_3 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("rt_v1_method3", length(dates[[Country]])),
        "rt" = 
        c(3.3960065,3.3960065,3.3960065,3.3960065,3.3960065,3.3960065,3.3960065,4.4971301,4.4971301,4.4971301,4.4971301,4.4971301,4.4971301,4.4971301,4.5010653,4.5010653,4.5010653,4.5010653,4.5010653,4.5010653,4.5010653,4.5495600,4.5495600,4.5495600
          ,4.5495600,4.5495600,4.5495600,4.5495600,2.0575008,2.0575008,2.0575008,2.0575008,2.0575008,2.0575008,2.0575008,3.0774328,3.0774328,3.0774328,3.0774328,3.0774328,3.0774328,3.0774328,0.5442687,0.5442687,0.5442687,0.5442687,0.5442687,0.5442687
          ,0.5442687,0.6186838,0.6186838,0.6186838,0.6186838,0.6186838,0.6186838,0.6186838,1.9128580,1.9128580,1.9128580,1.9128580,1.9128580,1.9128580,1.9128580,2.6690029,2.6690029,2.6690029,2.6690029,2.6690029,2.6690029,2.6690029,1.7752002,1.7752002
          ,1.7752002,1.7752002,1.7752002,1.7752002,1.7752002,2.0467355,2.0467355,2.0467355,2.0467355,2.0467355,2.0467355,2.0467355,2.6324411,2.6324411,2.6324411,2.6324411,2.6324411,2.6324411,2.6324411,2.4155615,2.4155615,2.4155615,2.4155615,2.4155615
          ,2.4155615,2.4155615,2.1643789,2.1643789,2.1643789,2.1643789,2.1643789,2.1643789,2.1643789,2.3735734,2.3735734,2.3735734,2.3735734,2.3735734,2.3735734,2.3735734,2.6263317,2.6263317,2.6263317,2.6263317,2.6263317,2.6263317,2.6263317,3.1303721
          ,3.1303721,3.1303721,3.1303721,3.1303721,3.1303721,3.1303721,3.4337481,3.4337481,3.4337481,3.4337481,3.4337481,3.4337481,3.4337481,3.3030747,3.3030747,3.3030747,3.3030747,3.3030747,3.3030747,3.3030747,2.8351726,2.8351726,2.8351726,2.8351726
          ,2.8351726,2.8351726,2.8351726,2.7736771,2.7736771,2.7736771,2.7736771,2.7736771,2.7736771,2.7736771,3.4302396,3.4302396,3.4302396,3.4302396,3.4302396,3.4302396,3.4302396,3.9770956,3.9770956,3.9770956,3.9770956,3.9770956,3.9770956,3.9770956
          ,3.6262617,3.6262617,3.6262617,3.6262617,3.6262617,3.6262617,3.6262617,3.0311078,3.0311078,3.0311078,3.0311078,3.0311078,3.0311078,3.0311078,3.2730004,3.2730004,3.2730004,3.2730004,3.2730004,3.2730004,3.2730004,3.7535772,3.7535772,3.7535772
          ,3.7535772,3.7535772,3.7535772,3.7535772,3.8940805,3.8940805,3.8940805,3.8940805,3.8940805,3.8940805,3.8940805,4.0981834,4.0981834,4.0981834,4.0981834,4.0981834,4.0981834,4.0981834,3.7872089,3.7872089,3.7872089,3.7872089,3.7872089,3.7872089
          ,3.7872089,3.4378670,3.4378670,3.4378670,3.4378670,3.4378670,3.4378670,3.4378670,3.7986477,3.7986477,3.7986477,3.7986477,3.7986477,3.7986477,3.7986477,3.8440765,3.8440765,3.8440765,3.8440765,3.8440765,3.8440765,3.8440765,3.6856863,3.6856863
          ,3.6856863,3.6856863,3.6856863,3.6856863,3.6856863,3.4831835,3.4831835,3.4831835,3.4831835,3.4831835,3.4831835,3.4831835,4.6444607,4.6444607,4.6444607,4.6444607,4.6444607,4.6444607,4.6444607,4.2965014,4.2965014,4.2965014,4.2965014,4.2965014
          ,4.2965014,4.2965014,4.8150719,4.8150719,4.8150719,4.8150719,4.8150719,4.8150719,4.8150719,4.9495381,4.9495381,4.9495381,4.9495381,4.9495381,4.9495381,4.9495381,4.7772613,4.7772613,4.7772613,4.7772613,4.7772613,4.7772613,4.7772613,3.7314170
          ,3.7314170,3.7314170,3.7314170,3.7314170,3.7314170,3.7314170,2.7739825,2.7739825,2.7739825,2.7739825,2.7739825,2.7739825,2.7739825,3.5141112,3.5141112,3.5141112,3.5141112,3.5141112,3.5141112,3.5141112,1.1874394,1.1874394,1.1874394,1.1874394
          ,1.1874394,1.1874394,1.1874394,2.9490029,2.9490029,2.9490029,2.9490029,2.9490029,2.9490029,2.9490029,3.9566905,3.9566905,3.9566905,3.9566905,3.9566905,3.9566905,3.9566905,3.2829147,3.2829147,3.2829147,3.2829147,3.2829147,3.2829147,3.2829147
          ,3.2667439,3.2667439,3.2667439,3.2667439,3.2667439,3.2667439)
    )
    rt_v2_3 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("rt_v2_method3", length(dates[[Country]])),
        "rt" = 
        c(0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
          ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
          ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
          ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
          ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
          ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
          ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
          ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
          ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
          ,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,6.302925,8.305559,8.305559,8.305559,8.305559,8.305559,8.305559,8.305559,7.707198,7.707198,7.707198,7.707198,7.707198,7.707198,7.707198,8.639674,8.639674,8.639674,8.639674
          ,8.639674,8.639674,8.639674,8.820092,8.820092,8.820092,8.820092,8.820092,8.820092,8.820092,8.486620,8.486620,8.486620,8.486620,8.486620,8.486620,8.486620,6.512778,6.512778,6.512778,6.512778,6.512778,6.512778,6.512778,4.957304,4.957304,4.957304
          ,4.957304,4.957304,4.957304,4.957304,6.314893,6.314893,6.314893,6.314893,6.314893,6.314893,6.314893,2.159420,2.159420,2.159420,2.159420,2.159420,2.159420,2.159420,5.403366,5.403366,5.403366,5.403366,5.403366,5.403366,5.403366,7.209335,7.209335
          ,7.209335,7.209335,7.209335,7.209335,7.209335,5.995512,5.995512,5.995512,5.995512,5.995512,5.995512,5.995512,5.941919,5.941919,5.941919,5.941919,5.941919,5.941919)
        
    )
    rt_df = rbind(rt_v1_1,rt_v1_2,rt_v1_3,rt_v2_1,rt_v2_2,rt_v2_3)
    ggplot() +
        geom_line(data=rt_df,aes(x = time, y = rt, color = type)) +
        scale_color_manual(name = "", labels = c("Non-P.1 (First scenario)","Non-P.1 (Second scenario)","Non-P.1 (Third scenario)","P.1 (First scenario)","P.1 (Second scenario)","P.1 (Third scenario)"),
                           values = c("red","blue","orange","chocolate","gold","green"))+
        guides(color=guide_legend(ncol=2)) +
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.25, 0.9)) +
        xlab("Date") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab(expression(Reproduction ~ number * "," ~ R[0]))+
        scale_x_date(date_breaks = "1 month", labels = date_format("%b/%y")) +
        theme(legend.text = element_text(size = 11))+
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))
}

p_F <- plt_rt(out)
p_F
ggsave(here(paste0("transmission_model/figures/fig_R0_Comparing_contour.pdf")),plot=p_F)








#################################Plot of Rt
#MMM<-read_excel("C://Users/林立新/Desktop/Rt_adj.xlsx",sheet=1,na="NA")
plt_rt <- function(outcheck)
{
    rt_v1_adj_3 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("rt_v1_adj_method3", length(dates[[Country]])),
        "rt"  = c(3.3960065,3.3960065,3.3960065,3.3960065,3.3960065,3.3960065,3.3952740,4.4960408,4.4956302,4.4951207,4.4944511,4.4935500,4.4923294,4.4906711,4.4941264,4.4919617,4.4891280,4.4853469,4.4802633,4.4734079,4.4641414,4.5053846,4.4911939,4.4727193
                  ,4.4481239,4.4151483,4.3709224,4.3117158,1.9318030,1.9126967,1.8894469,1.8620273,1.8303407,1.7941004,1.7528894,2.6038846,2.5158559,2.4238121,2.3221219,2.2074125,2.0785712,1.9360727,0.3317025,0.3247905,0.3177253,0.3107586,0.3040847,0.2978126
                  ,0.2919883,0.3362995,0.3307423,0.3260946,0.3220801,0.3185035,0.3152480,0.3122475,0.9812309,0.9669131,0.9548678,0.9436672,0.9324131,0.9206680,0.9082427,1.2823989,1.2666573,1.2519363,1.2370519,1.2212209,1.2040332,1.1852727,0.7913278,0.7838673
                  ,0.7767763,0.7697402,0.7625536,0.7551116,0.7473617,0.8660837,0.8587949,0.8521863,0.8457936,0.8393030,0.8325499,0.8254525,1.0712817,1.0634308,1.0562189,1.0491192,1.0417757,1.0339992,1.0256926,0.9423936,0.9367242,0.9313899,0.9260902,0.9206203
                  ,0.9148701,0.9087810,0.8203229,0.8166212,0.8132358,0.8099567,0.8066397,0.8032070,0.7996179,0.8766154,0.8735574,0.8708150,0.8681704,0.8654726,0.8626391,0.8596257,0.9570773,0.9546848,0.9525347,0.9504450,0.9482888,0.9459959,0.9435272,1.1264880
                  ,1.1239445,1.1215877,1.1192231,1.1167126,1.1139772,1.1109696,1.2194118,1.2165467,1.2137578,1.2108638,1.2077326,1.2042845,1.2004682,1.1615432,1.1586375,1.1557591,1.1527498,1.1494937,1.1459195,1.1419800,0.9797365,0.9775911,0.9755474,0.9734771
                  ,0.9712893,0.9689320,0.9663737,0.9445431,0.9428577,0.9413385,0.9398437,0.9382730,0.9365699,0.9347021,1.1649683,1.1629793,1.1611710,1.1593459,1.1573619,1.1551368,1.1526217,1.3400424,1.3372018,1.3343865,1.3313942,1.3280755,1.3243385,1.3201224
                  ,1.2031930,1.2000459,1.1968026,1.1933293,1.1895271,1.1853326,1.1807015,0.9863121,0.9838557,0.9814557,0.9789936,0.9763852,0.9735817,0.9705534,1.0535190,1.0509021,1.0484705,1.0460288,1.0434377,1.0406178,1.0375237,1.1977374,1.1945567,1.1915203
                  ,1.1883972,1.1850209,1.1812952,1.1771634,1.2164500,1.2126495,1.2088622,1.2048989,1.2006232,1.1959550,1.1908461,1.2580655,1.2527849,1.2474481,1.2418452,1.2358245,1.2292975,1.2222115,1.1253135,1.1197466,1.1141182,1.1082579,1.1020456,1.0954135
                  ,1.0883232,0.9938211,0.9890724,0.9844889,0.9798758,0.9750980,0.9700827,0.9647934,1.0678536,1.0627155,1.0578879,1.0531022,1.0481716,1.0429971,1.0375305,1.0517669,1.0469811,1.0423592,1.0377030,1.0328743,1.0277984,1.0224382,0.9838143,0.9798514
                  ,0.9760658,0.9722861,0.9683928,0.9643217,0.9600409,0.9050557,0.9025257,0.9002425,0.8980501,0.8958402,0.8935542,0.8911625,1.2011420,1.1970145,1.1931841,1.1893613,1.1853404,1.1810086,1.1763086,1.0757119,1.0722935,1.0688115,1.0651891,1.0613736
                  ,1.0573354,1.0530584,1.1870437,1.1818080,1.1766168,1.1712923,1.1657073,1.1597915,1.1535078,1.1864695,1.1797425,1.1728356,1.1656261,1.1580193,1.1499496,1.1413637,1.0962330,1.0881312,1.0796504,1.0706352,1.0609270,1.0503550,1.0387123,0.8073188
                  ,0.8002002,0.7923518,0.7835495,0.7735313,0.7619748,0.7484793,0.5584370,0.5508666,0.5424528,0.5327736,0.5214714,0.5082321,0.4927574,0.6211618,0.5954128,0.5684468,0.5384127,0.5043015,0.4659320,0.4237496,0.1396001,0.1339105,0.1283003,0.1228261
                  ,0.1175644,0.1125876,0.1079542,0.2677253,0.2458251,0.2286860,0.2138472,0.2001252,0.1871785,0.1750370,0.2648848,0.2484975,0.2346492,0.2224493,0.2114227,0.2013987,0.1923468,0.1698600,0.1647425,0.1605740,0.1571353,0.1542898,0.1519654,0.1501232
                  ,0.1721154,0.1709108,0.1703188,0.1701512,0.1702921,0.1706867)
    )
    rt_v2_adj_3 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("rt_v2_adj_method3", length(dates[[Country]])),
        "rt" = c(0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,4.4012036,4.3965145,4.3921614,4.3878171,4.3832492,4.3783301,4.3729950,4.0192133,4.0152101,4.0111329,4.0068886,4.0024118
                 ,3.9976657,3.9926289,4.5220117,4.5160153,4.5100459,4.5038988,4.4974238,4.4905316,4.4831659,4.6118288,4.6041412,4.5961491,4.5876776,4.5785682,4.5686737,4.5578321,4.3867966,4.3759963,4.3640893,4.3506360,4.3351177,4.3168939,4.2951347,3.3215736
                 ,3.3062113,3.2877845,3.2654478,3.2381679,3.2046376,3.1632113,2.3562882,2.3284063,2.2964568,2.2588802,2.2143123,2.1614962,2.0991873,2.6477491,2.5480610,2.4423600,2.3241142,2.1897558,2.0386551,1.8724328,0.6134278,0.5881681,0.5631583,0.5385885
                 ,0.5147396,0.4918879,0.4702628,1.1801426,1.0896686,1.0165729,0.9518125,0.8910614,0.8330924,0.7780426,1.0988503,1.0218521,0.9546971,0.8938254,0.8373624,0.7846548,0.7355995,0.6315143,0.5995280,0.5711667,0.5455498,0.5221315,0.5006313,0.4809046
                 ,0.5394509,0.5219469,0.5067810,0.4932070,0.4807719,0.4692680)
    )
    rt_v1_adj_2 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("rt_v1_adj_method2", length(dates[[Country]])),
        "rt" = c(3.3727426,3.3727426,3.3727426,3.3727426,3.3727426,3.3727426,3.3717011,3.9393348,3.9388597,3.9382755,3.9375323,3.9365722,3.9353256,3.9337030,3.9113563,3.9091039,3.9062335,3.9025439,3.8977830,3.8916274,3.8836548,3.7950738,3.7836040,3.7690703
                 ,3.7505073,3.7267380,3.6963039,3.6573802,2.5411797,2.5129952,2.4787843,2.4381014,2.3903856,2.3349021,2.2708341,1.9037739,1.8484335,1.7878571,1.7219187,1.6507513,1.5747218,1.4944180,0.4958716,0.4853055,0.4746157,0.4641728,0.4542492,0.4449881
                 ,0.4364409,0.4370394,0.4301661,0.4241122,0.4187222,0.4138758,0.4094893,0.4055041,0.6871766,0.6788552,0.6715485,0.6649103,0.6587021,0.6527913,0.6471115,0.9578225,0.9480404,0.9389759,0.9302575,0.9216322,0.9129697,0.9042131,0.8869827,0.8797304
                 ,0.8727745,0.8659651,0.8592009,0.8524298,0.8456283,0.8816001,0.8753162,0.8693367,0.8635211,0.8577722,0.8520388,0.8462973,0.9129490,0.9072780,0.9018627,0.8965726,0.8913163,0.8860451,0.8807357,0.8587400,0.8544422,0.8503384,0.8463393,0.8423825
                 ,0.8384342,0.8344772,0.7881170,0.7852212,0.7824884,0.7798626,0.7773036,0.7747894,0.7723079,0.8394433,0.8370704,0.8348696,0.8327666,0.8307079,0.8286638,0.8266186,0.9346936,0.9326190,0.9306773,0.9287973,0.9269280,0.9250407,0.9231199,1.0235074
                 ,1.0215241,1.0196296,1.0177554,1.0158514,1.0138890,1.0118520,1.1117225,1.1095535,1.1074379,1.1053027,1.1030943,1.1007815,1.0983462,1.0870214,1.0848862,1.0827555,1.0805809,1.0783269,1.0759725,1.0735043,0.9882574,0.9865242,0.9848126,0.9830883
                 ,0.9813269,0.9795139,0.9776402,0.9638498,0.9623971,0.9610088,0.9596381,0.9582511,0.9568281,0.9553580,1.0699638,1.0683825,1.0668821,1.0653878,1.0638450,1.0622221,1.0605014,1.1991888,1.1970626,1.1949489,1.1927634,1.1904433,1.1879511,1.1852637
                 ,1.1414572,1.1391486,1.1367811,1.1343092,1.1316991,1.1289299,1.1259874,1.0226586,1.0206025,1.0185367,1.0164270,1.0142487,1.0119874,1.0096341,1.0192920,1.0172730,1.0153140,1.0133546,1.0113510,1.0092785,1.0071233,1.1413861,1.1387543,1.1361811
                 ,1.1335664,1.1308371,1.1279512,1.1248855,1.1677749,1.1645582,1.1613008,1.1579254,1.1543753,1.1506174,1.1466324,1.1956699,1.1913854,1.1870226,1.1824950,1.1777394,1.1727194,1.1674145,1.0994412,1.0949882,1.0904609,1.0858127,1.0810119,1.0760420
                 ,1.0708951,1.0134886,1.0092783,1.0050993,1.0009019,0.9966524,0.9923342,0.9879403,1.0057963,1.0017006,0.9977140,0.9937595,0.9897827,0.9857553,0.9816643,0.9867549,0.9830638,0.9794875,0.9759510,0.9724013,0.9688105,0.9651653,0.9597671,0.9567470
                 ,0.9538458,0.9510009,0.9481681,0.9453242,0.9424579,0.9308227,0.9285944,0.9264975,0.9244768,0.9224935,0.9205263,0.9185647,1.1337328,1.1303555,1.1271341,1.1239231,1.1206172,1.1171583,1.1135184,1.0930307,1.0896919,1.0862937,1.0828000,1.0791858
                 ,1.0754381,1.0715506,1.1508753,1.1461966,1.1414766,1.1366373,1.1316214,1.1263968,1.1209461,1.1291153,1.1233597,1.1174527,1.1113362,1.1049614,1.0982893,1.0912819,1.0373643,1.0309076,1.0241369,1.0169606,1.0092715,1.0009396,0.9917991,0.8096757
                 ,0.8033620,0.7962615,0.7882313,0.7790778,0.7685393,0.7562839,0.6101585,0.6006985,0.5897656,0.5771319,0.5625676,0.5458350,0.5267047,0.4863050,0.4666825,0.4454498,0.4223464,0.3973260,0.3705547,0.3423713,0.1983250,0.1873124,0.1764429,0.1659068
                 ,0.1558761,0.1464815,0.1378128,0.2180784,0.1987247,0.1824388,0.1681570,0.1553351,0.1437851,0.1334707,0.1831993,0.1686963,0.1565675,0.1461871,0.1371827,0.1293712,0.1226603,0.1249161,0.1206974,0.1173947,0.1148547,0.1129683,0.1116634,0.1108886
                 ,0.1106026,0.1107697,0.1113569,0.1123338,0.1136719,0.1153446)
        
    )
    rt_v2_adj_2 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("rt_v2_adj_method2", length(dates[[Country]])),
        "rt" = c(0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,4.2205093,4.2156821,4.2110676,4.2064668,4.2017354,4.1967945,4.1916055,4.1230916,4.1183272,4.1134801,4.1085001,4.1033520
                 ,4.0980168,4.0924842,4.4063827,4.3997635,4.3930799,4.3862172,4.3790886,4.3716398,4.3638346,4.4054830,4.3973315,4.3888787,4.3800129,4.3706239,4.3605983,4.3498034,4.1397034,4.1292777,4.1178570,4.1051332,4.0907247,4.0741467,4.0547727,3.3132622
                 ,3.2977422,3.2791548,3.2568930,3.2301664,3.1979426,3.1589324,2.5503408,2.5167608,2.4772386,2.4309192,2.3769324,2.3143731,2.2423626,2.0758761,2.0014105,1.9204987,1.8322762,1.7366311,1.6342036,1.5262529,0.8971137,0.8532468,0.8097970,0.7674581
                 ,0.7268576,0.6884730,0.6526347,1.0217202,0.9440666,0.8774246,0.8179541,0.7636994,0.7139945,0.6687189,0.8850439,0.8216926,0.7669214,0.7184735,0.6750086,0.6358679,0.6007172,0.5782343,0.5523797,0.5298518,0.5100956,0.4927146,0.4774445,0.4640948
                 ,0.4525125,0.4425634,0.4341251,0.4270836,0.4213328,0.4167738)
    )
    rt_v1_adj_1 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("rt_v1_adj_method1", length(dates[[Country]])),
        "rt" = c(3.4562913,3.4562913,3.4562913,3.4562913,3.4562913,3.4562913,3.4549510,3.6631365,3.6626116,3.6619683,3.6611668,3.6601610,3.6588955,3.6573016,3.6417490,3.6393196,3.6362679,3.6324312,3.6276047,3.6215307,3.6138838,3.4201396,3.4095135,3.3962498
                 ,3.3797722,3.3593633,3.3341359,3.3030138,2.6776073,2.6467623,2.6095354,2.5655299,2.5142089,2.4548755,2.3867876,1.5249308,1.4875551,1.4459838,1.4013209,1.3544862,1.3061116,1.2566638,0.6577109,0.6429193,0.6280854,0.6136333,0.5998673,0.5869447
                 ,0.5749215,0.4987848,0.4909078,0.4836963,0.4771157,0.4711277,0.4656910,0.4607629,0.6069711,0.6000475,0.5938316,0.5881892,0.5830177,0.5782483,0.5738352,0.7726107,0.7658249,0.7595642,0.7537068,0.7481640,0.7428830,0.7378346,0.8775767,0.8711592
                 ,0.8650525,0.8591844,0.8535027,0.8479782,0.8425972,0.9015772,0.8958380,0.8902832,0.8848801,0.8796051,0.8744457,0.8693966,0.8685161,0.8638189,0.8592605,0.8548289,0.8505156,0.8463159,0.8422273,0.8360095,0.8323041,0.8287318,0.8252825,0.8219486
                 ,0.8187255,0.8156099,0.8063185,0.8035712,0.8009437,0.7984276,0.7960165,0.7937057,0.7914920,0.8273382,0.8251800,0.8231331,0.8211841,0.8193226,0.8175421,0.8158385,0.8980129,0.8961248,0.8943329,0.8926177,0.8909645,0.8893649,0.8878143,0.9804328
                 ,0.9786220,0.9768779,0.9751804,0.9735150,0.9718733,0.9702513,1.0316790,1.0298951,1.0281428,1.0264080,1.0246800,1.0229530,1.0212241,1.0294439,1.0277355,1.0260381,1.0243453,1.0226527,1.0209576,1.0192589,1.0084927,1.0069347,1.0053888,1.0038523
                 ,1.0023232,1.0008005,0.9992837,1.0176825,1.0161689,1.0146755,1.0131946,1.0117207,1.0102507,1.0087830,1.0482680,1.0466997,1.0451472,1.0435989,1.0420462,1.0404844,1.0389112,1.0859910,1.0842243,1.0824550,1.0806693,1.0788566,1.0770108,1.0751290
                 ,1.0702570,1.0684156,1.0665546,1.0646695,1.0627570,1.0608155,1.0588441,1.0485740,1.0466839,1.0447833,1.0428682,1.0409363,1.0389864,1.0370180,1.0637242,1.0616436,1.0595612,1.0574654,1.0553477,1.0532035,1.0510307,1.1008824,1.0983908,1.0958835
                 ,1.0933413,1.0907501,1.0881023,1.0853944,1.1274069,1.1243648,1.1212776,1.1181273,1.1149009,1.1115915,1.1081963,1.1374316,1.1336715,1.1298412,1.1259244,1.1219095,1.1177912,1.1135677,1.0955116,1.0913361,1.0870875,1.0827655,1.0783714,1.0739076
                 ,1.0693771,1.0291903,1.0250991,1.0209844,1.0168555,1.0127208,1.0085870,1.0044591,0.9772956,0.9735983,0.9699423,0.9663310,0.9627681,0.9592566,0.9557993,0.9427939,0.9396873,0.9366645,0.9337216,0.9308566,0.9280685,0.9253578,0.9401968,0.9375729
                 ,0.9350521,0.9326209,0.9302696,0.9279929,0.9257886,0.9808954,0.9784430,0.9760876,0.9738046,0.9715762,0.9693928,0.9672506,1.0687917,1.0658904,1.0630472,1.0602197,1.0573768,1.0545020,1.0515885,1.1239812,1.1202892,1.1165660,1.1127767,1.1088960
                 ,1.1049110,1.1008170,1.1386208,1.1339046,1.1290899,1.1241537,1.1190796,1.1138586,1.1084863,1.0973974,1.0919277,1.0863150,1.0805505,1.0746225,1.0685145,1.0622037,0.9819699,0.9765283,0.9708340,0.9648660,0.9585786,0.9518946,0.9447034,0.8071179
                 ,0.8014547,0.7950742,0.7878915,0.7797691,0.7705050,0.7598390,0.6104623,0.6017783,0.5915707,0.5798005,0.5663868,0.5511997,0.5340949,0.4158654,0.4024037,0.3875662,0.3715979,0.3547304,0.3371665,0.3191052,0.2591527,0.2445749,0.2301777,0.2161755
                 ,0.2027602,0.1900865,0.1782703,0.1856615,0.1724911,0.1608106,0.1504545,0.1413067,0.1332884,0.1263330,0.1410682,0.1339009,0.1279293,0.1229696,0.1188817,0.1155631,0.1129345,0.1250173,0.1230050,0.1215864,0.1206795,0.1202191,0.1201559,0.1204511
                 ,0.1210723,0.1219914,0.1231830,0.1246244,0.1262947,0.1281750)
    )
    rt_v2_adj_1 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("rt_v2_adj_method1", length(dates[[Country]])),
        "rt" = c(0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000
                 ,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,3.9928225,3.9886136,3.9844897,3.9803904,3.9762718,3.9721104,3.9678965,4.2458908,4.2406107,4.2352885,4.2298756,4.2243365
                 ,4.2186526,4.2128159,4.3663613,4.3596917,4.3528766,4.3458808,4.3386759,4.3312415,4.3235604,4.2915008,4.2836725,4.2755569,4.2671138,4.2582884,4.2490060,4.2391664,3.9259407,3.9169678,3.9071356,3.8962881,3.8841934,3.8705218,3.8548343,3.2994285
                 ,3.2848546,3.2674357,3.2467466,3.2221909,3.1929629,3.1580605,2.5413162,2.5098219,2.4723436,2.4286571,2.3784022,2.3210594,2.2560676,1.7656203,1.7131128,1.6552292,1.5928611,1.5268576,1.4579678,1.3869315,1.1437841,1.0866163,1.0298981,0.9744391
                 ,0.9209578,0.8700290,0.8220829,0.8664774,0.8129800,0.7647443,0.7211866,0.6818997,0.6466084,0.6150814,0.6777016,0.6442035,0.6149912,0.5893968,0.5669139,0.5471751,0.5298975,0.5616747,0.5458657,0.5322025,0.5203806,0.5101601,0.5013608,0.4938417
                 ,0.4874863,0.4821933,0.4778722,0.4744409,0.4718247,0.4699554)
        
    )
    
    rt_df = rbind(rt_v1_adj_1,rt_v1_adj_2,rt_v1_adj_3,rt_v2_adj_1,rt_v2_adj_2,rt_v2_adj_3)
    tail(rt_df,1)
    ggplot() +
        geom_line(data=rt_df,aes(x = time, y = rt, color = type)) +
        scale_color_manual(name = "", labels = c("Non-P.1 (First scenario)","Non-P.1 (Second scenario)","Non-P.1 (Third scenario)","P.1 (First scenario)","P.1 (Second scenario)","P.1 (Third scenario)"),
                           values = c("red","blue","orange","chocolate","gold","green"))+  
        guides(color=guide_legend(ncol=2)) +
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.25, 0.9)) +
        xlab("Date") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab(expression(Reproduction ~ number * "," ~ R[t]))+
        scale_x_date(date_breaks = "1 month", labels = date_format("%b/%y")) +
        theme(legend.text = element_text(size = 11))+
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))
}

p_F <- plt_rt(out)
p_F
ggsave(here(paste0("transmission_model/figures/fig_Rt_Comparing_contour.pdf")),plot=p_F)







#############################cummulative infection rate
#MMM<-read_excel("C://Users/林立新/Desktop/cum_IR.xlsx",sheet=1,na="NA")
plt_seroconv <- function(outcheck)
{
    seroconv_v1_1 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("seroconv_v1_method1", length(dates[[Country]])),
        "seroconv"  = c(0
                        ,
                        0.0000004585221
                        ,
                        0.0000015994790
                        ,
                        0.0000036863860
                        ,
                        0.0000070773530
                        ,
                        0.0000122048100
                        ,
                        0.0000195958000
                        ,
                        0.0000301019900
                        ,
                        0.0000445924200
                        ,
                        0.0000641312900
                        ,
                        0.0000900100900
                        ,
                        0.000123721
                        ,
                        0.000167104
                        ,
                        0.000222356
                        ,
                        0.000292193
                        ,
                        0.000379941
                        ,
                        0.000489817
                        ,
                        0.00062711
                        ,
                        0.000798479
                        ,
                        0.001012291
                        ,
                        0.001279126
                        ,
                        0.001612323
                        ,
                        0.002027204
                        ,
                        0.00254357
                        ,
                        0.003185799
                        ,
                        0.00398374
                        ,
                        0.004974031
                        ,
                        0.006201412
                        ,
                        0.007720749
                        ,
                        0.009581725
                        ,
                        0.01185236
                        ,
                        0.01461098
                        ,
                        0.01794526
                        ,
                        0.02195335
                        ,
                        0.02674311
                        ,
                        0.03243592
                        ,
                        0.03909231
                        ,
                        0.04680972
                        ,
                        0.05567978
                        ,
                        0.06577864
                        ,
                        0.0771642
                        ,
                        0.08986804
                        ,
                        0.1039089
                        ,
                        0.1191996
                        ,
                        0.1356637
                        ,
                        0.1531908
                        ,
                        0.1716388
                        ,
                        0.1908342
                        ,
                        0.2105708
                        ,
                        0.230636
                        ,
                        0.250804
                        ,
                        0.2708451
                        ,
                        0.2905366
                        ,
                        0.3096819
                        ,
                        0.3281114
                        ,
                        0.3456809
                        ,
                        0.3622758
                        ,
                        0.3778312
                        ,
                        0.392309
                        ,
                        0.4057047
                        ,
                        0.4180424
                        ,
                        0.4293691
                        ,
                        0.4397464
                        ,
                        0.4492394
                        ,
                        0.4579303
                        ,
                        0.4658968
                        ,
                        0.4732178
                        ,
                        0.4799684
                        ,
                        0.4862162
                        ,
                        0.4920282
                        ,
                        0.497461
                        ,
                        0.5025674
                        ,
                        0.5073929
                        ,
                        0.5119734
                        ,
                        0.5163401
                        ,
                        0.5205194
                        ,
                        0.5245354
                        ,
                        0.5284087
                        ,
                        0.5321565
                        ,
                        0.5357935
                        ,
                        0.5393311
                        ,
                        0.5427787
                        ,
                        0.5461439
                        ,
                        0.5494334
                        ,
                        0.552653
                        ,
                        0.5558068
                        ,
                        0.5588981
                        ,
                        0.5619299
                        ,
                        0.5649038
                        ,
                        0.5678213
                        ,
                        0.5706829
                        ,
                        0.573489
                        ,
                        0.5762397
                        ,
                        0.5789346
                        ,
                        0.5815734
                        ,
                        0.5841558
                        ,
                        0.5866812
                        ,
                        0.5891492
                        ,
                        0.5915594
                        ,
                        0.5939112
                        ,
                        0.5962045
                        ,
                        0.5984391
                        ,
                        0.6006152
                        ,
                        0.602733
                        ,
                        0.6047927
                        ,
                        0.6067946
                        ,
                        0.6087395
                        ,
                        0.6106282
                        ,
                        0.6124615
                        ,
                        0.6142405
                        ,
                        0.6159665
                        ,
                        0.6176408
                        ,
                        0.6192651
                        ,
                        0.6208415
                        ,
                        0.6223721
                        ,
                        0.6238591
                        ,
                        0.6253048
                        ,
                        0.6267115
                        ,
                        0.6280819
                        ,
                        0.6294181
                        ,
                        0.6307233
                        ,
                        0.632
                        ,
                        0.6332508
                        ,
                        0.6344783
                        ,
                        0.6356848
                        ,
                        0.6368728
                        ,
                        0.6380446
                        ,
                        0.6392025
                        ,
                        0.6403488
                        ,
                        0.6414855
                        ,
                        0.6426143
                        ,
                        0.643737
                        ,
                        0.6448552
                        ,
                        0.6459704
                        ,
                        0.6470838
                        ,
                        0.6481965
                        ,
                        0.6493097
                        ,
                        0.6504239
                        ,
                        0.65154
                        ,
                        0.6526583
                        ,
                        0.6537795
                        ,
                        0.6549034
                        ,
                        0.6560304
                        ,
                        0.6571603
                        ,
                        0.6582931
                        ,
                        0.6594285
                        ,
                        0.6605663
                        ,
                        0.6617061
                        ,
                        0.6628478
                        ,
                        0.6639909
                        ,
                        0.6651353
                        ,
                        0.6662806
                        ,
                        0.6674266
                        ,
                        0.6685733
                        ,
                        0.6697205
                        ,
                        0.6708684
                        ,
                        0.6720171
                        ,
                        0.673167
                        ,
                        0.6743184
                        ,
                        0.6754716
                        ,
                        0.6766272
                        ,
                        0.6777856
                        ,
                        0.6789477
                        ,
                        0.6801141
                        ,
                        0.6812856
                        ,
                        0.6824628
                        ,
                        0.6836467
                        ,
                        0.684838
                        ,
                        0.6860376
                        ,
                        0.687246
                        ,
                        0.6884642
                        ,
                        0.6896926
                        ,
                        0.6909318
                        ,
                        0.6921824
                        ,
                        0.6934445
                        ,
                        0.6947186
                        ,
                        0.6960045
                        ,
                        0.6973023
                        ,
                        0.6986117
                        ,
                        0.6999324
                        ,
                        0.7012641
                        ,
                        0.7026065
                        ,
                        0.703959
                        ,
                        0.7053215
                        ,
                        0.7066936
                        ,
                        0.7080751
                        ,
                        0.7094656
                        ,
                        0.7108652
                        ,
                        0.7122736
                        ,
                        0.7136911
                        ,
                        0.7151181
                        ,
                        0.7165549
                        ,
                        0.7180022
                        ,
                        0.7194605
                        ,
                        0.7209305
                        ,
                        0.7224133
                        ,
                        0.7239095
                        ,
                        0.7254204
                        ,
                        0.7269471
                        ,
                        0.7284907
                        ,
                        0.7300523
                        ,
                        0.731633
                        ,
                        0.733234
                        ,
                        0.7348563
                        ,
                        0.7365012
                        ,
                        0.7381698
                        ,
                        0.739863
                        ,
                        0.7415818
                        ,
                        0.743327
                        ,
                        0.7450997
                        ,
                        0.7469004
                        ,
                        0.7487298
                        ,
                        0.7505881
                        ,
                        0.7524758
                        ,
                        0.754393
                        ,
                        0.7563395
                        ,
                        0.7583153
                        ,
                        0.7603198
                        ,
                        0.7623518
                        ,
                        0.7644102
                        ,
                        0.7664936
                        ,
                        0.7686006
                        ,
                        0.7707292
                        ,
                        0.7728776
                        ,
                        0.7750435
                        ,
                        0.7772245
                        ,
                        0.779418
                        ,
                        0.7816215
                        ,
                        0.7838325
                        ,
                        0.7860485
                        ,
                        0.7882668
                        ,
                        0.7904852
                        ,
                        0.792701
                        ,
                        0.794912
                        ,
                        0.7971158
                        ,
                        0.7993106
                        ,
                        0.8014945
                        ,
                        0.8036655
                        ,
                        0.8058222
                        ,
                        0.8079632
                        ,
                        0.8100872
                        ,
                        0.8121931
                        ,
                        0.8142801
                        ,
                        0.8163474
                        ,
                        0.8183947
                        ,
                        0.8204215
                        ,
                        0.8224282
                        ,
                        0.824415
                        ,
                        0.8263824
                        ,
                        0.828331
                        ,
                        0.8302616
                        ,
                        0.8321753
                        ,
                        0.834073
                        ,
                        0.835957
                        ,
                        0.8378291
                        ,
                        0.8396912
                        ,
                        0.8415454
                        ,
                        0.8433939
                        ,
                        0.8452393
                        ,
                        0.8470837
                        ,
                        0.8489305
                        ,
                        0.8507825
                        ,
                        0.8526425
                        ,
                        0.8545131
                        ,
                        0.8563972
                        ,
                        0.8582973
                        ,
                        0.860216
                        ,
                        0.8621558
                        ,
                        0.8641192
                        ,
                        0.8661081
                        ,
                        0.8681244
                        ,
                        0.8701697
                        ,
                        0.8722456
                        ,
                        0.8743535
                        ,
                        0.8764941
                        ,
                        0.8786683
                        ,
                        0.8808766
                        ,
                        0.8831192
                        ,
                        0.885396
                        ,
                        0.8877066
                        ,
                        0.8900505
                        ,
                        0.8924255
                        ,
                        0.8948299
                        ,
                        0.8972613
                        ,
                        0.8997173
                        ,
                        0.9021947
                        ,
                        0.9046898
                        ,
                        0.907199
                        ,
                        0.9097159
                        ,
                        0.9122352
                        ,
                        0.9147511
                        ,
                        0.9172575
                        ,
                        0.9197477
                        ,
                        0.9222151
                        ,
                        0.9246527
                        ,
                        0.9270521
                        ,
                        0.9294057
                        ,
                        0.931706
                        ,
                        0.9339459
                        ,
                        0.9361185
                        ,
                        0.9382171
                        ,
                        0.9402358
                        ,
                        0.9421683
                        ,
                        0.9440096
                        ,
                        0.9457555
                        ,
                        0.9474027
                        ,
                        0.9489488
                        ,
                        0.9503923
                        ,
                        0.9517325
                        ,
                        0.9529696
                        ,
                        0.9541049
                        ,
                        0.9551405
                        ,
                        0.9560795
                        ,
                        0.9569256
                        ,
                        0.9576834
                        ,
                        0.9583576
                        ,
                        0.9589537
                        ,
                        0.9594773
                        ,
                        0.9599343
                        ,
                        0.9603307
                        ,
                        0.9606725
                        ,
                        0.9609655
                        ,
                        0.9612152
                        ,
                        0.9614269
                        ,
                        0.9616055
                        ,
                        0.9617556
                        ,
                        0.9618812
                        ,
                        0.9619859
                        ,
                        0.9620731
                        ,
                        0.9621457
                        ,
                        0.962206
                        ,
                        0.9622561
                        ,
                        0.9622979
                        ,
                        0.9623327
                        ,
                        0.9623618
                        ,
                        0.9623863
                        ,
                        0.9624069
                        ,
                        0.9624244
                        ,
                        0.9624394
                        ,
                        0.9624522
                        ,
                        0.9624633
                        ,
                        0.9624729)
    )
    seroconv_v2_1 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("seroconv_v2_method1", length(dates[[Country]])),
        "seroconv" = c(0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0.0000000073788
                       ,
                       0.0000000188941
                       ,
                       0.0000000368844
                       ,
                       0.0000000651520
                       ,
                       0.0000001079032
                       ,
                       0.0000001708470
                       ,
                       0.0000002592599
                       ,
                       0.0000003809913
                       ,
                       0.0000005475734
                       ,
                       0.0000007719053
                       ,
                       0.0000010694000
                       ,
                       0.0000014610730
                       ,
                       0.0000019738130
                       ,
                       0.0000026427240
                       ,
                       0.0000035129940
                       ,
                       0.0000046452190
                       ,
                       0.0000061182000
                       ,
                       0.0000080352050
                       ,
                       0.0000105311900
                       ,
                       0.0000137841700
                       ,
                       0.0000180267500
                       ,
                       0.0000235621500
                       ,
                       0.0000307785000
                       ,
                       0.0000401857600
                       ,
                       0.0000524481200
                       ,
                       0.0000684286100
                       ,
                       0.0000892480300
                       ,
                       0.000116364
                       ,
                       0.000151672
                       ,
                       0.000197362
                       ,
                       0.000256399
                       ,
                       0.000332553
                       ,
                       0.000430588
                       ,
                       0.000556542
                       ,
                       0.000718056
                       ,
                       0.000924861
                       ,
                       0.001186819
                       ,
                       0.001517515
                       ,
                       0.001933502
                       ,
                       0.002454667
                       ,
                       0.003105002
                       ,
                       0.003913257
                       ,
                       0.004914428
                       ,
                       0.006139182
                       ,
                       0.007629286
                       ,
                       0.009432084
                       ,
                       0.0115999
                       ,
                       0.01419043
                       ,
                       0.01726605
                       ,
                       0.02089617
                       ,
                       0.02513087
                       ,
                       0.03003372
                       ,
                       0.03566686
                       ,
                       0.04208788
                       ,
                       0.04934768
                       ,
                       0.05748667
                       ,
                       0.06653819
                       ,
                       0.07650553
                       ,
                       0.08738531
                       ,
                       0.09915684
                       ,
                       0.1117825
                       ,
                       0.1252062
                       ,
                       0.1393522
                       ,
                       0.1541301
                       ,
                       0.1694472
                       ,
                       0.1851903
                       ,
                       0.2012388
                       ,
                       0.2174712
                       ,
                       0.2337668
                       ,
                       0.2500079
                       ,
                       0.2660793
                       ,
                       0.2818836
                       ,
                       0.2973295
                       ,
                       0.3123384
                       ,
                       0.3268434
                       ,
                       0.3407902
                       ,
                       0.3541372
                       ,
                       0.3668513
                       ,
                       0.3789132
                       ,
                       0.3903111
                       ,
                       0.4010426
                       ,
                       0.4111124
                       ,
                       0.4205309
                       ,
                       0.4293167
                       ,
                       0.4374927
                       ,
                       0.4450843
                       ,
                       0.4521207
                       ,
                       0.4586311
                       ,
                       0.464646
                       ,
                       0.4701962)
    )
    seroconv_v1_2 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("seroconv_v1_method2", length(dates[[Country]])),
        "seroconv" = c(0
                       ,
                       0.0000002595068
                       ,
                       0.0000009052471
                       ,
                       0.0000020863610
                       ,
                       0.0000040055240
                       ,
                       0.0000069074800
                       ,
                       0.0000110905100
                       ,
                       0.0000170260700
                       ,
                       0.0000251805100
                       ,
                       0.0000361302400
                       ,
                       0.0000505647400
                       ,
                       0.0000692647800
                       ,
                       0.0000931772000
                       ,
                       0.000123415
                       ,
                       0.000161334
                       ,
                       0.000208585
                       ,
                       0.000267231
                       ,
                       0.000339838
                       ,
                       0.000429617
                       ,
                       0.000540574
                       ,
                       0.000677751
                       ,
                       0.000847463
                       ,
                       0.001056923
                       ,
                       0.001315403
                       ,
                       0.00163426
                       ,
                       0.002027331
                       ,
                       0.002511545
                       ,
                       0.003107517
                       ,
                       0.003840509
                       ,
                       0.004732795
                       ,
                       0.00581531
                       ,
                       0.007123694
                       ,
                       0.008697921
                       ,
                       0.01058298
                       ,
                       0.0128288
                       ,
                       0.01549258
                       ,
                       0.0185979
                       ,
                       0.02218895
                       ,
                       0.02630777
                       ,
                       0.03098936
                       ,
                       0.03626077
                       ,
                       0.04213718
                       ,
                       0.04863002
                       ,
                       0.05569674
                       ,
                       0.06330289
                       ,
                       0.07139898
                       ,
                       0.07992161
                       ,
                       0.08879302
                       ,
                       0.09792019
                       ,
                       0.1072081
                       ,
                       0.1165531
                       ,
                       0.125851
                       ,
                       0.1349998
                       ,
                       0.1439095
                       ,
                       0.1525021
                       ,
                       0.16071
                       ,
                       0.1684786
                       ,
                       0.1757752
                       ,
                       0.1825793
                       ,
                       0.1888863
                       ,
                       0.1947048
                       ,
                       0.2000541
                       ,
                       0.2049605
                       ,
                       0.2094518
                       ,
                       0.2135645
                       ,
                       0.2173331
                       ,
                       0.2207932
                       ,
                       0.2239795
                       ,
                       0.2269231
                       ,
                       0.2296556
                       ,
                       0.232204
                       ,
                       0.2345935
                       ,
                       0.2368463
                       ,
                       0.2389802
                       ,
                       0.2410109
                       ,
                       0.2429516
                       ,
                       0.2448147
                       ,
                       0.2466107
                       ,
                       0.2483482
                       ,
                       0.2500346
                       ,
                       0.2516757
                       ,
                       0.2532761
                       ,
                       0.2548396
                       ,
                       0.2563696
                       ,
                       0.2578688
                       ,
                       0.259339
                       ,
                       0.2607818
                       ,
                       0.2621985
                       ,
                       0.2635897
                       ,
                       0.2649561
                       ,
                       0.2662978
                       ,
                       0.2676151
                       ,
                       0.2689076
                       ,
                       0.270175
                       ,
                       0.2714172
                       ,
                       0.2726338
                       ,
                       0.2738245
                       ,
                       0.2749888
                       ,
                       0.2761264
                       ,
                       0.2772368
                       ,
                       0.2783194
                       ,
                       0.2793742
                       ,
                       0.2804007
                       ,
                       0.2813991
                       ,
                       0.282369
                       ,
                       0.2833105
                       ,
                       0.2842239
                       ,
                       0.2851094
                       ,
                       0.2859672
                       ,
                       0.286798
                       ,
                       0.2876024
                       ,
                       0.288381
                       ,
                       0.2891347
                       ,
                       0.2898647
                       ,
                       0.2905721
                       ,
                       0.291258
                       ,
                       0.2919236
                       ,
                       0.2925704
                       ,
                       0.2931996
                       ,
                       0.2938124
                       ,
                       0.2944105
                       ,
                       0.2949952
                       ,
                       0.2955677
                       ,
                       0.2961295
                       ,
                       0.2966817
                       ,
                       0.2972255
                       ,
                       0.2977622
                       ,
                       0.2982928
                       ,
                       0.2988186
                       ,
                       0.2993405
                       ,
                       0.2998594
                       ,
                       0.3003762
                       ,
                       0.3008918
                       ,
                       0.3014068
                       ,
                       0.301922
                       ,
                       0.3024379
                       ,
                       0.3029549
                       ,
                       0.3034735
                       ,
                       0.303994
                       ,
                       0.3045166
                       ,
                       0.3050417
                       ,
                       0.3055691
                       ,
                       0.3060989
                       ,
                       0.306631
                       ,
                       0.3071653
                       ,
                       0.3077016
                       ,
                       0.3082398
                       ,
                       0.3087797
                       ,
                       0.3093208
                       ,
                       0.3098631
                       ,
                       0.310406
                       ,
                       0.3109494
                       ,
                       0.3114931
                       ,
                       0.3120367
                       ,
                       0.31258
                       ,
                       0.313123
                       ,
                       0.3136656
                       ,
                       0.3142076
                       ,
                       0.3147491
                       ,
                       0.3152903
                       ,
                       0.3158313
                       ,
                       0.3163723
                       ,
                       0.3169138
                       ,
                       0.3174562
                       ,
                       0.318
                       ,
                       0.3185458
                       ,
                       0.3190941
                       ,
                       0.3196455
                       ,
                       0.3202008
                       ,
                       0.3207605
                       ,
                       0.321325
                       ,
                       0.321895
                       ,
                       0.322471
                       ,
                       0.3230532
                       ,
                       0.3236421
                       ,
                       0.324238
                       ,
                       0.3248408
                       ,
                       0.3254505
                       ,
                       0.326067
                       ,
                       0.3266903
                       ,
                       0.3273199
                       ,
                       0.3279557
                       ,
                       0.3285973
                       ,
                       0.3292444
                       ,
                       0.3298965
                       ,
                       0.3305533
                       ,
                       0.3312145
                       ,
                       0.3318798
                       ,
                       0.3325488
                       ,
                       0.3332215
                       ,
                       0.3338977
                       ,
                       0.3345776
                       ,
                       0.3352611
                       ,
                       0.3359485
                       ,
                       0.33664
                       ,
                       0.3373361
                       ,
                       0.338037
                       ,
                       0.3387435
                       ,
                       0.3394561
                       ,
                       0.3401755
                       ,
                       0.3409023
                       ,
                       0.3416373
                       ,
                       0.3423811
                       ,
                       0.3431347
                       ,
                       0.3438986
                       ,
                       0.3446736
                       ,
                       0.3454605
                       ,
                       0.3462598
                       ,
                       0.3470722
                       ,
                       0.3478984
                       ,
                       0.348739
                       ,
                       0.3495941
                       ,
                       0.350464
                       ,
                       0.3513491
                       ,
                       0.3522493
                       ,
                       0.3531646
                       ,
                       0.3540949
                       ,
                       0.35504
                       ,
                       0.3559989
                       ,
                       0.356971
                       ,
                       0.3579554
                       ,
                       0.3589511
                       ,
                       0.359957
                       ,
                       0.3609718
                       ,
                       0.3619944
                       ,
                       0.3630235
                       ,
                       0.3640576
                       ,
                       0.3650954
                       ,
                       0.3661359
                       ,
                       0.367178
                       ,
                       0.3682206
                       ,
                       0.369263
                       ,
                       0.370304
                       ,
                       0.3713429
                       ,
                       0.372379
                       ,
                       0.3734119
                       ,
                       0.3744408
                       ,
                       0.3754654
                       ,
                       0.3764849
                       ,
                       0.3774987
                       ,
                       0.3785061
                       ,
                       0.3795067
                       ,
                       0.3804996
                       ,
                       0.3814844
                       ,
                       0.3824605
                       ,
                       0.3834276
                       ,
                       0.3843852
                       ,
                       0.385333
                       ,
                       0.3862709
                       ,
                       0.3871986
                       ,
                       0.3881164
                       ,
                       0.3890241
                       ,
                       0.3899219
                       ,
                       0.3908106
                       ,
                       0.3916907
                       ,
                       0.3925628
                       ,
                       0.3934278
                       ,
                       0.3942865
                       ,
                       0.3951402
                       ,
                       0.3959899
                       ,
                       0.3968371
                       ,
                       0.3976832
                       ,
                       0.3985297
                       ,
                       0.3993779
                       ,
                       0.4002294
                       ,
                       0.4010855
                       ,
                       0.4019476
                       ,
                       0.4028172
                       ,
                       0.4036956
                       ,
                       0.4045841
                       ,
                       0.405484
                       ,
                       0.4063964
                       ,
                       0.4073227
                       ,
                       0.4082642
                       ,
                       0.409222
                       ,
                       0.4101975
                       ,
                       0.4111919
                       ,
                       0.4122064
                       ,
                       0.4132422
                       ,
                       0.4143003
                       ,
                       0.4153817
                       ,
                       0.4164868
                       ,
                       0.417616
                       ,
                       0.4187698
                       ,
                       0.4199482
                       ,
                       0.4211511
                       ,
                       0.4223782
                       ,
                       0.4236292
                       ,
                       0.4249019
                       ,
                       0.426195
                       ,
                       0.4275067
                       ,
                       0.4288349
                       ,
                       0.4301771
                       ,
                       0.4315303
                       ,
                       0.4328918
                       ,
                       0.4342568
                       ,
                       0.4356214
                       ,
                       0.4369813
                       ,
                       0.4383322
                       ,
                       0.4396697
                       ,
                       0.440989
                       ,
                       0.4422858
                       ,
                       0.4435546
                       ,
                       0.4447907
                       ,
                       0.4459897
                       ,
                       0.4471475
                       ,
                       0.4482603
                       ,
                       0.4493244
                       ,
                       0.4503369
                       ,
                       0.4512947
                       ,
                       0.4521953
                       ,
                       0.4530371
                       ,
                       0.4538191
                       ,
                       0.4545411
                       ,
                       0.4552032
                       ,
                       0.4558063
                       ,
                       0.456352
                       ,
                       0.4568421
                       ,
                       0.4572794
                       ,
                       0.4576669
                       ,
                       0.4580078
                       ,
                       0.458306
                       ,
                       0.458565
                       ,
                       0.4587887
                       ,
                       0.4589808
                       ,
                       0.459145
                       ,
                       0.4592847
                       ,
                       0.4594032
                       ,
                       0.4595035
                       ,
                       0.4595882
                       ,
                       0.4596597
                       ,
                       0.45972
                       ,
                       0.4597709
                       ,
                       0.459814
                       ,
                       0.4598503
                       ,
                       0.4598811
                       ,
                       0.4599073
                       ,
                       0.4599296
                       ,
                       0.4599486
                       ,
                       0.4599649
                       ,
                       0.4599789
                       ,
                       0.4599909)
    )
    seroconv_v2_2 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("seroconv_v2_method2", length(dates[[Country]])),
        "seroconv" = c(0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0.0000000072156
                       ,
                       0.0000000184257
                       ,
                       0.0000000357549
                       ,
                       0.0000000626901
                       ,
                       0.0000001030042
                       ,
                       0.0000001617587
                       ,
                       0.0000002433313
                       ,
                       0.0000003542643
                       ,
                       0.0000005037072
                       ,
                       0.0000007015854
                       ,
                       0.0000009591503
                       ,
                       0.0000012915510
                       ,
                       0.0000017176010
                       ,
                       0.0000022613200
                       ,
                       0.0000029529500
                       ,
                       0.0000038328240
                       ,
                       0.0000049522950
                       ,
                       0.0000063778730
                       ,
                       0.0000081954800
                       ,
                       0.0000105177500
                       ,
                       0.0000134903600
                       ,
                       0.0000173012500
                       ,
                       0.0000221875700
                       ,
                       0.0000284580200
                       ,
                       0.0000365121000
                       ,
                       0.0000468651600
                       ,
                       0.0000601817900
                       ,
                       0.0000773215900
                       ,
                       0.0000993983600
                       ,
                       0.0001277096000
                       ,
                       0.000163994
                       ,
                       0.000210463
                       ,
                       0.000269919
                       ,
                       0.000345919
                       ,
                       0.000442978
                       ,
                       0.00056687
                       ,
                       0.000723444
                       ,
                       0.000920788
                       ,
                       0.001168805
                       ,
                       0.001479468
                       ,
                       0.001867312
                       ,
                       0.002349902
                       ,
                       0.002948829
                       ,
                       0.003683487
                       ,
                       0.00458035
                       ,
                       0.005669947
                       ,
                       0.006986769
                       ,
                       0.008569813
                       ,
                       0.01046259
                       ,
                       0.01271508
                       ,
                       0.01535873
                       ,
                       0.0184393
                       ,
                       0.02200299
                       ,
                       0.02609362
                       ,
                       0.03075172
                       ,
                       0.03601175
                       ,
                       0.04190657
                       ,
                       0.04842879
                       ,
                       0.05557817
                       ,
                       0.06334227
                       ,
                       0.07169514
                       ,
                       0.08059639
                       ,
                       0.08998916
                       ,
                       0.0998093
                       ,
                       0.1099832
                       ,
                       0.1204263
                       ,
                       0.131049
                       ,
                       0.1417646
                       ,
                       0.152489
                       ,
                       0.1631419
                       ,
                       0.1736483
                       ,
                       0.1839492
                       ,
                       0.1939927
                       ,
                       0.2037384
                       ,
                       0.2131563
                       ,
                       0.2222271
                       ,
                       0.2309391
                       ,
                       0.2392846
                       ,
                       0.2472684
                       ,
                       0.2548961
                       ,
                       0.2621786
                       ,
                       0.2691299
                       ,
                       0.2757654
                       ,
                       0.2821051
                       ,
                       0.2881685
                       ,
                       0.2939751
                       ,
                       0.2995459
                       ,
                       0.3048997
                       ,
                       0.3100544
                       ,
                       0.315027)
    )
    seroconv_v1_3 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("seroconv_v1_method3", length(dates[[Country]])),
        "seroconv" = c(0
                       ,
                       0.0000001743120
                       ,
                       0.0000006080589
                       ,
                       0.0000014014190
                       ,
                       0.0000026905300
                       ,
                       0.0000046397880
                       ,
                       0.0000074495480
                       ,
                       0.0000114325000
                       ,
                       0.0000168917300
                       ,
                       0.0000242041900
                       ,
                       0.0000338166900
                       ,
                       0.0000462285300
                       ,
                       0.0000620389400
                       ,
                       0.0000819435600
                       ,
                       0.000106784
                       ,
                       0.000137564
                       ,
                       0.000175535
                       ,
                       0.000222238
                       ,
                       0.000279589
                       ,
                       0.000349959
                       ,
                       0.000436314
                       ,
                       0.000542347
                       ,
                       0.000672247
                       ,
                       0.000831366
                       ,
                       0.001026217
                       ,
                       0.001264694
                       ,
                       0.001556411
                       ,
                       0.001913019
                       ,
                       0.002348757
                       ,
                       0.002876218
                       ,
                       0.003512822
                       ,
                       0.004278671
                       ,
                       0.005196395
                       ,
                       0.006291608
                       ,
                       0.007592993
                       ,
                       0.009133768
                       ,
                       0.01092635
                       ,
                       0.01299616
                       ,
                       0.01536768
                       ,
                       0.01806151
                       ,
                       0.02109393
                       ,
                       0.02447464
                       ,
                       0.02821194
                       ,
                       0.0322806
                       ,
                       0.0366615
                       ,
                       0.04132694
                       ,
                       0.04624115
                       ,
                       0.05135991
                       ,
                       0.05662976
                       ,
                       0.06199629
                       ,
                       0.06739939
                       ,
                       0.07277838
                       ,
                       0.07807371
                       ,
                       0.08323266
                       ,
                       0.08820949
                       ,
                       0.09296423
                       ,
                       0.09746446
                       ,
                       0.1016907
                       ,
                       0.1056308
                       ,
                       0.1092817
                       ,
                       0.1126485
                       ,
                       0.1157429
                       ,
                       0.1185804
                       ,
                       0.1211775
                       ,
                       0.1235557
                       ,
                       0.1257355
                       ,
                       0.1277379
                       ,
                       0.1295833
                       ,
                       0.1312902
                       ,
                       0.1328769
                       ,
                       0.1343592
                       ,
                       0.1357513
                       ,
                       0.1370661
                       ,
                       0.1383132
                       ,
                       0.1395015
                       ,
                       0.1406382
                       ,
                       0.1417301
                       ,
                       0.1427829
                       ,
                       0.1438015
                       ,
                       0.1447898
                       ,
                       0.145751
                       ,
                       0.1466877
                       ,
                       0.1476022
                       ,
                       0.1484963
                       ,
                       0.1493716
                       ,
                       0.1502293
                       ,
                       0.1510704
                       ,
                       0.1518956
                       ,
                       0.1527056
                       ,
                       0.1535008
                       ,
                       0.1542815
                       ,
                       0.1550479
                       ,
                       0.1558
                       ,
                       0.1565377
                       ,
                       0.1572609
                       ,
                       0.1579697
                       ,
                       0.1586638
                       ,
                       0.1593429
                       ,
                       0.160007
                       ,
                       0.1606556
                       ,
                       0.1612884
                       ,
                       0.1619053
                       ,
                       0.162506
                       ,
                       0.1630904
                       ,
                       0.1636584
                       ,
                       0.1642098
                       ,
                       0.1647448
                       ,
                       0.1652635
                       ,
                       0.165766
                       ,
                       0.1662526
                       ,
                       0.1667238
                       ,
                       0.1671798
                       ,
                       0.1676213
                       ,
                       0.1680488
                       ,
                       0.168463
                       ,
                       0.1688646
                       ,
                       0.1692542
                       ,
                       0.1696327
                       ,
                       0.1700007
                       ,
                       0.170359
                       ,
                       0.1707086
                       ,
                       0.1710501
                       ,
                       0.1713843
                       ,
                       0.1717121
                       ,
                       0.1720341
                       ,
                       0.172351
                       ,
                       0.1726637
                       ,
                       0.1729728
                       ,
                       0.1732791
                       ,
                       0.173583
                       ,
                       0.1738853
                       ,
                       0.1741864
                       ,
                       0.1744869
                       ,
                       0.1747873
                       ,
                       0.1750879
                       ,
                       0.1753891
                       ,
                       0.1756911
                       ,
                       0.1759942
                       ,
                       0.1762985
                       ,
                       0.1766041
                       ,
                       0.1769112
                       ,
                       0.1772196
                       ,
                       0.1775294
                       ,
                       0.1778404
                       ,
                       0.1781525
                       ,
                       0.1784657
                       ,
                       0.1787797
                       ,
                       0.1790944
                       ,
                       0.1794096
                       ,
                       0.179725
                       ,
                       0.1800406
                       ,
                       0.1803561
                       ,
                       0.1806713
                       ,
                       0.1809863
                       ,
                       0.1813007
                       ,
                       0.1816145
                       ,
                       0.1819279
                       ,
                       0.1822407
                       ,
                       0.1825529
                       ,
                       0.1828649
                       ,
                       0.1831766
                       ,
                       0.1834883
                       ,
                       0.1838003
                       ,
                       0.184113
                       ,
                       0.1844268
                       ,
                       0.1847419
                       ,
                       0.1850588
                       ,
                       0.185378
                       ,
                       0.1856998
                       ,
                       0.1860246
                       ,
                       0.1863527
                       ,
                       0.1866844
                       ,
                       0.1870199
                       ,
                       0.1873595
                       ,
                       0.1877032
                       ,
                       0.1880512
                       ,
                       0.1884035
                       ,
                       0.1887599
                       ,
                       0.1891204
                       ,
                       0.1894848
                       ,
                       0.1898528
                       ,
                       0.1902244
                       ,
                       0.1905991
                       ,
                       0.1909768
                       ,
                       0.1913572
                       ,
                       0.19174
                       ,
                       0.192125
                       ,
                       0.1925121
                       ,
                       0.1929011
                       ,
                       0.1932919
                       ,
                       0.1936845
                       ,
                       0.1940789
                       ,
                       0.1944752
                       ,
                       0.1948735
                       ,
                       0.1952741
                       ,
                       0.1956772
                       ,
                       0.1960829
                       ,
                       0.1964918
                       ,
                       0.196904
                       ,
                       0.1973201
                       ,
                       0.1977403
                       ,
                       0.1981651
                       ,
                       0.198595
                       ,
                       0.1990305
                       ,
                       0.199472
                       ,
                       0.1999201
                       ,
                       0.2003751
                       ,
                       0.2008376
                       ,
                       0.201308
                       ,
                       0.2017869
                       ,
                       0.2022745
                       ,
                       0.202771
                       ,
                       0.2032767
                       ,
                       0.2037916
                       ,
                       0.2043158
                       ,
                       0.2048493
                       ,
                       0.2053918
                       ,
                       0.2059433
                       ,
                       0.2065032
                       ,
                       0.207071
                       ,
                       0.207646
                       ,
                       0.2082278
                       ,
                       0.2088157
                       ,
                       0.2094089
                       ,
                       0.2100067
                       ,
                       0.2106084
                       ,
                       0.2112131
                       ,
                       0.21182
                       ,
                       0.2124286
                       ,
                       0.2130382
                       ,
                       0.2136482
                       ,
                       0.214258
                       ,
                       0.2148671
                       ,
                       0.2154748
                       ,
                       0.2160807
                       ,
                       0.2166845
                       ,
                       0.2172857
                       ,
                       0.217884
                       ,
                       0.2184791
                       ,
                       0.2190706
                       ,
                       0.2196581
                       ,
                       0.2202414
                       ,
                       0.2208201
                       ,
                       0.2213939
                       ,
                       0.2219628
                       ,
                       0.2225263
                       ,
                       0.2230845
                       ,
                       0.223637
                       ,
                       0.2241839
                       ,
                       0.2247249
                       ,
                       0.2252602
                       ,
                       0.2257896
                       ,
                       0.2263133
                       ,
                       0.2268316
                       ,
                       0.2273449
                       ,
                       0.2278535
                       ,
                       0.2283579
                       ,
                       0.2288586
                       ,
                       0.2293566
                       ,
                       0.2298524
                       ,
                       0.2303468
                       ,
                       0.2308409
                       ,
                       0.2313354
                       ,
                       0.2318312
                       ,
                       0.2323292
                       ,
                       0.2328301
                       ,
                       0.2333346
                       ,
                       0.2338435
                       ,
                       0.2343573
                       ,
                       0.2348765
                       ,
                       0.2354016
                       ,
                       0.2359332
                       ,
                       0.2364717
                       ,
                       0.2370178
                       ,
                       0.237572
                       ,
                       0.2381349
                       ,
                       0.2387073
                       ,
                       0.2392898
                       ,
                       0.2398831
                       ,
                       0.2404879
                       ,
                       0.2411048
                       ,
                       0.2417342
                       ,
                       0.2423765
                       ,
                       0.2430322
                       ,
                       0.2437016
                       ,
                       0.2443848
                       ,
                       0.2450819
                       ,
                       0.2457931
                       ,
                       0.2465176
                       ,
                       0.2472551
                       ,
                       0.2480049
                       ,
                       0.2487664
                       ,
                       0.2495385
                       ,
                       0.2503202
                       ,
                       0.2511105
                       ,
                       0.2519069
                       ,
                       0.2527077
                       ,
                       0.2535109
                       ,
                       0.2543142
                       ,
                       0.2551155
                       ,
                       0.2559121
                       ,
                       0.2567016
                       ,
                       0.2574805
                       ,
                       0.2582457
                       ,
                       0.2589944
                       ,
                       0.2597236
                       ,
                       0.2604306
                       ,
                       0.2611123
                       ,
                       0.2617665
                       ,
                       0.2623902
                       ,
                       0.2629813
                       ,
                       0.263538
                       ,
                       0.2640589
                       ,
                       0.264543
                       ,
                       0.2649899
                       ,
                       0.2653994
                       ,
                       0.2657719
                       ,
                       0.2661082
                       ,
                       0.2664096
                       ,
                       0.2666778
                       ,
                       0.2669146
                       ,
                       0.2671223
                       ,
                       0.2673032
                       ,
                       0.2674597
                       ,
                       0.2675942
                       ,
                       0.2677092
                       ,
                       0.2678071
                       ,
                       0.26789
                       ,
                       0.2679602
                       ,
                       0.2680193
                       ,
                       0.2680691
                       ,
                       0.2681111
                       ,
                       0.2681464
                       ,
                       0.2681762
                       ,
                       0.2682012
                       ,
                       0.2682225
                       ,
                       0.2682405
                       ,
                       0.2682558
                       ,
                       0.2682689
                       ,
                       0.2682801
                       ,
                       0.2682897
                       ,
                       0.268298)
            
    )
    seroconv_v2_3 = data.frame(
        "time" = dates[[Country]],
        "type" = rep("seroconv_v2_method3", length(dates[[Country]])),
        "seroconv" = c(0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0
                       ,
                       0.0000000071969
                       ,
                       0.0000000183585
                       ,
                       0.0000000355405
                       ,
                       0.0000000621300
                       ,
                       0.0000001017582
                       ,
                       0.0000001592723
                       ,
                       0.0000002387377
                       ,
                       0.0000003462454
                       ,
                       0.0000004899919
                       ,
                       0.0000006787660
                       ,
                       0.0000009221917
                       ,
                       0.0000012330880
                       ,
                       0.0000016270160
                       ,
                       0.0000021234830
                       ,
                       0.0000027466030
                       ,
                       0.0000035274950
                       ,
                       0.0000045054310
                       ,
                       0.0000057303780
                       ,
                       0.0000072656640
                       ,
                       0.0000091931990
                       ,
                       0.0000116169700
                       ,
                       0.0000146687300
                       ,
                       0.0000185133600
                       ,
                       0.0000233609100
                       ,
                       0.0000294792200
                       ,
                       0.0000372088300
                       ,
                       0.0000469819800
                       ,
                       0.0000593496500
                       ,
                       0.0000750150000
                       ,
                       0.0000947897100
                       ,
                       0.0001197432000
                       ,
                       0.000151219
                       ,
                       0.000190896
                       ,
                       0.000240882
                       ,
                       0.000303819
                       ,
                       0.000383047
                       ,
                       0.000482017
                       ,
                       0.000605389
                       ,
                       0.000758849
                       ,
                       0.000949258
                       ,
                       0.001184933
                       ,
                       0.001475928
                       ,
                       0.001834592
                       ,
                       0.00227178
                       ,
                       0.002802502
                       ,
                       0.003444099
                       ,
                       0.004216184
                       ,
                       0.00514104
                       ,
                       0.006243709
                       ,
                       0.00755333
                       ,
                       0.009086092
                       ,
                       0.0108679
                       ,
                       0.01292513
                       ,
                       0.01528287
                       ,
                       0.01796461
                       ,
                       0.02099073
                       ,
                       0.02438191
                       ,
                       0.02812976
                       ,
                       0.0322339
                       ,
                       0.03668727
                       ,
                       0.04147503
                       ,
                       0.04657408
                       ,
                       0.05195152
                       ,
                       0.05757157
                       ,
                       0.06338834
                       ,
                       0.06935265
                       ,
                       0.07541243
                       ,
                       0.08151769
                       ,
                       0.08762033
                       ,
                       0.09367411
                       ,
                       0.09963706
                       ,
                       0.1054765
                       ,
                       0.1111643
                       ,
                       0.1166799
                       ,
                       0.1220094
                       ,
                       0.1271462
                       ,
                       0.1320877
                       ,
                       0.1368339
                       ,
                       0.1413933
                       ,
                       0.1457744
                       ,
                       0.1499886
                       ,
                       0.1540493
                       ,
                       0.1579702
                       ,
                       0.1617673
                       ,
                       0.1654552
                       ,
                       0.1690482
                       ,
                       0.1725611
                       ,
                       0.1760065
                       ,
                       0.1793962
                       ,
                       0.1827405)
            
        
    )
    
    seroconv_df1 = rbind(seroconv_v1_1,seroconv_v1_2,seroconv_v1_3)
    seroconv_df2 = rbind(seroconv_v2_1,seroconv_v2_2,seroconv_v2_3)
    ggplot() +
        geom_line(data=seroconv_v1_1,aes(x = time, y = seroconv),col="red") +
        geom_line(data=seroconv_v1_2,aes(x = time, y = seroconv),col="blue") +
        geom_line(data=seroconv_v1_3,aes(x = time, y = seroconv),col="orange") +
        geom_line(data=seroconv_df2,aes(x = time, y = seroconv, color = type)) +
        scale_color_manual(name = "", labels = c("First scenario","Forth scenario","Fifth scenario"),
                           values = c("red","blue","orange"))+ 
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.3, 0.9)) +
        geom_text(aes(label = 'non-P.1'), x=ymd("2021-01-07"), y=1.01,size=4)+
        geom_text(aes(label = 'P.1'), x=ymd("2021-01-07"), y=0.17,size=4)+
        xlab("Date") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab("Cumulative incidence per capita") +
        scale_x_date(date_breaks = "1 month", labels = date_format("%b/%y")) +
        theme(legend.text = element_text(size = 11))+
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))
}

p_F <- plt_seroconv(out)
p_F
ggsave(here(paste0("transmission_model/figures/fig_CIR_contour.pdf")),plot=p_F)



