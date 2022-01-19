data 
{
  int <lower=1> M; 
  int <lower=1> N0; 
  int<lower=1> N[M]; 
  int<lower=1> N2; 
  int deaths[N2, M]; 
  matrix[N2, M] f; 
  int EpidemicStart[M];
  real pop[M];
  int W; 
  int week_index[M,N2];
  real SI[N2]; 
  real AR_SD_MEAN;
  int CasesStart;
  int T2; 
  real par; 
  int phylo_N_len;
  int phylo_N[phylo_N_len]; 
  int phylo_PSamples[phylo_N_len];
  int phylo_NSamples[phylo_N_len];
  matrix[N2, M] PCR_pos_prob; 
  matrix[N2, M] seroconv_cdf; 
  matrix[N2, M] serorev_surv; 
}

parameters 
{
  real<lower=0> R_difference; 
  real<lower=1> y_v1[M];
  real<lower=1> y_v2[M];
  real<lower=0> phi;
  real<lower=0,upper=1> cross;
  real<lower=0> tau;
  real <lower=0, upper=100> ifr1[M];
  real <lower=0> RR[M];
  real<lower=0,upper=2> weekly_effect[W];
}

transformed parameters 
{
  matrix[N2,M] cumm_sum = rep_matrix(0,N2,M);
  matrix[N2, M] prediction = rep_matrix(0,N2,M);
  matrix[N2, M] E_deaths  = rep_matrix(0,N2,M);
  matrix[N2,M] immune = rep_matrix(0,N2,M);
  matrix[N2, M] prediction_v1 = rep_matrix(0,N2,M);
  matrix[N2, M] prediction_v2 = rep_matrix(0,N2,M);
  matrix[N2, M] E_deaths_v1  = rep_matrix(0,N2,M);
  matrix[N2, M] E_deaths_v2  = rep_matrix(0,N2,M);
  matrix[N2, M] Rt_v1 = rep_matrix(0,N2,M);
  matrix[N2, M] Rt_v2 = rep_matrix(0,N2,M);
  matrix[N2, M] Rt_adj_immune_v1 = rep_matrix(0,N2,M);
  matrix[N2, M] Rt_adj_immune_v2 = rep_matrix(0,N2,M);
  matrix[N2,M] cumm_sum_v1 = rep_matrix(0,N2,M);
  matrix[N2,M] cumm_sum_v2 = rep_matrix(0,N2,M);
  matrix[N2,M] immune_v1 = rep_matrix(0,N2,M);
  matrix[N2,M] immune_v2 = rep_matrix(0,N2,M);
  matrix[N2,M] alpha_sus1 = rep_matrix(0,N2,M);
  matrix[N2,M] alpha_sus2 = rep_matrix(0,N2,M);
  matrix[N2,M] n1 = rep_matrix(0,N2,M);
  matrix[N2,M] n2 = rep_matrix(0,N2,M);
  matrix[N2,M] seroconv_v1 = rep_matrix(0, N2, M);
  matrix[N2,M] seroconv_v2 = rep_matrix(0, N2, M);
  matrix[N2,M] seropos_v1 = rep_matrix(0, N2, M);
  matrix[N2,M] seropos_v2 = rep_matrix(0, N2, M);
  matrix[N2,M] pcr_pos_v1 = rep_matrix(0, N2, M);
  matrix[N2,M] pcr_pos_v2 = rep_matrix(0, N2, M);
  matrix[N2,M] E_fraction = rep_matrix(0,N2,M);
  real <lower=0> ifr2[M];
  
  for (m in 1:M)
  {
    for (i in 2:N0)
    {
      cumm_sum_v1[i,m] = cumm_sum_v1[i-1,m] + y_v1[m];
    }
    
    for ( i in ( T2+1 ):( T2+1 ) ) 
    {
      cumm_sum_v2[i,m] = cumm_sum_v2[i-1,m] + y_v2[m];
    }
    
    prediction_v1[1:N0,m] = rep_vector(y_v1[m],N0); 
    prediction_v2[T2:(T2+N0-1),m] = rep_vector(y_v2[m],N0); 
    Rt_v1[1:7,m] = rep_vector(3.28 * weekly_effect[1],7);
    Rt_v1[8:14,m] = rep_vector(3.28 * weekly_effect[2],7);
    Rt_v1[15:21,m] = rep_vector(3.28 * weekly_effect[3],7);
    Rt_v1[22:28,m] = rep_vector(3.28 * weekly_effect[4],7);
    Rt_v1[29:35,m] = rep_vector(3.28 * weekly_effect[5],7);
    Rt_v1[36:42,m] = rep_vector(3.28 * weekly_effect[6],7);
    Rt_v1[43:49,m] = rep_vector(3.28 * weekly_effect[7],7);
    Rt_v1[50:56,m] = rep_vector(3.28 * weekly_effect[8],7);
    Rt_v1[57:63,m] = rep_vector(3.28 * weekly_effect[9],7);
    Rt_v1[64:70,m] = rep_vector(3.28 * weekly_effect[10],7);
    Rt_v1[71:77,m] = rep_vector(3.28 * weekly_effect[11],7);
    Rt_v1[78:84,m] = rep_vector(3.28 * weekly_effect[12],7);
    Rt_v1[85:91,m] = rep_vector(3.28 * weekly_effect[13],7);
    Rt_v1[92:98,m] = rep_vector(3.28 * weekly_effect[14],7);
    Rt_v1[99:105,m] = rep_vector(3.28 * weekly_effect[15],7);
    Rt_v1[106:112,m] = rep_vector(3.28 * weekly_effect[16],7);
    Rt_v1[113:119,m] = rep_vector(3.28 * weekly_effect[17],7);
    Rt_v1[120:126,m] = rep_vector(3.28 * weekly_effect[18],7);
    Rt_v1[127:133,m] = rep_vector(3.28 * weekly_effect[19],7);
    Rt_v1[134:140,m] = rep_vector(3.28 * weekly_effect[20],7);
    Rt_v1[141:147,m] = rep_vector(3.28 * weekly_effect[21],7);
    Rt_v1[148:154,m] = rep_vector(3.28 * weekly_effect[22],7);
    Rt_v1[155:161,m] = rep_vector(3.28 * weekly_effect[23],7);
    Rt_v1[162:168,m] = rep_vector(3.28 * weekly_effect[24],7);
    Rt_v1[169:175,m] = rep_vector(3.28 * weekly_effect[25],7);
    Rt_v1[176:182,m] = rep_vector(3.28 * weekly_effect[26],7);
    Rt_v1[183:189,m] = rep_vector(3.28 * weekly_effect[27],7);
    Rt_v1[190:196,m] = rep_vector(3.28 * weekly_effect[28],7);
    Rt_v1[197:203,m] = rep_vector(3.28 * weekly_effect[29],7);
    Rt_v1[204:210,m] = rep_vector(3.28 * weekly_effect[30],7);
    Rt_v1[211:217,m] = rep_vector(3.28 * weekly_effect[31],7);
    Rt_v1[218:224,m] = rep_vector(3.28 * weekly_effect[32],7);
    Rt_v1[225:231,m] = rep_vector(3.28 * weekly_effect[33],7);
    Rt_v1[232:238,m] = rep_vector(3.28 * weekly_effect[34],7);
    Rt_v1[239:245,m] = rep_vector(3.28 * weekly_effect[35],7);
    Rt_v1[246:252,m] = rep_vector(3.28 * weekly_effect[36],7);
    Rt_v1[253:259,m] = rep_vector(3.28 * weekly_effect[37],7);
    Rt_v1[260:266,m] = rep_vector(3.28 * weekly_effect[38],7);
    Rt_v1[267:273,m] = rep_vector(3.28 * weekly_effect[39],7);
    Rt_v1[274:280,m] = rep_vector(3.28 * weekly_effect[40],7);
    Rt_v1[281:287,m] = rep_vector(3.28 * weekly_effect[41],7);
    Rt_v1[288:294,m] = rep_vector(3.28 * weekly_effect[42],7);
    Rt_v1[295:301,m] = rep_vector(3.28 * weekly_effect[43],7);
    Rt_v1[302:308,m] = rep_vector(3.28 * weekly_effect[44],7);
    Rt_v1[309:315,m] = rep_vector(3.28 * weekly_effect[45],7);
    Rt_v1[316:322,m] = rep_vector(3.28 * weekly_effect[46],7);
    Rt_v1[323:329,m] = rep_vector(3.28 * weekly_effect[47],7);
    Rt_v1[330:336,m] = rep_vector(3.28 * weekly_effect[48],7);
    Rt_v1[337:342,m] = rep_vector(3.28 * weekly_effect[49],6);
    Rt_v2[T2:N2,m] = Rt_v1[T2:N2,m] * (R_difference); 
    Rt_adj_immune_v1[1:N0,m] = Rt_v1[1:N0,m]; 
    for (i in (N0+1):N2) 
    {
      real convolution_v1 = 1e-15;
      for (j in 1:(i-1)) 
      {
        convolution_v1 += prediction_v1[j, m] * SI[i-j];
        immune_v1[i,m] += prediction_v1[j, m] * exp( - 0.5 * (i-j) * (i-j) / (par * par) );
      }
      if ( i > (T2) ) 
      {
        real convolution_v2 = 1e-15;
        for(j in T2:(i-1))  // start the v2 convolution at T2
        {
          convolution_v2 += prediction_v2[j, m] * SI[i-j];
          immune_v2[i,m] += prediction_v2[j, m] * exp( - 0.5 * (i-j) * (i-j) / ( par * par) );
        }
        alpha_sus2[i,m] = (1 - cross) * immune_v2[i,m] / pop[m];
        n2[i,m] = immune_v2[i,m] + cross * (immune_v1[i,m] * ( 1 - alpha_sus2[i,m] ));
        Rt_adj_immune_v2[i,m] = ( 1 - n2[i,m] / pop[m]) * Rt_v2[i,m]; 
        prediction_v2[i, m] = Rt_adj_immune_v2[i,m] * convolution_v2;
        cumm_sum_v2[i,m]  = cumm_sum_v2[i-1,m] +  prediction_v2[i,m];
      }
      alpha_sus1[i,m] = (1 - cross) * immune_v1[i,m] / pop[m];
      n1[i,m] = immune_v1[i,m] + cross * (immune_v2[i,m] * ( 1 - alpha_sus1[i,m] ));
      
      Rt_adj_immune_v1[i,m] = ( 1 - n1[i,m] / pop[m]) * Rt_v1[i,m];
      prediction_v1[i, m] = Rt_adj_immune_v1[i,m] * convolution_v1;
      cumm_sum_v1[i,m] = cumm_sum_v1[i-1,m] +  prediction_v1[i,m];
      
      cumm_sum[i,m] = cumm_sum_v1[i,m] + cumm_sum_v2[i,m];
      prediction[i, m] = prediction_v1[i, m] + prediction_v2[i, m];
      immune[i, m] = immune_v1[i, m] + immune_v2[i, m];
    }
    E_deaths_v1[1, m]= 1e-15 * prediction_v1[1,m];
    E_deaths_v2[1, m]= 1e-15 * prediction_v2[1,m];
    E_deaths[1, m]= 1e-15 * (prediction_v1[1,m] + prediction_v2[1,m]);
    for (i in 2:N2)
    {
      for(j in 1:(i-1))
      {	
        pcr_pos_v1[i,m] += prediction_v1[j,m] * PCR_pos_prob[i-j,m];
        seroconv_v1[i,m] += prediction_v1[j,m] * seroconv_cdf[i-j,m];
        E_deaths_v1[i,m] += prediction_v1[j,m] * f[i-j,m] * ifr1[m];
      }
      
      if (i > T2)  // start accumulating deaths after T2
      {
        for(j in T2:(i-1)) 
        {
          pcr_pos_v2[i,m] += prediction_v2[j,m] * PCR_pos_prob[i-j,m];
          seroconv_v2[i,m] += prediction_v2[j,m] * seroconv_cdf[i-j,m];
          ifr2[m] = ifr1[m] * RR[m];
          E_deaths_v2[i,m] += prediction_v2[j,m] * f[i-j,m] * ifr2[m];
        }
      }
      E_fraction[i,m] = pcr_pos_v2[i,m]/(pcr_pos_v1[i,m] + pcr_pos_v2[i,m]);	
      E_deaths[i,m] = E_deaths_v1[i,m] + E_deaths_v2[i,m];
    }
  }
}

model 
{
  ifr1 ~ normal(0.32,0.1);
  cross ~ beta(2,1);
  RR ~ lognormal(0,0.5);
  R_difference ~ normal(1,1);
  tau ~ exponential(0.03);
  for (m in 1:M)
  {
    y_v1[m] ~ exponential(1/tau);
    y_v2[m] ~ normal(1,1);
  }
  for (j in 1:W)
  {
    weekly_effect[j] ~ uniform(0,2);
  }
  
  phi ~ normal(0,5);
  for(m in 1:M)
  {
    deaths[EpidemicStart[m]:(N[m]), m] ~ neg_binomial_2(E_deaths[EpidemicStart[m]:(N[m]), m], phi);
    for ( i in 1:phylo_N_len )
    {
      phylo_PSamples[i] ~ binomial(phylo_NSamples[i]+phylo_PSamples[i],E_fraction[phylo_N[i],m]);
    }
  }
}

generated quantities 
{
  real RR_prior = lognormal_rng(0,0.5);
  real cross_prior = beta_rng(2,1);
  real R_difference_prior = normal_rng(1,1);
  real ifr1_prior = normal_rng(0.32,0.1);
}

