# has global warming accelerated relative to human ghg emissions?
library(tidyverse) 
library(lmtest)
# working directory
setwd("/PATH/TO/YOUR/DIRECTORY/")

lag <- dplyr::lag

# sum of two regression coefficients using vcov matrix
sum_of_coefs <- function(mdl,coefs_list){
  coef_ests <- as.vector(coef(mdl)[coefs_list])
  summed_coef <- sum(coef_ests)
  summed_se <- sqrt(sum(vcov(mdl)[coefs_list,coefs_list]))  
  t_stat <- summed_coef/summed_se
  p <- round(2*pt(-abs(t_stat),df=nobs(mdl)-1),4)
  ci2_5 <- summed_coef - 1.96*summed_se
  ci97_5 <- summed_coef + 1.96*summed_se
  return(c(summed_coef,summed_se,t_stat,p,ci2_5,ci97_5))
}

# utility df 
month_df <- data.frame("month_name"=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                       "month"=seq(1,12,by=1))

# global temperature, monthly
gmt <- read.csv("GLB.Ts+dSST.csv") %>%
  filter(Year<2024) %>% mutate_at(c(1:13),as.numeric) %>%
  select(Year,Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec) %>%
  pivot_longer(!Year,names_to="month_name",values_to="gmt") %>% rename(year=Year)

# southern oscillation index, monthly
soi <- read.table("soi.txt",header=TRUE) %>%
  filter(YEAR<2024) %>% 
  pivot_longer(!YEAR,names_to="month_name",values_to="soi") %>%
  rename(year=YEAR) %>% mutate(month_name=str_to_title(month_name))

# sunspots
sunspots <- read.table("SN_m_tot_V2.0.txt") %>%
  rename(year=V1,month=V2,sunspots=V4) %>% select(year,month,sunspots) %>%
  left_join(month_df,by="month")

# emissions
emissions <- read.csv("global_carbon_budget_extracted.csv") %>% rename(year=Year) %>%
  mutate(emissions_gtc = fossil.emissions.excluding.carbonation+land.use.change.emissions) %>%
  select(year,emissions_gtc) %>% 
  mutate(cumulative_emissions_gtc = cumsum(replace_na(emissions_gtc, 0)))

# make final data
# we'll use just constants/offsets for volcanic eruptions in 1982 (el chichon) and 1991 (pinatubo)
data <- gmt %>% left_join(soi,by=c("year","month_name")) %>%
  left_join(month_df,by="month_name") %>% 
  left_join(emissions,by="year") %>% 
  left_join(sunspots,by=c("year","month","month_name")) %>%
  mutate(volc = case_when(year %in% c(1982,1991) ~ 1, 
                          !year %in% c(1982,1991) ~ 0)) %>%
  mutate(soi_lag1 = lag(soi,1),soi_lag2 = lag(soi,2),soi_lag3 = lag(soi,3),
         soi_lag4 = lag(soi,4),volc_lag1 = lag(volc,12),volc_lag2 = lag(volc,24)) %>% 
  drop_na()

# adjusting global temperature data
# more or less following tamino's approach
# https://assets-eu.researchsquare.com/files/rs-6079807/v1_covered_778c44a8-77f0-4283-bf58-055aec3b288f.pdf?c=1740984327
loess_mdl <- loess(gmt ~ year,data=data)
data$gmt_loess_residual <- data$gmt - predict(loess_mdl,data$year)
mdl <- feols(gmt_loess_residual ~ soi + soi_lag1 + soi_lag2 + soi_lag3 + soi_lag4 + volc + volc_lag1 + volc_lag2 + sunspots,data=data)
summary(mdl)
data$gmt_adjusted <- data$gmt - predict(mdl,data)
data_ann <- data %>% group_by(year) %>% summarize(gmt_anomaly=mean(gmt_adjusted),
                                                  emissions=first(cumulative_emissions_gtc),
                                                  ann_emissions=first(emissions_gtc))


## look at first and second 25-year periods
# an arbitrary choice, but useful and easy to understand
data_periods <- data_ann %>% filter(year>=1974,year<=2023) %>% 
  mutate(period = case_when(year>=1974 & year<1999 ~ 0,
                            year>=1999 & year<=2023 ~ 1)) %>%
  mutate(period=factor(period))

yr_mdl_periods <- lm(gmt_anomaly ~ year*period,data=data_periods)
lmtest::coeftest(yr_mdl_periods,vcov=vcovHAC(yr_mdl_periods)) #adjust for autocorrelation
emis_mdl_periods <- lm(gmt_anomaly ~ emissions*period,data=data_periods)
lmtest::coeftest(emis_mdl_periods,vcov=vcovHAC(emis_mdl_periods)) # adjust for autocorrelation

# assemble data and write out some useful dataframes
yr_mdl_df <- data.frame("term"=c("estimate","se","t","p","ci_2_5","ci_97_5"),
                        "year_period0"=sum_of_coefs(yr_mdl_periods,c("year")),
                        "year_period1"=sum_of_coefs(yr_mdl_periods,c("year","year:period1")),
                        "difference"=sum_of_coefs(yr_mdl_periods,c("year:period1")))
emis_mdl_df <- data.frame("term"=c("estimate","se","t","p","ci_2_5","ci_97_5"),
                          "emissions_period0"=sum_of_coefs(emis_mdl_periods,c("emissions")),
                          "emissions_period1"=sum_of_coefs(emis_mdl_periods,c("emissions","emissions:period1")),
                          "difference"=sum_of_coefs(emis_mdl_periods,c("emissions:period1")))

write.csv(yr_mdl_df,"year_interaction_model_results.csv")
write.csv(emis_mdl_df,"emissions_interaction_model_results.csv")
write.csv(data_periods,"annual_temperature_emissions_1974-2023.csv")

