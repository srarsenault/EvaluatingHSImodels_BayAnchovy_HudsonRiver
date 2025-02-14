#==============================================================================
#
# Evaluating effects of data quality and variable weighting on habitat suitability modelling 
# Script 1) Environmental data spatial interpolation (Step 1) and cross-validation (Step 2)
#
# Authors: Stephanie Arsenault, Robyn Linner, and Yong Chen
# Date: 02.12.2025
#
## LOAD REQUIRED PACKAGES =====================================================
#
l_packages <- c("gstat", "automap", "ggplot2", "ggpubr", "sp", "dplyr", "Metrics", "gridExtra")
for (p in l_packages){
  if(! p %in% installed.packages()){
    install.packages(p, dependencies = TRUE)
  }
  library(p, character.only = T)
}
#

## Step 1- SPATIALLY INTERPOLATE ENVIRONMENTAL DATA ============================
#
#Set working directory
setwd() #add path

#Load grid for kriging
mesh2=read.csv() #add path
mesh2=mesh2[-c(1,4:8)]
coordinates(mesh2) = ~lon + lat
grid_data=mesh2

#For loop to spatially interpolate the environmental data, run separately for each variable 
WQSnames<- list.files()

for(i in 1:length(WQSnames)){
  print(i)
  idpdf<- unlist(strsplit(WQSnames[i], split='.', fixed=TRUE))[1]
  
  Data_Set=read.csv(WQSnames[i])
  
  #Specify the environmental variable to be interpolated
  colnames(Data_Set)[11]="lon"
  colnames(Data_Set)[12]="lat"
  Data_Set=Data_Set[,c(11,12,5)]
  Data_Set=Data_Set[!is.na(Data_Set$DISSOLVED.OXYGEN),]
  Data_Set=aggregate(DISSOLVED.OXYGEN ~lat + lon, Data_Set, mean)
  coordinates(Data_Set) = ~lon + lat
  
  #Fit variogram
  lzn.vgm <- variogram(DISSOLVED.OXYGEN ~ lat+lon, Data_Set, width=0.1)
  lzn.fit <- fit.variogram(lzn.vgm, vgm(c("Nug", "Sph","Mat", "Pen" , "Ste", "Cir" , "Bes", "Exc")))
  print(lzn.fit)
  test <- variogramLine(lzn.fit, maxdist = max(lzn.vgm$dist))
  
  #Plot variogram
  vario<- ggplot(lzn.vgm,aes(x=dist,y=gamma))+
    geom_point()+
    geom_line(data=test)+
    theme_gray()+
    labs(x="Distance", y=expression("Semivariance"~ (gamma)))+
    ggtitle("Variogram")+theme(plot.title = element_text(hjust = 0.5))

  #Perform kriging
  lzn.kriged <-krige(DISSOLVED.OXYGEN ~ lat+lon, Data_Set, grid_data, model = lzn.fit)
  
  k=data.frame(lzn.kriged)
  k=k[,c(1:4)]
  colnames(k)[3]="Sal"
  colnames(k)[4]="Variance"
  
  #Create new data frame with kriged data
  if (i==1){
    krige_all=k
    names(krige_all)[3]=idpdf
  } else {
    krige_all[,ncol(krige_all)+1]=k$Sal
    names(krige_all)[ncol(krige_all)]=idpdf
  }
}
#

## Step 2- CROSS-VALIDATE INTERPOLATED ENVIRONMENTAL DATA USING OBSERVED DATA=========
# 

# Load observed water quality data
obs=read.csv() #add path

# Load kriged modeled water quality data
mod=read.csv() #add path

# Join observed and modeled environmental data by latitude, longitude, and year
combo=left_join(averages,mod, by=c("lat", "lon", "Year"))

# Remove NA values
combo=na.omit(combo)

#
# COMPUTE THE QUANTITATIVE METRICS TO ASSESS THE MODELLED DATA
#

# Correlation coefficient
Temp_corr=cor(combo$Temperature.mod, combo$Temperature.obs)
Sal_corr=cor(combo$Salinity.mod, combo$Salinity.obs) 
DO_corr=cor(combo$Dissolved_Oxygen.mod, combo$Dissolved_Oxygen.obs)

# Average Relative error (ARE)
# Variable_Err is just modeled minus observed
combo$Temp_Err=combo$Temperature.mod - combo$Temperature.obs
combo$Sal_Err=combo$Salinity.mod - combo$Salinity.obs
combo$DO_Err=combo$Dissolved_Oxygen.mod - combo$Dissolved_Oxygen.obs

Temp_ARE= mean(combo$Temp_Err)
Sal_ARE= mean(combo$Sal_Err)
DO_ARE= mean(combo$DO_Err)

# Average absolute error (AAE)
Temp_AAE= mean(abs(combo$Temp_Err))
Sal_AAE= mean(abs(combo$Sal_Err))
DO_AAE= mean(abs(combo$DO_Err))

# Root Mean Square Error (RMSE)
Temp_RMSE= rmse(combo$Temperature.obs, combo$Temperature.mod)
Sal_RMSE= rmse(combo$Salinity.obs, combo$Salinity.mod)
DO_RMSE= rmse(combo$Dissolved_Oxygen.obs, combo$Dissolved_Oxygen.mod)

# Generate data frame to store all the metrics
error_metrics <- data.frame(
  Metric = c("DO_AAE", "DO_ARE", "DO_corr", "DO_RMSE", "Sal_AAE", "SAL_ARE","SAL_corr", "Sal_RMSE", "Temp_AAE", "Temp_ARE", "Temp_corr", "Temp_RMSE"),
  Value = c(DO_AAE, DO_ARE, DO_corr, DO_RMSE, Sal_AAE, Sal_ARE,Sal_corr, Sal_RMSE, Temp_AAE, Temp_ARE, Temp_corr, Temp_RMSE))

#
# PLOT LINEAR REGRESSION OF MODELED DATA TO OBSERVED DATA. 
#

# Temperature
lr_temp <- lm(Temperature.mod ~ Temperature.obs, data=combo)
m1=summary(lr_temp)

t <- ggplot(combo, aes(x=Temperature.obs, y=Temperature.mod)) +
  geom_point(color="lightgray") +
  theme_bw()+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth=0.75) +
  labs(x = "Observed Temperature (°C)", y = "Kriged Temperature (°C)") +
  ggtitle(paste0(#"dist ~ ",round(m1$coefficients[1],2)," + ", 
    # round(m1$coefficients[2],2),"x",
    "    R squared: ",round(m1$r.squared,2),
    "    P-value: ", format.pval(pf(m1$fstatistic[1], # F-statistic
                                    m1$fstatistic[2], # df
                                    m1$fstatistic[3], # df
                                    lower.tail = FALSE))))

# Salinity
lr_sal <- lm(Salinity.mod ~ Salinity.obs, data=combo)
m2=summary(lr_sal)

s <- ggplot(combo, aes(x=Salinity.obs, y=Salinity.mod)) +
  geom_point(color="lightgray") +
  theme_bw()+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method = "lm", se = FALSE, color = "black",linewidth=0.75) +
  labs(x = "Observed Conductivity (µS/cm)", y = "Kriged Conductivity (µS/cm)") +
  ggtitle(paste0(#"dist ~ ",round(m1$coefficients[1],2)," + ", 
    #round(m1$coefficients[2],2),"x",
    "    R squared: ",round(m2$r.squared,2),
    "    P-value: ", format.pval(pf(m2$fstatistic[1], # F-statistic
                                    m2$fstatistic[2], # df
                                    m2$fstatistic[3], # df
                                    lower.tail = FALSE))))

# Dissolved oxygen
lr_do <- lm(Dissolved_Oxygen.mod ~ Dissolved_Oxygen.obs, data=combo)
m3=summary(lr_do)

d <- ggplot(combo, aes(x=Dissolved_Oxygen.obs, y=Dissolved_Oxygen.mod)) +
  geom_point(color="lightgray") +
  theme_bw()+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth=0.75) +
  labs(x = "Observed Dissolved Oxygen (mg/l)", y = "Kriged Dissolved Oxygen (mg/l)") +
  ggtitle(paste0(#"dist ~ ",round(m1$coefficients[1],2)," + ", 
    #round(m1$coefficients[2],2),"x",
    "    R squared: ",round(m3$r.squared,2),
    "    P-value: ", format.pval(pf(m3$fstatistic[1], # F-statistic
                                    m3$fstatistic[2], # df
                                    m3$fstatistic[3], # df
                                    lower.tail = FALSE))))

# Arrange the plots
grid.arrange(d,s,t, ncol = 3)
#


