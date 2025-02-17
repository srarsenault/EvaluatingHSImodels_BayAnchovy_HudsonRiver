#=============================================================================
#
# Evaluating effects of data quality and variable weighting on habitat suitability modelling 
# Script 3) Habitat suitability index (HSI) modelling
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

## Step 0- DATA PREP  =========================================================
# 
# Load data
df=read.csv() #add path


# Variance Inflation Factor Test (VIF) to assess for multicollinearity among environmental variables
vif_func<-function(in_frame,thresh=10,trace=T,...){
  require(fmsb)
  
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]))
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2])))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}
vif <- df[, c("WATER.TEMPERATURE", "SALINITY", "RIVER.DEPTH", "DISSOLVED.OXYGEN")] #add names of environmental variables
vif_func(in_frame = vif, thresh=3, trace=T) #threshold of 3


# Outliers
## Check environmental data for outliers
Out_temp = identify_outliers(df, variable="WATER.TEMPERATURE")
Out_depth = identify_outliers(df, variable = "RIVER.DEPTH")
Out_do = identify_outliers(df, variable = "DISSOLVED.OXYGEN")
Out_sal = identify_outliers(df, variable="SALINITY")

## Check environmental data for extreme outliers
ext_DO = Out_do[Out_do$is.extreme == 'TRUE',]
ext_depth = Out_depth[Out_depth$is.extreme == 'TRUE',]
ext_temp = Out_temp[Out_temp$is.extreme =='TRUE',]
ext_sal = Out_sal[Out_sal$is.extreme == 'TRUE',]

## Remove extreme outliers
df$RIVER.DEPTH[df$RIVER.DEPTH >=29.3]<-NA # Remove extreme outliers


# Determine the number of bins for suitability index (SI) modelling using the Freedman-Diaconis method
## Remove NA values 
df_temp = df$WATER.TEMPERATURE
df_temp = df_temp[!is.na(df_temp)]
df_DO = df$DISSOLVED.OXYGEN
df_DO = df_DO[!is.na(df_DO)]
df_sal = df$SALINITY
df_sal = df_sal[!is.na(df_sal)]
df_dep = df$RIVER.DEPTH
df_dep = df_dep[!is.na(df_dep)]

## Determine number of bins 
nclass.FD(df_temp)
nclass.FD(df_DO)
nclass.FD(df_sal)
nclass.FD(df_dep)
#

## STEP 1- SUITABILITY INDICES (SI) FOR EACH ENVIRONMENTAL VARIABLE ===========
#
# Suitability index (SI) temperature
## Break up the variable into bins
temp = df$WATER.TEMPERATURE
temp = as.numeric(as.character(temp))
temperature_int = (classIntervals(temp, 20 ,style = "fisher"))
temperature_int[[2]][1] = temperature_int[[2]][1]-0.1
temperature_bins = (cut(temp, breaks = temperature_int$brks))

## Assign a bin to each survey station
df$temperature_bins = cut(temp, breaks = temperature_int$brks)

## Average catch within each bin
temperature = aggregate(cpue ~ temperature_bins, data = df, FUN="mean")
colnames(temperature)[1] = "temperature_bins"
colnames(temperature)[2] = "abundance"

## Find relative SI values
temperature$SI_temp = ((temperature$abundance - min(temperature[,2]))/(max(temperature[,2]) - min(temperature[,2])))
temperature_axis = temperature$temperature_bins
plot(temperature$temperature_bins, temperature$SI_temp, xlab = "temp bins", ylab = "SI temp")
temperature$bin = 1:length(temperature[,1])

## Run a GAM for final SI values
g = gam(abundance ~ s(bin), data = temperature)
temperature$abundance = predict(g, newdata = temperature, scale = "response")
temperature$SI_temp = ((temperature$abundance - min(temperature[,2]))/(max(temperature[,2]) - min(temperature[,2])))
plot(temperature$temperature_bins, temperature$SI_temp, xlab = "temp bins", ylab = "SI temp")


# SI Dissolved oxygen
DO = df$DISSOLVED.OXYGEN
DO = as.numeric(as.character(DO))
DO_int = (classIntervals(DO, 20, style = "fisher"))
DO_int[[2]][1] = DO_int[[2]][1]-0.1
DO_bins = (cut(DO, breaks = DO_int$brks))
df$DO_bins = cut(DO, breaks = DO_int$brks)
DO = aggregate(cpue ~ DO_bins, data = df, FUN="mean")
colnames(DO)[1] = "DO_bins"
colnames(DO)[2] = "abundance"
DO$SI_DO = ((DO$abundance - min(DO[,2]))/(max(DO[,2]) - min(DO[,2])))
DO_axis = DO$DO_bins
plot(DO$DO_bins, DO$SI_DO, xlab = "DO bins", ylab = "SI DO")
DO$bin = 1:length(DO[,1])
g = gam(abundance ~ s(bin), data = DO)
DO$abundance = predict(g, newdata = DO, scale = "response")
DO$SI_DO = ((DO$abundance - min(DO[,2]))/(max(DO[,2]) - min(DO[,2])))
plot(DO$DO_bins, DO$SI_DO, xlab = "DO bins", ylab = "SI DO")


# SI Depth
dep = df$RIVER.DEPTH
dep = as.numeric(as.character(dep))
dep_int = (classIntervals(dep, 20, style = "fisher"))
dep_bins = (cut(dep, breaks = dep_int$brks))
df$dep_bins = cut(dep, breaks = dep_int$brks)
dep = aggregate(cpue~ dep_bins, data = df, FUN = "mean")
colnames(dep)[1] = "dep_bins"
colnames(dep)[2] = "abundance"
dep$SI_dep = ((dep$abundance - min(dep[,2]))/(max(dep[,2])-min(dep[,2])))
dep_axis = dep$dep_bins
plot(dep$dep_bins, dep$SI_dep, xlab = "river depth bins", ylab = "SI river depth")
dep$bins = 1:length(dep[,1])
g = gam(abundance ~ s(bins), data = dep)
dep$abundance = predict(g, newdata = dep, scale = "response")
dep$SI_dep = ((dep$abundance - min(dep[,2]))/(max(dep[,2])-min(dep[,2])))
plot(dep$dep_bins, dep$SI_dep, xlab = "sample river bins", ylab = "SI river depth")


# SI Salinity
sal = df$SALINITY
sal = as.numeric(as.character(sal))
salinity_int = (classIntervals(sal, 20 ,style = "fisher"))
salinity_int[[2]][1] = salinity_int[[2]][1]-0.1
salinity_bins = (cut(sal, breaks = salinity_int$brks))
df$salinity_bins = cut(sal, breaks = salinity_int$brks)
salinity = aggregate(cpue ~ salinity_bins, data = df, FUN="mean")
colnames(salinity)[1] = "salinity_bins"
colnames(salinity)[2] = "abundance"
salinity$SI_sal = ((salinity$abundance - min(salinity[,2]))/(max(salinity[,2]) - min(salinity[,2])))
salinity_axis = salinity$salinity_bins
plot(salinity$salinity_bins, salinity$SI_sal, xlab = "sal bins", ylab = "SI sal")
salinity$bin = 1:length(salinity[,1])
g = gam(abundance ~ s(bin), data = salinity)
salinity$abundance = predict(g, newdata = salinity, scale = "response")
salinity$SI_sal = ((salinity$abundance - min(salinity[,2]))/(max(salinity[,2]) - min(salinity[,2])))
plot(salinity$salinity_bins, salinity$SI_sal, xlab = "sal bins", ylab = "SI sal")
#

## STEP 2- COMBINING INDIVIDUAL SUITABILITY INDICES FOR A COMPOSITE HSI =======
#
# Load annual environmental data files and create HSIs for each point #
## These MUST be saved in a different folder than your current working directory or they will be overwritten
setwd() #add path

wqs1988=read.csv("1988.csv")
wqs1989=read.csv("1989.csv")
wqs1990=read.csv("1990.csv")
wqs1991=read.csv("1991.csv")
wqs1992=read.csv("1992.csv")
wqs1993=read.csv("1993.csv")
wqs1994=read.csv("1994.csv")
wqs1995=read.csv("1995.csv")
wqs1996=read.csv("1996.csv")
wqs1997=read.csv("1997.csv")
wqs1998=read.csv("1998.csv")
wqs1999=read.csv("1999.csv")
wqs2000=read.csv("2000.csv")
wqs2001=read.csv("2001.csv")
wqs2002=read.csv("2002.csv")
wqs2003=read.csv("2003.csv")
wqs2004=read.csv("2004.csv")
wqs2005=read.csv("2005.csv")
wqs2006=read.csv("2006.csv")
wqs2007=read.csv("2007.csv")
wqs2008=read.csv("2008.csv")
wqs2009=read.csv("2009.csv")
wqs2010=read.csv("2010.csv")
wqs2011=read.csv("2011.csv")
wqs2012=read.csv("2012.csv")
wqs2013=read.csv("2013.csv")
wqs2014=read.csv("2014.csv")
wqs2015=read.csv("2015.csv")
wqs2016=read.csv("2016.csv")
wqs2017=read.csv("2017.csv")

x = list(wqs1988=wqs1988, wqs1989=wqs1989,
         wqs1990=wqs1990, wqs1991=wqs1991, wqs1992=wqs1992, wqs1993=wqs1993, wqs1994=wqs1994, wqs1995=wqs1995 ,wqs1996=wqs1996,wqs1997=wqs1997, wqs1998=wqs1998, wqs1999=wqs1999, 
         wqs2000=wqs2000, wqs2001=wqs2001, wqs2002=wqs2002, wqs2003=wqs2003, wqs2004=wqs2004, wqs2005=wqs2005, wqs2006=wqs2006, wqs2007=wqs2007, wqs2008=wqs2008, wqs2009=wqs2009,
         wqs2010=wqs2010, wqs2011=wqs2011, wqs2012=wqs2012, wqs2013=wqs2013, wqs2014=wqs2014, wqs2015=wqs2015, wqs2016=wqs2016, wqs2017=wqs2017)


# For each year, each geographic point is given an HSI value by combining SIs
setwd() #add path

# environmental variables given equal weights
for (i in 1:length(x)){
  
  tstemp = x[[i]]
  
  #SI_temperature
  temp = tstemp[, "Temperature"]
  temp = as.numeric(as.character(temp))
  tstemp$temperature_bins = cut(temp, breaks = temperature_int$brks)
  tstemp = merge(tstemp, temperature, by = "temperature_bins", all = TRUE)
  
  #SI_depth
  depth = tstemp[, "Depth"]
  depth = as.numeric(as.character(depth))
  tstemp$dep_bins = cut(depth, breaks = dep_int$brks)
  tstemp = merge(tstemp, dep, by = "dep_bins", all = TRUE)
  
  #SI_salinity
  salt = tstemp[, "Salinity"]
  salt = as.numeric(as.character(salt))
  tstemp$salinity_bins = cut(salt, breaks = salinity_int$brks)
  tstemp = merge(tstemp, salinity, by = "salinity_bins", all = TRUE)
  
  # SI_do
  doo = tstemp[, "Dissolved_Oxygen"]
  doo = as.numeric(as.character(doo))
  tstemp$DO_bins = cut(doo, breaks = DO_int$brks)
  tstemp = merge(tstemp, DO, by = "DO_bins", all = TRUE)
  
  #Create Composite HSIs at each location for each year
  tstemp$AMM_HSI = (1/4)(tstemp$SI_temp + tstemp$SI_dep + tstemp$SI_sal + tstemp$SI_DO) #weighting environmental variables equally
  tstemp = tstemp[,c("lon", "lat", "AMM_HSI")] 
  tstemp = tstemp[!is.na(tstemp$AMM_HSI),]
  write.csv(tstemp, paste(names(x)[i], ".csv", sep = ""),row.names=FALSE)
  
}


# Read in the files that we just created, median HSI is chosen across all years for each point & rasterized
setwd() # add path

wqs1988=read.csv("wqs1988.csv")
wqs1989=read.csv("wqs1989.csv")
wqs1990=read.csv("wqs1990.csv")
wqs1991=read.csv("wqs1991.csv")
wqs1992=read.csv("wqs1992.csv")
wqs1993=read.csv("wqs1993.csv")
wqs1994=read.csv("wqs1994.csv")
wqs1995=read.csv("wqs1995.csv")
wqs1996=read.csv("wqs1996.csv")
wqs1997=read.csv("wqs1997.csv")
wqs1998=read.csv("wqs1998.csv")
wqs1999=read.csv("wqs1999.csv")
wqs2000=read.csv("wqs2000.csv")
wqs2001=read.csv("wqs2001.csv")
wqs2002=read.csv("wqs2002.csv")
wqs2003=read.csv("wqs2003.csv")
wqs2004=read.csv("wqs2004.csv")
wqs2005=read.csv("wqs2005.csv")
wqs2006=read.csv("wqs2006.csv")
wqs2007=read.csv("wqs2007.csv")
wqs2008=read.csv("wqs2008.csv")
wqs2009=read.csv("wqs2009.csv")
wqs2010=read.csv("wqs2010.csv")
wqs2011=read.csv("wqs2011.csv")
wqs2012=read.csv("wqs2012.csv")
wqs2013=read.csv("wqs2013.csv")
wqs2014=read.csv("wqs2014.csv")
wqs2015=read.csv("wqs2015.csv")
wqs2016=read.csv("wqs2016.csv")
wqs2017=read.csv("wqs2017.csv")

x = list(wqs1988=wqs1988, wqs1989=wqs1989,
         wqs1990=wqs1990, wqs1991=wqs1991, wqs1992=wqs1992, wqs1993=wqs1993, wqs1994=wqs1994, wqs1995=wqs1995 ,wqs1996=wqs1996,wqs1997=wqs1997, wqs1998=wqs1998, wqs1999=wqs1999,
         wqs2000=wqs2000, wqs2001=wqs2001, wqs2002=wqs2002, wqs2003=wqs2003, wqs2004=wqs2004, wqs2005=wqs2005, wqs2006=wqs2006, wqs2007=wqs2007, wqs2008=wqs2008, wqs2009=wqs2009,
         wqs2010=wqs2010, wqs2011=wqs2011, wqs2012=wqs2012, wqs2013=wqs2013, wqs2014=wqs2014, wqs2015=wqs2015, wqs2016=wqs2016, wqs2017=wqs2017)


for (i in 1:length(x)){
  x[[i]] = x[[i]][, c("lon", "lat", "AMM_HSI")] #AMM_HSI
}

ts = plyr::join_all(x, by = c("lat", "lon"))

write.csv(ts, "compositive.unweighted.csv", row.names = F)
####

# To generate the HSI scores by weighting the environmental variables use this instead:
# environmental variables given equal weights
for (i in 1:length(x)){
  
  tstemp = x[[i]]
  
  #SI_temperature
  temp = tstemp[, "Temperature"]
  temp = as.numeric(as.character(temp))
  tstemp$temperature_bins = cut(temp, breaks = temperature_int$brks)
  tstemp = merge(tstemp, temperature, by = "temperature_bins", all = TRUE)
  
  #SI_depth
  depth = tstemp[, "Depth"]
  depth = as.numeric(as.character(depth))
  tstemp$dep_bins = cut(depth, breaks = dep_int$brks)
  tstemp = merge(tstemp, dep, by = "dep_bins", all = TRUE)
  
  #SI_salinity
  salt = tstemp[, "Salinity"]
  salt = as.numeric(as.character(salt))
  tstemp$salinity_bins = cut(salt, breaks = salinity_int$brks)
  tstemp = merge(tstemp, salinity, by = "salinity_bins", all = TRUE)
  
  # SI_do
  doo = tstemp[, "Dissolved_Oxygen"]
  doo = as.numeric(as.character(doo))
  tstemp$DO_bins = cut(doo, breaks = DO_int$brks)
  tstemp = merge(tstemp, DO, by = "DO_bins", all = TRUE)
  
  #Create Composite HSIs at each location for each year
  tstemp$AMM_HSI = (tstemp$SI_temp*0.321 + tstemp$SI_dep*0.107 + tstemp$SI_sal*0.433 + tstemp$SI_DO*0.138) #adding weights based on percent influence from the BRTs
  tstemp = tstemp[,c("lon", "lat", "AMM_HSI")] 
  tstemp = tstemp[!is.na(tstemp$AMM_HSI),]
  write.csv(tstemp, paste(names(x)[i], ".csv", sep = ""),row.names=FALSE)
  
}

# Read in the files that we just created, median HSI is chosen across all years for each point & rasterized
setwd() # add path

wqs1988=read.csv("wqs1988.csv")
wqs1989=read.csv("wqs1989.csv")
wqs1990=read.csv("wqs1990.csv")
wqs1991=read.csv("wqs1991.csv")
wqs1992=read.csv("wqs1992.csv")
wqs1993=read.csv("wqs1993.csv")
wqs1994=read.csv("wqs1994.csv")
wqs1995=read.csv("wqs1995.csv")
wqs1996=read.csv("wqs1996.csv")
wqs1997=read.csv("wqs1997.csv")
wqs1998=read.csv("wqs1998.csv")
wqs1999=read.csv("wqs1999.csv")
wqs2000=read.csv("wqs2000.csv")
wqs2001=read.csv("wqs2001.csv")
wqs2002=read.csv("wqs2002.csv")
wqs2003=read.csv("wqs2003.csv")
wqs2004=read.csv("wqs2004.csv")
wqs2005=read.csv("wqs2005.csv")
wqs2006=read.csv("wqs2006.csv")
wqs2007=read.csv("wqs2007.csv")
wqs2008=read.csv("wqs2008.csv")
wqs2009=read.csv("wqs2009.csv")
wqs2010=read.csv("wqs2010.csv")
wqs2011=read.csv("wqs2011.csv")
wqs2012=read.csv("wqs2012.csv")
wqs2013=read.csv("wqs2013.csv")
wqs2014=read.csv("wqs2014.csv")
wqs2015=read.csv("wqs2015.csv")
wqs2016=read.csv("wqs2016.csv")
wqs2017=read.csv("wqs2017.csv")

x = list(wqs1988=wqs1988, wqs1989=wqs1989,
         wqs1990=wqs1990, wqs1991=wqs1991, wqs1992=wqs1992, wqs1993=wqs1993, wqs1994=wqs1994, wqs1995=wqs1995 ,wqs1996=wqs1996,wqs1997=wqs1997, wqs1998=wqs1998, wqs1999=wqs1999,
         wqs2000=wqs2000, wqs2001=wqs2001, wqs2002=wqs2002, wqs2003=wqs2003, wqs2004=wqs2004, wqs2005=wqs2005, wqs2006=wqs2006, wqs2007=wqs2007, wqs2008=wqs2008, wqs2009=wqs2009,
         wqs2010=wqs2010, wqs2011=wqs2011, wqs2012=wqs2012, wqs2013=wqs2013, wqs2014=wqs2014, wqs2015=wqs2015, wqs2016=wqs2016, wqs2017=wqs2017)


for (i in 1:length(x)){
  x[[i]] = x[[i]][, c("lon", "lat", "AMM_HSI")] #AMM_HSI
}

ts = plyr::join_all(x, by = c("lat", "lon"))

write.csv(ts, "compositive.weighted.csv", row.names = F)

