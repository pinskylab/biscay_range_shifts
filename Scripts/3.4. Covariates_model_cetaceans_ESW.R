###########################################################################
#Density covariates >> this should be done after running a simple model
########################################################################

#https://climatedataguide.ucar.edu/climate-data/hurrell-north-atlantic-oscillation-nao-index-station-based


setwd("C:/use/PREDOC/Scripts/0ChapterIII/VAST/Output/New")

library(raster)
library(fields)
library(abind)

#Load the object
del<-readRDS("C:/use/PREDOC/Scripts/0ChapterIII/VAST/Output/final_del.RDS")
dim(del) #123945
del<-del[del$Year>=1992 & del$Year<2018,] #123596

#Our grid
grid<-readRDS("C:/use/PREDOC/Scripts/0ChapterIII/VAST/grid_01.RDS")
plot(grid$Lon,grid$Lat)
dim(grid) #10078
head(grid)

#Load our model
#Same resolution as the grid 
fit1<-readRDS("density_lognormal_grid_4444.RDS")
summary(fit1)
df<-fit1$data_frame #123596

optimize.args<-fit1[["optimize_args"]]
parameters.estimates<-fit1[["parameter_estimates"]]
AIC<-parameters.estimates[["AIC"]]
AIC # 20690.59

extrapolation_list<-fit1["extrapolation_list"]
data_extrap<-extrapolation_list$extrapolation_list["Data_Extrap"]
grid_extrap<-data_extrap$Data_Extrap    #10078
dim(grid_extrap)

model.args<-fit1[["model_args"]]
z_gm<-model.args[[1]]["Z_gm"]
z_gm<-as.data.frame(z_gm[[1]])
dim(z_gm) #10078

par(mfrow=c(1,2))
plot(grid_extrap$E_km,grid_extrap$N_km) #the whole BOB
plot(z_gm$E_km,z_gm$N_km)


#COVARIATES: required in the grid and in the observations
#------------------------------------------------------------
#############
#In our grid
#############

X_gtp<-model.args[[1]]["X_gtp"] #How it must be
str(X_gtp)
class(X_gtp[[1]])

xy<-grid[,1:2] #10078

#TEMP
temp.grid<-NULL
detach(package:tidyr)

yy<-c(1992:2017) #no data before 1992

for (i in 1:length(yy)){
  #i=2
  ras<-readRDS(paste("C:/use/PREDOC/Scripts/0ChapterIII/VAST/environmental_variables/Copernicus/resampled/SST10m_",yy[i],sep="","_.RDS_TRUE"))
  temp<-extract(ras, xy)
  temp.grid<-cbind(temp.grid,temp)
  colnames(temp.grid)[ncol(temp.grid)] <- paste(yy[i])
  
}

dim(temp.grid)
temp.grid<-as.data.frame(temp.grid)

#To use later
temp.grid$Lon<-grid$Lon
temp.grid$Lat<-grid$Lat

#Testing
print(yy[i])
plot(ras)
quilt.plot(temp.grid$Lon,temp.grid$Lat,temp.grid$`2017`)

#Na's
table(is.na(temp.grid)) #3744
temp.grid$sum<-rowSums(temp.grid[,1:26], na.rm=T)
table(temp.grid$sum==0) #144*26
temp.grid<-temp.grid[temp.grid$sum!=0,]
table(is.na(temp.grid))
temp.grid$sum<-NULL


#CHL
chl.grid<-NULL
detach(package:tidyr)

yy<-c(1992:2017) #no data before 1992

for (i in 1:length(yy)){
  #i=2
  ras<-readRDS(paste("C:/use/PREDOC/Scripts/0ChapterIII/VAST/environmental_variables/Copernicus/Chl_log_10m_",yy[i],sep="","_.RDS_TRUE"))
  chl<-extract(ras, xy)
  chl.grid<-cbind(chl.grid,chl)
  colnames(chl.grid)[ncol(chl.grid)] <- paste(yy[i])
  
}

plot(ras)
dim(chl.grid)
chl.grid<-as.data.frame(chl.grid)

#Checking
plot(ras)
quilt.plot(xy$Lon,xy$Lat,chl.grid$`2017`)

#Na's
table(is.na(chl.grid)) #3744
chl.grid<-as.data.frame(chl.grid)
chl.grid$sum<-rowSums(chl.grid, na.rm=T)
table(chl.grid$sum==0) #144*26
chl.grid<-chl.grid[chl.grid$sum!=0,]
table(is.na(chl.grid))
chl.grid$sum<-NULL

#NAO
library(dplyr)
library(tidyr)

nao<-read.csv("C:/use/PREDOC/Scripts/0ChapterIII/VAST/environmental_variables/NAO.csv",header=TRUE,sep=";")
colnames(nao)[1]<-c("NAO")
nao<-nao %>% separate("NAO", c("Year", "Value"),sep="  ")
nao<-nao[nao$Year>=1992 & nao$Year<2018,]
dim(nao)
head(nao)

rownames(nao)<-nao[,1]
nao$Year<-NULL
nao$Value<-as.numeric(nao$Value)

tnao<-t(nao)
tnao<-as.data.frame(tnao)
colnames(tnao)[1:length(tnao)]<-rownames(nao)
nao.grid<-tnao[rep(seq_len(nrow(tnao)), nrow(xy)), ]

#Removing from the grid
dim(temp.grid) #9934
dim(temp.grid)
dim(nao.grid)
dim(grid) #10078


plot(temp.grid$Lon,temp.grid$Lat)
plot(grid$Lon,grid$Lat)

new.grid<-semi_join(grid, temp.grid, by = c("Lon", "Lat")) #9934
plot(new.grid$Lon,new.grid$Lat)
dim(new.grid)

temp.grid$Lon<-NULL
temp.grid$Lat<-NULL

#Array
nao.grid<-nao.grid[1:9934,]
naoX_gtp<-array(as.matrix(nao.grid), dim=c(nrow(new.grid),length(yy),1))
n<-naoX_gtp[,,1] #repasar

cov<-as.matrix(abind(temp.grid,chl.grid))
dim(cov)
str(cov)

X_gtp<-array(cov, dim=c(nrow(cov),length(yy),2))
dimnames(X_gtp)[[2]]<-colnames(cov)[1:26]
dimnames(X_gtp)[[3]]<-c("Temp","Chl")
str(X_gtp)

a<-X_gtp[,,1] #compare to temp.grid
b<-X_gtp[,,2] #compare to chl.grid

#Quadratic
X_gtp<-abind(X_gtp, X_gtp[,,1,drop=FALSE]^2, along=3 ) #temp
X_gtp<-abind(X_gtp, X_gtp[,,2,drop=FALSE]^2, along=3 ) #chl

dim(X_gtp) #rows, years, cov
class(X_gtp)
str(X_gtp)

dimnames(X_gtp)[[3]][3:4]<-c("TempSq","ChlSq")
str(X_gtp)

t<-X_gtp[,,1]
tt<-X_gtp[,,3]

c<-X_gtp[,,2]
cc<-X_gtp[,,4]

#Years
y<-X_gtp[1,,1]
y

####################
#In our observations
####################

X_itp<-model.args[[1]]["X_itp"]
str(X_itp) #123596
dim(del) #123596


cov1<-NULL
yy<-c(1992:2017)
detach(package:tidyr)

#TEMP
for (i in 1:length(yy)){
  #i=2
  print(yy[i])
  subset<-del[del$Year==yy[i],]
  xy<-subset[,c("Lon","Lat")]
  ras<-readRDS(paste("C:/use/PREDOC/Scripts/0ChapterIII/VAST/environmental_variables/Copernicus/SST10m_",yy[i],sep="","_.RDS_TRUE"))
  temp<-extract(ras, xy)
  new<-cbind(subset,temp)
  cov1<-rbind(cov1, new)
}

dim(cov1) #123596
table(is.na(cov1$temp)) #2234


#CHL
new<-NULL
cov2<-NULL

for (i in 1:length(yy)){
  #i=1
  print(yy[i])
  subset<-cov1[cov1$Year==yy[i],]
  xy<-subset[,c("Lon","Lat")]
  ras<-readRDS(paste("C:/use/PREDOC/Scripts/0ChapterIII/VAST/environmental_variables/Copernicus/Chl_log_10m_",yy[i],sep="","_.RDS_TRUE"))
  chl<-extract(ras, xy)
  new<-cbind(subset,chl)
  cov2<-rbind(cov2, new)
  
}

plot(ras)
dim(cov2)
table(is.na(cov2$chl)) ##2234

#NAO
library(dplyr)
library(tidyr)

nao<-read.csv("C:/use/PREDOC/Scripts/0ChapterIII/VAST/environmental_variables/NAO.csv",header=TRUE,sep=";")
colnames(nao)[1]<-c("NAO")
nao<-nao %>% separate("NAO", c("Year", "Value"),sep="  ")
nao<-nao[nao$Year>=1992 & nao$Year<2018,]
dim(nao)


cov3<-NULL
for (i in 1:length(yy)){
  #i=1
  print(yy[i])
  subset<-cov2[cov2$Year==yy[i],]
  subset$nao<-paste(nao$Value[nao$Year==yy[i]])
  #ras<-readRDS(paste("C:/use/PREDOC/Scripts/0ChapterIII/VAST/environmental_variables/Copernicus/Chl10m_",yy[i],sep="","_.RDS_TRUE"))
  #chl<-extract(ras, xy)
  cov3<-rbind(cov3, subset)
  
}

dim(cov3)
str(cov3$nao)
cov3$nao<-as.numeric(as.character(cov3$nao))
length(unique(cov3$nao))
table(is.na(cov3$nao)) ##0

del<-cov3
del$sum<-rowSums(del[,c("temp", "chl")], na.rm=T)
table(del$sum==0)

plot(del$Lon[del$sum==0],del$Lat[del$sum==0])
points(del$Lon[del$sum==0 & del$Density>0], del$Lat[del$sum==0& del$Density>0], col="red", pch=19)

del<-del[del$sum!=0,] 
table(is.na(del$temp)) 
del$sum<-NULL
dim(del) #121362

quilt.plot(del$Lon,del$Lat,del$temp)
quilt.plot(del$Lon,del$Lat,del$chl)
quilt.plot(del$Lon,del$Lat,del$nao)

#Array

naoX_itp = aperm(outer(as.matrix(del$nao), rep(1,length(yy))), perm=c(1,3,2) )
nn<-naoX_itp[,,1]
  
cov<-del[,c("temp", "chl")]
head(cov)
class(cov)
cov.grid<-as.matrix(cov)

# Expand to expected size of input
X_itp = aperm( outer(cov.grid, rep(1,length(yy))), perm=c(1,3,2) )
str(X_itp)
class(X_itp)
dim(X_itp)

# Add quadratic temperature for demonstration purposes
X_itp = abind( X_itp, "Temp.sq"=X_itp[,,'temp',drop=FALSE]^2, along=3 )
X_itp = abind( X_itp, "chl.sq"=X_itp[,,'chl',drop=FALSE]^2, along=3 )

dimnames(X_itp)[[3]][3:4]<-c("TempSq","ChlSq")
str(X_itp)
dim(X_itp)

t<-X_itp[,,1]
head(t)
tt<-X_itp[,,3]
head(tt)

c<-X_itp[,,2]
head(c)
max(c); min(c)
cc<-X_itp[,,4]
head(cc)
max(cc);min(cc)




############################################################################################################################
#                                  Fitting the model: Density and covariates
############################################################################################################################


Sys.setenv(PATH="%PATH%;C:/Rtools/gcc-4.6.3/bin;c:/Rtools/bin")
Sys.getenv("PATH")

library(devtools)
library(VAST)
library(TMB)

# Make settings
settings = make_settings(n_x=25, Region="User", purpose="index", strata.limits=data.frame(STRATA="All_areas"),ObsModel=c(1,1))
#settings = make_settings(n_x=50, Region="Other", purpose="index", strata.limits=data.frame(STRATA="All_areas"),ObsModel=c(1,1))
#settings$FieldConfig=c(1,1,1,1)
#settings$RhoConfig=c(2,2,2,2)

settings$FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1)
settings$RhoConfig = c("Beta1"=4, "Beta2"=4, "Epsilon1"=4, "Epsilon2"=4)
#OverdispersionConfig = c("Eta1"=0, "Eta2"=0)

#RhoConfig[1] and RhoConfig[2] #spatio
#RhoConfig[3] and RhoConfig[4] #spatiotemporal

###0: fixed effects
###1: random following IID distribution
###2: random following a random walk
###3: constant
###4: each year as random following AR1 process

#-------------------------------------------------------------------------------------------------
#Lognormal:  Obsmodel(1,1) #the model is likely not converged (not anymore) >> AIC 21194.79
#Gamma:  Obsmodel(2,1): success >> AIC 21508.59

#Conventional delta model (0,2): error
#Error in nlminb(start = startpar, objective = fn, gradient = gr, control = nlminb.control,  : 
#gradient function must return a numeric vector of length 67
#In addition: Warning message:
#In nlminb(start = startpar, objective = fn, gradient = gr, control = nlminb.control,  :
#NA/NaN function evaluation

#Poisson link delta model (1,2): error
#Error in nlminb(start = startpar, objective = fn, gradient = gr, control = nlminb.control,  : 
#gradient function must return a numeric vector of length 67
#In addition: Warning message:
#In nlminb(start = startpar, objective = fn, gradient = gr, control = nlminb.control,  :
#NA/NaN function evaluation

#Tweedie:  Obsmodel(10,2): #Error in sparseHessianFun(env, skipFixedEffects = skipFixedEffects) : 
#Memory allocation fail in function 'MakeADHessObject2'

#Negative binomial (5,) #Error in nlminb(start = startpar, objective = fn, gradient = gr, control = nlminb.control,  : 
#gradient function must return a numeric vector of length 68

#Normal: Obsmodel(0,1): error
#---------------------------------------------------------------------------------------------------

# Renumber years to have same indexing as covariates
del$t_i = match(del$Year, sort(unique(del$Year)) )
str(t_i)

#New grid
dim(new.grid)
str(new.grid)
new.grid$Area_km2<-as.numeric(new.grid$Area_km2)

# Make Xconfig_zcp, using arbitrary decisions for illustration
# By default, linear effect for each variable
Xconfig_zcp = array( 1, dim=c(2,1,dim(X_gtp)[3]) )
dim(Xconfig_zcp) #2 rows, 1 column, 2 dim

# Turn off for all covariate effects on 1st linear predictor
#Xconfig_zcp[1,,] = 0


# Use spatially varying effect for "Nao" for 2nd predictor
Xconfig_zcp_nao = array( 3, dim=c(2,1,dim(naoX_gtp)[3]) )
#Xconfig_zcp_nao[2,1,1] = 3


#
# Fit model
density_esw = fit_model( "settings"=settings, 
                         "Lat_i"=del$Lat,  "Lon_i"=del$Lon, #not UTM
                         "t_i"=del$Year, 
                         "c_i"=rep(0,nrow(del)),
                         "b_i"=del$Density,
                         "a_i"=rep(1,nrow(del)),
                         #"Q_ik"=as.matrix(del$Month), 
                         "X_gtp"=array(X_gtp[,,c(2,4)],dim=c(nrow(X_gtp),ncol(X_gtp),2)),
                         #"X_itp"=array(X_itp[,,c(2,4)],dim=c(nrow(X_itp),ncol(X_itp),2)), 
                         #"Xconfig_zcp"=array(Xconfig_zcp[,,2],dim=c(2,1,2)),
                         #"X_gtp"=naoX_gtp,
                         #"X_itp"=naoX_itp,
                         #"X_config_zcp"=Xconfig_zcp_nao,
                         #extrapolation_args = list(observations_LL=del[,c("Lon","Lat")]), 
                         extrapolation_args = list(input_grid=as.data.frame(new.grid)),
                         bias.correct=TRUE) #FALSE when there is no space

setwd("C:/use/PREDOC/Scripts/0ChapterIII/VAST/Output/New")
getwd()
saveRDS(density_esw ,file="density_lognormal_grid_4444_nao.RDS")
saveRDS(settings,file="settings_lognormal_grid_4444_nao.new.RDS")

setwd("C:/use/PREDOC/Scripts/0ChapterIII/VAST/Output/New")
fit1<-readRDS("density_lognormal_grid_4444_nao.RDS") #from TOSHIBA
summary(fit1)
data.frame<-fit1["data_frame"] #121362
a<-data.frame$data_frame$a_i #area
v<-data.frame$data_frame$v_i #vessel effect
b<-data.frame$data_frame$b_i #density

extrapolation_list<-fit1["extrapolation_list"]
data_extrap<-extrapolation_list$extrapolation_list["Data_Extrap"]
grid_extrap<-data_extrap$Data_Extrap    #9934 rows
grid_extrap$id<-rownames(grid_extrap)
plot(grid_extrap$E_km,grid_extrap$N_km) #the whole BOB
plot(grid_extrap$Lon,grid_extrap$Lat)   

spatial_list<-fit1["spatial_list"]
extrapolation_args<-fit1[["extrapolation_args"]]

model.args<-fit1[["model_args"]]
Field<-model.args[[1]]["FieldConfig"]
Rho<-model.args[["RhoConfig"]]

obs_model<-model.args[[1]]["ObsModel_ez"] #Lognormal
x_gtp<-model.args[[1]]["X_gtp"]
x1<-x_gtp[[1]][,,1]
max(x1)
x2<-x_gtp[[1]][,,2]
max(x2)

x_itp<-model.args[[1]]["X_itp"]
z1<-x_itp[[1]][,,1]
max(z1)
z2<-x_itp[[1]][,,2]
max(z2)

n_g<-model.args[[1]]["n_g"] #9934; number of extrapolation grid cells when fine_scale=T
n_t<-model.args[[1]]["n_t"] #26
a_gl<-model.args[[1]]["a_gl"] #121 km2 ? #area of each extrapolation grid cell
n_s<-model.args[[1]]["n_s"][1] #?
z_gm<-model.args[[1]]["Z_gm"]
z_gm<-as.data.frame(z_gm[[1]])
z_gm$id<-rownames((z_gm))
plot(z_gm$E_km,z_gm$N_km)

optimize.args<-fit1[["optimize_args"]]
parameters.estimates<-fit1[["parameter_estimates"]]
AIC<-parameters.estimates[["AIC"]]
AIC 


# Plot results

plot_results(settings=settings, fit=fit1,
             working_dir = paste0("D:/"),use_biascorr = TRUE) # "C:/use/PREDOC/Scripts/0ChapterIII/VAST/Output/new"),

#plot( fit=fit1, plot_set=c(3,13,14) )

Extrapolation_List<-fit1["extrapolation_list"]

MapDetails_List<-make_map_info( Region=settings$Region,Extrapolation_List=fit1$extrapolation_list, 
                                spatial_list=fit1$spatial_list,
                                NN_Extrap=fit1$spatial_list$PolygonList$NN_Extrap,
                                fine_scale=fit1$spatial_list$fine_scale) 

MapDetails<-MapDetails_List[["MappingDetails"]]
#Year_Set = seq(min(del[,'Year']),max(del[,'Year']))
#Years2Include = which( Year_Set %in% sort(unique(del[,'Year'])))

#setwd("C:/use/PREDOC/Scripts/0ChapterIII/VAST/Output/new/Legend")
setwd("D:/")
plot_results(settings=settings, fit=fit1,working_dir = paste0("D:/"),
               #paste0("C:/use/PREDOC/Scripts/0ChapterIII/VAST/Output/new/Legend"),
             plot_maps(Legend=MapDetails_List[["Legend"]], 
                       Report=fit1$Report, PlotDF=MapDetails_List[["PlotDF"]],
                       MappingDetails = MapDetails, Xlim=MapDetails_List[["Xlim"]], 
                       Ylim = MapDetails_List[["Ylim"]]
                       #Year_Set = seq(min(del[,'Year']),max(del[,'Year']),
                       #Years2Include = which( Year_Set %in% sort(unique(del[,'Year'])))
             ))
