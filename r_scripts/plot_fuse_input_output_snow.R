rm(list=ls())

require(ncdf4)

### 1. SET NSE FUNCTION
compute_nse<-function(obs,sim){
  nse<-1-sum((sim-obs)^2)/sum((obs-mean(obs))^2)
  return(nse)
}

### 2. SET DIR
dir_home<-'/Volumes/hc1_home/'
dir_settings=paste(dir_home,'fuse/settings/',sep='')
dir_input=paste(dir_home,'fuse/settings/',sep='')
dir_output=paste(dir_home,'fuse/output/',sep='')

### 3. DEFINE CATCHMENT AND MODEL SETUP 
id<-'HETCHHETCHY' # catchment id
mod_setup<-paste('URS_MOPEX__',id,'__FUSE_1408__RFERR-2_ARCH1-3_ARCH2-1_QSURF-1_QPERC-1_ESOIL-2_QINTF-1_Q_TDH-1_SNOWM-2_NMETH-2_SSTEP-0__2-0__1.e-2-1.e-2__1.0000000000',sep='') # model setup - string part of the name of the output file

### 4. LOAD FUSE INPUT DATA
input_data<-read.table(paste(dir_input,'f_hetchhetchy.dat',sep=''),sep='',skip=4,header=TRUE)
t_input<-as.Date(apply(cbind(input_data$YYYY,sprintf('%02d',input_data$MM),sprintf('%02d',input_data$DD)),1,paste,collapse=''),'%Y%m%d') # convert string to date

### 5. LOAD FORCING INFO DATA
forcing_info_file<-paste(dir_settings,'forcinginfo.',id,'.txt',sep='')
sim_periods_details<-readLines(forcing_info_file)[3] # load third line
sim_periods_details_split<-strsplit(sim_periods_details,' ')[[1]] # split the line
index_start_sim<-as.numeric(sim_periods_details_split[2]) # retrive index of simulation start and convert it to numeric
index_end_warmup<-as.numeric(sim_periods_details_split[3])
index_end_sim<-as.numeric(sim_periods_details_split[4])

i_sim_input<-index_start_sim:index_end_sim # determine the indices correponding to the simulation period for the input data
i_sim_output<-(index_start_sim:index_end_sim)-index_start_sim+1 # same for the output data (shorter than the input data)

i_eval_input<-(index_end_warmup+1):index_end_sim # determine the indices correponding to the evaluation period (i.e. after spinup)
i_eval_output<-((index_end_warmup+1):index_end_sim)-index_start_sim+1 # same for the output data

t_sim<-t_input[i_sim_input] # dates used for the simulation
summary(t_sim)
t_eval<-t_input[i_eval_input] # dates used for the evaluation
summary(t_eval)

### 6. LOAD FUSE URS OUTPUT
fuse_urs<-paste(dir_output,mod_setup,'.nc',sep='') 

nc_urs<-nc_open(fuse_urs)               # open NedCDF file
q_obs<-ncvar_get(nc_urs,'obsq')         # load observed discharge
urs_q_sim<-ncvar_get(nc_urs,'q_routed') # load simulated discharge
urs_swe_sim<-ncvar_get(nc_urs,'swe_z01') # load simulated discharge
urs_d_sim<-ncatt_get(nc_urs,'tim') # 

t_thres<-ncvar_get(nc_urs,'PXTEMP') 
dd_min<-ncvar_get(nc_urs,'MFMIN') 
dd_max<-ncvar_get(nc_urs,'MFMAX') 

# check consistency between the input data, forcing info data and model ouput
if(length(t_sim)!=dim(urs_q_sim)[1]){stop('Unexpected number of time steps in FUSE output file')}
#if(any(input_data$OBS_Q[i_sim_input]!=q_obs[,1])){stop('Observed discharge from the input and output file do not correspond')} # problem here, there is 3-day shift between the two time series

# plot discharge simulations
plot(t_eval,q_obs[i_eval_output,1],type='n',xlab='',ylab='Discharge [mm/d]') # prepare plot
apply(urs_q_sim[i_eval_output,], 2, lines, x=t_eval,col='orange')
lines(t_eval,q_obs[i_eval_output,1])

apply(urs_q_sim[i_eval_output,], 2,compute_nse,obs=q_obs[i_eval_output,1]) # compute NSE
ncvar_get(nc_urs,'nash_sutt') 

# plot SWE simulations
plot(t_eval,urs_swe_sim[i_eval_output,1],type='n',xlab='',ylab='Snow pack [mm]') # prepare plot
apply(urs_swe_sim[i_eval_output,], 2, lines, x=t_eval,col='orange')

### 7. LOAD SCE PARAM FILE
sce_setup<-'SCE_5000_3_1.e-3' # SCE setup - string part of the name of the SCE output file

fuse_sce_params<-paste(dir_output,'sce_params_',mod_setup,'_',sce_setup,'.nc',sep='') 
nc_sce_params<-nc_open(fuse_sce_params)         # open NetCDF file
nse_sce<-ncvar_get(nc_sce_params,'nash_sutt')   # load NSE
rmse_sce<-ncvar_get(nc_sce_params,'raw_rmse')    # load RMSE

# plot(nse_sce,ylab='NSE') # check that overall increase of NSE with the number of model runs
                           # note that if more than one calibration is performed, the results
                           # are concatenated, leading to a "bridge-like" pattern.

### 8. LOAD CALIBRATED FUSE SIMULATIONS
fuse_sce_modrun<-paste(dir_output,'sce_modrun_',mod_setup,'_',sce_setup,'.nc',sep='')
nc_sce_modrun<-nc_open(fuse_sce_modrun)             # open NedCDF file
sce_q_sim<-ncvar_get(nc_sce_modrun,'q_routed')      # load simulated discharge 
lines(t_eval,sce_q_sim[i_eval_output],col='blue1')  # plot simulated discharge

compute_nse(obs=q_obs[i_eval_output,1],sim=sce_q_sim[i_eval_output]) # compute NSE

nse_sce[which.min(rmse_sce)]
nse_sce[length(nse_sce)]
