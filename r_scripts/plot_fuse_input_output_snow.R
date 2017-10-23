rm(list=ls())

require(ncdf4)

### 1. SET NSE FUNCTION
compute_nse<-function(obs,sim){
  nse<-1-sum((sim-obs)^2)/sum((obs-mean(obs))^2)
  return(nse)
}

### 2. SET DIR
dir_plots<-'/glade/u/home/naddor/plots/'
dir_home<-'/glade/scratch/naddor/fuse_case_study_1/'
dir_settings=paste(dir_home,'settings/',sep='')
dir_input=paste(dir_home,'settings/',sep='')
dir_output=paste(dir_home,'output/',sep='')

### 3. DEFINE CATCHMENT AND MODEL SETUP
dom_id<-'us_09066300' # catchment id
fuse_id<-'902'

### LOAD data

fuse_mode_list<-c('runs_def','runs_best')

fuse_qobs<-list()
fuse_qsim<-list()
fuse_swe<-list()
fuse_sm<-list()
fuse_et<-list()

fuse_d<-list()

for(fuse_mode in fuse_mode_list){

  fuse_runs_def<-paste0(dir_output,dom_id,'_',fuse_id,'_',fuse_mode,'.nc') # file name
  nc_runs<-nc_open(fuse_runs_def)                # open NedCDF file

  # load variables
  fuse_qobs[[fuse_mode]]<-ncvar_get(nc_runs,'obsq')      # load observed discharge
  fuse_qsim[[fuse_mode]]<-ncvar_get(nc_runs,'q_routed')  # load simulated discharge
  fuse_swe[[fuse_mode]]<-ncvar_get(nc_runs,'swe_tot')     # load SWE
  fuse_et[[fuse_mode]]<-ncvar_get(nc_runs,'evap_1')+ncvar_get(nc_runs,'evap_2')  # load ET
  fuse_sm[[fuse_mode]]<-ncvar_get(nc_runs,'watr_1')+ncvar_get(nc_runs,'watr_2')  # load soil moisture

  # load time
  d_origin_raw<-ncatt_get(nc_runs,'time')$units
  d_origin_split<-strsplit(d_origin_raw,' ')[[1]]
  if(d_origin_split[1]!='days'){stop('Unexpected time format')}
  if(d_origin_split[2]!='since'){stop('Unexpected time format')}
  d_origin<-as.Date(d_origin_split[3])
  fuse_d[[fuse_mode]]<-d_origin+ncvar_get(nc_runs,'time')

  nc_close(nc_runs)

}

# PLOT RESULT
i2plot<-fuse_d[['runs_def']]>=as.Date('2010-10-01')&fuse_d[['runs_def']]<=as.Date('2014-09-30') # indices to plot

col_obs<-'black' # colors
col_def<-'steelblue3'
col_sce<-'orange'

# COMPUTE NSE
nse_def<-compute_nse(fuse_qobs[['runs_def']][i2plot],fuse_qsim[['runs_def']][i2plot])
nse_sce<-compute_nse(fuse_qobs[['runs_best']][i2plot],fuse_qsim[['runs_best']][i2plot])

pdf(paste0(dir_plots,'fuse_test.pdf'),8,8)

par(mfrow=c(2,2),mar=c(4,4,3,0.5))

plot(fuse_d[['runs_def']][i2plot],fuse_qobs[['runs_def']][i2plot],xlab='',col=col_obs,'l',main='Discharge',ylab='[mm/day]',ylim=c(0,15))
lines(fuse_d[['runs_def']][i2plot],fuse_qsim[['runs_def']][i2plot],col=col_def)
lines(fuse_d[['runs_best']][i2plot],fuse_qsim[['runs_best']][i2plot],col=col_sce)

legend('top',col=c(col_obs,col_def,col_sce),c('OBS',
                                              paste0('SIM-DEF - NSE = ',round(nse_def,2)),
                                              paste0('SIM-SCE - NSE =',round(nse_sce,2))),lty=1,bty='n')

plot(fuse_d[['runs_def']][i2plot],fuse_swe[['runs_def']][i2plot],xlab='',col=col_def,'l',main='SWE',ylab='[mm]',
      ylim=c(min(unlist(fuse_swe)),max(unlist(fuse_swe))))
lines(fuse_d[['runs_best']][i2plot],fuse_swe[['runs_best']][i2plot],col=col_sce)

plot(fuse_d[['runs_def']][i2plot],fuse_et[['runs_def']][i2plot],xlab='',col=col_def,'l',main='ET',ylab='[mm/day]',
      ylim=c(min(unlist(fuse_et)),max(unlist(fuse_et))))
lines(fuse_d[['runs_best']][i2plot],fuse_et[['runs_best']][i2plot],col=col_sce)

plot(fuse_d[['runs_def']][i2plot],fuse_sm[['runs_def']][i2plot],xlab='',col=col_def,'l',main='Soil moisture',ylab='[mm]',
      ylim=c(min(unlist(fuse_sm)),max(unlist(fuse_sm))))
lines(fuse_d[['runs_best']][i2plot],fuse_sm[['runs_best']][i2plot],col=col_sce)

dev.off()
