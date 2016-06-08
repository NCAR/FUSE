rm(list=ls())

# load netcdf package
require(ncdf4) 

# set paths
dir_fuse<-'/Volumes/hydro-c1/fuse/'
dir_input<-paste(dir_fuse,'input/',sep='')
dir_output<-paste(dir_fuse,'output/',sep='')
dir_settings<-paste(dir_fuse,'settings/',sep='')

basin_id<-'08013000'

# load content of FUSE input file
fuse_input<-read.table(paste(dir_input,'input_',basin_id,'.txt',sep=''),sep='',skip=2,header=TRUE)
d_input<-as.Date(paste(fuse_input$YYYY,sprintf('%02s',fuse_input$MM),sprintf('%02s',fuse_input$DD),sep=''),'%Y%m%d') # date of each (daily) time step
summary(d_input)

# load FUSE output
fuse_output_ncfile<-paste(dir_output,'URS_MOPEX__',basin_id,'__FUSE_070__NMETH-2_SSTEP-0__2-0__1.e-2-1.e-2__1.0000000000.nc',sep='')
fuse_output<-nc_open(fuse_output_ncfile) 
q_obs<-ncvar_get(fuse_output,'obsq')
q_sim<-ncvar_get(fuse_output,'q_routed')
# d_sim<-ncvar_get(fuse_output,'tim') # FUSE does not seem to save a time vector...

# ... so load settings and input to retrive time
forcing_info_file<-paste(dir_settings,'forcinginfo.',basin_id,'.txt',sep='')
sim_periods_details<-as.numeric(strsplit(readLines(forcing_info_file)[3],' ')[[1]][2:4]) # start index, end of warm-up, end index
if((sim_periods_details[3]-sim_periods_details[1]+1)!=dim(q_sim)[1]){'The length of the simulation does not correspond to the settings'}

d_output<-d_input[(sim_periods_details[2]+1):sim_periods_details[3]] # time steps covered by the simulation after the spinup period
summary(d_output)

ind_sim<-((sim_periods_details[2]+1):sim_periods_details[3])-sim_periods_details[1] # time steps of the simulations after the spinup period
summary(d_output[ind_sim])

# plot results for the period after spinup
plot(d_output[ind_sim],q_obs[ind_sim,1],xlab='',ylab='Discharge [mm/day]',type='n')
apply(q_sim[ind_sim,],2,lines,col='orange',x=d_output[ind_sim])
lines(d_output[ind_sim],q_obs[ind_sim,1],type='l')

