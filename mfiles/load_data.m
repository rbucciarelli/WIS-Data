%% load_data.m
%- Load data from FRF and NCEP servers into Matlab arrays
%- 



%% Function: get_NCEP
function [wtime hs tp dp] = get_WIS(year,month,field)
 
%Field is the type of variable you want : geopotential_height,temperature,u,v,w
%Dimension of output value is (time,pressure_level,latitude,longitude)        
 
url = 'https://chlthredds.erdc.dren.mil/thredds/fileServer/wis/Pacific/ST46219/2000/WIS-ocean_waves_ST46219_200001.nc';
 
%List the files in selected folder
air_list=dir([folder,'air.',num2str(year,'%04i'),'*']);
hgt_list=dir([folder,'hgt.',num2str(year,'%04i'),'*']);
omega_list=dir([folder,'omega.',num2str(year,'%04i'),'*']);
uwnd_list=dir([folder,'uwnd.',num2str(year,'%04i'),'*']);
vwnd_list=dir([folder,'vwnd.',num2str(year,'%04i'),'*']);
 
if(strcmp(field,'temperature'))
filename=air_list.name;
filepath=[folder,filename];
ncid=netcdf.open(filepath,'NC_NOWRITE');
varid=netcdf.inqVarID(ncid,'air');                
end  
 
if(strcmp(field,'geopotential_height'))
filename=hgt_list.name;
filepath=[folder,filename];
ncid=netcdf.open(filepath,'NC_NOWRITE');
varid=netcdf.inqVarID(ncid,'hgt');                                               
end
 
if(strcmp(field,'w'))
filename=omega_list.name;
filepath=[folder,filename];
ncid=netcdf.open(filepath,'NC_NOWRITE');
varid=netcdf.inqVarID(ncid,'omega');
end
 
if(strcmp(field,'u'))
filename=uwnd_list.name;
filepath=[folder,filename];
ncid=netcdf.open(filepath,'NC_NOWRITE');
varid=netcdf.inqVarID(ncid,'uwnd');
end  
 
if(strcmp(field,'v'))
filename=vwnd_list.name;
filepath=[folder,filename];
ncid=netcdf.open(filepath,'NC_NOWRITE');
varid=netcdf.inqVarID(ncid,'vwnd');
end
 
varidp=netcdf.inqVarID(ncid,'level');
varidlat=netcdf.inqVarID(ncid,'lat');
varidlon=netcdf.inqVarID(ncid,'lon');
varidt=netcdf.inqVarID(ncid,'time');
 
value=netcdf.getVar(ncid,varid,'short');
scale = netcdf.getAtt(ncid,varid,'scale_factor');
offset= netcdf.getAtt(ncid,varid,'add_offset');
value=double(value)*scale+offset;
value=single(permute(value,[4 3 2 1]));
 
pressure=single(netcdf.getVar(ncid,varidp));
latitude=single(netcdf.getVar(ncid,varidlat));
longitude=single(netcdf.getVar(ncid,varidlon));
%to keep only the time asked for
time=netcdf.getVar(ncid,varidt);     
dayssince111=time/24;
datevalue=dayssince111+datenum(1,1,1)-2;	
mm=datevec(datevalue); mm=mm(:,2);
value=value(month==mm,:,:,:);                
netcdf.close(ncid);
 
value(value>10^14)=nan;
 
%pads the data for the month with NAN if not all month is available yet
theoreticaltime=datenum(year,month,1):(datenum(year,month+1,1)-1);
ts=length(theoreticaltime)*4;
dummy=zeros([ts,length(pressure),length(latitude),length(longitude)]);
dummy=dummy*0/0;
dummy(1:size(value,1),:,:,:)=value;
value=dummy;
 
end