%% load_WIS.m
%-------------------------------------------------------------------------
%- Load an individual WIS netcdf file into matlab struct and save to .mat
%- The WIS model data is stored in NetCDF format on following server:
%- https://chlthredds.erdc.dren.mil/thredds/catalog/wis/catalog.html
%- Files are stored as individual months in specific regions.
%-------------------------------------------------------------------------

function [ data ] = load_WIS(cdip_id,region,yyyymm)

data = {};

%% Initialize variables
%region = 'Atlantic';         %- Need to figure out which region buoy is in.
% cdip_id = '132';
% yyyymm = '201610';

%- Find ndbc_id using table: ../ndbc_id_table.csv
M = csvread('../ndbc_id_table.csv');
index = find(M(:,1) == str2num(cdip_id));
ndbc_id = num2str(M(index,2));       %'46219';

si = 1;     %- start_index
idx = 0;

url = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/wis/';
url = [url region '/ST' ndbc_id '/'];                  %Pacific/ST46219/

yyyy = yyyymm(1:4);                 %- '2009'
src_dir = [url,yyyy,'/'];
fname = ['WIS-ocean_waves_ST',ndbc_id,'_',yyyymm,'.nc'];
disp(['Loading ' fname]);

try
    finfo = ncinfo([src_dir fname]);
    ncid = netcdf.open([src_dir fname],'NC_NOWRITE');
    %[ndims,nvars,natts,unlimdimID] = netcdf.inq(ncid);
     
    %-- Read depth from global attributes
    depth = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'),'depth');
    
    %-- Create a list of variables from netcdf
    fvars = {};
    for vid = 1:length(finfo.Variables)
        fvars{vid} = finfo.Variables(vid).Name;
    end
    
    %- List of variables to extract from NetCDF file
    if (find(contains(fvars,'longitude')))
        var_list = {'time', 'waveFrequency', 'directionalWaveEnergyDensity', ...
        'waveDirectionBins', 'waveHs', 'waveTp', 'waveTm', ...
        'latitude', 'longitude'};
    else
        var_list = {'time', 'waveFrequency', 'directionalWaveEnergyDensity', ...
        'waveDirectionBins', 'waveHs', 'waveTp', 'waveTm', ...
        'lat', 'lon'};
    end

    %% Iterate over vars of interest and get data
    for vid = 1:length(var_list)
        the_var = var_list{vid};
        %-- Find the index of this variableda
        vindex = find(strcmp(fvars,the_var)) - 1;
        disp(['--> Loading ' the_var]);
        eval([the_var '=netcdf.getVar(ncid,' num2str(vindex) ');']);
    end
    netcdf.close(ncid);
    %% Aggregate monthly data based on start/end times
    mat_time = time_correct(double(time));
    %-- Figure out end index
    ei = si + length(time) - 2;
    if (find(contains(fvars,'longitude')))
        data.lat = latitude;
        data.lon = longitude;
    else
        data.lat = lat;
        data.lon = lon;
    end

    data.depth = depth;
    data.f = waveFrequency;
    data.dir = waveDirectionBins;
   
    %-- Calculate a1,b1,a2,b2
    [ND,NF,NT] = size(directionalWaveEnergyDensity);
    %-- Function requires 2d spectra to be [time x freq x dir]
    energy2d = directionalWaveEnergyDensity;    
    
    ab_data = calc_a_b_wis(waveFrequency,waveDirectionBins,energy2d);
    data.bw = ab_data.bw;
    data.time(si:ei) = mat_time(1:end-1);
    data.energy(:,si:ei) = ab_data.energy(:,1:end-1);
    data.a0(:,si:ei) = ab_data.a0(:,1:end-1);
    data.a1(:,si:ei) = ab_data.a1(:,1:end-1);
    data.b1(:,si:ei) = ab_data.b1(:,1:end-1);
    data.a2(:,si:ei) = ab_data.a2(:,1:end-1);
    data.b2(:,si:ei) = ab_data.b2(:,1:end-1);
    data.hs(si:ei) = waveHs(1:end-1)';
 
    %si = ei + 1;
    %-- Save data to .mat file
    eval(['A' cdip_id '=data;']);   
    out_dir = '../data/';
    savefile = ['A',cdip_id,'_',yyyymm,'.mat'];
    save([out_dir savefile],['A' cdip_id])
catch
    disp(['Could not open file: ' fname]);
    
end



end

%% Function to correct time epoch from Allie H cdipxww3 
function [ mat_time ] = time_correct(ts)
    toff = datenum(1970,1,1,0,0,0);
    mat_time = ts./(24*60*60) + toff;
end