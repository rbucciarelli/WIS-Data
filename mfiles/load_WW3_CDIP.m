%% load_WW3_CDIP.m
%-------------------------------------------------------------------------
%- WW3 data on CDIP THREDDS server
%- http://thredds.cdip.ucsd.edu/thredds/catalog/test/WW3/catalog.html
%- http://thredds.cdip.ucsd.edu/thredds/dodsC/test/WW3/41114_WW3_realtime.nc
%- Files are stored as complete time-series for historic stations
%- Usage:  data = load_WW3_CDIP(100,'201610','201612',1)    % all variables
%-         data = load_WW3_CDIP(100,'201610','201612','bulk',0)   %- bulk parmas
%-------------------------------------------------------------------------

function [data] = load_WW3_CDIP(cdip_id,start_time,end_time,params)

    data = {};

    %- Find ndbc_id using table: ../ndbc_id_table.csv
    M = csvread('../ndbc_id_table.csv');
    index = find(M(:,1) == str2num(cdip_id));
    ndbc_id = num2str(M(index,2));       %'46219';

    %- WW3 files are by month, create an array of dates
    date_list = get_dates(start_time,end_time); 

    if (params == 1)
        var_list = {'waveTime', 'waveFrequency', 'waveEnergyDensity', ...
            'waveBandwidth', 'waveA1Value', 'waveB1Value', 'waveA2Value', ...
            'waveB2Value', 'waveHs', 'waveTp', 'waveDp', 'waveTa', ...
            'metaLatitude', 'metaLongitude', 'metaWaterDepth', ...
            'metaSiteLabel'};
    else
        var_list = {'waveTime', 'waveFrequency', ...
            'waveBandwidth', 'waveHs', 'waveTp', 'waveDp', 'waveTa', ...
            'metaLatitude', 'metaLongitude',  ...
            'metaWaterDepth', 'metaSiteLabel'};
    end


    %% Load CDIP data from THREDDS Server
    %url = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/';
    url = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/test/WW3/';   
    %url = [url ndbc_id 'p1/'];
    fname = [ndbc_id '_WW3_realtime.nc']; 
    ncid = netcdf.open([url fname],'NC_NOWRITE');
    [ndims,nvars,natts,unlimdimID] = netcdf.inq(ncid);
    finfo = ncinfo([url fname]);
    %-- Create a list of variables from netcdf
    fvars = {};
    for vid = 1:length(finfo.Variables)
        fvars{vid} = finfo.Variables(vid).Name;
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

    %% Subset data based on start/end times
    mat_time = time_correct(double(waveTime));
    stime = datenum(str2num(start_time(1:4)),str2num(start_time(5:6)),1);
    eday = eomday(str2num(end_time(1:4)),str2num(end_time(5:6)));
    etime = datenum(str2num(end_time(1:4)),str2num(end_time(5:6)),eday);
    idx = find((mat_time >= stime) & (mat_time < etime));

    %% Output data to data structure
    data.time = mat_time(idx);
    data.name = string(cellstr(metaSiteLabel'));
    data.lat = metaLatitude;
    data.lon = metaLongitude;
    data.depth = metaWaterDepth;
    data.f = waveFrequency;
    data.bw = waveBandwidth;
    if (params == 1)
        data.energy = waveEnergyDensity(:,idx);
        data.a1 = waveA1Value(:,idx);
        data.a2 = waveA2Value(:,idx);
        data.b1 = waveB1Value(:,idx);
        data.b2 = waveB2Value(:,idx);
    end
    data.hs = waveHs(idx);
    

    %% Save data to .mat file
    eval(['W' cdip_id '=data;']);   
    out_dir = '../data/';
    savefile = ['W',cdip_id,'_',start_time,'.mat'];
    save([out_dir savefile],['W' cdip_id])

end

%% Function to correct time epoch from Allie H cdipxww3 
function [ mat_time ] = time_correct(ts)
    toff = datenum(1970,1,1,0,0,0);
    mat_time = ts./(24*60*60) + toff;
end








