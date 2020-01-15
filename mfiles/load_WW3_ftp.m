%% load_WW3.m
clear all;

%% Initialize variables
%url = 'https://data.nodc.noaa.gov/thredds/catalog/ncep/nww3/catalog.html';
region = 'Pacific';
buoy = '46219';
start_time = '200601';      % YYYYMM
end_time = '200612';

[data_WW3, N] = download_data(buoy,start_time,end_time);

%params = merge_data(data_WW3,N);

% time = time_correct(params.time);
% plot(time,params.hs)
% axis tight;
% datetick('x','YYYY-mm');
% ylabel('Hs (m)');
% xlabel ('UTC');
% title([buoy,' WIS Output ', start_time(1:4),'-',end_time(1:4),]);

%% Function to correct time epoch from Allie H cdipxww3 
function [ utc_time ] = time_correct(ts)
    toff = datenum(1970,1,1,0,0,0);
    utc_time = ts./(24*60*60) + toff;
end

function [ data, N ] = download_data(buoy, start_time, end_time)
    
%- FTP Server
src_url = 'ftp://polar.ncep.noaa.gov/pub/history/waves/multi_1/';

%- Get Year and Month
syear = start_time(1:4);
smonth = start_time(5:6);
eyear = end_time(1:4);
emonth = end_time(5:6);
YEARS = str2num(syear):str2num(eyear);
data = {};
N = 0;

    for i = 1:length(YEARS)
        yyyy = num2str(YEARS(i),'%04i');
        if (yyyy == syear)
            months = str2num(smonth):12;
        elseif (yyyy == eyear)
            months = 1:str2num(emonth);
        else
            months = 1:12;
        end
        
        %src_dir = [src_url,yyyy,'/'];
        data_year = struct();            %- Array with 12 months of data
        for j = 1:length(months)
            data_month = struct();
            mm = num2str(months(j),'%02i');
            ftime = [yyyy,mm];
            %- Data directory to download to
            data_dir = ['../WW3/data/', ftime]

            %- Spectral Filename
            spec_fname = ['multi_1_base.buoys_spec.', ftime, '.tar.gz'];
            url = [src_url, ftime, '/points/', spec_fname];
        
            %- WMO filename
            %- <year> <month> <day> <hour> <uabs> <udir> <Hs> <Tp>
            %- year  month  day  hour   wind-speed   wind-dir  Hs   Tpeak
            wmo_fname = ['multi_1_base.buoys_wmo.', ftime, '.tar.gz'];
            url = [src_url, ftime, '/points/', wmo_fname];

            %- Download and untar data (currently not on THREDDS server)
            %filenames = untar(url, data_dir );
            
            %- Access file of interest (e.g. multi_1.42036.HIND_WMO.200502)
            fname = ['multi_1.', buoy, '.HIND_WMO.', ftime];
            [ dn, wind_speed, wind_dir, Hs, Tp ] = load_wmo([data_dir,'/',fname]);
            data_month.dn = dn;
            data_month.Hs = Hs;
            data_month.month = mm;
            N = N + length(dn);
            data_year(j).data = data_month;       
            
            j = j + 1;
        end
        data(i).data = data_year;
        i = i + 1;       

    end
    
end

function [ dn, wind_speed, wind_dir, Hs, Tp ] = load_wmo(fname)
    fid = fopen(fname, 'r');
    formatSpec = '%d%d%d%d%f%d%f%f';
    data = textscan(fid,formatSpec,inf);
    mm = zeros(length(data{3}),1);
    ss = zeros(length(data{3}),1);
    dn = datenum(double(data{1}),double(data{2}),double(data{3}),double(data{4}),mm,ss);    %- datenum
    ts = [num2str(data{1}),num2str(data{2},'%02i'),num2str(data{3},'%02i'),num2str(data{4},'%02i')]; %- timestamp YYYYMMDDHH
    wind_speed = data{5};
    wind_dir = data{6};
    Hs = data{7};
    Tp = data{8};
    fclose(fid);
end

%% Function to merge years and months within struct to single struct with 
%- Fields:   hs and time
function [ data_params ] = merge_data(data_WW3,N)
data_params = struct;

% data_params.hs = zeros(N,1);
% data_params.time = zeros(N,1);
for i = 1:length(data_WW3)      %- iterate over years
    data_year = data_WW3(i);
    for j = 1:12            %- iterate over months
        ds = data_year{j};
        if (i==1) & (j == 1)
            count = ds.size('time');
            data_params.time(1:count) = ds.data('time');
            data_params.hs(1:count) = ds.data('waveHs');
        else        
            data_params.time(count+1:count+ds.size('time')) = ds.data('time');
            data_params.hs(count+1:count+ds.size('waveHs')) = ds.data('waveHs');
            count = count + ds.size('time');
        end
        %data_params.time(1:ds.size('time')) = 
        %data_params.time =  [data_params.time; ds.data('time')];
    end   
end

end
