%% load_WW3.m
%-------------------------------------------------------------------------
%- WW3 spectral data is pre- downloaded from NCEI as tarred ascii files.
%- This script reads and parses spectral data for individual buoys with
%- NDBC id
%-------------------------------------------------------------------------

clear all;

%% Initialize variables
cdip_id = '106';
start_time = '200502';      % YYYYMM
end_time = '201306';

%- Find ndbc_id using table: ../ndbc_id_table.csv
M = csvread('../ndbc_id_table.csv');
index = find(M(:,1) == str2num(cdip_id));
ndbc_id = num2str(M(index,2));       %'46219';

%- WW3 files are by month, create an array of dates
date_list = get_dates(start_time,end_time); 

ww3_dir = '\\d.cdip.ucsd.edu\data09_01\WW3\data\multi_1\';
data_dir = [ww3_dir,ndbc_id,'\'];

file_info = dir([data_dir,'*SPEC*']);
data_WW3 = {};
si = 1;     %- start_index
    
%- Iterate over dates and load spectral files
for i = 1:length(date_list)
    fdate = date_list(i);
    fname = ['multi_1.',ndbc_id,'.','HIND_SPEC.',num2str(fdate)];
    disp(['--> ' fname]);
    if (exist([data_dir,fname],'file'))
        data = ww3_to_wnc([data_dir,fname]);
    end
    %- The WW3 files have one extra record extending into next day,
    %- need to pop this off of relevant fields
    ei = si + length(data.time) - 2;
    data_WW3.time(si:ei) = data.time(1:end-1);
    data_WW3.lon = data.lon;
    data_WW3.lat = data.lat;
    data_WW3.depth = data.depth;
    data_WW3.f = data.f;
    data_WW3.dir = data.dir;
    data_WW3.bw = data.bw;
    data_WW3.energy(:,si:ei) = data.energy(:,1:end-1);
    data_WW3.a0(:,si:ei) = data.a0(:,1:end-1);
    data_WW3.a1(:,si:ei) = data.a1(:,1:end-1);
    data_WW3.b1(:,si:ei) = data.b1(:,1:end-1);
    data_WW3.a2(:,si:ei) = data.a2(:,1:end-1);
    data_WW3.b2(:,si:ei) = data.b2(:,1:end-1);
    data_WW3.hs(si:ei) = data.hs(1:end-1);
    si = ei + 1;
end

            
eval(['W' cdip_id '=data_WW3;']);                            %- W067

% %% Save data to .mat file
out_dir = '../data/';

savefile = [out_dir 'W' cdip_id '.mat'];
save(savefile,['W' cdip_id]);








