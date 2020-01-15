% GOAL: loading data from CDIP (THREDDS, ideally) 
%       create fxn
% CURRENT VERSION : V3
clear all; close all;


mydir = 'C:\project\cdip\WIS\'; %  --- > change this! < ---
addpath(genpath([mydir 'buoyxmodel/v3/'])) %  --- > change this! < ---

%% downloading and plotting yourself with dload_cdipvar and dload_ww3var 
ndbc_id = '46219';
cdip_id = '067';
varname = 'Hs';
start_date = '201601';
end_date = '201612';

tlims = [datenum(2016,1,1) datenum(2017,1,1)];
% -------------------------------------------------------------------------
di = 10;
t1 = dload_cdipvar(067, 'Time',tlims,di);
var1 = dload_cdipvar(067, varname,tlims,di);
disp(['--------'])
disp(['beg time: ' datestr(t1(1))])
disp(['end time: ' datestr(t1(end))])
disp(['      n = ' num2str(length(t1))])


di = ceil(di/5.7);
di = 1;
t2 = dload_ww3var(46219, 'Time',tlims,di);
var2 = dload_ww3var(46219, varname,tlims,di);
ww3_time = t2;
ww3_hs = var2;
disp(['--------'])
disp(['beg time: ' datestr(t2(1))])
disp(['end time: ' datestr(t2(end))])
disp(['      n = ' num2str(length(t2))])

%% Load in WIS data from .mat file
data_dir = '../data/';
fname_WIS = ['WIS_',ndbc_id,'_','198001','-','201612'];
load([data_dir,fname_WIS]);
wis_hs = data.waveHs;
wis_time = data.matTime;
%- Get start and end index corresponding to start and end dates
yr = str2num(start_date(1:4));
mo = str2num(start_date(5:6));
yr_end = str2num(end_date(1:4));
idx = find( (wis_time >= datenum(yr,mo,1) & (wis_time < datenum(yr_end+1,1,1))));
wis_hs = wis_hs(idx);
wis_time = wis_time(idx);

%% Make sure these datasets have same timestep and N


figure(1); clf; hold on; 
plot(wis_time,wis_hs,'b','DisplayName','WIS');
plot(ww3_time,ww3_hs,'r','DisplayName','WIS');
datetick('x','YYYY-mm');
ylabel('Hs (m)');
xlabel('UTC');

% 
% figure(1); clf; hold on; 
% plot(t1, var1, 'r.', 'DisplayName', 'CDIP')
% plot(t2, var2, 'b.', 'DisplayName', 'WW3')
% legend;
% datetick('x', 'mmm ''yy')
% ylabel(varname)
% 
% %% using comp_buoyxmodel
% 
% %%% OPTION 1: Plot and compare all data:
% % -------------------------------------------------------------------------
% figure(2)
% comp_buoyxmodel(067, 46219, varname, 'plot', 'on', 'InterpMode', 'fast');
% 
% 
% % OPTION 2: Plot and compare a subset and at a resolution:
% % -------------------------------------------------------------------------
% figure(3)
% tlims = [datenum(2017,3,1) datenum(2019,4,1)];
% comp_buoyxmodel(067, 46219, varname, 'plot', 'on', 'InterpMode', 'fast', ...
%     'DataRes', 2, 'TimeLimits', tlims);
% 
% % OPTION 3: ONLY get data:
% % -------------------------------------------------------------------------
% [A B] = comp_buoyxmodel(067, 46219, varname, 'plot', 'off', 'InterpMode', 'fast');
% 
% %% the five stations
% 
% % varname = 'Dp';o = 10;
% varname = 'Hs';o = 0;
% % varname = 'Tp';o = 20;
% 
% % -------------------------------------------------------------------------
% 
% % CALIFORNIA - SAN NICOLAS
% % PLOT
% figure(11+o);
% comp_buoyxmodel(067, 46219, varname, 'plot', 'on', 'InterpMode', 'fast');
% 
% % WASHINGTON - GRAYS HARBOR
% figure(12+o);
% comp_buoyxmodel(036, 46211, varname, 'plot', 'on', 'InterpMode', 'fast');
% 
% % HAWAII - WAIMEA BAY
% figure(13+o);
% comp_buoyxmodel(106, 51201, varname, 'plot', 'on', 'InterpMode', 'fast');
% 
% % GULF COAST - ST PETERSBURG OFFSHORE
% figure(14+o);
% comp_buoyxmodel(144, 42099, varname, 'plot', 'on', 'InterpMode', 'fast');
% 
% % EAST COAST - DUCK FRF 26M
% figure(15+o);
% comp_buoyxmodel(430, 44100, varname, 'plot', 'on', 'InterpMode', 'fast');
