function [var_cdip] = dload_cdipvar(cdipid,varname,tlims, tres)
% -------------------------------------------------------------------------
% DLOAD_CDIPVAR  Downloads model data from CDIP THREDDS (historic and rt)
%                Here we usually require a smaller subset of the data
%                because it's way too big to download all of it.
% -------------------------------------------------------------------------
%   Syntax:
%      [var_cdip] = dload_cdipvar(ndbcid,varname,tlims,tres)
%
%   Inputs:
%      CDIPID 
%      VARNAME      options: 'Time' or 'Hs' or 'Dp' or 'Tp' or 'Ta'
%      TLIMS        default = [Nt-100 to Nt];
%      TRES         default = 1;
% 
%   Output:
%      [var_cdip]
% 
%   Uses:
%      get_tinfo.m (subfunction)
% 
%   Sample: 
%       dload_cdipvar(067,'Hs',[735965 736330],1) 
%       NOTE: this is from 2015 to 2016 every index for the San Nic buoy
%
% Updated as of 07-18-2019 by Alli Ho
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% function set up 
if isnumeric(cdipid)
    cdipid = num2str(cdipid);
    if length(cdipid)<3
        cdipid = ['0' cdipid];
    end
end
cdipname = cdipid;

if ~exist('tlims')
    tlims = [];
end

if ~exist('tres')
    tres = 1;
end

%% LOAD HISTORIC CDIP DATA FROM CDIP THREDDS
% -------------------------------------------------------------------------
%%% observations (buoys)
mode = 'hist';
baseurl = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% FIRST: get tinfo and limit time series if need be
[~, ~, starti, endi, ~] = get_tinfo(baseurl, cdipname, tlims,mode);
% check if starti, endi, tres, are even/logical
tt = [starti:tres:endi];
endi = tt(end);
tstring = [ '[' num2str(starti) ':' num2str(tres) ':' num2str(endi) ']'];

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% SECOND: create URL name and load data
if starti~=endi
    nameurl = [cdipname 'p1/' cdipname 'p1_historic.nc.ascii?'];
    paramurl = ['wave' varname tstring]; % to do figure out subset!
    url = [baseurl nameurl paramurl];
    waveVar = dload(url);
    var_cdip_1 = waveVar;
else
    var_cdip_1 = [];
    disp('---> Not using CDIP historic');
end

%% LOAD REALTIME CDIP DATA FROM CDIP THREDDS
% -------------------------------------------------------------------------
%%% observations (buoys)
mode = 'rt';
baseurl = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% FIRST: get tinfo and limit time series if need be
[~, ~, starti, endi, ~] = get_tinfo(baseurl, cdipname, tlims,mode);
% check if starti, endi, tres, are even/logical
tt = [starti:tres:endi];
endi = tt(end);
tstring = [ '[' num2str(starti) ':' num2str(tres) ':' num2str(endi) ']'];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% SECOND: create URL name and load data
if starti~=endi
    nameurl = [cdipname 'p1_rt.nc.ascii?'];
    paramurl = ['wave' varname tstring]; % to do figure out subset!
    url = [baseurl nameurl paramurl];
    waveVar = dload(url);
    var_cdip_2 = waveVar;
else
    var_cdip_2 = [];
    disp('---> Not using CDIP realtime');
end

%% COMBINE REALTIME AND HISTORIC
% -------------------------------------------------------------------------
var_cdip = [var_cdip_1 var_cdip_2];

%% IF TIME TURN TO DATENUM FORMAT
% -------------------------------------------------------------------------
if strcmp(varname, 'Time')
    toff = datenum(1970,1,1,0,0,0);
    var_cdip = var_cdip./(24*60*60) + toff;
end

end

function [tstart tend starti endi Nt] = get_tinfo(baseurl, cdipname, tlims,mode)
    if strcmp(mode, 'hist')
        infourl = [baseurl cdipname 'p1/' cdipname 'p1_historic.nc.html'];
%         ilines = [122, 388,389];
    elseif strcmp(mode, 'rt')
        infourl = [baseurl cdipname 'p1_rt.nc.html'];
%         ilines = [120, 398,399];
    else
        error('Mode for CDIP data not compatible');
    end
    info = urlread(infourl); info = strsplit(info,'\n');

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% FIND VAR = Nt   
    idx = find(contains(info, 'Int32 waveTime[waveTime = ')); idx = idx(1);
    infot = strsplit(info{idx}, '= '); infot = infot{2};
    Nt = str2num(infot(1:end-2))-1;
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% FIND VAR = tstart, tend
    idx = find(contains(info, 'time_coverage_start'));
    tstr = strsplit(info{idx},' '); 
    tstr = tstr{2};
    yr = str2num(tstr(1:4)); mo = str2num(tstr(6:7)); day = str2num(tstr(9:10));
    tstart = datenum(yr,mo,day);
    
    idx = find(contains(info, 'time_coverage_end'));
    tstr = strsplit(info{idx},' '); 
    tstr = tstr{2};
    yr = str2num(tstr(1:4)); mo = str2num(tstr(6:7)); day = str2num(tstr(9:10));
    tend = datenum(yr,mo,day);
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%     disp(['     CDIP runs from ' datestr(tstart) ' to '  datestr(tend)])
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% FIND VAR = starti, endi
    if ~isempty(tlims)
        endi = Nt;
        starti = 1;
        % downloading coarse resolution time series data to get
        % approximate index location for tlims
        sdi = 100; % change this if you want to make the time cuttoff more accurate
        endi = floor(endi/sdi).*sdi;
        tstring = [ '[' num2str(starti) ':' num2str(sdi) ':' num2str(endi) ']']; 
        idx = [starti:sdi:endi];
        
        if strcmp(mode, 'hist')
            nameurl = [cdipname 'p1/' cdipname 'p1_historic.nc.ascii?'];
        elseif strcmp(mode, 'rt')
            nameurl = [cdipname 'p1_rt.nc.ascii?'];
        end        
        paramurl = ['waveTime' tstring]; % to do figure out subset!
        url = [baseurl nameurl paramurl];
        waveVar = dload(url);
        
        toff = datenum(1970,1,1,0,0,0);
        simpletime = waveVar./(24*60*60) + toff; % the course resolution data, corresponding indexes in VAR=idx

        [~, idxstart] = min(abs(simpletime-tlims(1)));
        [~, idxend] = min(abs(simpletime-tlims(2)));

        starti = idx(idxstart); endi = idx(idxend);
    else % or small subset (data way to big to pull entire thing)
        starti = Nt-1000;
        endi = Nt;
    end
end

function [var] = dload(url) % download thredds data from url
    data = urlread(url); data = strsplit(data,'\n');
    dstart = find(strcmp(data,'---------------------------------------------'))+1;

    var = strsplit(data{dstart+1}, ', ');
    var = cellfun(@str2double,var);
end