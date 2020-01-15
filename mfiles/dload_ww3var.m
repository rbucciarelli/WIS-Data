function [var_ww3] = dload_ww3var(ndbcid,varname,tlims,tres)
% -------------------------------------------------------------------------
% DLOAD_WW3VAR   Downloads model data from CDIP THREDDS 
%                Here we choose to download only realtime model, and all of
%                it, which goes from present back to June 2015.
% -------------------------------------------------------------------------
%   Syntax:
%      [var_ww3] = dload_ww3var(ndbcid,varname,tlims,tres)
%
%   Inputs:
%      NDBCID 
%      VARNAME      options: 'Time' or 'Hs' or 'Dp' or 'Tp' or 'Ta'
%      TLIMS        default = [Nt-1000 to Nt];
%      TRES         default = 1;
% 
%   Output:
%      [var_ww3]
% 
%   Uses:
%      get_tinfo.m (subfunction)
% 
%   Sample: 
%       dload_ww3var(46219,'Hs',[])
%       NOTE: this is for all indexes for the San Nic buoy
% 
% Updated as of 07-18-2019 by Alli Ho
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% function set up 
if isnumeric(ndbcid)
    ndbcid = num2str(ndbcid);
end

ww3name = ndbcid;
if ~exist('tlims')
    tlims = [];
end

if ~exist('tres')
    tres = 1;
end

dim = 1;
if strcmp(varname, 'EnergyDensity') | strcmp(varname, 'DirectionalSpectrum')
    dim = 2;
end

%% LOAD WW3 from CDIP THREDDS
% -------------------------------------------------------------------------

%%% ww3 spectral output
baseurl = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/test/WW3/';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% FIRST: get length of time variable
[tstart tend starti endi Nt] = get_tinfo(baseurl, ww3name, tlims);
% check if starti, endi, tres, are even/logical
tt = [starti:tres:endi];
endi = tt(end);
tstring = [ '[' num2str(starti) ':' num2str(tres) ':' num2str(endi) ']'];

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% SECOND: create URL name and load data
if starti~=endi
    nameurl = [ww3name '_WW3_realtime.nc.ascii?'];
    paramurl = ['wave' varname tstring];
    url = [baseurl nameurl paramurl]; %full url name
    waveVar = dload(url);
    var_ww3 = waveVar;
else
    var_ww3 = [];
    disp('---> Not using WW3 realtime');
end

%% IF TIME TURN TO DATENUM FORMAT
% -------------------------------------------------------------------------

if strcmp(varname, 'Time')
    toff = datenum(1970,1,1,0,0,0);
    var_ww3 = var_ww3./(24*60*60) + toff;
end

end


function [tstart tend starti endi Nt] = get_tinfo(baseurl, ww3name, tlims)
    infourl = [baseurl ww3name '_WW3_realtime.nc.html'];
    info = urlread(infourl); info = strsplit(info,'\n');

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% FIND VAR = Nt   
    idx = find(contains(info, 'Int32 waveTime[waveTime = ')); idx = idx(1);
    infot = strsplit(info{idx}, '= '); infot = infot{2};
    Nt = str2num(infot(1:end-2))-1;
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% FIND VAR = tstart, tend
    idx = find(contains(info, 'date_created'));
    tstr = strsplit(info{idx},' '); 
    tstr = tstr{2};
    yr = str2num(tstr(1:4)); mo = str2num(tstr(6:7)); day = str2num(tstr(9:10));
    tstart = datenum(yr,mo,day);
    
    idx = find(contains(info, 'date_modified'));
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
        sdi = 50; % change this if you want to make the time cuttoff more accurate
        endi = floor(endi/sdi).*sdi;
        tstring = [ '[' num2str(starti) ':' num2str(sdi) ':' num2str(endi) ']']; 
        idx = [starti:sdi:endi];
        
        nameurl = [ww3name '_WW3_realtime.nc.ascii?'];
        paramurl = ['waveTime' tstring]; % to do figure out subset!
        url = [baseurl nameurl paramurl];
        waveVar = dload(url);
        
        toff = datenum(1970,1,1,0,0,0);
        simpletime = waveVar./(24*60*60) + toff; % the course resolution data, corresponding indexes in VAR=idx

        [~, idxstart] = min(abs(simpletime-tlims(1)));
        [~, idxend] = min(abs(simpletime-tlims(2)));

        starti = idx(idxstart); endi = idx(idxend);
    else % get the entire thing
        starti = 1;
        endi = Nt;
    end
end


function [var] = dload(url) % download thredds data from url
    data = urlread(url); data = strsplit(data,'\n');
    dstart = find(strcmp(data,'---------------------------------------------'))+1;

    var = strsplit(data{dstart+1}, ', ');
    var = cellfun(@str2double,var);
end