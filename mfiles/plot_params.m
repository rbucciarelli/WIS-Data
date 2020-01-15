%% plot_params.m
clear all;
close all;

%-- Initialize
stn = '067';
stime = '2009';
etime = '2010';
var = 'hs';

figure;

%% Load datasets
dsets = ['C','A','W'];
data_dir = '../data/';
for i = 1:length(dsets)
    eval(['load ' data_dir dsets(i) stn '.mat']);
end

%% Plot var timeseries
for i = 1:length(dsets)
    h(i) = subplot(3,1,i);
    data = eval([dsets(i) stn]);
    plot(data.time,eval(['data.' var]));
    %-- Find data w/in lims
    axis tight;
    xlim([min([min(W067.time) min(C067.time) min(A067.time)]) ...
    max([max(W067.time) max(C067.time) max(A067.time)])]);
    autodatetick(gca, 'x');

end



%% Plot 2d at time = ti
figure
var2d = 'energy';
edate = datenum(2009,01,02);
%% Plot var timeseries
for i = 1:length(dsets)
    h(i) = subplot(3,1,i);
    data = eval([dsets(i) stn]);
    idx = find(data.time < edate);
    %plot(data.time,eval(['data.' var]));
    pcolor(data.b1(:,idx));
    shading flat;
    %caxis([min(min(C067.energy)) max(max(C067.energy))])
    %-- Find data w/in lims
    axis tight;
%     xlim([min([min(W067.time) min(C067.time) min(A067.time)]) ...
%     max([max(W067.time) max(C067.time) max(A067.time)])]);
%     autodatetick(gca, 'x');

end
% wi = W067.time(end);
% %-- Find closest cdip index to WW3 value ti
% [c cidx] = min(abs(C067.time-wi));
% [c cidx] = min(abs(wi - C067.time));


function [ output_args ] = nicedayticks( axeshandle, dt )
    %nicedayticks( axeshandle, dt ) 

    if ~exist('dt')
        dt = 1;
    end

    xlims = axeshandle.XLim;

    if isnumeric(dt)
        daylims = [round(xlims(1)):dt:round(xlims(2))];
    elseif strcmp(dt, 'month')
        [~, ~, firstday] = datevec(datestr(xlims(1)));
        startday = xlims(1)+(1-firstday);
        daylims = [startday:31:xlims(2)];
    end

    set(axeshandle, 'XTick', daylims);

end    
function [] = autodatetick(hax, dim,varargin)
    % AUTODATETICK 
    % 
    %   Syntax:
    %       [] = autodatetick(hax, dim)
    %
    % -------------------------------------------------------------------------
    %%% SETUP 
    % -------------------------------------------------------------------------
    if strcmp(dim, 'x') % only x axis
        tlims = hax.XLim;
    else
        error('We don''t want that axis')
    end

    %%% interval - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    Dt = abs(diff(tlims));
    dt = ceil(Dt/10);

    for i=1:length(varargin)
      vin = varargin{i};
      if isequal(vin,'dt')
        opt = varargin{i+1};
        dt = opt;
      end
    end

    %%% date style - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    if dt<30 
        dlabstyle = 'mm/dd/yy';
    else
        dt  = ceil(dt/30)*30;
        dlabstyle = 'mmm ''yy';
    end

    %%% information - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    disp(['Making date ticks with ' num2str(dt) '-day spacing'])

    % -------------------------------------------------------------------------
    %%% CHANGE AXIS DATE TICKS 
    % -------------------------------------------------------------------------
    nicedayticks(hax, dt);
    datetick(hax, 'x', dlabstyle,'keeplimits', 'keepticks')
end
%     legend(gca)
%     ylabel(varname)
%     box on; grid on;
%     set(gca, 'FontSize', 14)
%     title(['CDIPID=' cdipid ' (n=' num2str(length(t_cdip)) ') , NDBCID=' ndbcid  ' (n=' num2str(length(t_ww3)) ')']);
