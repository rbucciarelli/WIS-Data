%% run_stations.m
%-- Script to process WW3 and CDIP data for WIS comparisons

%% Initialize
clear all;

%-- Get a list of cdip stations to process
%% Initialize variables
cdip_id = '132';
start_time = '200602';      % YYYYMM
end_time = '201312';

%- Find ndbc_id using table: ../ndbc_id_table.csv
M = csvread('../cdip-westcoast.csv');
cdip_list = num2str(M(:,1),'%03g');

for i = 1:length(cdip_list)
    cdip_id = cdip_list(i,1:3);
    start_date = num2str(M(i,3));                % YYYYMM
    end_date = '201312';
    data = {};
    disp(['Processing CDIP ' cdip_id ' ---> ' start_date '-' end_date]);
    data = load_CDIP(cdip_id,start_date,end_date);
    disp('');
end

