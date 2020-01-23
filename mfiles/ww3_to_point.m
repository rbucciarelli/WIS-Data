%% ww3_to_point.m
%-- Load ww3 netcdf files remotely and interpolate data located at buoys
%-- THREDDS dataserver:
%-- https://data.nodc.noaa.gov/thredds/catalog/ncep/nww3/catalog.html
%---------------------------------------------------------------------

%% Initialize
clear all; close all;
%-- Use NCTOOLBOX to access grib data
run ./nctoolbox-master/setup_nctoolbox.m;
%-- Setup data source urls
src_url = 'ftp://polar.ncep.noaa.gov/pub/history/waves/';
src_model = 'multi_1';
src_file = 'multi_1/200502/partitions/multi_1.partition.wc_4m.200502.gz';
src_url = 'https://data.nodc.noaa.gov/thredds/catalog/ncep/nww3/2005/02/catalog.html'
src_grid = 'wc_4m';

src_dir = '../../WW3/data-grids/';
src_file = 'multi_1.partition.wc_10m.200502';
nco = ncdataset([src_dir src_file]);

