%- run_WIS.m
cdip_id = '132';
start_date = '200901';
end_date = '201012';
[region,buoy_start] = WIS_region(cdip_id);
date_list = get_dates(start_date,end_date);
for i = 1:length(date_list)
    the_date = num2str(date_list(i));
    WIS_to_mat(cdip_id,region,the_date);
end