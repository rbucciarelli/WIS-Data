%% WIS_region.m
%-- Given a cdip_id, determine WIS region in (i.e. 'Pacific' or 'Atlantic')
%-- Requires: ../cdip-Atlantic.csv and ../cdip-Pacific.csv
function [ region,start_date ] = WIS_region(cdip_id)

    %cdip_id = '132';
    start_date = '';
    region = '';
    
    reg_list = {'Atlantic';'Pacific'};
    for i=1:length(reg_list)
        rfile = ['../cdip-' reg_list{i} '.csv'];
        M = importdata(rfile,',');
        stn_list = num2str(M(:,1),'%03d');
        if(~ isempty(strcmp(cdip_id,cellstr(stn_list))))
            si = find(strcmp(cdip_id,cellstr(stn_list)));
            ndbc_id = M(si,2);
            region = reg_list{i};
            start_date = M(si,3);
            break;
        else
            disp('Cannot find station in any region');
            return
        end
        
    end

end