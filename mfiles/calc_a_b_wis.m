%-- Function: calc_a_b.m
%-- Inputs:  freq (vector NF x 1)
%--             dir (ND x 1)
%--             energy2d:  (time, freq, dir)
%-------------------------------------------------------------------------
function [ data ] = calc_a_b_wis(freq,dir,energy2d)

data = {};

rdir = deg2rad(dir);
dtheta = abs(rdir(2)-rdir(1));


%-- input energy is 3 Dims (time, freq, dir)
ds = energy2d;   % (ND x NF)
[ND, NF, NT ] = size(ds);   %-- # Time, # Freqs, # Dir
bands = NF;
dbins = ND;


deg_per_bin = 360/dbins;

%% Assign WW3 directions to correct bins (direction coming from, starting with 5 deg bin)
angle = round(rad2deg(rdir)) + 180;
idx = find(angle >= 360);
angle(idx) = angle(idx) - 360;
dir_idx = round(angle./deg_per_bin);
bin_dirs = zeros(dbins,1);
for i = 1:length(dir_idx)
    idx = dir_idx(i);
    bin_dirs(idx) = angle(i);
end

%% Calculate freq bandwidth
bw = zeros(bands,1);
for j = 2:bands-1
    bw(j) = (freq(j+1) - freq(j-1))/2;
end
bw(1) = freq(1) * (bw(2)/freq(2));
bw(bands) = freq(bands) * (bw(2)/freq(2));

data.bw = bw;

icnt=0;
pstnid = '';
a0 = zeros(NF,NT);
a1 = zeros(NF,NT);
b1 = zeros(NF,NT);
a2 = zeros(NF,NT);
b2 = zeros(NF,NT);
energy = zeros(NF,NT);
Hs = [];
Tp = [];
Dp = [];
Ta = [];

        
zero_energy = [];
energy = zeros(NF,1);

sp1d = zeros(NF,NT);
ds_set = zeros(dbins,bands,NT);

for icnt = 1:NT
    ds = squeeze(energy2d(:,:,icnt));       %-- (ND x NF)
              
    for i = 1:bands
        for j = 1:dbins
            ds_set(dir_idx(j),i,icnt) = ds(j,i);
        end
    end   
    
    %sp1d(icnt) = sum(ds')*dtheta;
    %b.sp1d{index+1} = sum(sp2d')*dtheta;
    
    %--   When calculating moments, rotate directions by pi so that the resulting 
    %-   coefficients are in true compass "arriving from" coordinates.
    zero_energy(icnt) = true;
    for i = 1:bands
        
        for j = 1:dbins
            

            %--   Little floating point exception check here to catch 
            %--   very small double-precision energies from ww3 model
            if(ds_set(j,i,icnt) < 1.0e-15) 
                ds_set(j,i,icnt) = 0.0;
            end
            
            ds_set(j,i,icnt) = ds_set(j,i,icnt)*(2*pi/dbins);
            a0(i,icnt) = a0(i,icnt)+ds_set(j,i,icnt);
            a1(i,icnt) = a1(i,icnt)+ds_set(j,i,icnt)*cos(deg2rad(bin_dirs(j)));
            b1(i,icnt) = b1(i,icnt)+ds_set(j,i,icnt)*sin(deg2rad(bin_dirs(j)));
            a2(i,icnt) = a2(i,icnt)+ds_set(j,i,icnt)*cos(deg2rad(2.*bin_dirs(j)));
            b2(i,icnt) = b2(i,icnt)+ds_set(j,i,icnt)*sin(deg2rad(2.*bin_dirs(j)));
            ds_set(j,i,icnt) = ds_set(j,i,icnt) / deg_per_bin;            
                 
        end     %-- End of 1st loop through directions (dbins)
        
        %--   Normalize fourier coeffiecients
        if(a0(i,icnt) > 0) 
            zero_energy(icnt) = false;
            a1(i,icnt) = a1(i,icnt)/a0(i,icnt);
            b1(i,icnt) = b1(i,icnt)/a0(i,icnt);
            a2(i,icnt) = a2(i,icnt)/a0(i,icnt);
            b2(i,icnt) = b2(i,icnt)/a0(i,icnt);
            dir(i,icnt) = rad2deg(atan2(b1(i,icnt),a1(i,icnt)));
            if (dir(i,icnt) < 0) 
                dir(i,icnt) = dir(i,icnt) + 360;
            end
        end     
        data.a0(i,icnt) = a0(i,icnt);
        data.a1(i,icnt) = a1(i,icnt);
        data.b1(i,icnt) = b1(i,icnt);
        data.a2(i,icnt) = a2(i,icnt);
        data.b2(i,icnt) = b2(i,icnt);
        data.energy(i,icnt) = sum(ds(:,i))*dtheta;
        
        
                

    end  %- End of freq loop
    
    %--   Calculate Hs, Tp, Dp, Ta

    if (~ zero_energy(icnt)) 
        total_energy = 0;
        M0 = 0;
        M1 = 0;
        for i = 1:bands
            energy(i) = a0(i,icnt) * bw(i);
            total_energy = total_energy + energy(i);
            M0 = M0 + energy(i);
            M1 = M1 + energy(i) * freq(i);
        end
        Hs(icnt) = 4 * sqrt(total_energy);
        [tmp,peak_band]=max(a0(:,icnt));
        Tp(icnt) = 1.0 / freq(peak_band);
        Dp(icnt) = dir(peak_band,icnt);
        Ta(icnt) = M0 / M1;
       
    end
    
end   
data.hs = Hs;
data.tp = Tp;
data.dp = Dp;
data.ta = Ta; 
end