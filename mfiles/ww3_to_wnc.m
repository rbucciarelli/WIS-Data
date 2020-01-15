%% www3_to_wnc
%- Read in WW3 spectral data and calculate parameters
%- Based on fortran code: ww3_to_wnc.f
%---------------------------------------------------------------------

function [ data ] = ww3_to_wnc(filename)

% data_dir = '../WW3/data/46219/';
% fname = 'multi_1.46219.HIND_SPEC.201307';
% filename = [data_dir,fname];

if (exist(filename,'file'))
    fid=fopen(filename);
else
    fprintf(1,'ERROR : File Not found !!');
    return;
end
% Scan basic properties of spectra from file
dum=fscanf(fid,'%c',[23]);
bands =fscanf(fid,'%g',[1]);
dbins =fscanf(fid,'%g',[1]);
deg_per_bin = 360/dbins;

NF = bands;
ND = dbins;

%% The array rdir contains the true 
%  compass head from which the waves are heading in radians
dum=fscanf(fid,'%s',[1]);
dum=fscanf(fid,'%c',[33]);
freq=fscanf(fid,'%g',NF);
rdir=fscanf(fid,'%g',ND);

%% Assign WW3 directions to correct bins (direction coming from, starting with 5 deg bin)
angle = round(rad2deg(rdir)) + 180;
idx = find(angle >= 360);
angle(idx) = angle(idx) - 360;
dir_idx = round(angle./deg_per_bin) + 1;    %-- Randy added
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

icnt=0;
pstnid = NaN;

%% Loop over spectral files within the WW3 output

dir = rdir;
f = freq;
        
fprop=f(2)/f(1);
fdum=[f(1)/fprop;f];

df=diff(fdum);
dtheta=abs(dir(2)-dir(1));

b.f = f;
b.dir = dir;
b.bw = bw;
dir(end+1) = dir(1); % Closing for circular integral

% Loop over times
index = 0;
b.time = [];
b.espt = [];
b.U10 = [];
b.Udir=[];
b.cU=[];

b.cdir=[];
b.sp1d = [];
b.hs = [];
b.dp = [];
b.fp = [];

%- Randy adding a1,b1,a2,b2
b.a1 = [];
b.a2 = [];
b.b1 = [];
b.b2 = [];

%b.sp2d = zeros(NF,ND);

while ~feof(fid)
    dumdate=fscanf(fid,'%g',2);
    if numel(dumdate)
        year = floor(dumdate(1)/10000);
        month = floor(dumdate(1)/100) - year*100;
        day = dumdate(1) - year*10000 - month*100;
        hr = floor(dumdate(2)/10000);
        mn = floor(dumdate(2)/100) - hr*100;
        sc = dumdate(2) - hr*10000 - mn*100;
        
        b.time(index+1) = datenum(year,month,day,hr,mn,sc);

        
        stn=fscanf(fid,'%c',13);
        poslat=fscanf(fid,'%g',1);
        poslon=fscanf(fid,'%g',1);
        depth=fscanf(fid,'%g',1);
        U10=fscanf(fid,'%g',1);
        Udir=fscanf(fid,'%g',1);
        cU = fscanf(fid,'%g',1);
        cdir = fscanf(fid,'%g',1);
        sp2d=fscanf(fid,'%g',[NF,ND]);

        if (index == 0)
            b.name = stn;
            b.lat = poslat;
            b.lon = poslon;
            b.depth = depth;
        end

        b.U10(index+1) = U10;
        b.Udir(index+1) = Udir;
        b.cU(index+1) = cU;
        b.cdir(index+1) = cdir;
        
        b.espt{index+1} = sp2d;
        b.sp1d{index+1} = sum(sp2d')*dtheta;
        b.hs(index+1) = 4*sqrt(trapz(f,b.sp1d{index+1}));
        
        
        % Computing peak direction
        sp2d(:,end+1)=sp2d(:,1);
        [tmp,loc]=max(b.sp1d{index+1});
        b1=trapz(dir,sin(dir)'.*sp2d(loc,:));
        a1=trapz(dir,cos(dir)'.*sp2d(loc,:));
        theta_m = atan2(b1,a1);
        if (theta_m < 0)
            theta_m = theta_m + 2*pi;
        end
        b.dp(index+1) = theta_m;
        b.fp(index+1) = f(loc);
        
        %- Randy adding to structure
%         b.a1(index+1) = a1;
%         b.b1(index+1) = b1;
%         a2=trapz(dir,cos(2.*(dir)'.*sp2d(loc,:)));   %- From Fortran code: ww3_to_sp.f
%         b2=trapz(dir,sin(2.*(dir)'.*sp2d(loc,:)));
%         b.a2(index+1) = a2;
%         b.b2(index+1) = b2;

        index = index + 1;
        
    end
    
end
fclose(fid);



%- Randy compute a1,b1,a2,b2
NT = length(b.time);
icnt=0;
pstnid = '';
a0 = zeros(NF,NT);
a1 = zeros(NF,NT);
b1 = zeros(NF,NT);
a2 = zeros(NF,NT);
b2 = zeros(NF,NT);

zero_energy = [];
energy = zeros(NF,1);
ds_set = zeros(dbins,bands,NT);

%%   Loop over spectal files within the WW3 output

for icnt = 1:NT
    times(icnt) = b.time(icnt);
    %% Read in directional spectrum (m^2/hz-rad)
    %- Convert from cell to matrix
    ds = cell2mat(b.espt(icnt));
    ds = ds';   % (ND x NF)
    for i = 1:bands
        for j = 1:dbins
            ds_set(dir_idx(j),i,icnt) = ds(j,i);
        end
    end

    %% Do the same for the 1D energy
    ds1d = cell2mat(b.sp1d(icnt));
    
%%   When calculating moments, rotate directions by pi so that the resulting 
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
        end
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
        end     %-- End dbins loop 
        b.a0(i,icnt) = a0(i,icnt);
        b.a1(i,icnt) = a1(i,icnt);
        b.b1(i,icnt) = b1(i,icnt);
        b.a2(i,icnt) = a2(i,icnt);
        b.b2(i,icnt) = b2(i,icnt);
        %-- Add in 1D energy
        b.energy(i,icnt) = ds1d(i);
    end     %-- End Freq loop
    
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
data = b;

end

