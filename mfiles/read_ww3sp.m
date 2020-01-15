% M-file to read in data from a 2D spectral file in WW3 V2.22 format
% J. Henrique Alves, March, 2000
% Modified into a function form by A Chawla, June 2011
% https://polar.ncep.noaa.gov/mmab/papers/tn302/MMAB_302.pdf


function b = read_ww3sp(filename)

if (exist(filename,'file'))
    fid=fopen(filename);
else
    fprintf(1,'ERROR : File Not found !!');
    return;
end
% Scan basic properties of spectra from file
dum=fscanf(fid,'%c',[23]);
NF =fscanf(fid,'%g',[1]);
ND =fscanf(fid,'%g',[1]);
dum=fscanf(fid,'%s',[1]);
dum=fscanf(fid,'%c',[33]);
f=fscanf(fid,'%g',NF);
dir=fscanf(fid,'%g',ND);
fprop=f(2)/f(1);
fdum=[f(1)/fprop;f];
df=diff(fdum);
dtheta=abs(dir(2)-dir(1));


b.f = f;
b.dir = dir;
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
        b.a1(index+1) = a1;
        b.b1(index+1) = b1;
        a2=trapz(dir,cos(2.*(dir)'.*sp2d(loc,:)));   %- From Fortran code: ww3_to_sp.f
        b2=trapz(dir,sin(2.*(dir)'.*sp2d(loc,:)));
        b.a2(index+1) = a2;
        b.b2(index+1) = b2;

        index = index + 1;
        
    end
    
end
fclose(fid);
return;
end