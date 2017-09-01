% This script reads in statistics related to the TKE budget outputted by diablo 
% when using HDF5 output

% User settings....
% Set the run directory
base_dir='/local/data/public/jrt51/testing';
NY=65; % Here, NY should match the value in grid_def.all
N_TH=1; % The number of scalars
Re = 1000; NU=1/Re; % Enter the Reynolds number or viscosity from input.dat
Pr=1;   kappa=NU/Pr; % Prandtl number
RI(1:N_TH)=0.15; % Enter the richardson number for each scalar

% Set the start and end time in code units for start of averaging
tstart=0;

% Set the filename
filename=[base_dir '/mean.h5'];

% Read in the number of samples (store in nk)
file_info=h5info(filename);
att_info=file_info.Groups.Attributes;
nk=att_info.Value

% Preallocate variables for speed
time=zeros(1,nk);
epsilon=zeros(NY,nk);

for k=1:nk
    if (k<10)
        timename=['000' int2str(k)];
    elseif (k<100)
        timename=['00' int2str(k)];
    elseif (k<1000)
        timename=['0' int2str(k)];
    else
        timename=[int2str(k)];
    end
    varname=['/time/' timename];            % TIME
    time(k)=h5read(filename,varname); 
    varname=['/gyf/' timename];             % Y-COORDINATE
    gyf(:)=h5read(filename,varname); 

% Here, read in the statistics calculated in save_stats
% If you add new output variables to save_stats, simply add
% lines here to read in the variables, following the examples below

    varname=['/epsilon/' timename];             % DISSIPATION RATE
    epsilon(:,k)=h5read(filename,varname);
end

% epsilon hasn't yet been multiplied by NU, do that here
epsilon=epsilon*NU;


% Get the time index based on start time
kstart=0;
for k=1:nk
  if (time(k) <= tstart)
     kstart=k;
  end
end
if (kstart == 0)
  kstart=1;
end
%'Start of time averaging window: ',time(kstart)

% Get the time index based on end time (if defined)
if exist('tend')
kend=0;
for k=1:nk
  if (time(k) <= tend)
     kend=k;
  end
end
if (kend == 0)
  kend=1;
end
else 
kend=nk;
end
%'End of time averaging window: ',time(kend)


for j=1:NY
  ume_mean(j)=trapz(time(kstart:kend),ume(j,kstart:kend))/(time(kend)-time(kstart));
  vme_mean(j)=trapz(time(kstart:kend),vme(j,kstart:kend))/(time(kend)-time(kstart));
  wme_mean(j)=trapz(time(kstart:kend),wme(j,kstart:kend))/(time(kend)-time(kstart));
  urms_mean(j)=trapz(time(kstart:kend),urms(j,kstart:kend))/(time(kend)-time(kstart));
  vrms_mean(j)=trapz(time(kstart:kend),vrms(j,kstart:kend))/(time(kend)-time(kstart));
  wrms_mean(j)=trapz(time(kstart:kend),wrms(j,kstart:kend))/(time(kend)-time(kstart));
  dudy_mean(j)=trapz(time(kstart:kend),dudy(j,kstart:kend))/(time(kend)-time(kstart));
  dwdy_mean(j)=trapz(time(kstart:kend),dwdy(j,kstart:kend))/(time(kend)-time(kstart));
  tke_mean(j)=trapz(time(kstart:kend),tke(j,kstart:kend))/(time(kend)-time(kstart));
  uv_mean(j)=trapz(time(kstart:kend),uv(j,kstart:kend))/(time(kend)-time(kstart));
  wv_mean(j)=trapz(time(kstart:kend),wv(j,kstart:kend))/(time(kend)-time(kstart));
  cp_mean(j)=trapz(time(kstart:kend),cp(j,kstart:kend))/(time(kend)-time(kstart));
  omega_x_mean(j)=trapz(time(kstart:kend),omega_x(j,kstart:kend))/(time(kend)-time(kstart));
  omega_y_mean(j)=trapz(time(kstart:kend),omega_y(j,kstart:kend))/(time(kend)-time(kstart));
  omega_z_mean(j)=trapz(time(kstart:kend),omega_z(j,kstart:kend))/(time(kend)-time(kstart));
  if (dudy_mean(j)~=0) 
    nu_t_mean(j)=-uv_mean(j)/dudy_mean(j);
  else
    nu_t_mean(j)=0;
  end
  for n=1:N_TH
    thv_mean(j,n)=trapz(time(kstart:kend),thv(j,kstart:kend,n))/(time(kend)-time(kstart));
    dthdy_mean(j,n)=trapz(time(kstart:kend),dthdy(j,kstart:kend,n))/(time(kend)-time(kstart));
    thrms_mean(j,n)=trapz(time(kstart:kend),thrms(j,kstart:kend,n))/(time(kend)-time(kstart));
    thme_mean(j,n)=trapz(time(kstart:kend),thme(j,kstart:kend,n))/(time(kend)-time(kstart));
    pe_diss_mean(j,n)=trapz(time(kstart:kend),pe_diss(j,kstart:kend,n))/(time(kend)-time(kstart));
    if (dthdy_mean(j,n)~=0) 
      kappa_t_mean(j,n)=-thv_mean(j,n)/dthdy_mean(j,n);
    else
      kappa_t_mean(j,n)=0;
    end 
  end
  shear_mean(j)=trapz(time(kstart:kend),shear(j,kstart:kend))/(time(kend)-time(kstart));
end

for j=2:NY
  gy(j)=(gyf(j)+gyf(j-1))/2;    
end
for j=2:NY-1
  dyf(j)=(gy(j+1)-gy(j));
end
for j=2:NY
  dy(j)=gyf(j)-gyf(j-1);
end

for j=2:NY-1
  if ((gyf(j)-gyf(j-1))~=0)
    ry(j)=(gyf(j+1)-gyf(j))/(gyf(j)-gyf(j-1));
  else
    ry(j)=1.0;
  end
end

% Calculate y-integrated quantities
urms_int=trapz(gyf,urms,1);
vrms_int=trapz(gy,vrms,1);
wrms_int=trapz(gyf,wrms,1);
tke_int=trapz(gyf,tke,1);
hke_int=trapz(gyf,(urms.^2+wrms.^2)/2,1);
vke_int=trapz(gy,vrms.^2/2,1);
for n=1:N_TH
  tpe_int(:,n)=RI(n)*trapz(gyf,thrms(:,:,n).^2,1)./(thme(end,:,n)-thme(1,:,n))/2;
end
thv_int=trapz(gyf,thv,1);
thrms_int=trapz(gyf,thrms,1);

