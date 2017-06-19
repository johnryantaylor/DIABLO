% This script reads in statistics outputted by diablo 
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
filename=[base_dir '/mean.h'];

% Read in the number of samples (store in nk)
file_info=h5info(filename);
att_info=file_info.Groups.Attributes;
nk=att_info.Value

% Preallocate variables for speed
time=zeros(1,nk);
ume=zeros(NY,nk); vme=zeros(NY,nk); wme=zeros(NY,nk);
urms=zeros(NY,nk); vrms=zeros(NY,nk); wrms=zeros(NY,nk);
uv=zeros(NY,nk); uw=zeros(NY,nk); wv=zeros(NY,nk);
dudy=zeros(NY,nk); dwdy=zeros(NY,nk); cp=zeros(NY,nk); shear=zeros(NY,nk);
omega_x=zeros(NY,nk); omega_y=zeros(NY,nk); omega_z=zeros(NY,nk);
thme=zeros(NY,nk,N_TH);dthdy=zeros(NY,nk,N_TH); thrms=zeros(NY,nk,N_TH);
thv=zeros(NY,nk,N_TH); pe_diss=zeros(NY,nk,N_TH);

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

    varname=['/ume/' timename];             % MEAN VELOCITIES
    ume(:,k)=h5read(filename,varname);
    varname=['/vme/' timename];
    vme(:,k)=h5read(filename,varname);
    varname=['/wme/' timename];
    wme(:,k)=h5read(filename,varname);
    varname=['/urms/' timename];            % RMS VELOCITIES
    urms(:,k)=h5read(filename,varname);
    varname=['/vrms/' timename];
    vrms(:,k)=h5read(filename,varname);
    varname=['/wrms/' timename];
    wrms(:,k)=h5read(filename,varname);
    varname=['/uv/' timename];              % REYNOLDS STRESSES
    uv(:,k)=h5read(filename,varname);
    varname=['/uw/' timename];
    uw(:,k)=h5read(filename,varname);
    varname=['/wv/' timename];
    wv(:,k)=h5read(filename,varname);
    varname=['/dudy/' timename];            % MEAN VELOCITY GRADIENTS
    dudy(:,k)=h5read(filename,varname);
    varname=['/dwdy/' timename];
    dwdy(:,k)=h5read(filename,varname);
    varname=['/cp/' timename];              % MEAN PRESSURE
    cp(:,k)=h5read(filename,varname);
    varname=['/shear/' timename];           % MEAN SQUARE SHEAR
    shear(:,k)=h5read(filename,varname);
    varname=['/omega_x/' timename];         % RMS VORTICITIES
    omega_x(:,k)=h5read(filename,varname);
    varname=['/omega_y/' timename];
    omega_y(:,k)=h5read(filename,varname);
    varname=['/omega_z/' timename];
    omega_z(:,k)=h5read(filename,varname);
    for n=1:N_TH
        varname=['/thme' num2str(n,'%2.2d') '/' timename];  % MEAN BUOYANCY
        thme(:,k)=h5read(filename,varname);
        varname=['/dthdy' num2str(n,'%2.2d') '/' timename]; % MEAN BUOYANCY GRADIENTS
        dthdy(:,k)=h5read(filename,varname);
        varname=['/thrms' num2str(n,'%2.2d') '/' timename]; % RMS BUOYANCY
        thrms(:,k)=h5read(filename,varname);
        varname=['/thv' num2str(n,'%2.2d') '/' timename];
        thv(:,k)=h5read(filename,varname);
        varname=['/pe_diss' num2str(n,'%2.2d') '/' timename]; % POTENTIAL ENERGY DISSIPATION
        pe_diss(:,k)=h5read(filename,varname);
    end
end

% Compute secondary quantities
% You might want to calculate other quantities, these are a few examples
for k=1:nk
  for j=1:NY
    tke(j,k)=0.5*(urms(j,k)^2.+vrms(j,k)^2.+wrms(j,k)^2.);
    if (dudy(j,k)~=0)
      nu_t(j,k)=-uv(j,k)/dudy(j,k);
    else
      nu_t(j,k)=0;
    end
% Calculate the vertical taylor scale
    if (shear(j,k)~=0) 
      taylor(j,k)=sqrt((ume(j,k)^2.+wme(j,k)^2.+urms(j,k)^2.+wrms(j,k)^2.)/shear(j,k));
    else
      taylor(j,k)=0;
    end
    if (N_TH > 0)
      for n=1:N_TH
% Brunt-Vaisalla (buoyancy) frequency
        brunt(j,k,n)=sqrt(RI(n)*(dthdy(j,k,n))); 
        if (shear(j,k)~=0) 
% Gradient richardson number
          grarich(j,k,n)=brunt(j,k,n)^2./shear(j,k); 
        else
          grarich(j,k,n)=0;
        end
% Resolved turbulent diffusivity
        if (dthdy(j,k,n)~=0)
          kappa_t(j,k,n)=-thv(j,k,n)/dthdy(j,k,n);
        else
          kappa_t(j,k,n)=0;
        end 
% Perturbation potential energy 
        tpe(j,k,n)=RI(n)*thrms(j,k,n).^2/dthdy(j,k,n);
      end

    end
  end
end

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

