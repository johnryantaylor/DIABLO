% This script reads in statistics outputted by diablo 
% when using HDF5 output

% User settings....
% Set the run directory
filename='/local/data/public/jrt51/testing/mean_les.h5';
NY=65; % Here, NY should match the value in grid_def.all
N_TH=1; % The number of scalars
Re = 1000; NU=1/Re; % Enter the Reynolds number or viscosity from input.dat
Pr=1;   kappa=NU/Pr; % Prandtl number
RI(1:N_TH)=0.15; % Enter the richardson number for each scalar

% Set the start and end time in code units for start of averaging
tstart=0;

% Read in the number of samples (store in nk)
file_info=h5info(filename);
att_info=file_info.Groups.Attributes;
nk_les=att_info.Value

% Preallocate variables for speed
time_les=zeros(1,nk_les);
nu_sgs=zeros(NY,nk_les); eps_sgs1=zeros(NY,nk_les);

for k=1:nk_les
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
    time_les(k)=h5read(filename,varname); 
    varname=['/gyf/' timename];             % Y-COORDINATE
    gyf(:)=h5read(filename,varname);
    varname=['/nu_sgs/' timename];             % MEAN VELOCITIES
    nu_sgs(:,k)=h5read(filename,varname);
    varname=['/eps_sgs1/' timename];
    eps_sgs1(:,k)=h5read(filename,varname);
end

% Get the time index based on start time
kstart=0;
for k=1:nk_les
  if (time_les(k) <= tstart)
     kstart=k;
  end
end
if (kstart == 0)
  kstart=1;
end
%'Start of time averaging window: ',time_les(kstart)

% Get the time index based on end time (if defined)
if exist('tend')
kend=0;
for k=1:nk_les
  if (time_les(k) <= tend)
     kend=k;
  end
end
if (kend == 0)
  kend=1;
end
else 
kend=nk_les;
end
%'End of time averaging window: ',time_les(kend)


for j=1:NY
  nu_sgs_mean(j)=trapz(time_les(kstart:kend),nu_sgs(j,kstart:kend))/(time_les(kend)-time_les(kstart));
  eps_sgs1_mean(j)=trapz(time_les(kstart:kend),eps_sgs1(j,kstart:kend))/(time_les(kend)-time_les(kstart));
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
