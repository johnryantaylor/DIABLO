% This script reads in information about the turbulent kinetic energy (TKE) 
% User settings....
% Set the grid and domain size in the y-direction
%NP=input('Enter the number of processes in the y-direction: ');
NP=2;
%NY_S=input('Insert the number of points per process in the y-direction: ');
NY_S=33;
% Enter the number of scalars
N_TH=1;
% Set the base directory where the mean*.txt files are located
base_dir='../../KH_test'
% Enter the viscosity from input.dat
NU=1e-3;
% Enter the Prandtl number
Pr=1;
kappa=NU/Pr;
% Enter the richardson number for each scalar
RI(1:N_TH)=0.15;
% Set the start and end time in code units for start of averaging
tstart=0; 
%tend=10033; % If tend isn't defined, tend will default to the final time

% Loop over all processes

for proc=1:NP

% Define input files
turbo=dlmread([base_dir 'tke' num2str(proc) '.txt']);

% Set the start and end indices corresponding to the full size array in y
if (proc==1) 
  jstart=2
  jend=NY_S-1
elseif (proc~=NP)
  jstart=(NY_S-2)*(proc-1)+1
  jend=(NY_S-2)*(proc-1)+(NY_S-2)
else
  jstart=(NY_S-2)*(proc-1)+1
  jend=(NY_S-2)*(proc-1)+(NY_S-2)
end

% Calculate the global size of NY
NY=(NY_S)*(NP-1)+(NY_S-1)

% Determine the number of records in the file based on its length
nk=floor(length(turbo(:,1))/(NY_S+1))

row=1;

for k=1:nk
  tii(k)=turbo(row,2);
  dt(k)=turbo(row,3);
  row=row+1;
  if (proc~=1) 
    row=row+1;
  end
  for j=jstart:jend
    gyf(j)=turbo(row,2);
    epsilon(j,k)=turbo(row,3);
    row=row+1;
  end
end

end  % End loop over NP procs

% Row 1 was a ghost cell, set to zero
epsilon(1,:)=0;

% Compute secondary quantities
for k=1:nk
  for j=1:NY-3
    eta(j,k)=abs(epsilon(j,k))^(-0.25d0)*NU^(3/4);
  end
end

% Buoyancy Reynolds number
Re_b=epsilon(:,1:nk)/NU./(squeeze(dthdy(:,1:nk,1)*RI(1)));
loz=sqrt(epsilon./(dthdy*RI(1)).^(3/2));

% Get the time index based on start time
kstart=0;
for k=1:nk
  if (tii(k) <= tstart)
     kstart=k;
  end
end
if (kstart == 0)
  kstart=1;
end
'Start of time averaging window: ',tii(kstart)

% Get the time index based on end time (if defined)
if exist('tend')
kend=0;
for k=1:nk
  if (tii(k) <= tend)
     kend=k;
  end
end
if (kend == 0)
  kend=1;
end
else 
kend=nk;
end
'End of time averaging window: ',tii(kend)

for j=1:NY
  epsilon_mean(j)=trapz(tii(kstart:kend),epsilon(j,kstart:kend))/(tii(kend)-tii(kstart));
  eta_mean(j)=trapz(tii(kstart:kend),eta(j,kstart:kend))/(tii(kend)-tii(kstart));
end

% Calculate y-integrated quantities
epsilon_int=trapz(gyf,epsilon,1);
eta_int=trapz(gyf,eta,1);

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



