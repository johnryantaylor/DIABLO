% Reads in statistics outputted by diablo 
% User settings....
% Set the run directory if it hasn't already been defined
if (~exist('base_dir'))
  base_dir='/data/oceanus/jrt51/Crowe/bz0p001';
end
run_dir='/data/oceanus/jrt51/Crowe';
% Set the grid and domain size in the y-direction
%NP=input('Enter the number of processes in the y-direction: ');
NP=1;
%NY_S=input('Insert the number of points per process in the y-direction: ');
NY_S=51;
% Enter the number of scalars
N_TH=1;
% Enter the viscosity from input.dat
NU=1e-6;
% Enter the Prandtl number
Pr=1;
kappa=NU/Pr;
% Enter the richardson number for each scalar
RI(1:N_TH)=0.15;
% Set the start and end time in code units for start of averaging
tstart=0; 
%tend=999; % If tend isn't defined, tend will default to the final time

% Loop over all processes

for proc=1:NP

% Define input files
turbo=dlmread([base_dir '/mean' num2str(proc) '.txt']);
if (N_TH>0)
  turbo_th=dlmread([base_dir '/mean_th' num2str(proc) '.txt']);
end

% Set the start and end indices corresponding to the full size array in y
if (proc==1) 
  jstart=1;
  jend=NY_S;
elseif (proc~=NP)
  jstart=(NY_S-1)*(proc-2)+NY_S+1;
  jend=(NY_S-1)*(proc-1)+NY_S;
else
  jstart=(NY_S-1)*(proc-2)+NY_S+1;
  jend=(NY_S-1)*(proc-1)+NY_S;
end

% Calculate the global size of NY
NY=(NY_S-1)*(NP-1)+(NY_S);

% Determine the number of records in the file based on its length
nk=floor(length(turbo(:,1))/(NY_S+2));


row=1;
for k=1:nk
  tii(k)=turbo(row,2);
  dt(k)=turbo(row,3);
  row=row+1;
  ubulk(k)=turbo(row,1);
  row=row+1;

  if (proc>1) row=row+1; end %Skip the bottom row which we already read in

  for j=jstart:jend
    gyf(j)=turbo(row,2);
    ume(j,k)=turbo(row,3);
    vme(j,k)=turbo(row,4);
    wme(j,k)=turbo(row,5);
    urms(j,k)=turbo(row,6);
    vrms(j,k)=turbo(row,7);
    wrms(j,k)=turbo(row,8);
    uv(j,k)=turbo(row,9);
    uw(j,k)=turbo(row,10);
    wv(j,k)=turbo(row,11);
    dudy(j,k)=turbo(row,12);
    dwdy(j,k)=turbo(row,13);
    cp(j,k)=turbo(row,14);
    shear(j,k)=turbo(row,15);
    omega_x(j,k)=turbo(row,16);
    omega_y(j,k)=turbo(row,17);
    omega_z(j,k)=turbo(row,18);
    row=row+1;
  end
end

% Now read in scalars
if (N_TH>0)
% Determine the number of records in the file based on its length
  nk_th=floor(length(turbo_th(:,1))/(N_TH*(NY_S)+2));
else
  nk_th=0;
end

% Sometimes if a run was halted during writing of the mean files the
% number of timesteps in mean.txt and mean_th.txt differ.  In this case,
% use the smallest number
nk_th=min(nk,nk_th);
nk=min(nk,nk_th);

row=1;
for k=1:nk_th
  tii(k)=turbo_th(row,2);
  dt(k)=turbo_th(row,3);
  row=row+1;
  ubulk(k)=turbo_th(row,1);
  row=row+1;
  for n=1:N_TH
  if (proc>1) row=row+1; end %Skip the bottom row which we already read in
  for j=jstart:jend
    thme(j,k,n)=turbo_th(row,3);
    dthdy(j,k,n)=turbo_th(row,4); % Add one here if a background was subtracted
    thrms(j,k,n)=turbo_th(row,5);
    thv(j,k,n)=turbo_th(row,6); 
    if (RI(n) ~= 0) 
      pe_diss(j,k,n)=turbo_th(row,7)*NU/RI(n);
    else
      pe_diss(j,k,n)=turbo_th(row,7)*NU;
    end
    row=row+1;
  end
  end
end

end  % End loop over NP procs


% Compute secondary quantities
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
        brunt(j,k,n)=sqrt(RI(n)*(dthdy(j,k,n))); 
        if (shear(j,k)~=0) 
          grarich(j,k,n)=brunt(j,k,n)^2./shear(j,k); 
        else
          grarich(j,k,n)=0;
        end
        if (dthdy(j,k,n)~=0)
          kappa_t(j,k,n)=-thv(j,k,n)/dthdy(j,k,n);
        else
          kappa_t(j,k,n)=0;
        end 
        tpe(j,k,n)=RI(n)*thrms(j,k,n).^2/dthdy(j,k,n);
      end

    end
  end
end

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
%'Start of time averaging window: ',tii(kstart)

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
%'End of time averaging window: ',tii(kend)


for j=1:NY
  ume_mean(j)=trapz(tii(kstart:kend),ume(j,kstart:kend))/(tii(kend)-tii(kstart));
  vme_mean(j)=trapz(tii(kstart:kend),vme(j,kstart:kend))/(tii(kend)-tii(kstart));
  wme_mean(j)=trapz(tii(kstart:kend),wme(j,kstart:kend))/(tii(kend)-tii(kstart));
  urms_mean(j)=trapz(tii(kstart:kend),urms(j,kstart:kend))/(tii(kend)-tii(kstart));
  vrms_mean(j)=trapz(tii(kstart:kend),vrms(j,kstart:kend))/(tii(kend)-tii(kstart));
  wrms_mean(j)=trapz(tii(kstart:kend),wrms(j,kstart:kend))/(tii(kend)-tii(kstart));
  dudy_mean(j)=trapz(tii(kstart:kend),dudy(j,kstart:kend))/(tii(kend)-tii(kstart));
  dwdy_mean(j)=trapz(tii(kstart:kend),dwdy(j,kstart:kend))/(tii(kend)-tii(kstart));
  tke_mean(j)=trapz(tii(kstart:kend),tke(j,kstart:kend))/(tii(kend)-tii(kstart));
  uv_mean(j)=trapz(tii(kstart:kend),uv(j,kstart:kend))/(tii(kend)-tii(kstart));
  wv_mean(j)=trapz(tii(kstart:kend),wv(j,kstart:kend))/(tii(kend)-tii(kstart));
  cp_mean(j)=trapz(tii(kstart:kend),cp(j,kstart:kend))/(tii(kend)-tii(kstart));
  omega_x_mean(j)=trapz(tii(kstart:kend),omega_x(j,kstart:kend))/(tii(kend)-tii(kstart));
  omega_y_mean(j)=trapz(tii(kstart:kend),omega_y(j,kstart:kend))/(tii(kend)-tii(kstart));
  omega_z_mean(j)=trapz(tii(kstart:kend),omega_z(j,kstart:kend))/(tii(kend)-tii(kstart));
  if (dudy_mean(j)~=0) 
    nu_t_mean(j)=-uv_mean(j)/dudy_mean(j);
  else
    nu_t_mean(j)=0;
  end
  for n=1:N_TH
    thv_mean(j,n)=trapz(tii(kstart:kend),thv(j,kstart:kend,n))/(tii(kend)-tii(kstart));
    dthdy_mean(j,n)=trapz(tii(kstart:kend),dthdy(j,kstart:kend,n))/(tii(kend)-tii(kstart));
    thrms_mean(j,n)=trapz(tii(kstart:kend),thrms(j,kstart:kend,n))/(tii(kend)-tii(kstart));
    thme_mean(j,n)=trapz(tii(kstart:kend),thme(j,kstart:kend,n))/(tii(kend)-tii(kstart));
    pe_diss_mean(j,n)=trapz(tii(kstart:kend),pe_diss(j,kstart:kend,n))/(tii(kend)-tii(kstart));
    if (dthdy_mean(j,n)~=0) 
      kappa_t_mean(j,n)=-thv_mean(j,n)/dthdy_mean(j,n);
    else
      kappa_t_mean(j,n)=0;
    end 
  end
  shear_mean(j)=trapz(tii(kstart:kend),shear(j,kstart:kend))/(tii(kend)-tii(kstart));
end

% Calculate y-integrated quantities
urms_int=trapz(gyf,urms,1);
vrms_int=trapz(gyf,vrms,1);
wrms_int=trapz(gyf,wrms,1);
tke_int=trapz(gyf,tke,1);
hke_int=trapz(gyf,(urms.^2+wrms.^2)/2,1);
vke_int=trapz(gyf,vrms.^2/2,1);
for n=1:N_TH
  tpe_int(1:nk,n)=RI(n)*trapz(gyf,thrms(:,1:nk,n).^2,1)./(thme(end,1:nk,n)-thme(1,1:nk,n))/2;
end
thv_int=trapz(gyf,thv,1);
thrms_int=trapz(gyf,thrms,1);
 
for j=2:NY
  gy(j)=(gyf(j)+gyf(j-1))/2;    
end

% Optionally, get GY from hdf5 grid file
%gy=h5read([run_dir '/grid.h5'],'/grids/y');

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


