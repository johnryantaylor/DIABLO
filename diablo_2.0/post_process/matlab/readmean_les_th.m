function readmean_les_th(base_dir)
% This script reads in diagnostics from the optional scalar LES model

% Reads in statistics outputted by diablo

turbo=dlmread([base_dir 'mean_les_th.txt']);

% Set the domain size
NY=96;
NYM=NY-1;

% Enter the viscosity
NU=1/960;

% Set the starting time in code units for start of averaging
tstart=0;

% Determine the number of records in the file based on its length
%nk=floor(length(turbo(:,1))/(NY+1))
% Don't reset nk, use value from readmean.m

row=1;
for k=1:nk
  tii(k)=turbo(row,2);
  dt(k)=turbo(row,3);
  row=row+1;
  for j=1:NY
%    gyf(j)=turbo(row,2);
    lambda1(j,k)=turbo(row,3); 
    lambda2(j,k)=turbo(row,4); 
    lambda3(j,k)=turbo(row,5); 
    c_dyn_th(j,k)=turbo(row,6);
    row=row+1;
  end
end

for k=1:nk
  for j=1:NY
    thu_sgs(j,k)=lambda1(j,k);
    thv_sgs(j,k)=lambda2(j,k);
    thw_sgs(j,k)=lambda3(j,k);
    kappa_t_sgs(j,k)=-thv_sgs(j,k)/dthdy(j,k); %readmean must be done first
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
'Start of time average: ',tii(kstart)

for j=1:NY
  thu_sgs_mean(j)=mean(thu_sgs(j,kstart:nk));
  thv_sgs_mean(j)=mean(thv_sgs(j,kstart:nk));
  thw_sgs_mean(j)=mean(thw_sgs(j,kstart:nk));
  kappa_t_sgs_mean(j)=-thv_sgs_mean(j)/dthdy_mean(j);
end

end
