function readmean_les(base_dir)
% This script reads in diagnostics from the LES model
% Reads in statistics outputted by diablo

turbo=dlmread([base_dir 'mean_les.txt']);

% Set the domain size
NY=32;
NYM=NY-1;

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
    tau11(j,k)=turbo(row,3); 
    tau22(j,k)=turbo(row,4); 
    tau33(j,k)=turbo(row,5); 
    tau12(j,k)=turbo(row,6); 
    tau13(j,k)=turbo(row,7); 
    tau23(j,k)=turbo(row,8); 
    nu_t_sgs(j,k)=turbo(row,9);
    nu_t_test(j,k)=turbo(row,10);
    row=row+1;
  end
end

% This part was treated implicitly and was not saved in tau_ij
for k=1:nk
  for j=1:NY
    tau12(j,k)=tau12(j,k)-nu_t_test(j,k)*dudy(j,k);
    tau23(j,k)=tau23(j,k)-nu_t_test(j,k)*dwdy(j,k);
  end
end

for k=1:nk
  for j=1:NY
    urms_sgs(j,k)=sqrt(abs(tau11(j,k)));
    vrms_sgs(j,k)=sqrt(abs(tau22(j,k)));
    wrms_sgs(j,k)=sqrt(abs(tau33(j,k)));
    uv_sgs(j,k)=tau12(j,k);
    uw_sgs(j,k)=tau13(j,k);
    vw_sgs(j,k)=tau23(j,k);
    if (dudy(j,k)~=0) 
      nu_t_sgs(j,k)=-uv_sgs(j,k)/dudy(j,k);  %readmean.m must be done first
    else
      nu_t_sgs(j,k)=0;
    end
    tke_sgs(j,k)=0.5*(tau11(j,k)+tau22(j,k)+tau33(j,k));
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
  urms_sgs_mean(j)=mean(urms_sgs(j,kstart:nk));
  vrms_sgs_mean(j)=mean(vrms_sgs(j,kstart:nk));
  wrms_sgs_mean(j)=mean(wrms_sgs(j,kstart:nk));
  uv_sgs_mean(j)=mean(uv_sgs(j,kstart:nk));
  uw_sgs_mean(j)=mean(uw_sgs(j,kstart:nk));
  vw_sgs_mean(j)=mean(vw_sgs(j,kstart:nk));
  tke_sgs_mean(j)=mean(tke_sgs(j,kstart:nk));
  nu_t_sgs_mean(j)=-uv_sgs_mean(j)/dudy_mean(j);
end

end
