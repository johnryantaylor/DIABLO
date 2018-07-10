%nruns=4;
%run_dir='./double_front_jet_thermocline';
nruns=2;
run_dir='/local/data/public/jrt51/T2016_slip';

tstart_save=4*3600*24;
tend_save=5*3600*24;

kindex=1;
for irun=1:nruns
  irun
  base_dir=[run_dir '/run' num2str(irun)];
  readmean;
  nk  
% Do the folowing for each variable you want to save
  tii_save(kindex:kindex+nk-1)=tii(1:nk);
  clear tii;

  tke_save(:,kindex:kindex+nk-1)=tke(:,1:nk);
  clear tke;

  wrms_save(:,kindex:kindex+nk-1)=wrms(:,1:nk);
  clear wrms;

  vrms_save(:,kindex:kindex+nk-1)=vrms(:,1:nk);
  clear wrms;

  urms_save(:,kindex:kindex+nk-1)=urms(:,1:nk);
  clear urms;

  ume_save(:,kindex:kindex+nk-1)=ume(:,1:nk);
  clear ume;

  wme_save(:,kindex:kindex+nk-1)=wme(:,1:nk);
  clear wme;

  thme_save(:,kindex:kindex+nk-1,1:N_TH)=thme(:,1:nk,1:N_TH);
  clear thme;

  thrms_save(:,kindex:kindex+nk-1,1:N_TH)=thrms(:,1:nk,1:N_TH);
  clear thrms;

  thv_save(:,kindex:kindex+nk-1,1:N_TH)=thv(:,1:nk,1:N_TH);
  clear thv;

  uw_save(:,kindex:kindex+nk-1)=uw(:,1:nk);
  clear uw;

  grarich_save(:,kindex:kindex+nk-1)=grarich(:,1:nk);
  clear grarich;

  dudy_save(:,kindex:kindex+nk-1)=dudy(:,1:nk);
  clear dudy;

  dthdy_save(:,kindex:kindex+nk-1,1:N_TH)=dthdy(:,1:nk,1:N_TH);
  clear dthdy;

  kappa_t_save(:,kindex:kindex+nk-1,1:N_TH)=kappa_t(:,1:nk,1:N_TH);
  clear kappa_t;

  kindex=kindex+nk

end

% Reset nk to the full number of timesteps
nk=length(tii_save);

% Get the time index based on start time
kstart=0;
for k=1:nk
  if (tii_save(k) <= tstart_save)
     kstart=k;
  end
end
if (kstart == 0)
  kstart=1;
end
'Start of time averaging window: ',tii_save(kstart)

% Get the time index based on end time (if defined)
if exist('tend_save')
kend=0;
for k=1:nk
  if (tii_save(k) <= tend_save)
     kend=k;
  end
end
if (kend == 0)
  kend=1;
end
else
kend=nk;
end
'End of time averaging window: ',tii_save(kend)

thv_mean(:,:)=squeeze(trapz(tii_save(kstart:kend),thv_save(:,kstart:kend,:),2));
dthdy_mean(:,:)=squeeze(trapz(tii_save(kstart:kend),dthdy_save(:,kstart:kend,:),2));
kappa_t_mean(:,:)=squeeze(trapz(tii_save(kstart:kend),kappa_t_save(:,kstart:kend,:),2));

