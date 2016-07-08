% This script illustrates how to combine several mean files
% Before running this, comment out the definition of "base_dir" in readmean.m

nruns=3;
run_dir='../../KH_test';

kindex=1;
for irun=1:nruns
  irun
  base_dir=[run_dir '/run' num2str(irun) '/'];
  readmean;

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

  uw_save(:,kindex:kindex+nk-1)=uw(:,1:nk);
  clear uw;

  grarich_save(:,kindex:kindex+nk-1)=grarich(:,1:nk);
  clear grarich;

  dudy_save(:,kindex:kindex+nk-1)=dudy(:,1:nk);
  clear dudy;

  dthdy_save(:,kindex:kindex+nk-1,1)=dthdy(:,1:nk,1);
  clear dthdy;

  kindex=kindex+nk
  nk

end
