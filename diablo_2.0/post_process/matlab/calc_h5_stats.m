% This script illustrates how to load in 3D model data and calculate 
% some basic diagnostics
% Run after readmean.m

filename='/data/proteus/jrt51/KH_test/KH_large/out.h5';

LX=400;
LZ=200;

U=h5read(filename,'/Timestep/U');
V=h5read(filename,'/Timestep/V');
W=h5read(filename,'/Timestep/W');
TH1=h5read(filename,'/Timestep/TH1');
%TH2=h5read(filename,'/Timestep/TH2');

NX=size(U,1); NZ=size(U,3); NY=size(U,2);

x=linspace(0,LX,NX); z=linspace(0,LZ,NZ);

drhodz1=0.0;
drhodz2=0.0;
% Add the background buoyancy gradient
for k=1:NZ
  TH1(:,:,k)=TH1(:,:,k)+drhodz1*z(k);
%  TH2(:,:,k)=TH2(:,:,k)+drhodz2*z(k);
end

% Calculate the x-average velocity
ume=squeeze(mean(U,1));
vme=squeeze(mean(V,1));
wme=squeeze(mean(W,1));
thme1=squeeze(mean(TH1,1));
%thme2=squeeze(mean(TH2,1));

% Calculate correlation terms
for k=1:NZ
  uw_mean(:,k)=mean((U(:,:,k)-U_BT(k)).*(W(:,:,k)-W_BT(k)),1);
  uu_mean(:,k)=mean((U(:,:,k)-U_BT(k)).*(U(:,:,k)-U_BT(k)),1);
  ww_mean(:,k)=mean((W(:,:,k)-W_BT(k)).*(W(:,:,k)-W_BT(k)),1);
  uw_BT(:,k)=trapz(gyf,uw_mean(:,k))/(gyf(end)-gyf(1));
  uu_BT(:,k)=trapz(gyf,uu_mean(:,k))/(gyf(end)-gyf(1));
  ww_BT(:,k)=trapz(gyf,ww_mean(:,k))/(gyf(end)-gyf(1));
end

% calculate the mean shear
dudy=zeros(size(ume));
dwdy=zeros(size(wme));
for j=2:NY-1
    dudy(j,:)=(ume(j+1,:)-ume(j-1,:))/(gyf(j+1)-gyf(j-1));
    dwdy(j,:)=(wme(j+1,:)-wme(j-1,:))/(gyf(j+1)-gyf(j-1));
end

% Calculate the local shear
dUdy=zeros(size(U));
for j=2:NY-1
dUdy(:,j,:)=(U(:,j+1,:)-U(:,j-1,:))/(gyf(j+1)-gyf(j-1));
end


