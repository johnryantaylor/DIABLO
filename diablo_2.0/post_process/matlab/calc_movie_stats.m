% Illustrates how to calculate diagnostics from 2D "movie" output files

% Run after readmean.m
% Set the following parameters
LX=10;
NX=128;
LZ=3.14;
NZ=64;

% Point to the movie file
filename='../../KH_test/movie.h5'

% Background density gradient
drhodz=0.0;

x=linspace(0,LX,NX);
z=linspace(0,LZ,NZ);

for k=1:1:nk
k
  if (k<10)
    timename=['000' int2str(k)];
  elseif (k<100) 
    timename=['00' int2str(k)];
  elseif (k<1000)
    timename=['0' int2str(k)];
  else
    timename=[int2str(k)];
  end

varname=['/th1_xz/' timename];
A_th=h5read(filename,varname);
varname=['/u_xz/' timename];
A_u=h5read(filename,varname);
varname=['/w_xz/' timename];
A_w=h5read(filename,varname);
varname=['/v_xz/' timename];
A_v=h5read(filename,varname);

for j=1:size(A_th,2)
  A_th(:,j)=A_th(:,j)+drhodz*z(j);
end

u(:,:,k)=A_u(:,:);
v(:,:,k)=A_v(:,:);
w(:,:,k)=A_w(:,:);
th(:,:,k)=A_th(:,:);

% Calculate statistics

% Mean quantities
ume_movie(:,k)=mean(A_u,1);
wme_movie(:,k)=mean(A_w,1);
vme_movie(:,k)=mean(A_v,1);
thme_movie(:,k)=mean(A_th,1);

% Derivatives
for j=2:NZ-1
dudz_movie(j,k)=(ume_movie(j+1,k)-ume_movie(j-1,k))/(z(j+1)-z(j-1));
dvdz_movie(j,k)=(vme_movie(j+1,k)-vme_movie(j-1,k))/(z(j+1)-z(j-1));
dwdz_movie(j,k)=(wme_movie(j+1,k)-wme_movie(j-1,k))/(z(j+1)-z(j-1));
dthdz_movie(j,k)=(thme_movie(j+1,k)-thme_movie(j-1,k))/(z(j+1)-z(j-1));
end

% RMS
for j=1:NZ
urms_movie(j,k)=sqrt(mean((A_u(:,j)-ume_movie(j,k)).^2));
vrms_movie(j,k)=sqrt(mean((A_v(:,j)-vme_movie(j,k)).^2));
wrms_movie(j,k)=sqrt(mean((A_w(:,j)-wme_movie(j,k)).^2));
thrms_movie(j,k)=sqrt(mean((A_th(:,j)-thme_movie(j,k)).^2));
end

% Reynolds stress
for j=1:NZ
  uw_movie(j,k)=mean((A_u(:,j)-ume_movie(j,k)).*(A_w(:,j)-wme_movie(j,k)));
  uv_movie(j,k)=mean((A_u(:,j)-ume_movie(j,k)).*(A_v(:,j)-vme_movie(j,k)));
  wv_movie(j,k)=mean((A_w(:,j)-wme_movie(j,k)).*(A_v(:,j)-vme_movie(j,k)));
  thv_movie(j,k)=mean((A_th(:,j)-thme_movie(j,k)).*(A_v(:,j)-vme_movie(j,k)));
  thu_movie(j,k)=mean((A_th(:,j)-thme_movie(j,k)).*(A_u(:,j)-ume_movie(j,k)));
  thw_movie(j,k)=mean((A_th(:,j)-thme_movie(j,k)).*(A_w(:,j)-wme_movie(j,k)));
end

end

% horizontal eddy viscosity/diffusivity
nu_t_movie(2:NZ-1,:)=-uw_movie(2:NZ-1,:)./dudz_movie(2:NZ-1,:);
kappa_t_movie(2:NZ-1,:)=-thw_movie(2:NZ-1,:)./dthdz_movie(2:NZ-1,:);



