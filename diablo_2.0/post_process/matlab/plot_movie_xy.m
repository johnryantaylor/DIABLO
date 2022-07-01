% This script shows how to load in 2D slices and make a movie of the simulation output
% Run after readmean.m
LX=30;
NX=128;

x=linspace(0,LX,NX);

filename=[base_dir '/movie.h5'];

drhodx=0;

for k=1:nk
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

varname_th=['/th1_xy/' timename];
varname_u=['/u_xy/' timename];
%varname=['/nu_t_xy/' timename];

A_th=h5read(filename,varname_th);
A_u=h5read(filename,varname_u);

for i=1:size(A_th,1)
   A_th(i,:)=A_th(i,:)+drhodx*x(i);
end

pcolor(x,gyf,A_u'); shading interp;
hold on
contour(x,gyf,A_th',linspace(0,3e-4,40),'w-');
caxis([-0.02 0.02]);
title(['time=' num2str(tii(k)/3600) ' hours']);

axis tight
shading interp
colormap(parula(256));
colorbar
M(k)=getframe(gcf);
clf;

end



