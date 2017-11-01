% Run after readmean.m
LX=1000;
NX=512;
LZ=1000;
NZ=512;

filename=[base_dir '/movie.h5'];

% Background density gradient
drhodx=0.00000003;
%drhodx=0.0;
drhodz=0.0;

x=linspace(0,LX,NX);
z=linspace(0,LZ,NZ);

for k=400:400
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
A_th1=h5read(filename,varname);
for j=1:size(A_th1,2)
   A_th1(:,j)=A_th1(:,j)+drhodz*z(j);
end
for i=1:size(A_th1,1)
   A_th1(i,:)=A_th1(i,:)+drhodx*x(i);
end
varname=['/th4_xz/' timename];
A_th4=h5read(filename,varname);
varname=['/th2_xz/' timename];
A_th2=h5read(filename,varname);

varname=['/v_xz/' timename];
A_v=h5read(filename,varname);
varname=['/w_xz/' timename];
A_w=h5read(filename,varname);
varname=['/u_xz/' timename];
A_u=h5read(filename,varname);


subplot(2,2,1)
pcolor(x,z,A_th1'); shading interp;
set(gca,'FontName','Times','FontSize',14);
xlabel('x'); ylabel('z'); title(['b, t=' num2str(floor(tii(k)/3600)) ' hours' ]);
axis equal; axis tight;
caxis([1.4e-4 1.8e-4]);
colormap(jet(256));
colorbar

subplot(2,2,2)
pcolor(x,z,A_v'); shading interp;
set(gca,'FontName','Times','FontSize',14);
xlabel('x'); ylabel('z'); title(['v']);
axis equal; axis tight;
caxis([-0.01 0.01]);
colormap(jet(256));
colorbar

subplot(2,2,3)
pcolor(x,z,A_th2'); shading interp;
set(gca,'FontName','Times','FontSize',14);
xlabel('x'); ylabel('z'); title(['th2']);
axis equal; axis tight;
caxis([0 1]);
colormap(jet(256));
colorbar

subplot(2,2,4)
pcolor(x,z,A_th4'); shading interp;
set(gca,'FontName','Times','FontSize',14);
xlabel('x'); ylabel('z'); title(['th4']);
axis equal; axis tight;
caxis([0 1]);
colormap(jet(256));
colorbar

colorbar
M(k)=getframe(gcf);

%clf;

end




