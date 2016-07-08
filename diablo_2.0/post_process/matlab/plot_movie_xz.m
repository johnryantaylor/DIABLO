% Run after readmean.m
LX=30;
NX=128;
LZ=4;
NZ=16;

filename=[base_dir '/movie.h5'];

% Background density gradient
drhodx=0.0;
drhodz=0.0;

x=linspace(0,LX,NX);
z=linspace(0,LZ,NZ);

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

varname=['/th1_xz/' timename];
A_th=h5read(filename,varname);
for j=1:size(A_th,2)
   A_th(:,j)=A_th(:,j)+drhodz*z(j);
 end
varname=['/v_xz/' timename];
A_v=h5read(filename,varname);
varname=['/w_xz/' timename];
A_w=h5read(filename,varname);


subplot(2,1,1)
pcolor(x,z,A_th'); shading interp;
set(gca,'FontName','Times','FontSize',14);
xlabel('x'); ylabel('z'); title(['b, t=' num2str(tii(k)) ]);
caxis([-1 1]);
colormap(jet(256));
colorbar

subplot(2,1,2)
pcolor(x,z,A_v'); shading interp;
set(gca,'FontName','Times','FontSize',14);
xlabel('x'); ylabel('z'); title(['v, t=' num2str(tii(k)) ]);
caxis([-0.15 0.15]);
colormap(jet(256));
colorbar

colorbar
M(k)=getframe(gcf);

clf;

end




