% Run after readmean.m
LZ=4;
NZ=16;

filename=[base_dir '/movie.h5'];

% Background density gradient
drhodz1=0.0;

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

varname=['/th1_zy/' timename];
A_th1=h5read(filename,varname);
varname=['/v_zy/' timename];
A_w=h5read(filename,varname);
varname=['/u_zy/' timename];
A_u=h5read(filename,varname);
varname=['/w_zy/' timename];
A_v=h5read(filename,varname);
varname=['/nu_t_zy/' timename];
A_nu_t=h5read(filename,varname);

% Add a background gradient
for j=1:size(A_th1,1)
  A_th1(j,:)=A_th1(j,:)+drhodz1*z(j);
end

A_u_save(:,:,k)=A_u;
A_th1_save(:,:,k)=A_th1;
A_w_save(:,:,k)=A_w;

surf(z,gyf,zeros(size(A_nu_t')),A_nu_t','EdgeColor','none'),view(0,90);
hold on
caxis([0 1e-4]);
%contour(z/1e3,gyf-gyf(end),A_v',10,'w-');

xlabel('y'); ylabel('z'); title(['t=' num2str(tii(k))]);

axis tight
shading interp
colorbar
M(k)=getframe(gcf);


clf;

end


