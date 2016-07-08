% Run after readmean.m
LX=10;
NX=30;
LZ=4;
NZ=16;

% Background density gradient
drhodz=0.0;
drhodz=0.0;

nruns=2;

x=linspace(0,LX,NX);
z=linspace(0,LZ,NZ);

movie_frame=1;

for run_num=1:nruns

filename=[base_dir '/run' num2str(run_num) '/movie.h5'];

file_info=h5info(filename);
att_info=file_info.Groups.Attributes;
nk=att_info.Value

for k=1:2:nk
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

varname=['/v_zy/' timename];
varname_th=['/th_zy/' timename];

A=h5read(filename,varname);
A_th=h5read(filename,varname_th);

for j=1:size(A_th,1)
  A_th(j,:)=A_th(j,:)+drhodz*z(j);
end

surf(x/1e3,gyf-gyf(end),zeros(size(A')),A','EdgeColor','none'),view(0,90); shading interp;
set(gca,'FontName','Times','FontSize',14);

colorbar
axis tight
M(movie_frame)=getframe(gcf);

movie_frame=movie_frame+1;
clf;

end

end

%movie2avi(M,'movie_xz_all.avi');

