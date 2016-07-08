% This script illustrates how to combine several runs to make a movie
% Run after readmean.m
LX=30;
NX=128;
LZ=4;
NZ=16;

% Background density gradient
drhodz=0.0;
drhodz=0.0;

% Set the number of individual runs
nruns=3;

x=linspace(0,LX,NX);
z=linspace(0,LZ,NZ);

movie_frame=1;

for run_num=1:nruns

filename=[base_dir '/run' num2str(run_num) '/movie.h5'];

file_info=h5info(filename);
att_info=file_info.Groups.Attributes;
nk=att_info.Value

dataset=file_info.Groups.Datasets;
for i=1:nk
  time(i)=dataset(i).Attributes.Value;
end

for k=1:nk
    
  tii_save(movie_frame)=time(k);
    
k
  if (k<10)
    tii_savename=['000' int2str(k)];
  elseif (k<100) 
    tii_savename=['00' int2str(k)];
  elseif (k<1000)
    tii_savename=['0' int2str(k)];
  else
    tii_savename=[int2str(k)];
  end

varname_u=['/v_xz/' tii_savename];
varname_th=['/th2_xz/' tii_savename];

A_u=h5read(filename,varname_u);
A_th=h5read(filename,varname_th);

for j=1:size(A_th,1)
  A_th(:,j)=A_th(:,j)+drhodz*z(j);
end

A_u_save(:,:,movie_frame)=A_u(:,:);
A_th_save(:,:,movie_frame)=A_th(:,:);


surf(x,z,zeros(size(A_th')),A_th','EdgeColor','none'),view(0,90); shading interp;
 
colorbar
axis tight
 
M(movie_frame)=getframe(gcf);
movie_frame=movie_frame+1;
clf;

end

end


