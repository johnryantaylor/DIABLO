fname_old='./Ret45_Ri0p05/Adjust_down_end/start.h5';
fname_new='./spot.h5';

NX_old=1024;
NZ_old=1024;
NY_old=65;

NX_new=512;
NZ_new=512;
NY_new=65;

istart=513; iend=1024;
jstart=1; jend=65;
kstart=513; kend=1024;

fname_new='spot.h5';

var_list={'P';'TH1';'U';'V';'W'};

time=h5readatt(fname_old,'/Timestep','Time');

for i=1:length(var_list)

  h5create(fname_new,char(strcat('/Timestep/',var_list(i))),[NX_new NY_new NZ_new],'Datatype','double','ChunkSize',[NX_new 1 NZ_new/8]);
end

h5writeatt(fname_new,'/','Resolution',floor([NX_new NY_new NZ_new]));

h5writeatt(fname_new,'/Timestep','Time',time);

disp('Done creating file structure');

for i=1:length(var_list)
  disp(['Reading variable, ' char(strcat('/Timestep/',var_list(i)))]);
  var_old=h5read(fname_old,char(strcat('/Timestep/',var_list(i))));
  var_new=var_old(istart:iend,jstart:jend,kstart:kend);
  h5write(fname_new,char(strcat('/Timestep/',var_list(i))),var_new);
end 


