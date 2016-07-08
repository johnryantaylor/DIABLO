function ah=h5plotxz(namefile,nVar,yc,varargin)
    
    gridfile='grid.h5';

    xl=201;zl=100;

    if length(varargin)>=1;xl=varargin{1};end
    if length(varargin)>=2;zl=varargin{2};end
    
    res=double(h5readatt(namefile,'/','Resolution'));
        
    NX=res(1);NZ=res(3);
    
    xgrid=[0:NX-1]*xl/NX;    
    zgrid=[0:NZ-1]*zl/NZ;    
    
    if exist(gridfile,'file')
        y=h5read(gridfile,'/grids/y');
        yf=0.5*(y(1:end-1)+y(2:end));
        [err,yc]=min(abs(yf-yc));
        yf(yc)
        disp([' Plotting xz-plane at y= ',num2str(yf(yc)), ... 
              ' ( error= ',num2str(err),')']);
    else
        disp([' Plotting xz-plane at position ' num2str(yf(yc))]);    
    end
    
    Data=h5read(namefile,['/Timestep/' nVar],[1 yc 1],[NX 1 NZ]);
%    fh=figure();
    ah=surf(xgrid,zgrid,double(squeeze(Data))');view(2);shading interp;set(gcf,'Renderer','Zbuffer');
 	colormap(bone(256));
    axis equal;
    axis([0 xl 0 zl]); 
%    axis equal; tight;
