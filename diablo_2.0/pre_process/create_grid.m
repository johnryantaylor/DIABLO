% This file creates a grid in the y-direction for input to DIABLO
% Written by John Taylor, 10/23/2005
% The grid definitions are produced as follows:
%   1 First, define the G grid based on a specified function with
%     loops over j=0,N+1 where the 0 and N+1 correspond to the bottom
%     and top wall values respectively.  Since matlab is unable to store
%     elements in the zero index, all arrays are indexed starting at 1.
%   2 Then, define the fractional grid GF halfway between neighboring
%     G grid locations
%   3 Stretch both GF and G so that GF(0) and GF(N) now correspond
%     to the upper and lower walls.

% Select the type of grid function to be used
disp('1) High resolution at both ends (tanh stretching function)');
disp('2) High resolution in center (tanh stretching function)');
disp('3) High resolution at both ends (starting at 0)');
GRID_TYPE=input('Select a grid type: ');

% Set the dimensions of the grid 
% This should match (NX,NY,NZ) in grid_def
N=input('Enter the grid size, N: ');

% Enter the domain size unless we are reading the grid
L=input('Enter the domain size: ');
% Select the stretching parameter
CS=input('Enter the stretching parameter, CS (e.g. 1.75): ');

if (GRID_TYPE==1) 
  % Closed Channel
  for J=1:N+1
    G(J+1)=(L/2.0)*tanh(CS*((2.0*(J-1))/(N)-1.0))/tanh(CS);
  end
elseif (GRID_TYPE==2)
% Grid for a shear layer, with higher resolution in the center
  for J=1:N+1
    G(J+1)=(L/2.0)*(CS*((2.0*(J-1))/N-1.0)^3+(2.0*(J-1))/N-1.0)/(CS+1);
  end
elseif (GRID_TYPE==3)
  % Closed Channel
  for J=1:N+1
    G(J+1)=(L/2.0)*tanh(CS*((2.0*(J-1))/(N)-1.0))/tanh(CS)+L/2.0;
  end
else
   disp('Error, entered grid type unknown');
end

% The following lines are done for all cases
% First, define the half (fractional) grid points
for J=1:N
  GF(J+1)=(G(J+1)+G(J+2))/2.0;
end

% Scale both grids to place GF(2) and GF(NY+1) at the walls
gf_lower=GF(2);
gf_upper=GF(N+1);
g_lower=G(2);
g_upper=G(N+2);
for J=1:N
  GF(J+1)=GF(J+1)*(g_upper-g_lower)/(gf_upper-gf_lower);
end
for J=1:N+1
  G(J+1)=G(J+1)*(g_upper-g_lower)/(gf_upper-gf_lower);
end
gf_lower=GF(2);
gf_upper=GF(N+1);
% And shift the grids if necessary
shift=g_lower-gf_lower;
for J=1:N
  GF(J+1)=GF(J+1)+shift;
end
for J=1:N+1
  G(J+1)=G(J+1)+shift;
end
gf_lower=GF(2);
gf_upper=GF(N+1);

% Now, write the grid to file
count=1;
matwrite(1)=N;
for j=1:N+1
  count=count+1;
  matwrite(count)=G(j+1);
end
for j=1:N
  count=count+1;
  matwrite(count)=GF(j+1);
end

% For testing the grid stretching
for j=2:N-1
  r(j)=(GF(j+1)-GF(j))/(GF(j)-GF(j-1));
end
disp('The maximum grid-stretching ratio is:'),max(r)
disp(' ');
dlmwrite('grid.txt',matwrite','\t');
disp('grid.txt has been written to the current directory');
disp('  re-name this file to either xgrid.txt, ygrid.txt, or zgrid.txt');
disp('  and move to the desired case directory');





