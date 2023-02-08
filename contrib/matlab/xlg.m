%
% Creates initial positions for NP Lag. particles to be read
% by GHOST.
%

N       = 64;
Np      = 10;
fpref   = 'xlginit';
kmin    = 11;    % kmin of inertial range; should be >= 3
kmax    = 30;    % kmax of inertial range
bdrfixed= 0;     % used fixed separation distance, with random orientation
dr      = 1/12;  % if bdrfixed, the sep. distance in terms of box size
append  = 0;     % append to existing file?

% Create '1st' points with 1/2 total number of particles:

xm   = (N)*[1/kmax 1/kmin];
del  = xm(2) - xm(1);

% Uniform on [0+del,N-del-1]: Avoid boundaries, so 
% we don't have to worry about periodicity:
xl = random('Uniform',del,N-del,Np,3); 
inan = find(isnan(xl));
if length(inan) > 0
  error('xlg: bounds must be increasing');
end

%x:
if bdrfixed == 1
  th  = random('Uniform',0.0,pi  ,Np,1);
  phi = random('Uniform',0.0,2*pi,Np,1);
  dxp = dr*sin(th).*cos(phi)*(N-1);
else
  dxp = random('Uniform',xm(1),xm(2),Np,1)/sqrt(3.0);
end
xl(2:2:Np,1) = xl(1:2:Np,1) + dxp(1:2:Np);
xl(:,1) = mod(xl(:,1),N-1);

%y:
if bdrfixed == 1
  th  = random('Uniform',0.0,pi  ,Np,1);
  phi = random('Uniform',0.0,2*pi,Np,1);
  dxp = dr*sin(th).*sin(phi)*(N-1);
else
  dxp = random('Uniform',xm(1),xm(2),Np,1)/sqrt(3.0);
end
xl(2:2:Np,2) = xl(1:2:Np,2) + dxp(1:2:Np);
xl(:,2) = mod(xl(:,2),N-1);

%z:
if bdrfixed == 1
  th  = random('Uniform',0.0,pi  ,Np,1);
  phi = random('Uniform',0.0,2*pi,Np,1);
  dxp = dr*cos(th)*(N-1);
else
  dxp = random('Uniform',xm(1),xm(2),Np,1)/sqrt(3.0);
end
xl(2:2:Np,3) = xl(1:2:Np,3) + dxp(1:2:Np);
xl(:,3) = mod(xl(:,3),N-1);

% Write initial data file (no header)
if bdrfixed == 1
  fileout = sprintf('%s_dr%f_n%d_np%d.dat',fpref,dr,N,Np); 
else
  fileout = sprintf('%s_k%d_%d_n%d_np%d.dat',fpref,kmin,kmax,N,Np); 
end
dlmwrite(fileout,xl,'delimiter',' ','precision','%.6f');

% Write initial data file (w/ header)
if bdrfixed == 1
  fileout = sprintf('%s_dr%f_n%d_np%d_h.dat',fpref,dr,N,Np); 
else
  fileout = sprintf('%s_k%d_%d_n%d_np%d_h.dat',fpref,kmin,kmax,N,Np); 
end

if append==1
  wun = fopen(fileout,'ta');
else
  wun = fopen(fileout,'wt');
end
if  wun == -1
  error(['File ' fileout ' cannot be opened for writing']);
end
fprintf(wun,'%f\n',Np);
fprintf(wun,'%f\n',0.0);
fprintf(wun,'%f %f %f\n',xl(:,1:3)');
fclose(wun);
warning(['xlg: data written to file ' fileout '.']);


