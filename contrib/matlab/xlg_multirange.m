%
% Creates initial positions for NP Lag. particles to be read
% by GHOST. Produces a set of ICs with per-ordered-pair separations
% that fall within specified k-ranges.
%

N       = 512;
Np      = 1.0e6;
fpref   = 'xlginit';
%kmin    = [1:2:40];   % kmin of inertial range; should be >= 3
for j = 3:7
  kmin(j-2) = 2^(j-1)
end
kmax    = kmin+1; ;    % kmax of inertial range
nk      = length(kmin);% no. k-sub-ranges or fixed sep distances
bdrfixed= 0;     % used fixed separation distance, with random orientation; uses nk ranges tho
dr      = 1/15;  % if bdrfixed, the sep. distance in terms of box size
 
% Create '1st' points with 1/2 total number of particles;
% In the following, we only care about every other particle location.
% Uniform base points on [0+del,N-del-1]: Avoid boundaries, so 
% we don't have to worry about periodicity:
del = 0.0;
for j = 1:nk
  xm   = (N)*[1/kmax(j) 1/kmin(j)]; % put in terms of requred particle coords.
  del  = max(del,xm(2) - xm(1));
end
xl = random('Uniform',del,N-del,Np,3); 


% Now, modify each pair according to initial relative position type:
is = 1;
nn = 0;
for j = 1:nk
  
  xm   = (N-1)*[1/kmax(j) 1/kmin(j)]; % put in terms of requred particle coords.

  Npi = floor(Np/(nk));
  if j <= rem(Np,nk)
    Npi = Npi + 1;
  end

  nn = nn + Npi
 
  inan = find(isnan(xl));
  if length(inan) > 0
    error('xlg: bounds must be increasing');
  end
  
  %x:
  if bdrfixed == 1
    th  = myrandom(Npi,1,[0.0,pi/2]);
    phi = myrandom(Npi,1,[0.0,pi/2]);
    dxp = dr*sin(th).*cos(phi)*(N-1);
  else
    th  = myrandom(Npi,1,[0.0,pi/2]);
    phi = myrandom(Npi,1,[0.0,pi/2]);
%   dxp = myrandom(Npi,1,[xm(1),xm(2)],)/sqrt(3.0);
    dxp = myrandom(Npi,1,[xm(1),xm(2)]);
  end
% xl(is+1:2:is+Npi-1,1) = xl(is:2:is+Npi-2,1) + dxp(1:2:Npi-1);
% xl(is:Npi,1) = mod(xl(is:Npi,1),N-1);
  for k=1:3
    if k == 1
      dx = dxp.* sin(th).*cos(phi);
     elseif k == 2
      dx = dxp.* sin(th).*sin(phi);
     else
      dx = dxp.* cos(th);
     end
    xl(is+1:2:is+Npi-1,k) = xl(is:2:is+Npi-2,k) + dx(1:2:Npi-1);
    xl(is+1:2:is+Npi-1,k) = mod(xl(is+1:2:is+Npi-1,k),N-1);
%   xl(is:Npi,k) = mod(xl(is:Npi,k),N-1);
  end
  
% %y:
% if bdrfixed == 1
%   th  = random('Uniform',0.0,pi  ,Npi/2,1);
%   phi = random('Uniform',0.0,2*pi,Npi/2,1);
%   dxp = dr*sin(th).*sin(phi)*(N-1);
% else
%   dxp = random('Uniform',xm(1),xm(2),Npi,1)/sqrt(3.0);
% end
% xl(is+1:2:is+Npi-1,2) = xl(is:2:is+Npi-2,2) + dxp(1:2:Npi-1);
% xl(is:Npi,2) = mod(xl(is:Npi,2),N-1);

% %z:
% if bdrfixed == 1
%   th  = random('Uniform',0.0,pi  ,Npi/2,1);
%   phi = random('Uniform',0.0,2*pi,Npi/2,1);
%   dxp = dr*cos(th)*(N-1);
% else
%   dxp = random('Uniform',xm(1),xm(2),Npi,1)/sqrt(3.0);
% end
% xl(is+1:2:is+Npi-1,3) = xl(is:2:is+Npi-2,3) + dxp(1:2:Npi-1);
% xl(is:Npi,3) = mod(xl(is:Npi,3),N-1);

  is = is + Npi;

end

% Write initial data file (no header)
if bdrfixed == 1
  fileout = sprintf('%s_dr%f_n%d_np%d.dat',fpref,dr,N,Np); 
else

  ftmp = fpref
  for j = 1:nk
  ftmp = sprintf('%s_k%dr%d_%d',ftmp,j,kmin(j),kmax(j))
  end
  fileout = sprintf('%s_n%d_np%d.dat',ftmp,N,Np); 
end
dlmwrite(fileout,xl,'delimiter',' ','precision','%.6f');

% Write initial data file (w/ header)
if bdrfixed == 1
  fileout = sprintf('%s_dr%f_n%d_np%d_h.dat',fpref,dr,N,Np); 
else
  ftmp = fpref;
  for j = 1:nk
  ftmp = sprintf('%s_k%dr%d_%d',ftmp,j,kmin(j),kmax(j)); 
  end
  fileout = sprintf('%s_n%d_np%d_h.dat',ftmp,N,Np); 
end

wun = fopen(fileout,'wt');
if  wun == -1
  error(['File ' fileout ' cannot be opened for writing']);
end
fprintf(wun,'%f\n',Np);
fprintf(wun,'%f\n',0.0);
fprintf(wun,'%f %f %f\n',xl(:,1:3)');
fclose(wun);
warning(['xlg: data written to file ' fileout '.']);


