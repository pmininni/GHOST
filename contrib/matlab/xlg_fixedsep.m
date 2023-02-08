%
% Creates initial positions for NP Lag. particles to be read
% by GHOST. Produces a set of ICs with per-ordered-pair separations
% that are a range of fixed distances, as opposed to being within
% a range of wavenumber bands.
%

N       = 512;
Np      = 1.005e6;
fpref   = 'xlginit';
dr(1)   = 1/30;

%dr(1)   = 1/15;  % if bdrfixed, the sep. distance in terms of box size
%dr(1)   = 1/20;  % if bdrfixed, the sep. distance in terms of box size
nk      = 1;     % no. k-sub-ranges or fixed sep distances
 
% Create '1st' points with 1/2 total number of particles;
% In the following, we only care about every other particle location.
% Uniform base points on [0+del,N-del-1]: Avoid boundaries, so 
% we don't have to worry about periodicity:
del = 0.0;
for j = 1:nk
  xm   = (N-1)*dr(j); % put in terms of requred particle coords.
  del  = max(del,xm);
end
xl = random('Uniform',del,N-del,Np,3); 


% Now, modify each pair according to initial relative position type:
is = 1;
for j = 1:nk
  
  Npi = Np/(nk);
  if j <= rem(Np,nk)
    Npi = Npi + 1;
  end
 
  inan = find(isnan(xl));
  if length(inan) > 0
    error('xlg: bounds must be increasing');
  end
  
  %x:
  th  = myrandom(Npi,1,[0.0,pi/2]);
  phi = myrandom(Npi,1,[0.0,pi/2]);
  dxp = dr(j)*(N-1);
  for k=1:3
    if k == 1
      dx = dxp.* sin(th).*cos(phi);
     elseif k == 2
      dx = dxp.* sin(th).*sin(phi);
     else
      dx = dxp.* cos(th);
     end
    xl(is+1:2:is+Npi-1,k) = xl(is:2:is+Npi-2,k) + dx(1:2:Npi-1);
    xl(is:Npi,k) = mod(xl(is:Npi,k),N-1);
  end
  
  is = is + Npi;

end

% Write initial data file (no header)
ftmp = fpref
for j = 1:nk
  ftmp = sprintf('%s_k%dr%f',ftmp,j,dr(j))
end
fileout = sprintf('%s_n%d_np%d.dat',ftmp,N,Np); 
dlmwrite(fileout,xl,'delimiter',' ','precision','%.6f');

% Write initial data file (w/ header)
ftmp = fpref;
for j = 1:nk
 ftmp = sprintf('%s_k%dr%f',ftmp,j,dr(j)); 
end
fileout = sprintf('%s_n%d_np%d_h.dat',ftmp,N,Np); 
wun = fopen(fileout,'wt');
if  wun == -1
  error(['File ' fileout ' cannot be opened for writing']);
end
fprintf(wun,'%f\n',Np);
fprintf(wun,'%f\n',0.0);
fprintf(wun,'%f %f %f\n',xl(:,1:3)');
fclose(wun);
warning(['xlg: data written to file ' fileout '.']);


