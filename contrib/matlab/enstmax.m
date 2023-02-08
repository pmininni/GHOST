 function [zmax t0 i0 trng irng] = enstmax(fn,eps,dosym)
% 
% function [zmax t0 i0 trng irng] = enstmax(fn,eps,dosym)
%
% Function enstmat reads GHOST 'balance.txt' file, computes enstrophy max, zmax,
% and the time range that encompasses it, as well as the corresponding index range.
% 
% Usage:
%
%    [zmax t0 i0 trng irng] = enstmax('balance.txt')
%
% Arguments:
%  fn   : filename containing enstrophy. Default is 'balance.txt' in current directory.
%  eps  : relative value of constancy at peak for which enstrophy
%         is considered constant, as a percentage of zmax. Default is 0.05.
%  dosym: If 1, find time interval symmetrically located about time peak; else,
%         find time interval only _after_ peak. Default is 0.
%
% Output: 
%  zmax : max enstrophy
%  t0   : time at enstrophy max
%  i0   : time index at enstrophy max
%  trng : 1d array of rank 2 containing bounding times for peak 
%  irng : 1d array of rank 2 containing bounding indices for peak into balance.txt
%

if ( nargin < 1 )
  fn = 'balance.txt';
  eps = 0.05;
  dosym = 0;
end

if ( nargin < 2 )
  eps =0.05;
  dosym = 0;
end

if ( nargin < 3 )
  dosym = 0;
end

bal = gfclean(fn,5,[2 3]);

z = bal(:,3);
t = bal(:,1);

%plot(t,z)
%pause(2)

[zmax J] = max(abs(z));

zmax  = zmax(1);

i0    = J(1);
t0    = t(i0);

if dosym > 0 
  % Find time range symmetrically about peak:
  tmin  = t0-1;
  tmax  = t0+1;
  ii    = find(t >= tmin & t <= tmax);
  imin  = min(ii);
    imax  = max(ii) ;
  dz    = abs(z(imin:imax)-zmax);
    rat   = dz  ./ abs(zmax);
  itrng = find(abs(rat) < eps) + imin - 1;
  irng  = [ min(itrng) max(itrng) ];
  trng  = [t(irng(1)) t(irng(2))];

else 
  % Check variation only _after_ peak of enstrophy:
  dz    = abs(z(i0:end)-zmax);
  rat   = dz  ./ abs(zmax);

  itrng = find(abs(rat) < eps) + i0 - 1;

  irng  = [ min(itrng) max(itrng) ];
  trng  = [t(irng(1)) t(irng(2))];
end


end
