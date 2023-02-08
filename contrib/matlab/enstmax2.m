 function [zmax t0 i0 trng irng] = enstmax(fn,nwave,bvfreq,omega)
% 
% function [zmax t0 i0 trng irng] = enstmax(fn,eps)
%
% Function enstmat reads GHOST 'balance.txt' file, computes enstrophy max, zmax,
% and the time range that encompasses it, as well as the corresponding index range.
% 
% Usage:
%
%    [zmax t0 i0 trng irng] = enstmax('balance.txt',10,1.3,30)
%
% Arguments:
%  fn   : filename containing enstrophy. Default is 'balance.txt' in current directory.
%  nwave: how many wave periods to retain at peak
% bvfreq: Brunt-Vaissala freq
%  omega: rotation rate

%
% Output: 
%  zmax : max enstrophy
%  t0   : time at enstrophy max
%  i0   : time index at enstrophy max
%  trng : 1d array of rank 1 containing bounding times for peak 
%  irng : 1d array of rank 1 containing bounding indices for peak
%

if ( nargin < 1 )
  fn = 'balance.txt';
  eps = 1.0e-4;
end

if ( nargin < 2 )
  eps = 1.0e-4;
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

wfreq  = sqrt(bvfreq^2); %+ 4*omega^2)
abs(t(i0:end)-t0);
rat   = abs(t(i0:end)-t0)*wfreq;
max(rat)
II    = find(rat >= nwave);

II

irng  = [ i0 max(II) ];
trng  = [t(irng(1)) t(irng(2))];

end
