function [nehavg neavg nhavg neht net nht t] = ehindex(sdir, dt, sstep, krange, srange, tstar, tfin, iplot)
%
% Computes average spectral e+h slope over time for tstar < t < tfin,
% using spectral range krange, and time index range, srange.
% The time is computed from the timestep, and the spectral output
% time cycle, sstep.
%
% Usage:
%  [avg slope t] = ehindex('ABCEk7R13', 1.25e-4, 1000, '15:45', '2:110', 5.0);
%
% sdir  : directory name (in)
% dt    : time step (in)
% sstep : spectral output cycle interval (ij)
% krange: string of form 'kbeg:kend' giving beg. and end k for fit to spectrum (in)
% srange: cycles between spectral output (in)
% tstar : time s.t. average spectral indexn is computed for tstar < t < tfin (in)
% tfin  : time s.t. average spectral indexn is computed for tstar < t < tfin (in)
% iplot : (optional) if 1, do plot; else don't. Default is to do no plot.
%
% Returns:
%  nehavg : average spectral index for E*H: E*H is computed as an average over time 
%           tfin > t > tstar; then fit is done on this avg spectrum over krange to get nehavg
%  neavg  : average spectral index for E: E is computed as an average over time 
%           tfin > t > tstar; then fit is done on this avg spectrum over krange to get neavg
%  nhavg  : average spectral index for H: H is computed as an average over time 
%           tfin > t > tstar; then fit is done on this avg spectrum over krange to get nhavg
%  neht   : spetral index of E*H  as a function of t
%  net    : spectral index of E as function of t; not computed from a time avg of spectra
%  nht    : spectral index of H as function of t; not computed from a time avg of spectra
%  t      : times corresponging to 'neht', etc. 

iplot = 1;
if nargin < 6
  error('ehindex: Not enough input arguments. Do a "help ehindex"');
end
if nargin < 7 
  iplot = 0;
end

icolon = strfind(krange,':');
if isempty(icolon)
  error('ehindex: invalid "krange"');
end
km(1) = str2num(krange(1:icolon-1));
km(2) = str2num(krange(icolon+1:length(krange)));


[eh eht et ht tindex] = ehavg(srange);
k  = [1:numel(eh)]';
nt = numel(tindex);
t = [];
neht= [];

% compute time-varying quantities:
for i = 1:nt
  t  (i) = dt*sstep*tindex(i);
  kx = log10(k(km(1):km(2)));
  ey = log10(eht{i}(km(1):km(2)));
  p = polyfit(kx,ey,1);
  neht (i) =  p(1);

  ey = log10(et{i}(km(1):km(2)));
  p = polyfit(kx,ey,1);
  net(i) =  p(1);

  ey = log10(ht{i}(km(1):km(2)));
  p = polyfit(kx,ey,1);
  nht(i) =  p(1);
end

% compute indices based on time-averaged spectra:
e  = zeros(numel(eht{1}),1);
h  = zeros(numel(eht{1}),1);
eh = zeros(numel(eht{1}),1);
I = find(t >= tstar & t <= tfin);
for i = 1:numel(I)
  j = I(i);
  e  =  e  + et {j};
  h  =  h  + et {j};
  eh =  eh + eht{j};
end
e  = e  / nt;
h  = h  / nt;
eh = eh / nt;

kx     = log10(k(km(1):km(2)));
ey     = log10(eh(km(1):km(2)));
p      = polyfit(kx,ey,1);
nehavg =  p(1);

ey    = log10(e(km(1):km(2)));
p     = polyfit(kx,ey,1);
neavg =  p(1);

hy    = log10(h(km(1):km(2)));
p     = polyfit(kx,hy,1);
nhavg =  p(1);

%I = find(t>tstar);
%nehavg = mean(ns(I));

%if iplot > 0
%  plot(t,ns); 
%  stext = sprintf('Time-averaged slope: %f\n',mean(ns(I)) );
%  title('e+h slope');
%  xlabel('t');
%  ylabel('EH slope');
%end
