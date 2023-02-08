function [favg ftk tindex kw] = gflux(sprefix,ivers,ncol,icol,trange,dt,sstep,tstar,tfin)
%  function [favg ftk tindex] = gflux(sprefix,ncol,icol,trange,dt,sstep,tstar,tfin)
%
% Reads ASCII data files produced by GHOST, and computes fluxes, providing
% both a time average over the specified time range, and the individual
% k-dependent fluxes for each time found. If time range contains fluxes with different
% sizes (e.g., due to bootstrapping ), the average and time-dependent fluxes
% will be sized to the maximum found.
%
%           Usage:
% 
%           [favg ftk tindex] = gflux('ktransfer',2,'',dt,sstep,9.0,10.0)
%
%           Input:
%  spref  : prefix of file set (e.g., 'ktransfer', 'hktransfer' from which to read transfers (required)
%  ivers  : GHOST spectral output version. If 0, this is the 'old' version, 
%           for which wavenumber is implicit, and doesn't exist in the file; 
%           if 1, this is the quantity, if requested is read from column 1 of 
%           the file. The quantities ncol and icol should be provided just
%           as they are for ivers=0, as though there is no column of wavenumbers
%  ncol   : no. columns in file. Must be specified correctly (required).
%           This number should not include the column of wave numbers in the file,
%           that exists when iver=1.
%  icol   : if specified, use this column of the file; reqired. The column
%           of wavenumbers (column 0) that appears in the file when ivers=1, 
%           should be negelected when specifying icol.
%  trange : string of form 'start:stop' indices; if  stop is 'end', all files from 
%           start until the end of file list are averaged. If 'trange' is specified, then
%           this range of time indices is used in computing averages. If this string is '', then
%           'tstar' and 'tfin' are used for computing averages, and are then required. 
%  dt     : timestep (required)
%  sstep  : time output index interval for spectra (required)
%  tstar  : start time for averaging, if averaging over time interval
%  tfin   : end time for averaging, if averaging over time interval
%
%           Output:
%  favg   : average of spectra in range, trange
%  ftk    : cell array of the individual flux spectra that went into computing favg,
%           of form ftk{time_index}. Each cell is a function of k.
%  tindex : array of the output index of the spectrum; when multiplied by dt * sstep, it
%           yields the time  at which the spectum is taken
%  kw     : if specified, contains the wave numbers (not required). 
%
favg  = [];
ftk   = {};
tindex = [];

if nargin < 7
  error('Not enough parameters. Do a "help gflux".')
end
if ivers ~= 0 && ivers ~= 1
  error('Incorrect version number');
end

if length(trange) <= 0 
  if nargin < 8
    error('If trange = \'\', then tstar and tend are required. Do a "help gflux".')
  end
  i1     = uint16(tstar/(dt*sstep)+0.50001);
  i2     = uint16(tfin /(dt*sstep)+0.50001);
  trange = sprintf('%d:%d',i1,i2);
end

efiles = sprintf('%s.*',sprefix);
de=dir(efiles); 
ltad=length(de);                    % length of time series
if ltad <= 0 
  error('Files not found');
end
bCheck = true;
if  strcmp(trange,'1:end')
  ib = 1;
  ie = ltad;
  bCheck = false;
end

for i = 1:ltad
  efilename{i} = de(i).name;
end

%display(trange)

[ibe iee] = filerange(efilename, trange);

sform = '';
for n=1:ncol+ivers
  sform  = [sform '%f'];
end

nmax = 0;
nmin = 1e7; 
ftk   = {};
tindex = [];
time   = [];
n = 1;
for i = 1:(iee-ibe+1)
  ie = ibe+i-1;
  fn = efilename{ie}(1:strfind(efilename{ie},'.txt')-1);
  [t r]  = strtok(fn,'.');
% ssp    = textread(efilename{ie});
  fid    = fopen(efilename{ie});
% ssp    = textscan(fid,sform, 'delimiter','\n','CommentStyle',{'#'});
  xdat   = textscan(fid,sform, 'delimiter','\n');
  fclose(fid);
  ssp    = xdat{icol+ivers};
  if ivers == 1
    ktmp = xdat{1};
  end
  nn     = int32(2*(length(ssp)-1));
  ndeal  = nn/2 + 1;
  sie    = ssp(1:ndeal);
  sie    = reshape(sie,numel(sie),1);
  ind    = find(isfinite(sie)==0 ); 
  if nargin <= 6 % make sure time is within set time interval
    ttime = str2num(r(2:end)) * sstep * dt;
  end
  if isempty(ind) 
    tindex(n) = str2num(r(2:end));
    time  (n) = tindex(n) * sstep * dt;
%   for j = 1:length(sie)
%     ftk{n}(j) = -sum(sie(1:j),1);
%   end
    ftk{n} = -1*cumsum(sie,1);
    nmax = max(length(ftk{n}),nmax);
    nmin = min(length(ftk{n}),nmin);
    n = n + 1;
  end
end

favg = zeros(nmax,1);

if nargin <= 6
  tstar = time(1);
  tfin  = time(end);
end
I = find(time >= tstar & time <= tfin);

if isempty(I)
  max(time)
  min(time)
  sstep
  dt
  tstar
  tfin
  error('Invalid averaging interval: time out of range:');
end
NT = numel(I);

nmax = 0;
nmin = 1e7; 
for i = 1:NT
  ie   = I(i);
  si   = deal(ftk{ie});
  nfa  = length(favg); 
  nsi  = length(si); 
  nmax = max(max(nsi),nmax);
  nmin = min(min(nsi),nmin);
  if nsi > nfa
    tmp  = favg;
    favg = zeros(nsi,1);
    favg(1:nfa) = tmp;
    favg = favg + si;
  elseif nsi < nfa
    favg(1:nsi) = favg(1:nsi) + si;
  else
    favg = favg + si;
  end
end
favg = favg ./ NT;

% Get wavenumbers, if requested:
if nargout >= 4
  Dkk = 1.0;
  if ivers == 1
    Dkk = ktmp(2)-ktmp(1); % get (constant) shell width
  end
  kw = [1:length(favg)]*Dkk;
end

favg = Dkk*favg;

% Synch up time-dependent spectral sizes:
if nmax ~= nmin
  for i = 1:NT
    ftk{I(i)} = Dkk*ftk{I(i)};
    si   = deal(ftk{I(i)});
    nsi  = length(si);
    if nsi < nmax
      ftk{I(i)} = zeros(nmax,1);
      ftk{I(i)}(1:nsi) = tmp;
      
    end
  end
end


