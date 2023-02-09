function [st spect tindex kw] = spectavgn(sprefix, ivers, trange, ncol, icol, nfield)
%
% function [st spect tindex] = spectavgn(sprefix, ivers, trange, ncol, icol, nfield)
%
% Reads ASCII data files produced by GHOST, and averages over specified
% index range, if there is one. If time range contains spectra with different
% sizes (e.g., due to bootstrapping ), the average and time-dependent spectra
% will be sized to the maximum found.
%
%	Usage:
%
%           [spec_avg spec_t tindex k] = spectavgn('khelicity',0,'1:end');
%
%           Input:
%  sprefix: file prefix
%  ivers  : GHOST spectral output version. If 0, this is the 'old' version, 
%           for which wavenumber is implicit, and doesn't exist in the file; 
%           if 1, this is the quantity, if requested is read from column 1 of 
%           the file. The quantities ncol and icol should be provided just
%           as they are for ivers=0, as though there is no column of wavenumbers
%  trange : string of indices to read, of form '2 3 4 50:100  200 300'
%           listing individual time indices, or ranges. The spaces may be 
%           replaced with either ',' or ';'. Alternatively, '1:end' may
%           be specified to read all indices in the directory.
%           Spectra corresponding to specifed indices will be averaged.
%  ncol   : no. columns in file. Must be specified correctly; default is 1. 
%           This number should not include the column of wave numbers in the file,
%           that exists when ivers=1.
%         
%  icol   : if specified, use this column of the file. Default is 1. The column
%           of wavenumbers (column 0) that appears in the file when ivers=1, 
%           should be negelected when specifying icol.
%  nfield : length of time index field (default is 4).
%
%           Output:
%  st     : average of spectra in range, trange
%  spec   : cell array of the individual spectra that went into computing st
%  tindex : array of the output index of the spectrum; when multiplied by dt * sstep, it
%           yields the time  at which the spectrum is taken
%  kw     : if specified, contains the wave numbers. 
%
st     = [];
spect  = {};
tindex = [];

if nargin<2
  error('Must provide a filename prefix and version number');
end
if nargin < 3
  trange = '1:end';
  ncol   = 1;
  icol   = 1;
  nfield = 4;
end
if nargin < 4
  ncol   = 1;
  icol   = 1;
  nfield = 4;
end
if nargin < 5
  icol   = 1;
  nfield = 4;
end
if nargin < 6
  nfield = 4;
end
sprec = sprintf('%s%d%s','%0',nfield,'d');

if ivers ~= 0 && ivers ~= 1
  error('Incorrect version number');
end

if nargout < 1
  error('At least one output argument is required')
end

sfiles = sprintf('%s.*', sprefix);
d=dir(sfiles);
ltad=length(d) ;                    % length of time series
if ltad <= 0 
  error('Files not found');
end

ir = strfind(trange,'1:end');
if  ~isempty(ir)
  fnrange = [1:ltad];
else
  fnrange = iparse(trange);
end

sform = '';
for n=1:ncol+ivers
  sform  = [sform '%f']; 
end

for i = 1:ltad
  filename{i} = d(i).name;
end

fid = fopen(filename{1});
sc = textscan(fid,sform, 'delimiter','\n','CommentStyle',{'#'});
fclose(fid);
s  = cell2mat(sc);
st = zeros([ numel(s(:,1)) 1]);
if nargout >= 2 
  spect = {};
end
if nargout >= 2 
  tindex = [];
end
if nargout > 3 
  kw = [];
end

ktmp = [];

nmax = 0;
nspect = 0;
nmin = 1e7;
NT = 0;
k = 0;
for j = 1:length(fnrange)
    i = fnrange(j); 
%   si = textread(filename{i});
    sfile = sprintf(['%s.' sprec '.txt'],sprefix,i);
    fid= fopen(sfile);
    if fid < 0
      warning(['spectavgn: file not found: ' sfile]);
      continue;
    end
    sc = textscan(fid,sform, 'delimiter','\n','CommentStyle',{'#'});
    fclose(fid);
    NT = NT + 1;
    si = cell2mat(sc);
    nspect = nspect + 1;
    ktmp = [1:length(si)];
    if ivers == 1
      ktmp = si(:,1);
    end
    if nargout >= 2 
      spect{nspect} = si(:,icol+ivers) ;
      nmax = max(length(si(:,icol+ivers)),nmax);
      nmin = min(length(si(:,icol+ivers)),nmin);
    end
    if nargout == 3 
      fn = sfile(1:strfind(sfile,'.txt')-1);
      [t r] = strtok(fn,'.');
      ttmp = str2num(r(2:end));
      tindex(NT) = ttmp;
    end
    nsi = length(si(:,icol+ivers));
    nst = length(st);
    if nsi > nst
      tmp = st;
      st = zeros(nsi,1);
      st(1:nst) = tmp;
      st = st + si(:,icol+ivers);
    elseif nsi < nst
      st(1:nsi) = st(1:nsi) + si(:,icol+ivers);
    else
      st = st + si(:,icol+ivers);
    end
    if nargout == 3 
      fn = sfile(1:strfind(sfile,'.txt')-1);
      [t r] = strtok(fn,'.');
      ttmp = str2num(r(2:end));
      tindex(NT) = ttmp;
    end
    % Get wavenumbers, if requested:
    if nargout >= 4 
      Dkk = 1.0;
      if ivers == 1
        Dkk = ktmp(2)-ktmp(1); % get (constant) shell width
      end
      kw = [1:length(st)]*Dkk;
    end

end
st = st ./ NT;

% Synch up time-dependent spectral sizes:
if nargout >= 2 & nmax ~= nmin

  II = find(tindex > 0);
  tt = tindex(II);

  for i = 1:length(spect)
    si   = deal(spect{i});
    nsi  = length(si);
    if nsi < nmax
      tmp  = spect{i};
      spect{i} = zeros(nmax,1);
      spect{i}(1:nsi) = tmp;
    end
  end
end

