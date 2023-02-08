function [st spect tindex] = spectavg(sprefix, trange, ncol, icol)
%
% function [st spect tindex] = spectavg(sprefix, trange, ncol, icol)
%
% Reads ASCII data files produced by GHOST, and averages over specified
% index range, if there is one. If time range contains spectra with different
% sizes (e.g., due to bootstrapping ), the average and time-dependent spectra
% will be sized to the maximum found.
%
%	Usage:
%
%           [spec_avg spec_t tindex] = spectavg('khelicity','1:end');
%
%           Input:
%  sprefix: file prefix
%  trange : string of form 'start:stop' indices; if  stop is 'end', all files from 
%           start until the end of file list are averaged
%  ncol   : no. columns in file. Must be specified correctly; default is 1.
%  icol   : if specified, use this column of the file. Default is 1.
%
%           Output:
%  st     : average of spectra in range, trange
%  spec   : cell array of the individual spectra that went into computing st
%  tindex : array of the output index of the spectrum; when multiplied by dt * sstep, it
%           yields the time  at which the spectrum is taken
%
st     = [];
spect  = {};
tindex = [];

if nargin<1
  error('Must provide a filename prefix');
end
if nargin < 2
  trange = '1:end';
  ncol   = 1;
  icol   = 1;
end
if nargin < 3
  ncol   = 1;
  icol   = 1;
end
if nargin < 4
  icol   = 1;
end

sfiles = sprintf('%s.*', sprefix);
d=dir(sfiles);
ltad=length(d) ;                    % length of time series
if ltad <= 0 
  error('Files not found');
end
bCheck = true;
if  strcmp(trange,'1:end')
  fnrange = [1, ltad];
  bCheck = false;
end

sform = '';
for n=1:ncol
  sform  = [sform '%f']; 
end

for i = 1:ltad
  filename{i} = d(i).name;
end

if bCheck 
  [t r] = strtok(trange,':');
  ibeg = str2num(t);
  if isempty(ibeg) 
    error('invalid range initialization');
  end
  ib = 0;
  bfound = false;
  while ib<ltad && ~bfound 
    ib = ib + 1;
    fn = filename{ib}(1:strfind(filename{ib},'.txt')-1);
    [t r] = strtok(fn,'.');
    index = str2num(r(2:end));
    bfound = index >= ibeg;
  end
  if ~bfound
    error('unable to find start');
  end
  
  [t r] = strtok(trange,':');
  if strcmp(strtrim(r(2:end)),'end')  
    iend = ltad;
  else
    iend = str2num(r(2:end));
    if isempty(iend) 
      error('invalid range termination');
    end
  end
  ie = ltad+1;
  bfound = false;
  while ie>=1 && ~bfound 
    ie = ie - 1;
    fn = filename{ie}(1:strfind(filename{ie},'.txt')-1);
    [t r] = strtok(fn,'.');
    index = str2num(r(2:end));
    bfound = index <= iend;
  end
  if ~bfound
    error('unable to find end');
  end
  fnrange = [ib, ie];
end % bCheck 

if ( fnrange(2) < fnrange(1) ) 
  error('invalid file range');
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

nmax = 0;
nmin = 1e7;
NT   = 0;
for i = fnrange(1):fnrange(2)
%   si = textread(filename{i});
    fid= fopen(filename{i});
    sc = textscan(fid,sform, 'delimiter','\n','CommentStyle',{'#'});
    fclose(fid);
    si = cell2mat(sc);
    NT = NT + 1;
    if nargout >= 2 
      spect{i-fnrange(1)+1} = si(:,icol) ;
      nmax = max(length(si(:,icol)),nmax);
      nmin = min(length(si(:,icol)),nmin);
    end
    if nargout == 3 
      fn = filename{i}(1:strfind(filename{i},'.txt')-1);
      [t r] = strtok(fn,'.');
      tindex(i-fnrange(1)+1) = str2num(r(2:end));
    end
    nsi = length(si(:,icol));
    nst = length(st);
    if nsi > nst
      tmp = st;
      st = zeros(nsi,1);
      st(1:nst) = tmp;
      st = st + si(:,icol);
    elseif nsi < nst
      st(1:nsi) = st(1:nsi) + si(:,icol);
    else
      st = st + si(:,icol);
    end
end
st = st ./ NT;

% Synch up time-dependent spectral sizes:
if nargout >= 2 & nmax ~= nmin
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


