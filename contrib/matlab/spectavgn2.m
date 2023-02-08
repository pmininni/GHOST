function [spavg tindex] = spectavgn2(sprefix, nr, trange, isz, sformat, nfield)
%
% function [spavg tindex] = spectavgn2(sprefix, [nx nz], trange, isz, sform, nfield)
%
% Reads 2D spectra data files produced by GHOST, and averages over specified
% index range, if there is one. If time range contains spectra with different
% sizes (e.g., due to bootstrapping ), the average and time-dependent spectra
% will be sized to the maximum found in that range.
%
%	Usage:
%
%           [spec_avg tindex] = spectavgn2('sspec2D','1:end');
%
%           Input:
%  sprefix: file prefix
%  nr     : vector of size 2 containing full x, z dimensions: [nx nz];
%  trange : string of indices to read, of form '2 3 4 50:100  200 300'
%           listing individual time indices, or ranges. The spaces may be 
%           replaced with either ',' or ';'. Alternatively, '1:end' may
%           be specified to read all indices in the directory.
%           Spectra corresponding to specifed indices will be averaged.
%           '1:end' is the default.
%  isz    : data size (in bytes: either 4 or 8, e.g.); isz=4 is default.
%  sformat: data format of file: 'ieee-be' or 'ieee-le' for big-endian or 
%           little-endian; 'ieee-le' is the default.
%  nfield : length of time index field (default is 3).
%
%           Output:
%  spavg  : average of spectra in range, trange
%  tindex : (optional) array of the output index of the spectrum; when multiplied by dt * sstep, it
%           yields the time  at which the spectrum is taken
%

if nargin<2
  error('Must provide a filename prefix, rank');
end

if length(nr) < 2
  error('Rank of nr must be >= 2.');
end
% function [spavg tindex] = spectavgn2(sprefix, [nx nz], trange, isz, sform, nfield)

if nargin < 3
  trange = '1:end';
  isz    = 4;
  sformat= 'ieee-le';
  nfield = 3;
end
if nargin < 4
  isz    = 4;
  sformat= 'ieee-le';
  nfield = 3;
end
if nargin < 5
  sformat= 'ieee-le';
  nfield = 3;
end
if nargin < 6
  nfield = 3;
end

sprec = sprintf('%s%d%s','%0',nfield,'d');

sfiles = sprintf('%s.*', sprefix);
d=dir(sfiles);
ltad=length(d) ;                    % length of time series
if ltad <= 0 
  spavg = zeros(nr(1),nr(2)/2+1);
  return;
% error('Files not found');
end

ir = strfind(trange,'1:end');

if  ~isempty(ir)
  fnrange = [1:ltad];
else
  fnrange = iparse(trange);
end

for i = 1:ltad
  filename{i} = d(i).name;
end

if nargout <= 1 
  error('Must have at least one return parameter');
end

if nargout >= 2 
  tindex = [];
end

nspect = 0;
NT     = 0;
nx     = nr(1)/2+1;
nz     = nr(2)/2+1;
accum   = zeros(nx, nz);
for j = 1:length(fnrange)
    i = fnrange(j);
    sfile    = sprintf(['%s.' sprec '.out'],sprefix,i);
    sspec    = r2dspec(sfile,[nx nz],'b',isz,sformat);
    sspec    = sspec';
    accum    = accum + sspec;
    NT = NT + 1;
    nspect = nspect + 1;
    if nargout == 2 
      fn = sfile(1:strfind(sfile,'.out')-1);
      [t r] = strtok(fn,'.');
      ttmp = str2num(r(2:end));
      if nargout >= 2 
        tindex(NT) = ttmp;
      end
    end
end

if NT <= 0
  spavg = zeros(nx,nz);
  return;
% error('No spectra read!');
end
spavg = accum ./ NT;
