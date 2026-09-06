function data = rghost(filein, N, rank, isz, sformat)
%
% Reads binary GHOST data, and stores in local variable data.
% This function assumes data in a cubic box with N^3 grid points.
%
%  Usage:
%      data = rghost(filename, N, isz, 'ieee-be');
%
%  filein  : input file to read
%  N       : data cube dimension
%  rank    : problem rank (2 or 3)
%  isz     : data size (in bytes: either 4 or 8)
%  sformat : data format of file: 'ieee-be' or 'ieee-le' for big-endian or little
%            endian if isz=4, or 'ieee-be.l64', 'ieee-le.l64' if isz=8.
%
if nargin < 2
  error('Input file name, and cube dimension and must be specified');
end
if nargin < 3
  rank = 3;
  isz = 4;
  sformat = 'ieee-be';
end
if nargin < 4
  isz = 4;
  sformat = 'ieee-be';
end
if nargin < 5
  sformat = 'ieee-be';
end
if rank ~= 2 && rank ~= 3 
  error('rank must be 2 or 3');
end

%sformat
lun =fopen(filein,'r',sformat);
if  lun == -1
  error(['File ' filein ' cannot be opened for reading']);
end
[fn permission thismachineformat] = fopen(lun); %machine format is for machine that reads, not that wrote
if ~strcmp(permission,'r')
   error('Invalid file')
end
ssize = sprintf('real*%d',isz);
if strcmp(ssize,'real*4' )
  zsize = 'single';
elseif strcmp(ssize,'real*8')
  zsize = 'double';
else
  error('Type must be "real*4" or "real*8"');
end

isize = N^rank

data = zeros(isize,1,zsize);
data = fread(lun, isize, ssize);
fclose(lun);
if isempty(data)
  error(['File ' filein ' corrupted']);
end
