function [data np time] = rglag_binary(sfile,irank,nbytes,machform,hdonly)
%
% Reads GHOST particle data from binary formatted file
%
%  Usage:
%      [data np t] = rglag_binary('glag.010.txt');
%
%  Inputs:
%
%  sfile    : input file (required)
%  irank    : rank of each record. Default is 3.
%  nbytes   : byte size (4, 8) of real data; default is 4
%  machform : machine format (endian-ness) of real data: 'ieee-be' for
%             big-endian, or 'ieee-le' for little-endian; default is 'ieee-le'.
%  hdonly   : = 1: fill header info only (np and time; data==null); if 0, neglect; default is 0.

%  Outputs:
%
%  data    : particle data in (irank,np) array with x,y,z in each record
%  np      : no. particles required, as found in first file
%  time    : time stamp of file
%
if nargin < 1
  error('Input file name prefix at least! Do a "help rglag".');
end
if nargin < 2
  irank    = 3;
  nbytes = 4;
  machform = 'ieee-le'
  hdonly   = 0;
end
if nargin < 3
  nbytes = 4;
  machform = 'ieee-le'
  hdonly   = 0;
end
if nargin < 4
  machform = 'ieee-le'
  hdonly   = 0;
end
if nargin < 5
  hdonly = 0;
end

if      nbytes == 4
   prec = 'float32';
elseif nbytes == 8
   prec = 'float64';
else
   error('nbytes must be 4 or 8');
end

fp = fopen(sfile,'r');
if  fp  == -1
  error(sprintf('File %s cannot be opened for reading',sfile));
end
np   = fread(fp,1   ,prec,0,machform);
time = fread(fp,1   ,prec,0,machform);
if hdonly > 0
 data = 0;
 fclose(fp);
 return;
end

%data = zeros(np,irank);
data = fread(fp,irank*np,prec,0,machform);
fclose(fp);
data = reshape(data,irank,np);
data = data';

