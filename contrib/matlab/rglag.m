function [data np time] = rglag(sfile, trange, stype, irank, nbytes, machform,hdonly)
%
% Reads GHOST particle data from ASCII or binary formatted files,
% fiven a time index range, trange.
%
%  Usage:
%      [data np t] = rglag('xlg',[10,20],'a');
%
%  Inputs:
%
%  sfile    : input file prefix (required)
%  trange   : array [tstart tstop] of integer start and stop time indices.
%             If == [] (empty), then no name mangling is done. 
%  stype    : 'b','B' for binary file; 'a','A' for ASCII file
%  irank    : rank of each record. Default is 3.
%  nbytes   : byte size (4, 8) of real data
%  machform : machine format (endian-ness) of real data: 'ieee-be' for
%             big-endian, or 'ieee-le' for little-endian.
%  hdonly   : if > 0, read in header: np and time only, with data=null

%  Outputs:
%
%  data    : particle data in (nt,np,3) array where nt is the number of
%            time stamps recorded, np the number of particles, and 3
%            for x, y, z positions.
%  np      : no. particles required, as found in first file
%  time    : time stamp of file
%
if nargin < 1
  error('Input file name prefix at least! Do a "help rglag".');
end
if nargin < 2
  trange = [-1 1];
  stype  = 'a';
  irank  = 3;
  nbytes = 4;
  machform = 'ieee-le';
  hdonly = 0;
end
if nargin < 3
  stype  = 'a';
  irank  = 3;
  nbytes = 4;
  machform = 'ieee-le';
  hdonly = 0;
end
if nargin < 4
  irank  = 3;
  nbytes = 4;
  machform = 'ieee-le';
  hdonly = 0;
end
if nargin < 5
  nbytes = 4;
  machform = 'ieee-le';
  hdonly = 0;
end
if nargin < 6
  machform = 'ieee-le';
  hdonly = 0;
end
if nargin < 7
  hdonly = 0;
end

if length(trange) > 1 
  if lower(stype) == 'a'
    efiles = fast_filerange(sfile,'txt',trange)
  else
    efiles = fast_filerange(sfile,'lag',trange)
  end
else
  efiles{1} = sfile;
end
nt     = length(efiles);


if nt <= 0 
  error(sprintf('Files not found for %s', sfile));
end
if lower(stype) ~= 'a' & lower(stype) ~= 'b' 
  error('stype must be either "a","A","b", or "B"');
end

time = zeros(nt,1);
if lower(stype) == 'a'
  [tmp np time(1)] = rglag_ascii(efiles{1},irank,hdonly);
  if hdonly <=0 
    data = zeros(nt,np,irank);
    data(1,1:np,1:irank) = tmp(1:np,1:irank);
  else
    data = 0;
  end
  for itime = 2:length(efiles)
    [tmp nn time(itime)] = rglag_ascii(efiles{itime},irank,hdonly);
    if nn ~= np
      error(sprintf('File %s has bad data count: required: %d; found: %d',efiles{itime},np,nn));
    end
    if hdonly <=0 
      data(itime,1:np,1:irank) = tmp;
    end
  end % time loop

  return;
end 

  [tmp np time(1)] = rglag_binary(efiles{1},irank,nbytes,machform,hdonly);
  if hdonly <=0 
    data = zeros(nt,np,irank);
    data(1,1:np,1:irank) = tmp(1:np,1:irank);
  else
    data = 0;
  end
  for itime = 2:length(efiles)
    [tmp nn time(itime)] = rglag_binary(efiles{itime},irank,nbytes,machform,hdonly);
    if hdonly <=0 
      data(itime,1:np,1:irank) = tmp;
    end
    if nn ~= np
      error(sprintf('File %s has bad data count: required: %d; found: %d',efiles{itime},np,nn));
    end
  end % time loop

