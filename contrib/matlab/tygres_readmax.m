function [tm M] = tygres_readmax(sdir, sfile, itype, trange, npl)
%
% reads the all quantities from the raw TYGRES diagnostic extremum output file,
% sfile. Only unique time values are returned. The quantities in the file are
% returned; these are usually the square of the value we're interested in.
%
% Usage:
%  [t M] = tygres_readmax(',', 'maximum.txt', 1);
%
% Arguments:
%
% sdir  : directory name (in), required
% sfile : name of file (e.g., 'spectrum.txt'), required
% itype : type of file: 1==HD; 2==MHD, optional; default is 2.
% trange: size 2 array giving start and stop time. If not specified, all 
%         data will be read. If either array value == -1, then the min/max 
%         value found in the file will be used.
% npl   : no. places to compare times for uniquess. Default is 15.
% 
% Returns:
%  tm()   : vector of times in maximum.txt, not repeated; final time in file
%           is the one provided
%  M()    : Maxima:
%            if itype == 1: 
%                 index 1: v2max
%                 index 2: om2max
%            if itype == 2: 
%                 index 1: v2max
%                 index 2: om2max
%                 index 3: b2max
%                 index 4: j2max
if nargin < 2
  error('tygres_read: Not enough input arguments. Do a "help tygres_read"');
end
if nargin < 3
  itype  = 2;
  trange = [-1 -1];
  npl    = 15;
end
if nargin < 4
  trange = [-1 -1];
  npl    = 15;
end
if nargin < 5
  npl    = 15;
end
if nargin < 6
end
if nargout < 2 
  error('tygres_read: Not enough output variables. Do a "help tygres_read"');
end

if ( itype == 3  ) 
  itype = 2;
end
if ( itype ~= 1 && itype ~= 2 )
  error('tygres_read: itype must be 1 (MHD) or 2 (HD)');
end



filen = [sdir '/' sfile];
fid = fopen(filen);
if ( itype == 1 )
   C    = textscan(fid,'%f%f%f','commentstyle','shell');
   t    = C{1};
   v2m  = C{2};
   om2m = C{3};
end
if ( itype == 2 )
   C    = textscan(fid,'%f%f%f%f%f','commentstyle','shell');
   t    = C{1};
   v2m  = C{2};
   om2m = C{3};
   b2m  = C{4};
   j2m  = C{5};
end



[tm mt nt] = iunique(t, npl, 'last');

if trange(1) < 0 
  istart = 1;
else
  istart = find(tm >= trange(1));
  if isempty(istart) 
    error('tygres_read: start range value invalid');
  end
  istart = istart(1);
end

if trange(2) < 0 
  iend = numel(mt);
else
  iend = find(tm <= trange(2));
  if isempty(iend) 
    error('tygres_read: end range value invalid');
  end
  iend = iend(numel(iend));
end

tm = t(mt(istart:iend));
NN = iend-istart+1;

if ( itype == 1 )
  %  [t v2m om2m ] = textread('maximum.txt','%f %f %f','commentstyle','shell');
  M = zeros(2,NN);
  M(1,1:NN) = v2m (mt(istart:iend));
  M(2,1:NN) = om2m(mt(istart:iend));
end

if ( itype == 2 )
%  [t v2m om2m b2m j2m] = textread('maximum.txt','%f %f %f %f %f','commentstyle','shell');
  M = zeros(4,NN);
  v2m (mt(istart:iend));
  M(1,1:NN) = v2m (mt(istart:iend));
  M(2,1:NN) = om2m(mt(istart:iend));
  M(3,1:NN) = b2m (mt(istart:iend));
  M(4,1:NN) = j2m (mt(istart:iend));
end

