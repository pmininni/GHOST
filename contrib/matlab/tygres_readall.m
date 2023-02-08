function [t G S] = tygres_readall(sdir, sfile, N, itype, trange, npl )
%
% reads the all quantities from the raw TYGRES output file,
% sfile.
%
% Usage:
%  [t G S] = tygres_readall(',', 'spectrum.txt', 3072);
%
% Arguments:
%
% sdir  : directory name (in), required
% sfile : name of file (e.g., 'spectrum.txt'), required
% N     : Actual resolution (that is, the 'equivalent resolution), required
% itype : type of file: 1==HD; 2==MHD, optional. Default is 2.
% trange: size 2 array giving start and stop time. If not specified, all 
%         data will be read. If either array value == -1, then the min/max 
%         value found in the file will be used.
% npl   : no. places to compare times in 'sfile' when determining uniqueness. Default is 15.
% 
% Returns:
%  t()    : vector of times
%  G()    : matrix array (one array for each time, t(i)) of global quantities:
%           G(i,j,ko) = global cell array at time t(i), for globals, where j = 
%             1-omega   : global kinetic energy, enstrophy, palinstrphy, next order of ko
%             2-domega  : 
%             3-omeb    : global magnetic energy, enstrophy. palinstophy, next order of ko
%             4-domeb   : 
%             5-omub    : <u.b> cross hel, enstrophy, palinstropy, next order of ko
%             6-domub   :
%            So, G(1:ntimes, 1, 1) is the kinetic energy for all times, 
%                G(1:ntimes, 3, 1) is the magnetic energy, 
%                G(1:ntimes, 1, 2) is the kinetic enstophy, 
%                G(1:ntimes, 3, 2) is the magnetic enstophy, <j^2>, etc....
%   S()   : matrix array (one row for each t(i) containing spectra. Each spectrum is
%           of size nspec = N/2 + 1: S(t(i),j,1:nspec), where j = 
%             1-spTot   : kinetic spectra cell array (one for each t(i))
%             2-dspTot  : 
%             3-sbTot   : magnetic spectra cell array (one for each t(i))
%             4-dsbTot  : 
%             5-spubTot : spectrum of u.b
%             6-dspubTot: 
%            So, S(10, 1, 1:nspec) is the kinetic energy spectrum at time index 10, 
%                S(10, 3, 1:nspec) is the magnetic energy spectrum at time index 10, etc.

if nargin < 3
  error('tygres_readall: Not enough input arguments. Do a "help tygres_readall"');
end
if nargin < 4
  itype  = 2;
  trange = [-1 -1];
  npl    = 15;
end
if nargin < 5
  trange = [-1 -1];
  npl    = 15;
end
if nargin < 6
  npl    = 15;
end
if nargin < 7
end
if nargout < 3 
  error('tygres_readall: Not enough output variables. Do a "help tygres_readall"');
end
if ( itype ~= 1 && itype ~= 2 && itype ~= 3)
  error('tygres_readall: itype must be 1 (MHD) or 2 (HD)');
end


if ( itype == 3 ) 
  ng = 7;
  ns = 6;
end
if ( itype == 2 ) 
  ng = 6;
  ns = 6;
end
if ( itype == 1 ) 
  ng = 2;
  ns = 2;
end

  sfmt = '%f';
  s4fmt = '%f%f%f%f';

fid = fopen([sdir '/' sfile]);
if fid == -1
  error (['file cannot be opened for reading: ' sdir '/' sfile]);
end

nspec   = N/2 + 1;

ttot = tygres_times(sdir, sfile, N, itype);
[tb mt nt] = iunique(ttot, npl, 'last');
G = zeros(numel(tb),ng,4);
S = zeros(numel(tb),ns,nspec);
tmp =- zeros(4);

if trange(1) < 0 
  istart = 1;
else
  istart = find(tb >= trange(1))
  istart = istart(1);
  if isempty(istart) 
    error('tygres_readall: start range value invalid');
  end
end

if trange(2) < 0 
  iend = mt(numel(tb));
else
  iend = find(tb <= trange(2));
  iend = mt(numel(iend));
  if isempty(iend) 
    error('tygres_readall: end range value invalid');
  end
end

t = [];

i = 1; % index into ttot

n = 0; % current number of values in trange
while feof(fid) == 0  
  inib = find(mt == i); % check if index i is among unique indices into ttot
  bbad = 0;
  if ( i > iend ) 
    break;
  end
  if ~isempty(inib) && i >= istart && i <= iend
    %  read and store data...
    CT      = textscan(fid, '%s', 1);
    for kk = 1:ng
      sC       = textscan(fid, '%s', 4, 'Bufsize',1000000);
      C        = str2double(sC{1});
      nana     = isnan(C);
      infa     = isinf(C);
      inana    = find(nana~=0);
      iinfa    = find(infa~=0);
      bbad     =  bbad && (~isempty(iinfa) || ~isempty(iinfa));
if bbad
  sprintf('readall: global: bbad=%d',bbad)
  inana
  iinfa
end
      if ~bbad
        tmp(1:4) = C(:);
        G(n+1,kk,1:4) = tmp(1:4);
      end
    end

    for kk = 1:ns
      sC       = textscan(fid, '%s', nspec, 'Bufsize',1000000);
      C        = str2double(sC{1});
      nana     = isnan(C);
      infa     = isinf(C);
      inana    = find(nana~=0);
      iinfa    = find(infa~=0);
      bbad     =  bbad && (~isempty(iinfa) || ~isempty(iinfa));
if bbad
  sprintf('readall: spectrum: bbad=%d',bbad)
  inana
  iinfa
end
      if ~bbad
        S(n+1,kk,1:nspec) = C;
      end
    end
    if ~bbad
      t(n+1)    = str2double(CT{1});
      n = n + 1;
    else
      warning('tygres_readall: time[%d]=%d is rejected',i,CT{1})
    end
  else
    %  just fast-forward...
    C       = textscan(fid, '%s', 1);
    for kk = 1:ng
      C        = textscan(fid, '%s', 4,'Bufsize',1000000);
    end
    for kk = 1:ns
      C        = textscan(fid, '%s', nspec, 'Bufsize',1000000);
    end
    
  end
  i       = i + 1;
end

fclose(fid);


