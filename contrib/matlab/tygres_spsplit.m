function S = tygres_spsplit(sindir, soutdir, spref, N, sfile, sconf, itype, trange, npl, fsize)
%
% Reads TYGRES 'spectrum.txt' concatentated spectra file, and splits the 
% spectra into separate files, with prefix spref, indexed by place number.
% The run config file, 'parameter.txt' must reside in the directory 'sindir'.
%
% Usage:
%  = tygres_spsplit(',', 'spectrum.txt', 3072);
%
% Arguments:
%
% sindir  : input directory name (in), required
% soutdir : output directory name (in), required
% spref : prefix of output, required
% N     : Actual resolution (that is, the 'equivalent resolution), required
% sfile : name of file (e.g., 'spectrum.txt'). Default is 'spectrum.txt'.
% sconf : name of config file. Default is 'parameter.txt'
% itype : type of file: 1==HD; 2==MHD, optional. Default is 2.
% trange: size 2 array giving start and stop time. If not specified, all 
%         data will be read. If either array value == -1, then the min/max 
%         value found in the file will be used.
% npl   : no. places to compare times in 'sfile' when determining uniqueness. Default is 15.
% fsize : 32 or 64. Default is assumed to be 32.
% 
% Returns:
%   S()   : matrix array (one row for each t(i) containing spectra. Each spectrum is
%           of size nspec = N/2 + 1: S(t(i),j,1:nspec), where j = 
%             1-spTot   : kinetic spectra cell array (one for each t(i))
%             2-dspTot  : enstrophy spectra
%             3-sbTot   : magnetic spectra cell array (one for each t(i))
%             4-dsbTot  : magnetic enstrophy (j^2) spectra
%             5-spubTot : spectrum of u.b
%             6-dspubTot: 
%            So, S(10, 1, 1:nspec) is the kinetic energy spectrum at time index 10, 
%                S(10, 3, 1:nspec) is the magnetic energy spectrum at time index 10, etc.

if nargin < 4
  error('tygres_spsplit: Not enough input arguments. Do a "help tygres_spsplit"');
end
if nargin < 5
  sfile  = 'spectrum.txt';
  sconf  = 'parameter.txt';
  itype  = 2;
  trange = [-1 -1];
  npl    = 15;
  fsize  = 32;
end
if nargin < 6
  sconf  = 'parameter.txt';
  itype  = 2;
  trange = [-1 -1];
  npl    = 15;
  fsize  = 32;
end
if nargin < 7
  itype  = 2;
  trange = [-1 -1];
  npl    = 15;
  fsize  = 32;
end
if nargin < 8
  trange = [-1 -1];
  npl    = 15;
  fsize  = 32;
end
if nargin < 9
  npl    = 15;
  fsize  = 32;
end
if nargin < 10
  fsize  = 32;
end
if ( itype ~= 1 && itype ~= 2 && itype ~= 3)
  error('tygres_read: itype must be 1 (MHD) or 2 (HD)');
end


if ( itype == 3 ) 
  ng = 7;
  ns = 6;
  sspect = {'en','enst','ben','benst','ub','dub'};
end
if ( itype == 2 ) 
  ng = 6;
  ns = 6;
  sspect = {'en','enst','ben','benst','ub','dub'};
end
if ( itype == 1 ) 
  ng = 2;
  ns = 2;
  sspect = {'en','enst'};
end

sparam = sprintf('%s/parameter.txt',sindir);
[dt nouts]       = getgparam(sparam,'dt','nouts');

[trng Grng Srng] = tygres_readall(sindir, sfile, N, itype, trange, npl, fsize);
disp(sprintf('tygres_spsplit: %s file read...',[sindir '/' sfile]));
for i = 1:numel(trng)
  for kk = 1:ns
    sp = reshape(Srng(i,kk,1:nspec));
    n = round(trng(i)/(dt*nouts));
    outfile=sprintf('%s/%s_%s.%03d.txt',soutdir,spref,sspect{kk},n);
%   sp = squeeze(sp);
    dlmwrite(outfile, sp);
    sout = sprintf('tygres_spsplit: file %s written',outfile);
    disp(sout);
  end
end

fclose(fid);


