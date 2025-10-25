function [kw favg] = plotgflux(sprefix,ncol,icol,trange,nwidth,do_cumsum)
%  function [kw favg] = plotgflux(sprefix,ncol,icol,trange)
%
% Reads ASCII data files produced by GHOST, computes fluxes, 
% and plots them. Either a single time index or a range may be
% specified. If the latter, averaging is done. This function
% may also be used to plot spectra, for which do_cumsum = 0.
%
%           Usage:
% 
%           [k spect] = plotgflux('ktransfer',2,2,[1:10]);
%
%           Input:
%  spref  : prefix of file set (e.g., 'ktransfer', 'hktransfer' from which to read transfers; required. See however 'do_cumsum' description.
%  ncol   : no. columns in file; required 
%  icol   : use this column of the file; reqired. 
%  trange : array of indices over which to average. May be a single index; required.
%  nwidth : time index sector width. Default is 4.
%  do_cumsum
%         : do cumulative sum. Default is 1. Setting to 0 
%           means we can apply averaging framework to any
%           spectral files. Default behavior means that files 
%           specified must be transfer functions for result to
%           be meaningful.
%
%           Output:
%  kw     : contains the wave numbers (not required). 
%  favg   : average of spectra in range, trange
%           yields the time  at which the spectum is taken
%
favg  = [];
tindex = [];

if nargin < 4
  error('Not enough parameters. Do a "help plotgflux".')
end

if nargin < 5
  nwidth = 4;
  do_cumsum = 1;
end
if nargin < 6
   do_cumsum = 1;
end

if do_cumsum > 0
  do_cumsum = 1;
end

if length(trange) == 0 
  error('trange not specified properly');
end

lwidth   = 2;
szfont   = 18;

tiny = 1e-8;
symb     = {'k-' 'g-' 'r-' 'b-' 'm-' 'c-' 'y-','k--','g--','r--','b--','m--','c--','y--','k:','g:','r:','b:','m:','c:','y:'};
symbd    = {'k--' 'g--' 'r--' 'b--' 'm--' 'c--' 'y-'};
nsymb    = numel(symb);

scurrdir = pwd;

sform = ''; % file format
for n=1:ncol
  sform  = [sform '%f'];
end

tindex = [];
n = 0;
for i = 1:numel(trange)
  tindex(n+1) = trange(i);
  sindex = pad(num2str(tindex(n+1)), nwidth, 'left', '0');
  fn = [sprefix '.' sindex '.txt'];
  fid    = fopen(fn);
  xdat   = textscan(fid, sform, 'delimiter','\n');
  fclose(fid);
  ksp    = xdat{1};    % wavenumbers
  csp    = xdat{icol}; % current spectral coeffs
csp
  csplen = numel(csp); % current length
  if n == 0
    flen  = csplen; % average spectrum's length
    favg  = csp;    % average spectrum
    kw    = ksp;
  end
  if n > 0 & csplen ~= flen % make same lengths
    nmax = max(csplen,flen);
    if csplen >= flen
      favg(flen+1:nmax) = 0; % expand
      kw   = ksp;       % store new wavenumbers
      flen = csplen;    % reset avg spect len
    else
      csp(csplen+1:nmax) = 0; % expand new spect
    end
    favg = favg + csp;
  end
  n      = n + 1;
end % file reads

favg = favg / n;

if do_cumsum == 1 % if actually doing fluxes:
  favg = -1*cumsum(favg,1);
favg
end


figure;
if do_cumsum == 0
  h = loglog(kw,favg,'k-', 'LineWidth',lwidth);
  ylabel('Spectral density','FontSize',szfont,'FontWeight','Bold');
else
  h = semilogx(kw,favg,'k-', 'LineWidth',lwidth);
  ylabel('Flux','FontSize',szfont,'FontWeight','Bold');
end

hold on;
xlabel('k','FontSize',szfont,'FontWeight','Bold');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',lwidth,'FontSize',szfont,'FontWeight','Bold');

%set(hleg,'FontSize',szfont,'FontWeight','Bold');
%sfile = sprintf('%s/figs/%s_%d.pdf',sdir,sout,N);

