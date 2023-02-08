function [st] = tygres_ideal(sdir, sfile, runid, N, trange, kfitrng, itype, npl, btitle, binfonly)

%
% Reads ASCII data files produced by TYGRES, and computes a variety of
% relevant quantities
%
% Call sequence:
%
%    tygres_ideal('.', sfile, 'PVEC',N, trange)
%
%  sdir   : directory where 'sfile' resides
%  sfile  : 'spectrum.txt' file produced by TYGRES
%  runid  : string desc. of run for display on output plots
%  N      : spatial size of run (equivalent resolution)
%  trange : (optional) array of form [tstart tstop] times ; if this array is [-1 -1],
%           all times available will be plotted. [-1 -1] is the default.
%  kfitrng: wavenumber fitting range. Default is [0, N/2+1].
%  itype  : type of file: 1==HD; 2==MHD, optional. Default is 2.
%  npl    : number of places to compare times for uniqueness. Default is 15.
%  btitle : set titles in plots (1 to set, 0 not: e.g., for publication plots). Default is 1;
%  binfonly: if == 1, only info page is to be printed to screen; if == 2, then we compute return 
%            only, and return; if == 0, then we print everything to screen. Default is 0.
%
if nargin<4
  error('Do "help tygres_ideal" to set arguments');
end
ksta = 0;
if nargin<5
  trange   = [-1 -1];
  kfitrng  = [ksta N/2+1];
  itype    = 2;
  npl      = 15;
  btitle = 1;
  binfonly = 0;
end
if nargin<6
  kfitrng  = [ksta N/2+1];
  itype    = 2;
  npl      = 15;
  btitle = 1;
  binfonly = 0;
end
if nargin<7
  itype    = 2;
  npl      = 15;
  btitle = 1;
  binfonly = 0;
end
if nargin<8
  npl      = 15;
  btitle = 1;
  binfonly = 0;
end
if nargin<9
  btitle = 1;
  binfonly = 0;
end
if nargin<10
  btitle = 1;
  binfonly = 0;
end
if nargin<11
  binfonly = 0;
end
btitle = 0;
TINY   = 1.0e-30;
if ( itype ~= 1 && itype ~= 2 && itype ~= 3 )
  error('tygres_ideal: itype must be 1 (HD) or 2 (MHD) or 3 (modified MHD)');
end

% read data required for all times:
[tall Gall Sall] = tygres_readall(sdir, sfile, N, itype, [-1 -1], npl);

% read data over rng trange:
kmax = N/3;
x0   = [1; 1; 1];
opt   = optimset('MaxFunEvals',1000,'Maxiter',1000,'TolX',1e-8,'TolCon',1e-6,'TolFun',1.0e-15); %,'Algorithm','levenberg-marquardt');

[C n delt tfit]  = file_expdec(sdir, sfile, N, itype, trange, ksta, kfitrng, 1, x0, 2, npl, opt);
if numel(n) < numel(tfit)
  numel(n)
  numel(tfit)
  sfile
  N
  itype
  trange
  ksta
  kfitrng
  error('main: fit did not include all time values');
end

% read 'maximum.txt' file:
[tm S2m] = tygres_readmax(sdir, 'maximum.txt', itype, trange, npl);


[trng Grng Srng] = tygres_readall(sdir, sfile, N, itype, trange, npl);
nt = numel(tall);
ntrng = numel(trng);
nspec = N/2 + 1;
k = [ksta:nspec+ksta-1];


symb  = {'b-' 'g-' 'r-' 'k-' 'm-' 'c-'};
symbd = {'b--' 'g--' 'r--' 'k--' 'm--' 'c--'};
nsymb = numel(symb);


if binfonly == 2 
  return;
end

% check for existence of 'mystuff.ps':
a = dir('mystuff.ps');
if exist('a','var')
  eval('!rm mystuff.ps');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT STUFF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nf = 1;

%------------------------------------------------------------------
%figure(nf);
%nc = round(max(numel(trng)/nsymb,1));
%j  = 1;
%for i=1:numel(trng)
%% j = mod(i,nsymb)+1;
%  sp = reshape(Srng(i,1,1:nspec),nspec,1);
%  loglog(k, sp, symb{j});
%  if mod(i,nc) == 0
%    j = min(nsymb,j + 1);
%  end
%  hold on;
%end;
%title([runid ': Kinetic Energy spectra']);
%xlabel('k');
%ylabel('E_V(k)');
%print -dpsc -append 'mystuff.ps'
%nf = nf + 1;

%------------------------------------------------------------------
%
if itype == 2 || itype == 3 
%  figure(nf);
%  nc = round(max(numel(trng)/nsymb,1));
%  j  = 1;
%  for i=1:numel(trng)
%    sp = reshape(Srng(i,3,1:nspec),nspec,1);
%    loglog(k, sp, symb{j});
%    if mod(i,nc) == 0 
%      j = min(nsymb,j + 1);
%    end
%    hold on;
%  end;
%  title([runid ': Magnetic Energy spectra']);
%  xlabel('k');
%  ylabel('E_M(k)');
%  print -dpsc -append 'mystuff.ps'
%  nf = nf + 1;

  figure(nf);
% nc = round(max(numel(trng)/nsymb,1));
  display(trng);
  display(nsymb);
  j  = 1;
  ikeep = 1;
  iskip = 16;
  idel = int32(numel(trng)/iskip/nsymb);
  for i=1:iskip:numel(trng)
    sp0 = reshape(Srng(i,1,1:nspec)+Srng(i,3,1:nspec),nspec,1);
    sp = sp0;
    sp(2:nspec-1) = 0.25 * ( sp0(1:nspec-2) + sp0(3:nspec) + 2.0*sp0(2:nspec-1) );
    y  = C(i) * k.^(-n(i)) .* exp(-2*delt(i).*k);
    y  = max(y,TINY);
    if i == 1
      st{1} = sprintf('> %s starts at t=%f',symb{j},trng(i));
    end
    if (i-ikeep)==idel
      j = min(nsymb,j+1);
      st{j} = sprintf('> %s ends and %s begins at t=%f',symb{j-1},symb{j},trng(i));
      display(st{j})
      ikeep = i;
    end
    if i == numel(trng)
      st{numel(st)+1} = sprintf('> %s ends at t=%f',symb{j},trng(i));
    end
    loglog(k, sp, symb{j}, k, y, symbd{j});
%   loglog(k, sp, symb{j});
    set(gca,'YLim',[1.0e-4 10]);
    hold on;
  end;
  set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5, 'FontSize',11,'FontWeight','Bold','Box','on');
% ht = text(1,1e-20,st);
  if btitle == 1
  title([runid ': Total Energy spectra']);
  end
  xlabel('k','FontSize',11,'FontWeight','Bold');
  ylabel('E_V(k)+E_M(k)','FontSize',11,'FontWeight','Bold');
  print -dpsc -append 'mystuff.ps'
  nf = nf + 1;

  figure(nf);
  nc = round(max(numel(trng)/nsymb,1));
  j  = 1;
  ikeep = 1;
  iskip = 4;
  idel = int32(numel(trng)/iskip/nsymb)
  st{1} = sprintf('> %s starts at t=%f',symb{1},trng(1));
  disp(st{1})
  for i=1:iskip:numel(trng)
    sp0 = reshape(Srng(i,3,1:nspec)./(Srng(i,1,1:nspec)+1e-30),nspec,1);
    sp  = sp0;
    sp(2:nspec-1) = 0.25 * ( sp0(1:nspec-2) + sp0(3:nspec) + 2.0*sp0(2:nspec-1) );
    if (i-ikeep)==idel
      j = min(nsymb,j+1);
      st{j} = sprintf('> %s ends and %s begins at t=%f',symb{j-1},symb{j},trng(i));
      disp(st{j})
      ikeep = i;
    end
    loglog(k, sp, symb{j});
    set(gca,'YLim',[1.0e-4 10]);
    hold on;
%   disp(sprintf('t=%f; Em(k=1)=%f; Em(k=2)=%f; Em(k=3)=%f\n',trng(i),Srng(i,3,1),Srng(i,3,2),Srng(i,3,3)));
  end;
  st{numel(st)+1} = sprintf('> %s ends at t=%f',symb{j},trng(i-1));
  disp(st{numel(st)})
  set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5, 'FontSize',11,'FontWeight','Bold','Box','on');
  if btitle == 1 
  title([runid ': E_M(k) / E_V(k)'])
  end

  xlabel('k','FontSize',11,'FontWeight','Bold');
  ylabel('E_M(k)/E_V(k)','FontSize',11,'FontWeight','Bold');
  print -dpsc -append 'mystuff.ps'
  nf = nf + 1;
end

%------------------------------------------------------------------

figure(nf);
h(1)=subplot(2,1,1);
y   = [];
y(1:numel(tfit)) = 2/kmax;
y1 = 0.5 .* y;
y2 = 0.5 .* y1;
%semilogy(tfit, abs(delt), symb{1}, tfit, y, symbd{1}, tfit, y1, symbd{2}, tfit, y2, symbd{3});
semilogy(tfit, abs(delt), symb{1}, tfit, y, symbd{1});
set(h(1),'YMinorTick','on','XMinorTick','on','LineWidth',1.5, 'FontSize',11,'FontWeight','Bold','Box','on');
%set(h(1),'YLim',[1e-3 7]);
if btitle == 1
title([runid ': Logarithmic decrement'],'FontWeight','Bold');
end
set(h(1),'XLim',[trange(1) trange(2)]);
xlabel('t','FontSize',11,'FontWeight','Bold');
ylabel('\delta(t)','FontSize',11,'FontWeight','Bold');
hold on;

h(2)=subplot(2,1,2);
plot(tfit, n , 'r-');
set(h(2),'YLim',[0 7]);
set(h(2),'XLim',[trange(1) trange(2)]);
set(h(1),'YMinorTick','on','XMinorTick','on','LineWidth',1.5, 'FontSize',11,'FontWeight','Bold','Box','on');
if btitle == 1
title([runid ': Spectral index from fit'],'FontWeight','Bold');
end
xlabel('t','FontSize',11,'FontWeight','Bold');
ylabel('n(t)','FontSize',11,'FontWeight','Bold');

set(h(2),'XMinorTick','on','YMinorTick','on','LineWidth',1.5, 'FontSize',11,'FontWeight','Bold','Box','on');

print -dpsc -append 'mystuff.ps'
hold off;
nf = nf + 1;

%------------------------------------------------------------------


if itype == 1
  figure(nf);
  plot(trng, Grng(1:ntrng,1,1), 'k-');
  set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5, 'FontSize',11,'FontWeight','Bold','Box','on');
  if btitle == 1
  title([runid ': E_{V}(t)'],'FontSize',11,'FontWeight','Bold');
  end
  xlabel('t','FontSize',11,'FontWeight','Bold');
  ylabel('E_V','FontSize',11,'FontWeight','Bold');
  print -dpsc -append 'mystuff.ps'
  nf = nf + 1;

  figure(nf);
  plot(trng, Grng(1:ntrng,1,2), 'k-');
  set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5, 'FontSize',11,'FontWeight','Bold','Box','on');
  if btitle == 1
  title([runid ': <\omega^2>_{V}(t)'],'FontSize',11,'FontWeight','Bold');
  end
  xlabel('t','FontSize',11,'FontWeight','Bold');
  ylabel('<\omega^2>','FontSize',11,'FontWeight','Bold');
  print -dpsc -append 'mystuff.ps'
  nf = nf + 1;
end

if itype == 2 || itype == 3 
  figure(nf);
  semilogy(trng, Grng(1:ntrng,1,1), 'k-', trng, Grng(1:ntrng,3,1), 'r-', trng,  Grng(1:ntrng,1,1)+Grng(1:ntrng,3,1), 'b-');
  set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5, 'FontSize',11,'FontWeight','Bold','Box','on');
  if btitle == 1
  title([runid ': E_{V}(t) and E_{M}(t)'],'FontSize',11,'FontWeight','Bold');
  end
  xlabel('t','FontSize',11,'FontWeight','Bold');
  ylabel('E_V/E_M','FontSize',11,'FontWeight','Bold');
  legend('E_{V}', 'E_{M}', 'E_{tot}','Location','NorthWest');
  print -dpsc -append 'mystuff.ps'
  nf = nf + 1;

  figure(nf);
  semilogy(trng, Grng(1:ntrng,1,2), 'k-', trng, Grng(1:ntrng,3,2), 'r-', trng, Grng(1:ntrng,1,2)+Grng(1:ntrng,3,2), 'b-');
  set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5, 'FontSize',11,'FontWeight','Bold','Box','on');
  if btitle == 1
  title([runid ': <\omega^2>(t), <j^2>(t)'],'FontSize',11,'FontWeight','Bold');
  end
  xlabel('t','FontSize',11,'FontWeight','Bold');
  ylabel('<\omega^2>/<j^2>','FontSize',11,'FontWeight','Bold');
  legend('<\omega^2>', '<j^2>', '<\omega^2>+<j^2>','Location','NorthWest');
  print -dpsc -append 'mystuff.ps'
  nf = nf + 1;


  figure(nf);
  plot(trng, Grng(1:ntrng,3,1)./Grng(1:ntrng,1,1), 'k-');
  set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5, 'FontSize',11,'FontWeight','Bold','Box','on');
  if btitle == 1
  title([runid ': E_M / E_V'],'FontSize',11,'FontWeight','Bold');
  end
  xlabel('t','FontSize',11,'FontWeight','Bold');
  ylabel('E_M / E_V','FontSize',11,'FontWeight','Bold');
  print -dpsc -append 'mystuff.ps'
  nf = nf + 1;
end


if itype == 1 
  figure(nf);
  semilogy(tm, sqrt(S2m(2,:)), 'k-');
  set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5, 'FontSize',11,'FontWeight','Bold','Box','on');
  if btitle == 1
  title([runid ': \Omega_{max}'],'FontSize',11,'FontWeight','Bold');
  end
  xlabel('t','FontSize',11,'FontWeight','Bold');
  ylabel('\Omega_{max}','FontSize',11,'FontWeight','Bold');
  print -dpsc -append 'mystuff.ps'
  nf = nf + 1;
end 
if itype == 2 || itype == 3 
  figure(nf);
  semilogy(tm, sqrt(S2m(2,:)), 'k-', tm, sqrt(S2m(4,:)),'r-');
  set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5, 'FontSize',11,'FontWeight','Bold','Box','on');
  if btitle == 1
  title([runid ': \Omega_{max}, J_{max}'],'FontSize',11,'FontWeight','Bold');
  end
  xlabel('t','FontSize',11,'FontWeight','Bold');
  ylabel('\Omega_{max}, J_{max}','FontSize',11,'FontWeight','Bold');
  legend('\omega_{max}', 'j_{max}','Location','NorthWest');
  print -dpsc -append 'mystuff.ps'
  nf = nf + 1;
end
