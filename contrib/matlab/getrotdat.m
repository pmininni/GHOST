function getrotdat(sdir, tavg, srundesc, soutfile, ikeep, stype)

% function getrotdat(sdir, tavg, srundesc, soutfile, ikeep, stype)
%
% Reads ASCII data files produced by GHOST ROTH solver, and computes a variety of
% relevant quantities. All relevant files must reside in the specified directories.
% If more than one run is specified, then the quantities for the runs will be plotted 
% on the same plot for comparison.
%
% Call sequence:
%
%    getrotdat({'Run1Dir', 'Run2Dir'},[[0 10];[5 15]], {'ABCK5R14','ABCK5R15'},'compare_1_2.ps') 
%
%      		This call will cd to run directories 'Run1Dir' and 'Run2Dir' (relative to launch 
%               directory), and plot the same quantities on the same plots, placing 
%               all results in the file 'compare_1_2.ps' in the launch directory. Run1Dir will
%               average over the interval t =[0 10], and Run2Dir will average over the 
%               interval t = [5 15] for quantities that require time averaging. The plots will
%               contain the names 'ABCK5R14' and 'ABCK5R15' to distinguish the plots.
%
%    getrotdat({'MyRunDir'},[[0 10]], {'TGrun'},'mystuff.ps') 
%
%      		This call will cd to directory MyRunDir (relative to launch directory), and
%               place the plots in file 'mystuff.ps' in the launch directory. The
%               name 'TGrun' will be placed in the title in each plot.
%    
%               If file 'refparams.dat' exists, reference slope parameters 
%               will be plotted for selected plots. Search for refparams below.
%               Specification of these variables must be in the same way
%               as for the GHOST 'parameter.txt' file:
%                   var1 = 0.0000;
%                   var2 = 10; ...
%   
%
%  sdir    : cell array of strings representing run directories. If not specified, 
%            launch directory is used. All relevant files must reside in this directory.
%  tavg    : 2-element array list specifying [ start_i stop_i ] bounds of evolutionary time
%            over which to average spectra for each run specified in 'sdir'. If not specified,
%            the entire time range will be used for each run. The number of rows of this ndir-by-2
%            matrix must equal the number of runs specified in 'sdir'.
%  srundesc: cell array of length numel(sdir) with a description of each run. If not specified,
%            each description will be the NULL string.
%  soutfile: If specified, this will be the file name to which all image output is put. The 
%            default will be to put it in the launch directory under 'mystuff.ps'. This is
%            always a postscript file. 
%  ikeep   : flag tellilng whether items specified in 'stype' (below) are to be discarded (0) or the
%            only items retained in the plots (1). The stypes are given here:
%  stype   : cell array of strings spcifying plot types. If specified, and  'ikeep' flag (above) is 0 or 
%            not set, then these quantities will not be plotted. If stype is specified and ikeep=1,
%            then only these quantities will be plotted. The stypes are given below. A check on these
%            valid types _is_ case sensitive.
%     
%         Valid stype values:
%         'EvT'      -- Tot. energy v. time
%         'E23DvT'   -- 2D / 3D energy v. time
%         'EnstvT'   -- Tot. enstrophy v. time
%         'Hvt'      -- Tot. helicity v. time
%         'rhovT'    -- cos(u omega)
%         'Lint'     -- Integral scale v. time
%         'Ltay'     -- Taylor scale v. time
%         'Tint'     -- Integral time scale
%         'Reint'    -- Integral scale Re v. time
%         'Retay'    -- Taylor scale Re v. time
%         'Roint'    -- Taylor scale Ro v. time
%         'Rbint'    -- Buoyancy Re with N/f given by file, or hardcoded
%         'Dtr'      -- Turbulent diss. rate  (~U^3/L_int) v. time
%         'Kdr'      -- Kinetic diss. rate  (nu \int omega^2) v. time
%         'kOm'      -- Zeman scale v. time
%         'spTau'    -- tau_NL spectrum
%         'spE'      -- Energy  spectrum
%         'spEpar'   -- Parallel Energy spectrum
%         'spEperp'  -- Perpendicular Energy spectrum
%         'spH'      -- Helicity spectrum
%         'spH2kE'   -- H(k) / k E(k) spectrum
%         'PiE'      -- Energy flux spectrum
%         'PiH'      -- Helicity flux spectrum  
%
shome    = pwd;
ndir     = numel(sdir);
% Cell array of valid plot types (type 'database'). The order in the array
% is the page or plot  number; there should be one sdbtype element
% for each type of plot. Some types may have more than one plot
% associated with it. The array ndbtyp contains the max. number of
% plots associated with a given type:
sdbtype = {...
'EvT'     , ... 
'E23DvT'  , ...
'EnstvT'  , ...
'HvT'     , ...
'rhovT'   , ...
'Lint'    , ...
'Ltay'    , ...
'Tint'    , ...
'Reint'   , ...
'Retay'   , ...
'Roint'   , ...
'Rbint'   , ...
'Dtr'     , ...
'Kdr'     , ...
'kOm'     , ...
'spTau'   , ...
'spE'     , ...
'spEpar'  , ...
'spEperp' , ...
'spH'     , ...
'spH2kE'  , ...
'PiE'     , ...
'PiH'     , ...
}
ntype = numel(sdbtype);

% ndbtype describes no. figure windows (nw), and 
% no. curves in each figure window (nc), s.t.
% ndbtype = { [nw_1,nc1_w1,nc1_w2...], [nw_2,nc1_w1, nc1_w2,...], ...}
% where nw_1 is the no. windows for plot item 1, and nc1_w1 is the number 
%| of curves in the first plot window for item 1, and nc1_w2 is 
% the number of curves in the second plot window for item 1, etc..
ndbtype = {...
[1,1]  ,...
[1,2]  ,...
[1,1]  ,...
[1,1]  ,...
[1,1]  ,...
[2,1,2]  ,...
[1,1]  ,...
[1,1]  ,...
[1,1]  ,...
[1,1]  ,...
[1,1]  ,...
[1,1]  ,...
[1,1]  ,...
[1,1]  ,...
[2,1,1],...
[1,1]  ,...
[1,2]  ,...
[1,1]  ,...
[2,2,2]  ,...
[1,1]  ,...
[1,1]  ,...
[1,1]  ,...
[1,1]  ,...
};

nmaxplt= 0;
for j = 1:numel(ndbtype)
  nmaxplt = nmaxplt + ndbtype{j}(1);
end

% check ndbtype structure:
for j = 1:ntype
  nw = ndbtype{j}(1);
  if numel(ndbtype{j}) ~= (nw+1) 
    error('getrotdat: incorrect ndbtype structure!');
  end
end

nf = 0;
for j = 1:ntype
  stitle {j}  = cell(ndbtype{j}(1),1);
  sxlabel{j}  = cell(ndbtype{j}(1),1);
  sylabel{j}  = cell(ndbtype{j}(1),1);
  nw = ndbtype{j}(1);                     % no. windows for type j
  nc = sum(ndbtype{j}(2:nw+1));           % total no. curves for 
  ifig{j} = cell(nw,1);                   % fig. window number
  for k = 1:nw
    nf = nf + 1;
    stitle {j}{k,1} = '';
    sxlabel{j}{k,1} = '';
    sylabel{j}{k,1} = '';
    ifig{j}{k} = nf;
    nc = ndbtype{j}(k+1);                 % total no. curves for one run on window k
    sclabel{j}  = cell(ndir,1);           % curve label for each curve
    for m = 1:ndir
      sclabel{j}{m}  = cell(nw,nc);       % curve label for each curve
      for nk = 1:nw
        for nl = 1:nc
          sclabel{j}{m}{nk,nl} = '';
        end
      end
    end
  end
  hh     {j}  = cell(nw,1);
end


% Create global indirection reference into sdbtype:
itypestart(1) = 1;
for j = 2:ntype
% itypestart(j) = itypestart(j-1) + ndbtype{j-1}(1);
  itypestart(j) = j;
% sprintf('stype=%s, num=%d, itypestart=%d \n', sdbtype{j}, ndbtype(j), itypestart(j))
end

if nargin<1
  sdir{1} = pwd;
  srundesc{1} = '';
  tavg    = [0.0 1e7];
  soutfile= 'mystuff.ps';
  ikeep   = 1;
  stype   = sdbtype;
end
if nargin<2
  srundesc{1} = '';
  tavg    = [0.0 1e7];
  soutfile= 'mystuff.ps';
  ikeep   = 1;
  stype   = sdbtype;
end
if nargin<3
  tavg    = [0.0 1e7];
  soutfile= 'mystuff.ps';
  ikeep   = 1;
  stype   = sdbtype;
end
if nargin<4
  soutfile= 'mystuff.ps';
  ikeep   = 1;
  stype   = sdbtype;
end
if nargin<5
  ikeep   = 1;
  stype   = sdbtype;
end
if nargin<6
  if ikeep == 0
    error('getrotdat: Nothing to do!');
  end
  stype   = sdbtype;
end

if ikeep > 0 
  ikeep = 1;
end
if ikeep < 0 
  ikeep = 0;
end

if ikeep == 1
    bkeep = zeros(ntype,1);
else
    bkeep = ones(ntype,1);
end

for i=1:numel(stype)
% stype{i}
  ic = strcmp(stype{i},sdbtype);
  ip = find(ic == 1);

  jstart = itypestart(ip(1));
  bkeep(jstart) = ~bkeep(jstart);
%   sprintf('stype=%s; keep_1=%d; index=%d',stype{i},bkeep(j),j)
end

ns       = size(tavg);
if ns(1) ~= ndir
  error('getrotdat: number of tavg rows must equal number of runs, ndir');
end
if numel(srundesc) ~= ndir && ndir > 1
  error('getrotdat: number of srundesc entries must equal number of runs, ndir');
end

%
color= {'k' 'b' 'r' 'g' 'm' 'c'};
symb  = {'k-' 'b-' 'r-' 'g-' 'm-' 'c-'};
symb2  = {'k-' 'k--' 'k-.' '-' '--' '-.'};
symbd = {'k--' 'b--' 'r--' 'g--' 'm--' 'c--'};
nsymb = numel(symb);
TINY = 1.0e-30;


% check for existence of output image file:
cd(shome)
a = dir(soutfile);
if exist('a','var')
  eval(['!rm ' soutfile]);
end


for n=1:ndir
  sdir{n}
  % Time-dependent comparisons:
  cd(sdir{n});
  balp    = textread('balance.txt'   ,'','commentstyle','shell');
  help    = textread('helicity.txt'  ,'','commentstyle','shell');

  [bunique IIb JJb] = iunique(balp(:,1),7);
  [hunique IIh JJh] = iunique(help(:,1),7);

  clear balpp helpp bal hel
  balpp(:,1) = balp(IIb,1);
  balpp(:,2) = balp(IIb,2);
  balpp(:,3) = balp(IIb,3);
  balpp(:,4) = balp(IIb,4);
  helpp(:,1) = help(IIh,1);
  helpp(:,2) = help(IIh,2);
  [CC Ibal Ihel] = iintersect(balpp(:,1),helpp(:,1),5);
  bal(:,1) = balpp(Ibal,1);
  bal(:,2) = balpp(Ibal,2);
  bal(:,3) = balpp(Ibal,3);
  bal(:,4) = balpp(Ibal,4);
  hel(:,1) = helpp(Ihel,1);
  hel(:,2) = helpp(Ihel,2);
  clear bunique hunique IIb JJb IIh JJh balpp balp helpp help CC Ibal Ihel 

  Initial_hel(n)=hel(1,2);
  

  % Get relevant parameters from run config file:
  [dt sstep, nu, omega bvfreq] =getgparam('parameter.txt','sstep','dt','nu','omega','bvfreq');
  if omega ~= 0.0
     N2f = bvfreq/(2*omega)
  end
  if N2f == 0.0
     N2f = 4;
  end
  eps_kin = (2.0*nu)*bal(:,3);
  komega_kin = (omega^3 ./eps_kin).^(0.5);

  % Time-averaged quantities:
  irange(1) = uint16(tavg(n,1)/(dt*sstep)+0.50001);
  irange(2) = uint16(tavg(n,2)/(dt*sstep)+0.50001);
  trange = sprintf('%d:%d',irange(1),irange(2));
  [Hkavg Hk thindex]   = spectavg('khelicity', trange);
  [Ekavg Ek teindex]   = spectavg('kspectrum', trange);
  [Ekpar_avg Ekpar tepar_index]   = spectavg('kspecpara', trange);
  [Ekperp_avg Ekperp teperp_index]   = spectavg('kspecperp', trange, 3, 1);
  [Ek2D_avg Ek2D teperp_index]   = spectavg('kspecperp', trange, 3, 2);
% [Ekw_avg Ekw teperp_index]   = spectavg('kspecperp', trange,3,2);
  Ek3D_avg = Ekperp_avg - Ek2D_avg;

% [Ekpw_avg Ekpw teperp_index]   = spectavg('kspecperp', trange, 3, 3);

  ntavg = numel(thindex);
  ip = strmatch('PiE',sdbtype);
  if bkeep(itypestart(ip(1)))
    [efavg eftk etindex] = gflux('ktransfer' ,'', dt, sstep, tavg(n,1), tavg(n,2));
    ntavg = min(ntavg,numel(teindex));
  end
  ip = strmatch('PiH',sdbtype);
  if bkeep(itypestart(ip(1)))
     [hfavg hftk htindex] = gflux('hktransfer','', dt, sstep, tavg(n,1), tavg(n,2));
    ntavg = min(ntavg,numel(htindex));
  end
 
  [dum Hk thindex]    = spectavg('khelicity', '1:end');
  [dum Ek teindex]    = spectavg('kspectrum', '1:end');
  nt   = min(numel(thindex),numel(teindex));
% nn   = find(Ekavg < TINY);
  NN   = numel(Ekavg);
% if ~isempty(nn)
%   NN   = nn(1)-2;
% end
  Ekavg    = 0.5.*Ekavg(1:NN);
  Hkavg    = Hkavg(1:NN);
  ip = strmatch('PiE',sdbtype);
  if bkeep(itypestart(ip(1)))
  efavg    = efavg(1:NN);
  end
  ip = strmatch('PiH',sdbtype);
  if bkeep(itypestart(ip(1)))
    hfavg    = hfavg(1:NN);
  end

  for i = 1:nt
    Ek    {i} = Ek    {i}(1:NN);
    Hk    {i} = Hk    {i}(1:NN);
  end
  for i = 1:ntavg
    ip = strmatch('PiE',sdbtype);
    if bkeep(itypestart(ip(1)))
      eflux {i} = eftk  {i}(1:NN);
    end
    ip = strmatch('PiH',sdbtype);
    if bkeep(itypestart(ip(1)))
      hflux {i} = hftk  {i}(1:NN);
    end
  end
  E_tot    = [];
  E2D      = [];
  E3D      = [];
  H_tot    = [];
  L_k      = {};
  urms     = [];
  Enst     = {};
  Enst_tot = [];
  L_int    = [];
  L_tay    = [];
  tau_int  = [];
  eps_dtr  = [];
  Re_int   = [];
  Re_tay   = [];
  Ro_int   = []
  Rb_int   = [];

  N    = numel(Ekavg);
  k    = [1:N] ;
  tt   = teindex .* (sstep*dt);
  tau_nla = ( sqrt( 2.0*Ekavg .* ((k.^3)') ) ) .^ (-1);

  % ...Computed from time-dependent spectra:
  m = 0;
  for i = 1:nt
    tmp = sum(Ek{i});
    if isnan(tmp) | isinf(tmp) 
      break;
    end
    m = m + 1
    Ek        {i} = 0.5*Ek{i};
    E_tot     (i) = sum(Ek{i});
    H_tot     (i) = sum(Hk{i});
    L_k       {i} = Ek{i}./k';
    urms      (i) = sqrt(2.0*E_tot(i));
    Enst      {i} = Ek{i} .* k' .* k';
    Enst_tot  (i) = sum(Enst{i},1);
    L_int     (i) = 2*pi*sum(L_k{i},1)/E_tot(i);
    tau_int   (i) = L_int(i) / urms(i);
    eps_dtr   (i) = urms(i)^3/L_int(i); % turb. diss. rate
    Re_int    (i) = urms(i) * L_int(i) / nu;
    Ro_int    (i) = 0.0;
    Rb_int    (i) = 0.0;
    if omega ~= 0.0
      Ro_int  (i) = urms(i) / (2.0*L_int(i)*omega);
      Rb_int  (i) = Re_int(i)*(Ro_int(i)/N2f)^2;
    end
    L_tay     (i) = 2*pi*sqrt(E_tot(i)/(Enst_tot(i)));
    Re_tay    (i) = urms(i) * L_tay(i) / nu;
    if isnan(E_tot(i)) | isinf(E_tot(i))
      break;
    end 
  end
  nt = m;
  t  = tt(1:nt);
  komega_dtr = (omega^3./eps_dtr).^(0.5);
  

  nf = 0;  % figure index
  %  ****************************time-dependent quantities ***************************
 
  ip = strmatch('EvT',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (1) ...Energy:
    nf = nf + 1;     
    figure(ifig{itypestart(ip(1))}{1});
    hh     {itypestart(ip(1))}{1,1} = plot(bal(:,1),bal(:,2),symb{n});
    stitle {itypestart(ip(1))}{1,1} = 'Energy vs time';
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = 'Energy';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = srundesc{n};
    end
    hold on;
%   n_lim=4;
%   axis([0 max(bal(:,1)) 0 max(bal(:,2))])   
%   if n==ndir
%     % (1) ...Energy: 
%     nf = nf + 1;
%     [h_main, h_inset]=inset(fig1,fig1,0.5,ifig{itypestart(ip(1))}{2});
%     set(h_inset,'xlim',[0 bal(n_lim,1)],'ylim', [min(bal(1:n_lim,2)) max(bal(1:n_lim,2))])
%   end
  end

  ip = strmatch('E23DvT',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (2) ...2D, 3D Energy:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    for i = 1:numel(teperp_index)
      Ek3D      {i} = Ekperp{i} - Ek2D{i};
      E2D       (i) = 0.5*sum(Ek2D{i});
      E3D       (i) = 0.5*sum(Ek3D{i});
    end
    t2d = teperp_index*dt*sstep
    hh     {itypestart(ip(1))}{1,1} = plot(t2d,E2D,symb{n},t2d,E3D,symbd{n});
    stitle {itypestart(ip(1))}{1,1} = 'Energy vs time';
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = 'Energy';
    if  ndir > 1
      sclabel{itypestart(ip(1))}{n}{1} = sprintf('%s:%s',srundesc{n},'E_{3D}');
      sclabel{itypestart(ip(1))}{n}{2} = sprintf('%s:%s',srundesc{n},'E_{2D}');
    else
      sclabel{itypestart(ip(1))}{n}{1} = sprintf('%s','E_{3D}');
      sclabel{itypestart(ip(1))}{n}{2} = sprintf('%s','E_{2D}');
    end
    hold on;
  end


  ip = strmatch('EnstvT',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (3) ...Enstrophy:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh     {itypestart(ip(1))}{1,1} = plot(bal(:,1),bal(:,3),symb{n});
    stitle {itypestart(ip(1))}{1,1} = 'Enstrophy vs time';
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = '\int \omega^2';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = srundesc{n};
    end
    hold on;
  end


  ip = strmatch('HvT',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (4) ...Helicity:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh     {itypestart(ip(1))}{1,1} = plot(hel(:,1),hel(:,2),symb{n});
    stitle {itypestart(ip(1))}{1,1} = ['Helicity vs time,' 10 ' H(t=0) = ' num2str(Initial_hel)];
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = '<\omega.u>';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = srundesc{n};
    end
    hold on
    line([0 max(hel(:,1))],[mean(hel(:,2)) mean(hel(:,2))],'Color',color{n})
    hold on;
  end
 
 
  ip = strmatch('rhovT',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (5) ...Helicity:
    nf = nf + 1;
whos bal
whos hel
rm = hel(:,2)./sqrt(bal(:,2).*bal(:,3));
whos rm
    figure(ifig{itypestart(ip(1))}{1});
    hh     {itypestart(ip(1))}{1,1} = plot(hel(:,1),hel(:,2)./sqrt(bal(:,2).*bal(:,3)),symb{n});
    stitle {itypestart(ip(1))}{1,1} = 'Rel. Hel. vs time';
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = '\rho';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = srundesc{n};
    end
    hold on;
  end

  ip = strmatch('Lint',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (6) ...Integral scale: lin-lin
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh     {itypestart(ip(1))}{1,1} = plot(t,L_int,symb{n});
    stitle {itypestart(ip(1))}{1,1} = 'L_{int} vs time';
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = 'L_{int}';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = srundesc{n};
    end
    hold on;

  end

  ip = strmatch('Lint',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (6) ...Integral scale: log-log
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{2});
    hh{itypestart(ip(1))}{2,1} = loglog(t,L_int,symb{n});
    hold on;
    stitle {itypestart(ip(1))}{2,1} = 'L_{int} vs time';
    sxlabel{itypestart(ip(1))}{2,1} = 'time';
    sylabel{itypestart(ip(1))}{2,1} = 'L_{int}';
    if exist('refparams.dat','file')
      [teimin teimax teiindex teimult] =getgparam('refparams.dat','teimin','teimax','teiindex','teimult') 
      x  = teimin:1:teimax
      y  = teimult*(x.^(teiindex))./(teimin^teiindex)
      if ndir == 1 | n == ndir
      loglog(x, y, symbd{n});
      sclabel{itypestart(ip(1))}{n}{2,1} = 'L_{int}';
      sclabel{itypestart(ip(1))}{n}{2,2} = sprintf('t**(%f)',teiindex);
sprintf('.............................................here.')
      end
    else
%     error('No refparams: E(k)');
      sclabel{itypestart(ip(1))}{n}{2,1} = 'L_{int}';
      sclabel{itypestart(ip(1))}{n}{2,2} = ''
    end
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{2,1} = srundesc{n};
    end
    hold on;
  end


  ip = strmatch('Ltay',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (7) ...Taylor scale:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh     {itypestart(ip(1))}{1,1} = plot(t,L_tay,symb{n});
    stitle {itypestart(ip(1))}{1,1}= 'L_{Tay} vs time';
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = 'L_{Tay}';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = srundesc{n};
    end
    hold on;
  end

  ip = strmatch('Tint',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (7) ...Integral scale timescale:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh     {itypestart(ip(1))}{1,1} = plot(t,tau_int,symb{n});
    stitle {itypestart(ip(1))}{1,1}= '\tau_{int} vs time';
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = '\tau_{int}';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = srundesc{n};
    end
    hold on;
  end

  ip = strmatch('Reint',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (8) ...Integral-scale Re:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh     {itypestart(ip(1))}{1,1} = plot(t,Re_int,symb{n});
    stitle {itypestart(ip(1))}{1,1} = 'Re_{int} vs time';
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = 'Re_{int}';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = srundesc{n};
    end
    hold on;
  end

  ip = strmatch('Roint',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (8) ...Integral-scale Re:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh     {itypestart(ip(1))}{1,1} = plot(t,Ro_int,symb{n});
    stitle {itypestart(ip(1))}{1,1} = 'Ro_{int} vs time';
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = 'Ro_{int}';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = srundesc{n};
    end
    hold on;
  end

  ip = strmatch('Rbint',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (8) ...Integral-scale Rb:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh     {itypestart(ip(1))}{1,1} = plot(t,Rb_int,symb{n});
    stitle {itypestart(ip(1))}{1,1} = 'Rb_{int} vs time';
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = 'Rb_{int}';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = srundesc{n};
    end
    hold on;
  end


  ip = strmatch('Retay',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (9) ...Taylor scale Re:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh     {itypestart(ip(1))}{1,1} = plot(t,Re_tay,symb{n});
    stitle {itypestart(ip(1))}{1,1} = 'Re_{Tay} vs time';
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = 'Re_{Tay}';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = srundesc{n};
    end
    hold on;
  end
  
  ip = strmatch('Dtr',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (12) ...Turb. dissipatn rate: 
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh     {itypestart(ip(1))}{1,1} = plot(t,eps_dtr,symb{n});
    stitle {itypestart(ip(1))}{1,1} = 'Turb. dissipation rate (u^3/L)  vs time';
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = 'u_{rms}^3/L_{int}';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = srundesc{n};
    end
    hold on;
  end

  ip = strmatch('Kdr',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (13) ...Energy dissipatn rate from <\omega^2>i ('kinetic dissipatn rate'): 
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh     {itypestart(ip(1))}{1,1} = plot(bal(:,1),eps_kin,symb{n});
    stitle {itypestart(ip(1))}{1,1} = 'Kinetic dissipation rate (nu \int \omega^2 )  vs time';
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = '2 nu \int \omega^2';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = srundesc{n};
    end
    hold on;
  end
  
  ip = strmatch('kOm',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (14) ...k_N from turb. diss. rate:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh     {itypestart(ip(1))}{1,1} = semilogy(t,komega_dtr,symb{n});
    stitle {itypestart(ip(1))}{1,1} = 'k_{\Omega,turb}  vs time';
    sxlabel{itypestart(ip(1))}{1,1} = 'time';
    sylabel{itypestart(ip(1))}{1,1} = 'k_{\Omega,turb}';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = srundesc{n};
    end
    hold on;
  end
  
  ip = strmatch('kOm',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (15) ...k_N from kinetic diss. rate:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{2});
    hh     {itypestart(ip(1))}{2,1} = semilogy(bal(:,1),komega_kin,symbd{n});
    stitle {itypestart(ip(1))}{2,1} = 'k_{\Omega,kinetic}  vs time';
    sxlabel{itypestart(ip(1))}{2,1} = 'time';
    sylabel{itypestart(ip(1))}{2,1} = 'k_{\Omega,kin}';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{2,1} = srundesc{n};
    end
    hold on;
  end

  %  **************************** k-dependent quantities ***************************

  ip = strmatch('spTau',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (19) ...spectra of tau_{NL} avgd over time interval: 
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    stitle {itypestart(ip(1))}{1,1} = 'Spectra of <\tau_{NL}>_t';
    sxlabel{itypestart(ip(1))}{1,1} = 'k';
    sylabel{itypestart(ip(1))}{1,1} = '\tau(k)';
    sclabel{itypestart(ip(1))}{n}{1,1} = sprintf('%s:%s',srundesc{n},'\tau_k');
    y = ones(1,numel(k));
    hh{itypestart(ip(1))} = loglog(k,tau_nla,symb{n});
    hold on;
    axis([0 max(k) 0 max(tau_nla)])
  end

  ip = strmatch('spE',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (20) ...<E(k)>_t spectrum:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh{itypestart(ip(1))} = loglog(k,Ekavg,symb{n});
    hold on;
%   x  = 6:1:50;
%   y  = 0.1*x.^(-5.0/3)
%   loglog(x, y, symbd{n});
    hold on;
    if exist('refparams.dat','file')
      [keimin keimax keiindex keimult] =getgparam('refparams.dat','keimin','keimax','keiindex','keimult') 
      x  = keimin:1:keimax;
      y  = keimult*x.^(keiindex)
      if ndir == 1 | n == ndir
      loglog(x, y, symbd{n});
      sclabel{itypestart(ip(1))}{n}{1,1} = 'E(k)';
      sclabel{itypestart(ip(1))}{n}{1,2} = sprintf('k**(%f)',keiindex);
      end
    else
%     error('No refparams: E(k)');
      sclabel{itypestart(ip(1))}{n}{1,1} = 'E(k)';
      sclabel{itypestart(ip(1))}{n}{1,2} = ''
    end
    stitle {itypestart(ip(1))}{1,1} = '<E(k)>_t' 
    sxlabel{itypestart(ip(1))}{1,1} = 'k';
    sylabel{itypestart(ip(1))}{1,1} = 'E(k)';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = sprintf('%s:%s',srundesc{n},'E{k}')
    end
    hold on;
    max(k)
    axis([0 max(k) 0 max(Ekavg)])
  end
  
  ip = strmatch('spEpar',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (21) ...<E(k)_{||}>_t spectrum:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh{itypestart(ip(1))} = loglog([1:numel(Ekpar_avg)],Ekpar_avg,symb{n});
    stitle {itypestart(ip(1))}{1,1} = '<E(k)_{||}>_t';
    sxlabel{itypestart(ip(1))}{1,1} = 'k';
    sylabel{itypestart(ip(1))}{1,1} = 'E_{||}(k)';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = sprintf('%s:%s',srundesc{n},'E_{||}');
    end
    hold on;
    max(k)
  end
  
  ip = strmatch('spEperp',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (22) ...<E(k_{\perp})>_t spectrum:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh{itypestart(ip(1))}{1,1} = loglog([1:numel(Ekperp_avg)],Ekperp_avg,symb{n});
    hold on;
    if exist('refparams.dat','file')
      [kpmin kpmax kpindex kpmult] =getgparam('refparams.dat','kpmin','kpmax','kpindex','kpmult') 
      x  = kpmin:1:kpmax;
      y  = kpmult*x.^(kpindex)
      if ndir == 1 | n == ndir
      loglog(x, y, symbd{n});
      sclabel{itypestart(ip(1))}{n}{1,1} = 'E(k_\perp)';
      sclabel{itypestart(ip(1))}{n}{1,2} = sprintf('k**(%f)',kpindex);
      end
    else
%     error('No refparams: E(k_\perp)');
      sclabel{itypestart(ip(1))}{n}{1,1} = 'E(k_\perp)';
      sclabel{itypestart(ip(1))}{n}{1,2} = ''
    end
    stitle {itypestart(ip(1))}{1,1} = ['<E(k_\perp)>_t'];
    sxlabel{itypestart(ip(1))}{1,1} = 'k_\perp';
    sylabel{itypestart(ip(1))}{1,1} = 'E(k_\perp)';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = sprintf('%s:%s',srundesc{n},'E(k\perp)');
    end
    hold on;
  end

%  [Ek2D_avg Ek2D teperp_index]   = spectavg('kspecperp', trange, 3,2);
  ip = strmatch('spEperp',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (23) ...<e(k_\perp,kz=0)>_t spectrum:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{2});
    hh{itypestart(ip(1))}{2,1} = loglog(k,Ek2D_avg,symb{n});
    hold on;
    if exist('refparams.dat','file')
      [kemin kemax keindex kemult] =getgparam('refparams.dat','kemin','kemax','keindex','kemult') 
      x  = kemin:1:kemax;
      y  = kemult*x.^(keindex)
      if ndir == 1 | n == ndir
      loglog(x, y, symbd{n});
      sclabel{itypestart(ip(1))}{n}{2,1} = 'e(k_\perp,kz=0)';
      sclabel{itypestart(ip(1))}{n}{2,2} = sprintf('k**(%f)',keindex);
      end
    else
%     error('No refparams: E(k_\perp,kz=0)');
      sclabel{itypestart(ip(1))}{n}{1,1} = 'e(k_\perp,kz=0)';
      sclabel{itypestart(ip(1))}{n}{1,2} = ''
    end
    stitle {itypestart(ip(1))}{2,1} = ['<e(k_\perp,kz=0)>_t' ];
    sxlabel{itypestart(ip(1))}{2,1} = 'k';
    sylabel{itypestart(ip(1))}{2,1} = 'e(k_\perp,kz=0)';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{2,1} = sprintf('%s:%s',srundesc{n},'e(k\perp,kz=0)');
    end
    hold on;
  end

  
  ip = strmatch('spH',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (25) ...<H(k)>_t spectrum:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh{itypestart(ip(1))} = loglog(k,abs(Hkavg),symb{n});
    stitle {itypestart(ip(1))}{1,1} = '<H(k)>_t';
    sxlabel{itypestart(ip(1))}{1,1} = 'k';
    sylabel{itypestart(ip(1))}{1,1} = '<H(k)>_t';
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = sprintf('%s:%s',srundesc{n},'H_{k}');
    end
    hold on;
  end
  

  ip = strmatch('spH2kE',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (25) ...<H(k)>_t / <k E(k)>_t  spectrum:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    hh{itypestart(ip(1))} = loglog(k,abs(Hkavg)./abs(Ekavg.*k'),symb{n});
    stitle {itypestart(ip(1))}{1,1} = '<H(k)>_t / <k E(k)>_t ';
    sxlabel{itypestart(ip(1))}{1,1} = 'k';
    sylabel{itypestart(ip(1))}{1,1} = '<H(k)>_t / <k E(k)>_t' ;
    if ndir > 1
      sclabel{itypestart(ip(1))}{n}{1,1} = sprintf('%s:%s',srundesc{n},'H_k/ k E_k');
    end
    hold on;
  end
  

  ip = strmatch('PiE',sdbtype);
  if bkeep(itypestart(ip(1)))
    % (26) ...<Pi(k)>_t spectrum: energy flux spectrum:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    %hh{itypestart(ip(1))} = loglog(k,efavg,symb{n});
    hh{itypestart(ip(1))} = semilogx(k,efavg,symb{n});
    stitle {itypestart(ip(1))}{1,1} = '<\Pi_E(k)>_t';
    sxlabel{itypestart(ip(1))}{1,1} = 'k';
    sylabel{itypestart(ip(1))}{1,1} = '<\Pi_E(k)>_t';
    if ndir > 1 
      sclabel{itypestart(ip(1))}{n}{1,1} = sprintf('%s:%s',srundesc{n},'\Pi_E');
    end
    hold on;
  end

  ip = strmatch('PiH',sdbtype);
  if bkeep(itypestart(ip(1)))
   % (27) ...<Pi_H(k)>_t spectrum: helicity flux spectrum:
    nf = nf + 1;
    figure(ifig{itypestart(ip(1))}{1});
    %hh{itypestart(ip(1))} = loglog(k,hfavg,symb{n});
    hh{itypestart(ip(1))} = semilogx(k,hfavg,symb{n});
    stitle {itypestart(ip(1))}{1,1} = '<\Pi_H(k)>_t';
    sxlabel{itypestart(ip(1))}{1,1} = 'k';
    sylabel{itypestart(ip(1))}{1,1} = '<\Pi_H(k)>_t';
    if ndir > 1 
      sclabel{itypestart(ip(1))}{n}{1,1} = sprintf('%s:%s',srundesc{n},'\Pi_H');
    end
    hold on;
  end

end % end of file-loop

if nf > nmaxplt
  error('getrotdat: Invalid number of plots');
end

cd(shome);
for j = 1:ntype
% i = itypestart(j);
  if ( ~bkeep(j) ) 
    continue;
  end
  
  nw = ndbtype{j}(1);
  for k = 1:nw
    nc = ndbtype{j}(k+1);
    figure(ifig{j}{k});
    hold on;
    snew = stitle{j}{k};
    % build legend for this window:
    sleg = cell(nc*ndir,1);
    ic = 1;
    doleg = 0;
    for n=1:ndir
      for i=1:nc
        sleg{ic} = sclabel{j}{n}{k,i};
        doleg = doleg || ~isempty(sleg{ic});
        ic = ic + 1;
      end
    end
    title(snew);
    xlabel(sxlabel{j}{k});
    ylabel(sylabel{j}{k});
    if doleg
      legend(sleg);
    end
    print( '-dpsc','-append', soutfile);
  end
end

