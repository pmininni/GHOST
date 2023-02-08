function t = tygres_times(sdir, sfile, N, itype)
%
% reads the times from the TYGRES output text file, sfile
%
% Usage:
%  t = tygres_times(',', 'spectrum.txt', 3072);
%
% sdir  : directory name (in), required
% sfile : name of file (e.g., 'spectrum.txt'), required
% N     : Actual resolution (that is, the 'equivalent resolution), required
% itype : type of file: 1==HD; 2==MHD, optional
% 
% Returns:
%  t      : vector of times

if nargin < 3
  error('tygres_times: Not enough input arguments. Do a "help tygres_times"');
end
if nargin < 4
  itype = 2;
end
if nargout < 1 
  error('tygres_times: Not enough output variables. Do a "help tygres_times"');
end

fid = fopen([sdir '/' sfile]);
if fid == -1
  error (['file cannot be opened for reading: ' sdir '/' sfile]);
end

nspec   = N/2 + 1;
t = zeros(0);

if ( itype ~= 1 && itype ~= 2 && itype ~= 3 )
  error('tygres_times: itype must be 1 (MHD) or 2 (HD)');
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


i     = 1;
while feof(fid) == 0 
  T     = textscan(fid, '%s', 1, 'Bufsize',1000000);
  isbad = 0;
  for kk = 1:ng
    G       = textscan(fid, '%s', 4, 'Bufsize',1000000);
    [mm nn] = size(G{1});
    if mm < 4
      swarn = sprintf('tygres_times: globals format error: t=%f\n',str2double(T{1}));
      warning(swarn);
      isbad = isbad + 1;
    end
  end


  for kk = 1:ns
    S  = textscan(fid, '%s', nspec, 'Bufsize',1000000);
    [mm nn] = size(S{1});
    if mm < nspec
      swarn = sprintf('tygres_times: spectrum format error: t=%f\n',str2double(T{1}));
      warning(swarn);
      isbad = isbad + 1;
    end
  end

  if isbad == 0
    t(i)    = str2double(T{1});
    i       = i + 1;
  end

end

fclose(fid);


