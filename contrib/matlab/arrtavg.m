function tavg  = arrtavg(tarr,trng,varargin)
%
% function tavg  = arrunique(tarr,trng,varargin)
%
% Assumes rows of tarr advance in time, given by first column.
% Time averages specified rows over specified time interval, and returns.
%
% Usage:
%     tarr = arrtavg(tarr,[1.0 1.25],3,4) 
%
% Arguments:
%    tarr     : array we want to average specified columns in time for
%    trng     : len 2 array of start and stop time. If [], then use entire
%               time interval.
%    varagrin : list of column ids to find averages for. If not specified, 
%               compute for all columns.
%
%  Output:
%    tavg : array of size num(varargin) with averages for each column over
%           specified time interval
%    
  if nargin < 2
    error('Must provide time-based array, trng, and trng (which may be empty)');
  end

  % Try to 'clean' data a bit:
  [tt0 JJ, KK] = iunique(tarr(:,1),4,'last'); % find unique times to 4 digits
  sz = size(tarr);
  tarrnew = zeros(numel(JJ),sz(2));
  for j = 1:sz(2)
    tarrnew(:,j) = tarr(JJ,j);
  end

  sz = size(tarrnew);
  icol = 2:sz(2);
  ncol = numel(icol);
  if isempty(trng) 
    trng = [tarrnew(1,1) tarrnew(end,1)];
  else
    if numel(trng) < 2
      error('trng must be 2-element array (but may be empty)');
    end
  end
  if nargin > 2
    ncol = size(varargin,2);
    icol = cell2mat(varargin);
  end
  
  MM  = find(tarrnew(:,1)>=trng(1));
  NN  = find(MM==MM(1));

  II   = find(tarrnew(NN(end):end,1)>=trng(1) & tarrnew(NN(end):end,1)<=trng(2) );
  if length(II) < 1
    [~, idx] = min(abs(tarrnew(:,1) - trng(1)));
    if idx<0 | idx>length(tarrnew(:,1))
      error('Time interval not found in data');
    end
    II = idx;
  end

  if ( numel(II) >= 2 )
    I1   = II(1:end-1);
    I2   = II(2:end);
    t0   = tarrnew(II,1);
%   dt   = diff(t0); 
    ttot = abs(max(t0)-min(t0));
    dt   = ttot/(numel(II)-1); 
    for j = 1:ncol
      ic      = icol(j);
%     q       = 0.5*(tarrnew(I1,ic) + tarrnew(I2,ic));   % find quantity at midpoint of interval
      q       = tarrnew(II,ic);
%%%   tavg(j) = sum(q*dt)/ttot;   % Integral q dt / T
      tavg(j) = mean(tarrnew(II,ic));
    end
  else
    if numel(II) == 1
      for j = 1:ncol
        ic      = icol(j);
        tavg(j) = tarrnew(II(1),ic);
      end
    end
  end

end
