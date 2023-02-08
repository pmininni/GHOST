 function C = gfclean(fn, nunq, icols, bdoz)
% 
% function C  = gfclean(fn,4)
%
% Function gfclean takes as input a GHOST time-dep text file, 
% (e.g., bance) and (1) removes rows that have the specified columns 
% 0 or NaN, (2) removes repeated time stamps ( time in col 1), retaining data 
% that corresponds to the final repeated time.
%
% Usage:
% function C  = gfclean('balance.txt',4, [2 3])
% 
% Arguments:
%  X    : 2D array with rows consisting of (Time , Q, Q2, ..., Qn). Required.
%  nunq : number of digits to compare when determining uniqueness (optionsl)
%         of time variable. Default is 4 digits.
%  icols: vector of table columns to consider that are checked for either 
%         0's or Nan's. If any NaNs are found in _any_ columns, or if 
%         0's are found in all specified columns, that row is removed in 
%         the 'clean' return array. Default is '[2 3]'. 
%  bdoz : do check for zeros in all columns. Default is 1.
%
% Output: 
%  C   : 'cleaned' array.
%

if nargin < 1 
    error('Must specify GHOST time-file');
end

if nargin < 2 
  nunq  = 4;
  icols = [2 3];
  bdoz  = 1;
end

if nargin < 3
  icols = [2 3];
  bdoz  = 1;
end

if nargin < 4
  bdoz  = 1;
end

if length(icols) <= 1
  error('Must specify column(s) to examine');
end

nicols  = length(icols);
%b0      = textread(fn,'','commentstyle','shell','TreatAsEmpty','"Infinity"','TreatAsEmpty','"Inf"');
b0      = textread(fn,'','commentstyle','shell');
ncols   = size(b0,2);

if bdoz > 0
  % Check for 0's in _all specified_ columns:
  express = sprintf(' find(~(b0(:,%d) == 0',icols(1));
  for j = 2:nicols
    i = icols(j);
    express = strcat(express,sprintf(' & b0(:,%d) == 0',i)) ;  
  end
  express = strcat(express,'));');
  [I J] = eval(express);

  b1   = zeros(length(I),ncols);
  for i = 1:ncols
    b1(:,i) = b0(I,i);
  end

else
  b1 = b0;
end
clear b0;

% Check for NaN in _any_ column:
express = sprintf(' find(~isnan(b1(:,%d))',2);
for i = 3:ncols
  express = strcat(express,sprintf(' & ~isnan(b1(:,%d))',i)) ;  
end
express = strcat(express,');');
[I J] = eval(express);

b2   = zeros(length(I),ncols);
for i = 1:ncols
  b2(:,i) = b1(I,i);
end


% With first col as time, find unique times, choosing last:
t0      = b2(:,1); 
[t I J] = iunique(t0,nunq,'last'); % find unique indices of valid data
C   = zeros(length(I),ncols);
C(:,1) = t;
for i = 2:ncols
  C(:,i) = b2(I,i);
end

end
