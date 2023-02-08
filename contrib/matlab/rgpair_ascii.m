function [data time] = rgpair_ascii(sfile,irank)
%
% Reads GHOST particle pair data from ASCII formatted file 
% with records of the form
% 
%   time   x_component  y_component  z_component
%
%  Usage:
%      [data t] = rgpair_ascii('glag.010.txt',3);
%
%  Inputs:
%
%  sfile    : input file (required)
%  irank    : Rank of each record. Default is 3.

%  Outputs:
%
%  data    : field components arra of size (irank,ntime)
%  time    : time array of records, of size ntime
%
if nargin < 1
  error('Input file name prefix at least! Do a "help rglag".');
end
if nargin < 2
  irank = 3;
end

fp = fopen(sfile,'r');
if  fp  == -1
  error(sprintf('File %s cannot be opened for reading',sfile));
end
ntime = 0;
while ~feof(fp) 
  sdat =  textscan(fp,'%f',irank+1,'MultipleDelimsAsOne',1);
  ndat = size(sdat{1});
  if ndat(1) == irank+1
    ntime = ntime+1;
  end
end
if ntime == 0 
  error(sprintf('No valid records in file %s. Check irank parameter.',sfile));
end

frewind(fp);
data = zeros(irank,ntime);
time = zeros(ntime,1);
i = 0;
while ~feof(fp) & i < ntime 
  i = i + 1;
  sdat =  textscan(fp,'%f',irank+1,'MultipleDelimsAsOne',1);
  time(i) = sdat{1}(1);
  data(1:irank,i) = sdat{1}(2:irank+1);;
end
fclose(fp);

if i ~= ntime
  error('re-read of file failed');
end
