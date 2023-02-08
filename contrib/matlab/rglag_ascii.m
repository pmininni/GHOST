function [data np time] = rglag_ascii(sfile,irank,hdonly)
%
% Reads GHOST particle data from ASCII formatted file
%
%  Usage:
%      [data np t] = rglag_ascii('glag.010.txt');
%
%  Inputs:
%
%  sfile    : input file (required)
%  irank    : Rank of each record. Default is 3.
%  hdonly   : if >0, fill only np, and t, with data=null

%  Outputs:
%
%  data    : particle data in (irank,np) array with x,y,z in each record
%  np      : no. particles required, as found in first file
%  time    : time stamp of file
%
if nargin < 1
  error('Input file name prefix at least! Do a "help rglag".');
end
if nargin < 2
  irank = 3;
  hdonly = 0;
end
if nargin < 3
  hdonly = 0;
end


fp = fopen(sfile,'r');
if  fp  == -1
  error(sprintf('File %s cannot be opened for reading',sfile));
end
snp   = textscan(fp,'%d\n',1,'MultipleDelimsAsOne',1);
np = snp{1};
stime = textscan(fp,'%f\n',1,'MultipleDelimsAsOne',1);
time = stime{1};
data = [];
if ~isnumeric(np)
  serr = sprintf('Bad data in file %s,: np',sfile);
  error(serr);
end
if hdonly > 0
 data = 0;
 fclose(fp);
 return;
end

data = zeros(np,irank);
i = 0;
while ~feof(fp) & i < np
  i = i+1;
  sdat =  textscan(fp,'%f %f %f\n',1,'MultipleDelimsAsOne',1);
  data(i,1:irank) = cell2mat(sdat);
end
fclose(fp);

