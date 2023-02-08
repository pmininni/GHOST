function [newarr I] = arrunique(oldarr,n)
%
% function arr = arrunique(oldarr,4)
%
% Assumes first column of oldarr is time, and removes all but last version
% of that time stamp frmo array for all columns.
%
% Usage:
%     newarr = arrunique(tbous, 4)
%
% Arguments:
%    oldarr  : array we want to clean up
%    n       :  number of digits to compare when determining uniqueness (optionsl; default is 15)
%
%  Output:
%    newarr : cleaned up array
%    I      : indices s.t. newarr(:,:) = oldarr(I,:)
%    

  if nargin < 1
    error('Must provide filename and time index');
  end
  if nargin < 2
    n = 15;
  end

  sz = size(oldarr);

  t0     = oldarr(:,1);
  [tnewarr I J] = iunique(t0,n,'last');
  newarr     = zeros(length(I),sz(2));
  
  for i = 1:sz(2)
    newarr(:,i)=oldarr(I,i);
  end


end
