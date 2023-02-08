function [mom lx xpdf] = getpdfmom(fn, imom)
%
% function mom = getpdfmom(fn, imom)
%
% getpdfmom reads the 1d pdf file, fn, and computes the imom-moment of the
% distribution, returning in mom
%
% Usage:
%     mom2 = findkr('dtdzpdf.025.txt', 2)
%
% Arguments:
%    fn     : 1d pdf filei (required).
%    imom   : time stamp of file. Default is 1 (average computed).
%
%  Output:
%    mom    : moment returned 
%    

  if nargin < 1
    error('Must provide filename');
  end
  if nargin < 2
    imom = 1
  end

  xfp = fopen(fn,'r');
  if  xfp  == -1
    error(sprintf('File %s cannot be opened for reading',fn));
  end


  xdat  =  textscan(xfp,'%s',1,'delimiter','\n');
  xhead = char(xdat{1});
  E     = getpdfheader(xhead);
  renst = [E{1}(1) E{1}(2)];
  qavg  = E{2}(1);
  qsig  = E{3}(1);
  nbins = E{4}(1);
  dolog = E{5};
  nkeep = E{6};
  xdat =  textscan(xfp,'%f','delimiter','\n','CommentStyle',{'#'});
  fclose(xfp);

imom
dolog
nbins
qavg
qsig
renst
nkeep

  xpdf  = cell2mat(xdat);
  xpdf  = xpdf';

  if dolog <= 0
    del(1) = ( renst(2) - renst(1) ) / double(nbins);
    lx = (renst(1)+double([0:nbins-1]+0.5)*del(1));
  else
    del(1) = ( log10(renst(2)) - log10(renst(1)) ) / double(nbins);
    lx = log10(renst(1))+double([0:nbins-1]+0.5)*del(1);
  end

% xpdf = xpdf *double(nbins) / (double(nkeep) * (renst(2)-renst(1)));
  mom = sum((lx.^imom).*xpdf) / sum(xpdf);

  return;
