function varargout = getgparam(fn, varargin)

% function varargout = getgparam(fn, varargin)
%
% getgparam reads the GHOST 'parameter' datafile ( usually 'parameter.txt'), fn, and
% provides the data in the file corresponding to the parameter names
% specified as 'varargin' arguments. The output order is the same as the
% order of the parameter names given as input.
%
% Usage: 
%         
%      [sstep dt] = getgparam('parameter.txt','sstep','dt');
%
% Arguments:
%	fn: name of parameter file (required)
%       varargin: comma-separated list of strings that provide the parameter name 
%                 in the file
% Output:
%        values for the parameter names in the order in which the  names were
%        provided on input.
%
if ( nargin < 2 )
   error('getgparam: Must provide at least 2 arguments!');
end
nparams = size(varargin,2);

[fp message] = fopen(fn,'rt');
if fp < 0 
  error(sprintf('getgparam: Error opening file %s: %s',fn, message));
end
C = textscan(fp,'%s','Delimiter','\n','endOfLine','\n');
maxparams = numel(C{1});
fclose(fp);

for n = 1:nparams
  bfound = 0;
  varargin{n};
  lvar  = numel(varargin{n});
  for j = 1:maxparams
    str  = char(C{1}(j));
    ii = strfind(str,'=');
    if isempty(ii)
      continue;
    end
    [tok rem] = strtok(str,'=');
    
    svar = '';
    if isempty(tok) | ~strcmp(strtrim(tok), strtrim(varargin{n}))
      continue;
    end
    svar = strtrim(tok); % variable name in file
    
     % Now, get value for variable:
     ss = rem(2:end);
     i2   = strfind(ss,'"');
    if ~isempty(i2)  % is string data
      if length(i2) < 2
        error(sprintf('getgparam: Syntax error for string variable: %s', svar));
      end
      val = ' ';
      if ( i2(2) > i2(1) )
        val = ss(i2(1)+1:i2(2)-1);
      end
      bfound = 1;
      break;
    else              % numeric data

      % find comment
      i2 = strfind(ss,'!'); % check for comments
      if ~isempty(i2)  % is comment:
        imax = i2-1;
      else             % no comment
        imax = length(ss);
      end
      sval = strtrim(ss(1:imax));
      val  = str2num(sval);
      if exist('val','var')
        bfound = 1;
        break;
      else
        error(sprintf('getgparam: Problem with variable %s assignment.', svar));
      end
    end % end, data type check
  end % end j-loop
  if bfound == 0
    error(sprintf('getgparam: Parameter %s not found in %s',varargin{n},fn));
  end
  varargout{n} = val;
end

