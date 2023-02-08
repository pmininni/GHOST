function [ib ie] = filerange(files, trange)
%
% Determines file index range, ib-ie, for GHOST files, given desired time
% index range, trange
% 
%           [ib ie] = filerange(files, trange)
%
%           Input:
%  files  : list of files; array of strings
%  trange : string of form 'start:stop' time indices; if  stop is 'end', all files from 
%           start until the end of file list are averaged
%
%           Output:
%  ib     : beginning index into 'files' array
%  ie     : ending index into 'files' array
%

  fnum = numel(files);
  [t r] = strtok(trange,':');
  ibeg = str2num(t);
  if isempty(ibeg)
    error('invalid range initialization');
  end
  ib = 0;
  bfound = false;
  while ib<fnum && ~bfound
    ib = ib + 1;
    fn = files{ib}(1:strfind(files{ib},'.txt')-1);
    [t r] = strtok(fn,'.');
    index = str2num(r(2:end));
    if ~isempty(index)
      bfound = index >= ibeg;
    else
      bfound = false;
    end
  end
  if ~bfound
    error('unable to find start');
  end

  [t r] = strtok(trange,':');
  if strcmp(strtrim(r(2:end)),'end')
    iend = fnum;
  else
    iend = str2num(r(2:end));
    if isempty(iend)
      error('invalid range termination');
    end
  end
  ie = fnum+1;
  bfound = false;
  while ie>=1 && ~bfound
    ie = ie - 1;
    fn = files{ie}(1:strfind(files{ie},'.txt')-1);
    [t r] = strtok(fn,'.');
    index = str2num(r(2:end));
    if ~isempty(index)
      bfound = index <= iend;
    else
      bfound = false;
    end
  end
  if ~bfound
    error('unable to find end');
  end

