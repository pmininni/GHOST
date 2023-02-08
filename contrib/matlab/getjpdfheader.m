function C = getjpdfheader(sheader)

%
%function C = getjpdfheader(sheader)
%
% parsepdfheader parses the GHOST joint PDF  header of
% the form:
%  #enst_rng=[  0.83757297E-03,  0.59941225E+06]; enst_avg=  0.42427500E+04; \
%   enst_sig= 0.04675998e+03; diss_rng=[  0.32331958E-01,  0.14248777E+06]; diss_avg=  0.  24903667E+04; diss_sig= 0.4675988e+03; nbin=[ 100, 100]; blog=[1,1]; nkeep=[100000,120000]
%
% returning the range arrays, the averages, the number of bins, and the 
% flag for whether data was binned logarithmically
%
% Usage: 
%         
%      C = getpdfheader(sheader);
%
% Arguments:
%	sheader: header string (or cell)
%
% Output:
%        Cell array with values for the dynamic range, average, no. bins, and log flag
%        C{1} = x var. [rng0 rng1];
%        C{2} = x var. avg.
%        C{3} = x var. sig.
%        C{4} = y var. :[rng0 rng1];
%        C{5} = y var. avg.
%        C{6} = y var. sig.
%        C{7} = [nbin_x nbin_y]
%        C{8} = ilog;
%        C{9} = nkeep;
%

if ( nargin < 1 )
  error('parsejpdfheader: Must provide the header!');
end

%  #enst_rng=[  0.83757297E-03,  0.59941225E+06]; enst_avg=  0.42427500E+04; \
%   diss_rng=[  0.32331958E-01,  0.14248777E+06]; diss_avg=  0.24903667E+04; nbin=[ 100, 100]; blog=1

C      = cell(6,1);
i0     = 0 ;
i1     = i0+strfind(sheader(i0+1:end),'['); 
i2     = i0+strfind(sheader(i0+1:end),','); 
i3     = i0+strfind(sheader(i0+1:end),']'); 
rng(1) = str2num(sheader(i1(1)+1:i2(1)-1)); 
rng(2) = str2num(sheader(i2(1)+1:i3(1)-1)); 
C{1,1} = [rng(1) rng(2)]; 

i0     = i3(1)+2; 
i4     = i0+strfind(sheader(i0+1:end),'='); 
i5     = i0+strfind(sheader(i0+1:end),';'); 
avg    = str2num(sheader(i4(1)+1:i5(1)-1));
C{2,1} = [avg];
i6     = i5(1)+strfind(sheader(i5(1)+1:end),'=') ;
i7     = i5(1)+strfind(sheader(i5(1)+1:end),';') ;
sig    = str2num(sheader(i6(1)+1:i7(1)-1)) ;
C{3,1} = [sig];

i0     = i7(1)+1;
i6     = i0+strfind(sheader(i0+1:end),'[') ;
i7     = i0+strfind(sheader(i0+1:end),',') ;
i8     = i0+strfind(sheader(i0+1:end),']') ;
rng(1) = str2num(sheader(i6(1)+1:i7(1)-1)) ;
rng(2) = str2num(sheader(i7(1)+1:i8(1)-1)) ;
C{4,1} = [rng(1) rng(2)];

i0     = i8(1)+1;
i9     = i0+strfind(sheader(i0+1:end),'=');
i10    = i0+strfind(sheader(i0+1:end),';');
avg    = str2num(sheader(i9(1)+1:i10(1)-1));
C{5,1} = [avg];
i11    = i10(1)+strfind(sheader(i10(1)+1:end),'=');
i12    = i10(1)+strfind(sheader(i10(1)+1:end),';');
sig    = str2num(sheader(i11(1)+1:i12(1)-1));
C{6,1} = [sig];

i0     = i12(1)+1;
i13    = i0+strfind(sheader(i0+1:end),'[');
i14    = i0+strfind(sheader(i0+1:end),',');
i15    = i0+strfind(sheader(i0+1:end),']');
rng(1) = str2num(sheader(i13(1)+1:i14(1)-1));
rng(2) = str2num(sheader(i14(1)+1:i15(1)-1));
C{7,1} = [int32(rng(1)) int32(rng(2))];


i0     = i15(1)+1;
i16    = i0+strfind(sheader(i0+1:end),'[');
i17    = i0+strfind(sheader(i0+1:end),',');
i18    = i0+strfind(sheader(i0+1:end),']');
ilog1  = str2num(sheader(i16(1)+1:i17(1)-1));
ilog2  = str2num(sheader(i17(1)+1:i18(1)-1));
C{8,1} = [int32(ilog1), int32(ilog2)];

i0     = i18(1)+1;
i19    = i0+strfind(sheader(i0+1:end),'=');
ilog1  = str2num(sheader(i19(1)+1:end));
%i20    = i0+strfind(sheader(i0+1:end),',');
%ilog2  = str2num(sheader(i20(1)+1:end));
%C{9,1} = [int32(ilog1), int32(ilog2)];
C{9,1} = [ilog1];

