function cintrp = interpspect2d(spec2d, nnx, nny, dangles)
%
% function c = interpspect2d(spec2d, 512, 512, [10 20 40]);
%
% Interpolates 2D spectrum, spec2d, along a different lines at angles,
% dangles to produce 1d line spectra.
%
%           Input:
%  spec2d : 2d spectrum, obtained, e.g., from spect2davgn
%  nnx    : full problem size in perp direction, s.t. spectrum in x is nnx/2+1
%  nny    : full problem size in parallel direction, s.t. spectrum in y is nny/2+1
%  dangles: array of angles defining lines through 2d spectrum, in degrees
%
%           Returns:
%  cintrp : cell array of spectra at each angle
%
c     = [];

if nargin<4
  error('Incorrect number of arguments');
end

  sang  = dangles * pi/180;

  nx      =    nnx/2 + 1;
  ny      =    nny/2 + 1;
  k_perp  = [1:nx];
  k_para  = [1:ny];
  [kk_perp, kk_para] = meshgrid(k_perp, k_para);

  cintrp = zeros(numel(dangles),nx);
  for i=1:numel(sang)
      ang      = sang(i);
      yq       = tan(0.5*pi-ang) .* (k_perp-1) + 1;
%     k_vector = sqrt(k_perp.^2 + yq.^2);
      c        = interp2(kk_perp,kk_para,spec2d/sin(ang),k_perp,yq,'cubic');
      ind      = find(isnan(c));
      c(ind)   = 0.0;
      cintrp(i,:) = c;
  end

