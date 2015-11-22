% Reads binary files in a directory and does
% postprocessing (e.g., to visualize later with
% VAPOR). Note that for big runs, GHOST can do some
% automatic postprocessing (e.g., compute vorticity)
% at run time.
% Type 'postprocess.m' in MATLAB to execute

% Path to the binary data
path = '../../3D/bin/outs/';

% Spatial resolution
N = 128;
dx = 2*pi/N;

% Reads binary files, computes vertical vorticity, 
% and saves in a new binary file named 'wz.NNNN.out'
filelist = dir(strcat(path,'vx.*.out'));
nfiles = numel(filelist);
for i = 1:nfiles,
  fid = fopen(strcat(path,filelist(i).name));
  vx = fread(fid,N*N*N,'float32');
  fclose(fid);
  vx = reshape(vx,[N N N]);
  ind = filelist(i).name(4:7);
  str = ['vy.' ind '.out'];
  fid = fopen(strcat(path,str));
  vy = fread(fid,N*N*N,'float32');
  fclose(fid);
  vy = reshape(vy,[N N N]);
  adv = circshift(vy,[-1 0 0]);
  ret = circshift(vy,[1 0 0]);
  vy = (adv-ret)/(2*dx); % dv_y/dx
  adv = circshift(vx,[0 -1 0]);
  ret = circshift(vx,[0 1 0]);
  vx = (adv-ret)/(2*dx); % dv_x/dy 
  wz = vy-vx;
  str = ['wz.' ind '.out'];
  fid = fopen(strcat(path,str),'w');
  fwrite(fid,wz,'float32');
  fclose(fid);
end
