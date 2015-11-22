% Reads a binary file and plots a cut in the x-y plane
% Type 'plot_bindata' in MATLAB to execute

% Path to the data
path = '../../3D/bin/outs/';

% Spatial resolution
N = 128;

% Reads binary files
fid = fopen(strcat(path,'vx.0001.out'));
vx = fread(fid,N*N*N,'float32');
fclose(fid);
vx = reshape(vx,[N N N]);

% Show a horizontal cut of the field in the middle of the box
figure(1)
imagesc(vx(:,:,N/2))
xlabel('x')
ylabel('y')
