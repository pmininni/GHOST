% Plots energy spectra as a function of time
% Type 'plot_spectrum' in MATLAB to execute

% Path to the data
path = '../../3D/bin/';

% Spatial resolution
N = 128;

% Array with wavenumbers starting at 1
k = [1:N/2+1];

% Reads and plots all spectra in the directory.
% We only plot one every five spectra, starting
% from the second.
filelist = dir(strcat(path,'kspectrum.*.txt'));
nfiles = numel(filelist);
figure(1)
for i = 2:5:nfiles,
  ene = load(strcat(path,filelist(i).name));
  loglog(k,ene)
  hold all
end
xlim([1 N/3])
xlabel('k')
ylabel('E(k)')
hold off
