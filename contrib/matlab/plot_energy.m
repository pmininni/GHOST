% Plots energy as a function of time
% Assumes balance.txt is the output of an HD run
% Type 'plot_energy' in MATLAB to execute

% Path to the data
path = '../../3D/bin/';

% Reads balance.txt
%  balance(:,1) = time
%  balance(:,2) = energy (v^2)
%  balance(:,3) = enstrophy (w^2)
%  balance(:,4) = energy injection rate
balance = load(strcat(path,'balance.txt'));

% Plots energy vs. time in a window
figure(1)
plot(balance(:,1),balance(:,2))
xlabel('time')
ylabel('Energy')

% Saves plot to an EPS file
print('-f1','figure','-deps')
