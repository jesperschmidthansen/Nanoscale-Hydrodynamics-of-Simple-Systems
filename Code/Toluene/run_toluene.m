%%%%
%% Simulation of model toluene (featuring the molsim task manager) 
%% (Further Explorations  4.2)
%%
%% Calculates the transverse momentum autocorrelation
%%
%% Out-put data file: mgh-trans-momentum-acf.dat
%%                    mgh-wavevector.dat
%%                    msacf.dat
%%
%% Toluene model from Mol. Sim., 47:17, 1391-1401
%%
%% Tested with molsim 0.9.5 under GNU Octave 7.2.0 and Matlab R2020b
%%
%% NOTE: Under some GNU Octave versions the work space must be restarted after
%% each simulation
%%%%%

clear all;

%% Simulation parameters
temp0 = 4.0;
dens0 = 1.95;
dt  = 0.001;
nloops = 10000;  % <- Small value for testing purpose, change when simulating
nk = 10;

%% Interaction specifications
cutoff = 2.5; epsilon = 1.0; sigma = 1.0; %% aw=1 per default

bondlength_0 = 0.4;
bondlength_1 = 0.38;
bondconstant = 48910;

bondangle = 2.09;
angleconstant = 1173;

torsionparam_0 = [0.0, 133.0, 0.0 0.0 0.0];
torsionparam_1 = [0.0, -133.0, 0.0 0.0 0.0];

%% Set temp, remove intra-molecular pair-interaction etc
molsim('set','timestep', dt);
molsim('set', 'temperature', temp0);
molsim('set', 'exclusion', 'molecule');

%% Load positions and top file 
molsim('load', 'xyz', 'start_toluene.xyz');
molsim('load', 'top', 'start_toluene.top');

%% Sampler
molsim('sample', 'mhydrocorrelations', 500, 10.0, nk); 

%% Task manager
molsim('task', 'lj', 'CC', cutoff, sigma, epsilon, 1);
molsim('task', 'bond', 0, bondlength_0, bondconstant, 2);
molsim('task', 'bond', 1, bondlength_1, bondconstant, 2);
molsim('task', 'angle', 0, bondangle, angleconstant, 2);
molsim('task', 'torsion', 0, torsionparam_0, 2);
molsim('task', 'torsion', 1, torsionparam_1, 2);

molsim('set', 'omp', 2);

%% Main loop
for n=1:nloops

  molsim('reset')

  molsim('task', 'do', 2);

  molsim('thermostat', 'nosehoover', 'C', temp0, 10.0);
  molsim('integrate', 'leapfrog');

  molsim('sample', 'do');

  if rem(n,1000) == 0
    fprintf("\r Loop no. %d  ", n);
  end

end

fprintf('...done\n');

molsim('save', 'C', 'final_toluene.xyz');
molsim('clear');

%% Post-run processing
data = load('mgh-trans-momentum-acf.dat');
k = load('mgh-wavevector.dat');

omega= linspace(0,1);
for n=1:nk
  Cw = fltrans(data(:,1), hann(data(:,2*n))./dens0^2, omega);
  eta(n) = temp0/(k(n)^2*real(Cw(1)));
end

plot(k, eta, 'o-');

