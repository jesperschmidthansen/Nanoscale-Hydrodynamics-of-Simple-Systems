%%%%
%%
%% Simulation of the modified Ryckaert-Bellemann butane model
%% (Further Explorations  4.1)
%%
%% Molecular shear pressure correlation function calculated
%%
%% Out-put data file: msacf.dat
%%
%% Tested with molsim 0.9.5 under GNU Octave 7.2.0 and Matlab R2020b
%%
%% NOTE: Under some GNU Octave versions the work space must be restarted after
%% each simulation
%%%%%

clear all;

%% Simulation parameters
temp0 = 4.0;
dt  = 0.002;
nloops = 1000000;

%% Interaction specifications
cutoff = 2.5; epsCC = 1.0; sigmaCC = 1.0; aw = 1.0;

bondlength = 0.4;
bondconstant = 2000.0;
bondangle = 1.9;
angleconstant = 400.0;
torsionparam = [15.5000,  20.3050, -21.9170, -5.1150,  43.8340, -52.6070];

%% Set stuff
molsim('set', 'exclusion', 'molecule');
molsim('set','timestep', dt);
molsim('set', 'temperature', temp0);

%% Load positions, set temp, remove intra-molecular pair-interaction etc
molsim('load', 'xyz', 'start_butane.xyz');
molsim('load', 'top', 'start_butane.top');

%% Sampler
molsim('sample', 'msacf', 300, 2.0); 

molsim('set', 'omp', 4);

%% Main loop
for n=1:nloops

  molsim('reset')

  molsim('calcforce', 'lj', 'CC', cutoff, sigmaCC, epsCC, aw);
  molsim('calcforce', 'bond', 0, bondlength, bondconstant);
  molsim('calcforce', 'angle', 0, bondangle, angleconstant);
  molsim('calcforce', 'torsion', 0, torsionparam);

  molsim('thermostat', 'nosehoover', 'C', temp0, 10.0);
  molsim('integrate', 'leapfrog');

  molsim('sample', 'do');

  %% Print misc. info 
  if rem(n,1000)==0
    molsim('print');
  end
  
end

fprintf('\n');

x=load('msacf.dat');
plot(x(:,1), x(:,2), '-o');


molsim('save', 'C', 'final_butane.xyz');
molsim('clear');



