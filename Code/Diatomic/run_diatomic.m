%%%%
%%
%% Simulation of the generic diatomic molecule
%% (Further Explorations  4.3 & 4.4)
%%
%% Out-put data files:
%%   mgh-wavevector.dat, mgh-long-angmomentum-acf.dat (4.3)
%%
%% Tested with molsim 0.9.5 under GNU Octave 7.2.0 and Matlab R2020b
%%
%% NOTE: Under some GNU Octave versions the work space must be restarted after
%% each simulation
%%%%%

clear all;

%% Simulation parameters
temp0 = 2.0;
dens0 = 0.85;

nloops = 10000;  % <- Small value for testing purpose, change when simulating
charge = false;
q = 1.0;

%% Interaction specifications
cutoff = 2.^(1/6);
eps = 1.0; sigma = 1.0; aw=1.0;

bondlength = 1.0;
bondconstant = 2500.0;

%% Load positions, set temp, remove intra-molecular pair-interaction etc
molsim('load', 'xyz', 'start_diatomic.xyz');
molsim('load', 'top', 'start_diatomic.top');

molsim('set', 'temperature', temp0);
molsim('set', 'exclusion', 'molecule');

npart = molsim('get', 'numbpart');

if charge
  z = (-1).^(1:npart);
  molsim('set', 'charges', q*z);
end

molsim('sample', 'msacf', 100, 10.0);
molsim('sample', 'mhydrocorrelations', 200, 20.0, 10);

molsim('set', 'omp', 4);

%% Main loop
for n=1:nloops

  molsim('reset')

  molsim('calcforce', 'lj', 'CC', cutoff, sigma, eps, aw);
  molsim('calcforce', 'bond', 0, bondlength, bondconstant);

  if charge
    molsim('calcforce', 'coulomb', 2.5);
  end
  
  molsim('thermostat', 'nosehoover', 'C', temp0, 10.0);
  molsim('integrate', 'leapfrog');
  
  molsim('sample', 'do');

  %% Print misc. info 
  if rem(n,1000)==0
    molsim('print');
  end
  
end

fprintf('\n');

molsim('save', 'C', 'final_diatomic.xyz');
molsim('clear');

%% Post-run data analysis
if ~charge
  
  data = load("mgh-long-angmomentum-acf.dat");
  k = load("mgh-wavevector.dat");
  omega = logspace(0, 2,1000);

else

  data = load("mgh-dipole-acf.dat");
  k = load("mgh-wavevector.dat");
  omega = linspace(0, 20);

end


for n=1:length(k)
  Cw = fltrans(data(:,1), hann(data(:,2*n)), omega);
  
  iCw = -imag(Cw);
  i = find(iCw==max(iCw));
  omega_peak(n) = omega(i);
end

c = polyfit(k.^2, omega_peak', 1);
plot(k.^2, omega_peak, 'o', k.^2, c(1).*k.^2 + c(2), '--'); 
