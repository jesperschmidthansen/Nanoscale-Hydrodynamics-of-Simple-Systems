%%%%
%% Simulation of the Kob-Andersen binary system.
%% (Further Explorations  3.5)
%%
%% Hydrodynamic corr. functions and stress autocorrelation function
%% calculated.
%%
%% Change the 'temp' variable to change the temperature; ensure the
%% system has relaxed before sampling (load final_KA.xyz)
%% 
%% Data files: sacf.dat gh-rho-acf.dat gh-trans-momentum-acf.dat
%%             gh-wavevector.dat
%%
%% You can use the m-file fltrans.m if ou wish to work in the freq. domain
%% 
%% The simulation is long - patience :)
%%
%% Tested with molsim 0.9.5 under GNU Octave 7.2.0 and Matlab R2020b
%% 
%% NOTE: Under some GNU Octave versions the work space must be restarted after
%% each simulation
%%%%%

clear all;

%% Simulation parameters
temp = 1.5;
dens = 1.2;
nloops = 1000000;

%% Interaction specifications
cutoff = 2.5; aw = 1.0;

epsAA = 1.0; sigmaAA = 1.0;
epsBB = 0.5; sigmaBB = 0.88;
epsAB = 1.5; sigmaAB = 0.8; 

%% Load 
molsim('load', 'xyz', 'start_KA.xyz');

molsim('set', 'temperature', temp);

%% Samplers
molsim('sample', 'sacf', 200, 20.0);
molsim('sample', 'hydrocorrelations', 200, 10.0, 10);

molsim('set', 'omp', 4);

%% Main loop
for n=1:nloops

  %% Reset forces etc
  molsim('reset');

  %% Calc. conservative interaction forces
  molsim('calcforce', 'lj', 'AA', cutoff, sigmaAA, epsAA, aw);
  molsim('calcforce', 'lj', 'BB', cutoff, sigmaBB, epsBB, aw);
  molsim('calcforce', 'lj', 'AB', cutoff, sigmaAB, epsAB, aw);

  %% Integrate forward in time
  molsim('integrate', 'leapfrog');

  %% Simple relaxation surfice
  molsim('thermostat', 'relax', 'A', temp, 0.01);
  molsim('thermostat', 'relax', 'B', temp, 0.01);
	 
  %% Sample
  molsim('sample', 'do');

  %% Print misc. info 
  if rem(n,1000)==0
    molsim('print');
  end
  
end

fprintf("\n Sim. done ... post processing\n");

%% Saving and freeing memory
molsim('save', 'AB', 'final_KA.xyz');
molsim('clear');

%% Post run analysis

%% 1. Calculate the viscosity
data = load('sacf.dat');
figure(1);
plot(data(:,1), data(:,2), data(:,1), hann(data(:,2)));
eta0 = trapz(data(:,1), hann(data(:,2)))./temp;


%% 2. Compare data for the transverse relaxation with predictions
%% Prefactor different from book due to conversion from momentum
%% current to velocity
k = load('gh-wavevector.dat');

data = load('gh-trans-momentum-acf.dat'); 
Cuu_pred = temp*dens*exp(-eta0*k(1).^2.*data(:,1)./dens); 

figure(2);
semilogx(data(2:200,1), data(2:200,2:2:10), 'o-', ...
	 data(2:200,1), Cuu_pred(2:200));
xlabel('t'); ylabel('Trans. vel. autocorr.');
 
%% 3. Plot the density autocorrelation funnction
%% Notice the Brillouin process 
data = load('gh-rho-acf.dat');
figure(3);
plot(data(:,1), data(:,2:2:10), '-o');
xlabel('t'); ylabel('Dens. autocorr');


