%%%%
%% Simulation of a Couette flow.
%% (Further Explorations  5.1, 5.2, & 6.2)
%%
%% Calculates the density, streaming velocity, and temperature profiles.
%%
%% Out-put data file: profs.dat
%%                    Column 1 spatial coordinate
%%                    Column 2 momentum density profile
%%                    Column 3 mass density profile
%%                    Column 4 temprature profile
%%                    Column 5 streaming velocity profile
%%
%% Tested with molsim 0.9.5 under GNU Octave 7.2.0 and Matlab R2020b
%%
%% NOTE: Under some GNU Octave versions the work space must be restarted after
%% each simulation
%%%%%


clear all;

%% Simulation parameters
temp0 = 1.0;
nloops = 1000000;

%% Interaction parameters (Change accordingly) 
cutoffFF = 2.5; epsilonFF = 1.0; sigmaFF = 1.0; awFF = 1.0; 
cutoffWW = 2^(1/6); epsilonWW = 1.0; sigmaWW = 1.0; awWW = 1.0; 
cutoffww = 2^(1/6); epsilonww = 1.0; sigmaww = 1.0; awww = 1.0; 

cutoffFW = 2.5; epsilonFW = 0.5; sigmaFW = 1.0; awFW = 1.0; 
cutoffFw = 2.5; epsilonFw = 0.5; sigmaFw = 1.0; awFw = 1.0;

cutoffwW = 2^(1/6); epsilonwW = 1.0; sigmawW = 1.0; awwW = 1.0; 

%% Wall restoring spring constant, wall speed, and thermostat relax. time
kspring = 500.0; Vwall = 0.2; tau=0.01;

%% Set misc. 
molsim('set','timestep', 0.001);
cutoffs=[cutoffFF, cutoffWW, cutoffww, cutoffFW, cutoffFw, cutoffwW];
molsim('set', 'cutoff', max(cutoffs));

%% Load conf. 
molsim('load', 'xyz', 'start_Cflow.xyz');

%% Set sample - sample only for fluid particle labeled F
molsim('sample', 'profiles', 'F', 300, 20);

%% Prepare the wall increment per time step
dx = 0.001*Vwall;

types = molsim('get', 'types');
indexw = find(types=='w');
npart = molsim('get', 'numbpart');
increment = zeros(1, npart);
increment(indexw) = dx;

%% Set the wall lattice sites from the current positions
molsim('set', 'virtualsites');

%% OMP enabled
molsim('set', 'omp', 4);

%% Main loop
for n=1:nloops

  molsim('reset');

  molsim('calcforce', 'lj', 'FF', cutoffFF, sigmaFF, epsilonFF, awFF);
  molsim('calcforce', 'lj', 'WW', cutoffWW, sigmaWW, epsilonWW, awWW);
  molsim('calcforce', 'lj', 'ww', cutoffww, sigmaww, epsilonww, awww);
  
  molsim('calcforce', 'lj', 'FW', cutoffFW, sigmaFW, epsilonFW, awFW);
  molsim('calcforce', 'lj', 'Fw', cutoffFw, sigmaFw, epsilonFw, awFw);

  molsim('calcforce', 'lj', 'wW', cutoffwW, sigmawW, epsilonwW, awwW);
  
  molsim('calcforce', 'lattice', 'W', kspring);
  molsim('calcforce', 'lattice', 'w', kspring);
   
  molsim('integrate', 'leapfrog');

  molsim('thermostat', 'relax', 'W', temp0, tau);
  molsim('thermostat', 'relax', 'w', temp0, tau);

  molsim('add', 'tolattice', increment, 0);

  if n>100000 %% Need equilibration
    molsim('sample', 'do');
  end
  
  if rem(n,1000)==0 
    molsim('print');
  end

end

molsim('save', 'FwW', 'final_Cflow.xyz');
molsim('clear');

%% Post-run analysis
data = load("profs.dat");
a = min(find(data(:,3)>0.01));
b = max(find(data(:,3)>0.01));

z = data(a:b,1)-data(a,1);
u = data(a:b, 5);

coeff=polyfit(z,u,1);

Ls = coeff(2)/coeff(1);
h = z(end);
upred = (Vwall*(z+Ls))./(h + 2*Ls);

zfit = linspace(-2, z(end)+2);
plot(z, u, 'o', z, upred, 'r--', 'linewidth', 2, ...
     zfit, coeff(1).*zfit + coeff(2));
