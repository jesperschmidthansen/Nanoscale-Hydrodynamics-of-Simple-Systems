%%%%
%% Simulation of the Lennard-Jones fluid system. 
%% (Further Explorations: 3.4)
%%
%% Calculates the shear stress correlation function
%% from the non-advective balance equation
%%
%% Change the wavevector from the wavenumber variable
%% k = 2*pi*wavenumber/L
%%
%% Always do a simulation ensemble and perform simple
%% descriptive statistics
%%
%% Tested with molsim 0.9.5 under GNU Octave 7.2.0 and Matlab R2020b
%%%%

clear all

octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if ( octave )
  mkoctfile evcorr.cpp
else
  mex evcorr.c
end


%% Simulation parameters
temp = 1.2;
dens = 0.8;
nloops = 200000;

wavenumber = 5;  %% <- change
sample_lvec = 200;

%% Interaction specifications
cutoff = 2.5; eps = 1.0; sigma = 1.0; aw = 1.0;

%% Load positions and set kin. temperature
molsim('load', 'xyz', 'start_LJ.xyz');
molsim('set', 'temperature', temp);

%% Initialisation of correlation array and counters
PP_corr = zeros(sample_lvec,1);
i_counter = 1;
s_counter = 0;

%% Need the number of particles
npart = molsim('get', 'numbpart');

%% Compatibility GNU Octave/Matlab
I = i;
%% Main loop
for n=1:nloops
  
  %% Reset forces etc
  molsim('reset');

  %% Calc. pair interaction forces
  molsim('calcforce', 'lj', 'AA', cutoff, sigma, eps, aw);

  %% Integrate forward in time
  molsim('integrate', 'leapfrog');

  %% Simple relxation method suffice
  molsim('thermostat', 'relax', 'A', temp, 0.01);
  
  %% Print misc. info to screen
  if rem(n,1000)==0
    molsim('print');
  end

  %% Sample 
  if rem(n,5)==0
    [P(i_counter,:) Pm(i_counter,:)]=molsim('get', 'shearpressure', wavenumber);

    i_counter = i_counter + 1;
    if ( i_counter > sample_lvec )
      a = P(:,1) + I*P(:,2); b = Pm(:,1) + I*Pm(:,2);

      if octave
	PP_corr = PP_corr + real(evcorr(a, b));
      else
	rc = evcorr(a,b);
	PP_corr = PP_corr + rc(:,1);
      end
      
      s_counter = s_counter + 1;
      i_counter = 1;
    end
  end
  
end

fprintf('\n');

%% Saving configuration
molsim('save', 'A', 'final_LJ.xyz');

volume = npart/dens;
PP_corr = PP_corr./(s_counter*volume);
tsample = linspace(0, 5*0.005*sample_lvec, sample_lvec);
plot(tsample, PP_corr);


molsim('clear');
