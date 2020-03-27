% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: March 11th, 2020

% Produces the baseline parameters, runs baseline simulation, and saves
% results
tic
clear;

%Choose which FcR: FcR3aV = 1 FcR3aF = 2 FcR2aH = 3 FcR2aR = 4
fcr_id = 1;

% Read in parameters and ICs
file_to_load = Parameters(fcr_id); 
load(file_to_load)

% Run Baseline simulation
[ybase steadystate] = Simulate(params, paramnames, [],[],fcr_ttl);

% Display whether or not simulation reached steady state as defined within
% "simulate" function
if steadystate>0
              disp(["Baseline did not run to steady state"])
        else
            disp(["Baseline ran to steady state"])
end
        
% Save the results
save(strcat(datestr(today()),'_Base_Simulation_'+fcr_ttl+'.mat'),'ybase','params',...
    'paramnames','complexes');

toc