% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: March 11th, 2020

tic
clear;

%Choose which FcR: FcR3aV = 1 FcR3aF = 2 FcR2aH = 3 FcR2aR = 4
fcr_id = 1;

% Read in parameters and ICs
file_to_load = Parameters(fcr_id);
load(file_to_load)

dParams = logspace(log10(0.004),log10(20),50); %makes vector of 50
    % logarithmically spaced values between 0.004X and 20X to 
    % multiply our baseline parameters by

% Run unaltered simulation
[ybase steadystate] = Simulate(params, paramnames, [],[]);
if steadystate>0
    disp(["Baseline did not run to steady state"])
else
    disp(["Baseline ran to steady state"])
end
save(strcat(datestr(today()),'_Base_Simulation_3av.mat'),'ybase','params',...
    'paramnames','complexes');

% Run full 2D sensitivity analysis
delete(gcp('nocreate')) %shutdown any parallel pools open
p1 = 17;%IgG1 conc index
p2 = 19;%IgG3 conc index
[yfull, filettl] = Sensitivity_2D(params, dParams,...
    paramnames,length(ybase(1,:)),p1,p2);

yfull = reshape(yfull, [length(dParams),length(dParams),length(yfull(1,:))]);

save(filettl)

toc