% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: March 11th, 2020

% Wrapper code to input each individuals IgG subclass concentration and FcR
% concentration: The outputs to save in a mat file are the following:
%   ybase params patient id paramnames complex names

% Read the data from a file

function Run_Me_personal_simulations()
clear();
tic()
%% choose which strain to simulate, A244 or BAL
strain = 'A244';
% strain = 'BAL';

file_name = ['RV144_',strain,'_Personal_Input_v1-v105.xlsx'];

% Read in personal data
indiv_data = readtable(file_name);

%Descibes column index of each desired parameter
IgG_FcR_data = indiv_data{:,2:end};

patient_id = indiv_data{:,1};

IgG_FcR_data(isnan(IgG_FcR_data))=0;%turn nans to zero in data
g1_idx = 1;
g2_idx = 2;
g3_idx = 3;
g4_idx = 4;

complexname=['env-IgG1', 'env-IgG2', 'env-IgG3', 'env-IgG4',...
    'env-IgG1-IgG1', 'env-IgG1-IgG2','env-IgG1-IgG3', 'env-IgG1-IgG4',...
    'env-IgG2-IgG2', 'env-IgG2-IgG3','env-IgG2-IgG4','env-IgG3-IgG3',...
    'env-IgG3-IgG4','env-IgG4-IgG4','FcR-env-IgG1-IgG1',...
    'FcR-env-IgG1-IgG2','FcR-env-IgG1-IgG3', 'FcR-env-IgG1-IgG4',...
    'FcR-env-IgG2-IgG2', 'FcR-env-IgG2-IgG3','FcR-env-IgG2-IgG4',...
    'FcR-env-IgG3-IgG3','FcR-env-IgG3-IgG4','FcR-env-IgG4-IgG4',...
    "IgG1", "IgG2", "IgG3", "IgG4", "env", "FcR",...
    "FcR complexes with IgG1 and IgG3 only", ...
    "FcR complexes including IgG1 and IgG3","All FcR complexes"];

% Different FcR affinity values(mM-1s-1), rows are IgG1-IgG4 Fc affinity,
% columns are the FcR (FcR3aV FcR3aF FcR2aH FcR2aR  respectively)
FcR_kon = [20	11.7	52  35; % IgG1-Fc kon
    0.7	0.3	4.5     1; % IgG2-Fc kon
    98	77	8.9     9.1; % IgG3-Fc kon
    2.5	2	1.7     2.1]; % IgG4-Fc kon
FcR_names = ["FcR3a-V158" "FcR3a-F158" "FcR2a-H131" "FcR2a-R131"];

% Run simulations
clear all_run
clear param_idv
fcr_id = 0;
for fcr_idx = 1:length(FcR_kon(1,:)) % columns in the spreadsheet
    fcr_id = fcr_id + 1;
    for p_num = 1:length(patient_id)
        
        % get personal paramters for the given FcR
        [params,paramnames,~] = ...
            Parameters_indiv_FcR(FcR_kon(1,fcr_id), FcR_kon(2, fcr_id), ...
            FcR_kon(3, fcr_id), FcR_kon(4, fcr_id), ...
            IgG_FcR_data(p_num,g1_idx), IgG_FcR_data(p_num,g2_idx), ...
            IgG_FcR_data(p_num,g3_idx), IgG_FcR_data(p_num,g4_idx));
        
        % run the simulation
        [ybase, steadystate] = Simulate(params,paramnames,[],[]);
        
        % display if steady state was reached or not (defined in
        % "simulate.m"
        if steadystate>0
            disp(["Patient "+p_num+" did not reach steady state"])
        else
            disp(["Patient "+p_num+" at steady state"])
        end
        % save data in the larger matrix
        all_run(fcr_id, p_num, :) = ybase;
        param_idv(fcr_id, p_num,:) = params;
    end
end

save([strain,'_personal_baseline_all_fcrs_IgG_v1-v105_',datestr(today()),'.mat'],...
    'all_run', 'param_idv','complexname','paramnames', 'IgG_FcR_data',...
    'patient_id','FcR_names');
toc()
end
