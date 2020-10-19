% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: March 11th, 2020

function file_to_load = Parameters(fcr_id)
%% define parameters
complexes = {'env-IgG1', 'env-IgG2', 'env-IgG3', 'env-IgG4','env-IgG1-IgG1', 'env-IgG1-IgG2','env-IgG1-IgG3', 'env-IgG1-IgG4', 'env-IgG2-IgG2', 'env-IgG2-IgG3'...
    'env-IgG2-IgG4','env-IgG3-IgG3','env-IgG3-IgG4','env-IgG4-IgG4','FcR-env-IgG1-IgG1', 'FcR-env-IgG1-IgG2','FcR-env-IgG1-IgG3', 'FcR-env-IgG1-IgG4', 'FcR-env-IgG2-IgG2', 'FcR-env-IgG2-IgG3'...
    'FcR-env-IgG2-IgG4','FcR-env-IgG3-IgG3','FcR-env-IgG3-IgG4','FcR-env-IgG4-IgG4'};

% Different FcR affinity values(nM-1s-1), rows are IgG1-IgG4 Fc affinity, columns
% are the FcR (FcR3aV FcR3aF FcR2aH FcR2aR respectively)
all_FcR_kon = [20	11.7	52  35; % IgG1-Fc kon
    0.7	0.3	4.5     1 ; % IgG2-Fc kon
    98	77	8.9     9.1; % IgG3-Fc kon
    2.5	2	1.7     2.1]*10^-6; % IgG4-Fc kon
fcr_names = ["FcR3a-V158" "FcR3a-F158" "FcR2a-H131" "FcR2a-R131"];

fcr_kons = all_FcR_kon(:,fcr_id);
fcr_ttl = fcr_names(fcr_id);

f1 = 10*10^-6;%env-IgG1 from IgG vs IgA model
r1 = 2e-4;%reduced OM by 2
f2 = 10*10^-6;%env-IgG2
r2 = 2e-4;
f3 = 10*10^-6;%env-IgG3
r3 = 2e-4;
f4 = 10*10^-6;%env-IgG4
r4 = 2e-4;
f5 = fcr_kons(1);%IgG1-FcR (mM-1)
r5 = 1e-2;
f6 = fcr_kons(2);%IgG2-FcR
r6 = 1e-2;
f7 = fcr_kons(3);%IgG3-FcR
r7 = 1e-2;
f8 = fcr_kons(4);%IgG4-FcR
r8 = 1e-2;

g1tot    = (8.54195E-05)*10^6;%nM These are the averages from all patients
g2tot    = (1.32837E-07)*10^6;%nM 
g3tot    = (3.06998E-06)*10^6;%nM 
g4tot    = (1.02194E-07)*10^6;%nM 
etot    = 25;%nM 
ftot    = 20;%nM 

params = [f1,r1,f2,r2,f3,r3,f4,r4,f5,r5,f6,r6,f7,r7,f8,r8,...
    g1tot,g2tot, g3tot,g4tot, etot, ftot];
 
paramnames = ["k_{on} IgG1-env","k_{off} IgG1-env","k_{on} IgG2-env","k_{off} IgG2-env",...
     "k_{on} IgG3-env","k_{off} IgG3-env","k_{on} IgG4-env",...
     "k_{off} IgG4-env","k_{on} IgG1-FcR","k_{off} IgG1-FcR","k_{on} IgG2-FcR",...
     "k_{off} IgG2-FcR","k_{on} IgG3-FcR","k_{off} IgG3-FcR","k_{on} IgG4-FcR",...
     "k_{off} IgG4-FcR","IgG1 conc","IgG2 conc", "IgG3 conc",...
     "IgG4 conc", "env conc", "FcR conc"];
 
file_to_load = 'params_'+fcr_ttl+'.mat';
 save(file_to_load,'params', 'paramnames', 'complexes','fcr_ttl');
end