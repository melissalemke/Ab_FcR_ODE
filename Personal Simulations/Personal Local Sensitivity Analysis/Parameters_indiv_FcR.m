% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: March 11th, 2020

function [params, paramnames, complexes] = Parameters_indiv_FcR(g1_F, ...
    g2_F, g3_F, g4_F, g1tot_ind, g2tot_ind, ...
    g3tot_ind, g4tot_ind)

%% define parameters
colors = [[66, 134, 244]/255; [219, 48, 48]/255; [137, 196, 74]/255; [178, 91, 175]/255; [99, 210, 216]/255; [255, 248, 58]/255; [0 0 0]; [1 1 1]];

complexes = {'env-IgG1', 'env-IgG2', 'env-IgG3', 'env-IgG4','env-IgG1-IgG1', 'env-IgG1-IgG2','env-IgG1-IgG3', 'env-IgG1-IgG4', 'env-IgG2-IgG2', 'env-IgG2-IgG3'...
    'env-IgG2-IgG4','env-IgG3-IgG3','env-IgG3-IgG4','env-IgG4-IgG4','FcR-env-IgG1-IgG1', 'FcR-env-IgG1-IgG2','FcR-env-IgG1-IgG3', 'FcR-env-IgG1-IgG4', 'FcR-env-IgG2-IgG2', 'FcR-env-IgG2-IgG3'...
    'FcR-env-IgG2-IgG4','FcR-env-IgG3-IgG3','FcR-env-IgG3-IgG4','FcR-env-IgG4-IgG4'};

f1 = 10*10^-6;%env-IgG1 from IgG vs IgA model
r1 = 2e-4;
f2 = 10*10^-6;%env-IgG2
r2 = 2e-4;
f3 = 10*10^-6;%env-IgG3
r3 = 2e-4;
f4 = 10*10^-6;%env-IgG4
r4 = 2e-4;
f5 = g1_F*10^-6;%IgG1-FcR (nM-1)
r5 = 1e-2;
f6 = g2_F*10^-6;%IgG2-FcR
r6 = 1e-2;
f7 = g3_F*10^-6;%IgG3-FcR
r7 = 1e-2;
f8 = g4_F*10^-6;%IgG4-FcR
r8 = 1e-2;

g1tot    = g1tot_ind*10^6;%nM 
g2tot    = g2tot_ind*10^6;%nM 
g3tot    = g3tot_ind*10^6;%nM 
g4tot    = g4tot_ind*10^6;%nM 
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
end