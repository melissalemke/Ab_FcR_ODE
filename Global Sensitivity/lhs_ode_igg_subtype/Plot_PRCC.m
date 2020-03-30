% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: March 11th, 2020

clear;
try
load('Model_LHS_n_50.mat')

paramnames = ["k_{on} IgG1-env","k_{off} IgG1-env","k_{on} IgG2-env","k_{off} IgG2-env",...
     "k_{on} IgG3-env","k_{off} IgG3-env","k_{on} IgG4-env",...
     "k_{off} IgG4-env","k_{on} IgG1-FcR","k_{off} IgG1-FcR","k_{on} IgG2-FcR",...
     "k_{off} IgG2-FcR","k_{on} IgG3-FcR","k_{off} IgG3-FcR","k_{on} IgG4-FcR",...
     "k_{off} IgG4-FcR","IgG1 conc","IgG2 conc", "IgG3 conc",...
     "IgG4 conc", "env conc", "FcR conc"];
N = 25; %25 is the complex formation output
 
figure()
superbar(totalPrcc(end,:,N),'P',totalSign(end,:,N))
xticklabels(paramnames)
ylabel('PRCC')
xticks(1:length(PRCC_var))
xtickangle(45)
title(['Time = ',num2str(odeSettings.tspan(end)),' NS = ',num2str(length(LHSmatrix(:,1)))])
set(gca,'Fontsize',16)
catch
    disp('Make sure to replace the load(XX) function inline 3 with your data file name')
end

