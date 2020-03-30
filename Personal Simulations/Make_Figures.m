% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: March 11th, 2020

clear
%choose the strain
strain = 'A244';
% strain = 'BAL';

d = dir([strain,'_personal_baseline_all_fcrs_IgG_v1-v105*']);
file_name = d.name;
load(file_name);

% choosing the complexes we want as "outputs" in these figures
complex = [33]; %this is the sum of all fcr complexes
com_name = complexname(complex);
disp(com_name) %displaying the chosen complex

FcR_model = squeeze(all_run(:,:,complex))';

%This is the validation data saved into the spreadsheet containing personal
%IgG concentrations (only available for FcR3a-V158 FcR3a-F158 FcR2a-H131)
FcR_val = IgG_FcR_data(:,5:7);

all_com_formation_data(FcR_model,FcR_names,patient_id,strain)

validationcor(FcR_model(1:30,:),FcR_val(1:30,:),1,com_name,...
    patient_id(1:30),["Predicted Complex";"Formation (zscored mM)"],...
    FcR_names,strain)

load("IgG1_v_IgG3_surface_data_fc_all_wo_IgGn.mat")
surfacefig(17,19,FcR_model,FcR_val,3,patient_id, com_name,...
FcR_names,param_idv,x,y,z,paramnames)

%% Histogram to show all complex formation v1-v105 all FcRS
function [] = all_com_formation_data(FcR_model,FcR_names,patient_id,strain)
figure()
for i = 1:length(FcR_names)
   subplot(length(FcR_names),1,i)
   bar(FcR_model(:,i))
   set(gca,'xtick',1:105,'xticklabels',patient_id)
   xtickangle(90)
   title([strain+" "+FcR_names(i)])
   xlabel("Vaccinee ID")
   ylabel(["Complex"; "Formation (mM)"])
end
end

function [] = validationcor(FcR_model,FcR_val,fcr_id,com_name,patient_id,xlab,FcR_names,strain)
figure()
x = zscore(FcR_model(:,fcr_id));
y = zscore(FcR_val(:,fcr_id));

names = patient_id;

scatter(x,y);
text(x*1.1,y,names)

% Spearman correlation alpha = 0.05 - > final stats done in GraphPrism
test = 'Spearman';
[corrl pval] = corr(x,y,'type',test);
str = num2str(round(corrl,3));
if pval <0.05
    corleg = [test+" R = "+str+"*"];
else
    corleg =[test+" R = "+str+"ns"];
end
text(min(x)*0.5,max(y)*.9,corleg)

set(gca,'FontSize',12)
xtickangle(45)
title([strain+" "+FcR_names(fcr_id)])
xlabel(xlab)
ylabel(["Measured Complex";"Formation (zscored mM)"])
end