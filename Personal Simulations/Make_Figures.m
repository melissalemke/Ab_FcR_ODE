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

% %This is the validation data in terms of MFIs
input_file_name = 'RV144_A244_v1-30_MFI_2020_04_06.xlsx';
raw_val = readtable(input_file_name);

patient_idval = table2array(raw_val(:,1));

%FcR3aV, FcR3aF
FcR_val = table2array(raw_val(:,6:8));

remove = isnan(FcR_val(:,1));

FcR_val = FcR_val(~remove,:);
patient_idval = string(patient_idval(~remove));

if sum(patient_id(1:30)~=['v'+patient_idval(1:30)])
    disp("Patients not aligned!")
    pause()
end

%making figures
all_com_formation_data(FcR_model,FcR_names,patient_id,strain)

dozscore = false;
validationcor(FcR_model(1:30,:),FcR_val(1:30,:),1,com_name,...
    patient_id(1:30),FcR_names,strain,dozscore)

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
   ylabel(["Complex"; "Formation (nM)"])
end
end

function [] = validationcor(FcR_model,FcR_val,fcr_id,com_name,patient_id,FcR_names,strain,dozscore)
figure()
x = FcR_model(:,fcr_id);
y = FcR_val(:,fcr_id);
ttl="";

if dozscore
    x = zscore(x);
y = zscore(y);
ttl = "zscored"
end

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
% Pearson correlation alpha = 0.05 - > final stats done in GraphPrism
test = 'Pearson';
[corrl pval] = corr(x,y,'type',test);
str = num2str(round(corrl,3));
if pval <0.05
    corleg = [corleg;[test+" R = "+str+"*"]];
else
    corleg =[corleg;[test+" R = "+str+"ns"]];
end
text(min(x)*1.1,max(y)*.9,corleg)

set(gca,'FontSize',12)
xtickangle(45)
title([strain+" "+FcR_names(fcr_id)])
xlabel(["Predicted Complex";["Formation ("+ttl+" nM)"]])
ylabel(["Measured Complex";["Formation ("+ttl+" MFI)"]])
end