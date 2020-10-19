% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: March 11th, 2020

%% Individual Sensitivity further analysis and figure
clear
%choose the strain
strain = 'A244';
% strain = 'BAL';

%load data
d = dir([strain,'_personal_local_sensitivity_all_fcrs_IgG_v1-v105*']);
file_name = d.name;
load(file_name);

%choose complex of interest
com_id=33;
%choose FcR of interest
fcr_id = 1;
%display choice in command window
disp([complexname(com_id)+" "+FcR_names(fcr_id)])

%get appropriate data
patient_sense = squeeze(all_ysense(fcr_id,:,:,com_id)); 

%Reshape
patient_sense_3D = reshape(patient_sense',[length(paramnames) length(dParams) length(patient_id)]);

%Add baseline (1X) value
patient_sense_3D(:,5:7,:) = patient_sense_3D(:,4:6,:);
patient_sense_3D(:,4,:) = repmat(all_run(fcr_id,:,com_id), [length(paramnames),1]);
dParams(5:7) = dParams(4:6);
dParams(4) = 1;


%Sensitivity metric calculation
sense_metric = (squeeze(max(patient_sense_3D,[],2))-squeeze(min(patient_sense_3D,[],2)))./(max(dParams)-min(dParams))';

% Heatmap
figure()
imagesc(sense_metric)
xlabel("Vaccinee IDs")
set(gca,'ytick',[1:22],'yticklabel',paramnames,'xtick',[1:105],'xticklabel',patient_id)
xtickangle(270)
title([strain+" "+FcR_names(fcr_id)+" Personal Local Sensitivity"])
c = colorbar;
ylabel(c, "Sensitivity Metric")
clim = caxis();


%% Max complex formation level acheived between 0.004X-20X for each parameter
com_form = squeeze(max(patient_sense_3D,[],2))';
labels = paramnames;
med_thresh = 8e-1; %in nM
high_thresh = 1.0646; %in nM

medium = sum(com_form>=med_thresh  & com_form<high_thresh);
high = sum(com_form>=high_thresh);

figure()
y = [high; medium; 105-medium-high]';
[highmed ind] = sort(y(:,1)+y(:,2));
b=bar(y(ind,:),'stacked');
b(3).FaceColor = 'w';
b(1).FaceColor = [0 158 115]/255;
b(2).FaceColor = [204 121 167]/255;
set(gca, 'xtick',1:length(paramnames))
xticklabels([labels(ind)])
xtickangle(45);
ylabel("Number of vaccinees")
legend(["High complex formation" "Medium complex formation" "Low complex formation" ])
title([strain+" "+FcR_names(fcr_id)])
