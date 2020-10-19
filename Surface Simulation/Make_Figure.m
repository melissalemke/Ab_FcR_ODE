% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: March 11th, 2020

clear
%choose the output complex 
com_id = 33;

%Choose which FcR: FcR3aV = 1 FcR3aF = 2 FcR2aH = 3 FcR2aR = 4
fcr_id = 1;

%choose the strain
strain = 'A244';
% strain = 'BAL';

%load the surface data
d = dir(['*2D_IgG1_conc*']);
file_name = d.name;
load(file_name);

disp([paramnames(p1)+" "+paramnames(p2)])

x = dParams*params(p1);
y = dParams*params(p2);
z = yfull(:,:,com_id);



%% Load the personal baseline data!!! Move it to this folder
try % will only run if you've moved the personal baseline data into this folder
d = dir([strain,'_personal_baseline_all_fcrs_IgG_v1-v105*']);
file_name = d.name;
load(file_name);
com_name = complexname(com_id);

IgG_FcR_data = IgG_FcR_data*10^6;%mM -> nM

% Prep data
FcR_model = squeeze(all_run(:,:,com_id))';

%plot the surface
surfacefig(p1,p2,FcR_model,fcr_id,patient_id, com_name,...
FcR_names,param_idv,x,y,z,paramnames)
catch
    disp('Be sure to copy the "strain"_personal_baseline_all_fcrs...mat file generated in the "Personal Simulations" folder into this folder!')
end

log_x = log10(x);
log_y = log10(y);
z = z;

% Plot and calculate gradient on non-log scale for z-axis
[fx, fy] = gradient(z, abs(log_x(2) - log_x(1)), abs(log_y(2) - log_y(1)));

%% Grab gradient value based on personalized IgG3 and IgG1 concentrations

Vx = log10(IgG_FcR_data(:,1));
Vy = log10(IgG_FcR_data(:,3));
Nx = log_x;
Ny = log_y;

clear patient_gradient closeX closeY idx_xy
for i = 1:length(Vx)
    [minValueX,closestIndexX] = min(abs(Nx-Vx(i)));
    [minValueX,closestIndexY] = min(abs(Ny-Vy(i)));
    closeX(i) = Nx(closestIndexX);
    closeY(i) = Ny(closestIndexY);
    idx_xy(i,:) = [closestIndexX closestIndexY];
    patient_gradient(i,:) = [fx(idx_xy(i,2),idx_xy(i,1)) fy(idx_xy(i,2),idx_xy(i,1))];
end

%% Gradient Dot Plots
figure
r1 = (0.5 - rand(length(patient_gradient(:,1)),1))/8;
 plot(1*ones(length(patient_gradient(:,1))) + r1, patient_gradient(:,1),'o', ...
     'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k')
 hold on
 plot(1.2*ones(length(patient_gradient(:,1))) + r1, patient_gradient(:,2),'o', ...
     'MarkerFaceColor', 'r', 'MarkerEdgeColor','k')
 set(gca, 'xtick', [1, 1.2], 'xticklabels', {'IgG1', 'IgG3'})
 ylabel('Gradient')
 set(gca, 'fontsize', 14)


%% Surface fig
function [] = surfacefig(p1,p2,FcR_model,fcr_id,patient_id, com_name, FcR_names,param_idv,x,y,z,paramnames)
figure()
%plot personal baseline

rest_ind = 1:105;
%to mark responders and non-resdoners uncomment lines 92,93,98-102
% imp = [27,12;24,3;18,13;15,29;9,22;30,1;26,16;19,4]; %responders, non-responders Bal
% rest_ind = ~ismember(1:105, imp);

        rest = scatter3(param_idv(fcr_id,rest_ind,p1),param_idv(fcr_id,rest_ind,p2),...
            FcR_model(rest_ind,fcr_id)+1e-7,40,'filled','MarkerEdgeColor','k','MarkerFaceColor','none','linewidth',2);%[128 128 128]/255);
        hold on
%         scatter3(param_idv(fcr_id,imp(:,1),p1),param_idv(fcr_id,imp(:,1),p2),...
%             FcR_model(imp(:,1),fcr_id),200,'filled','MarkerEdgeColor','k','MarkerFaceColor',[0 114 178]/255,'linewidth',1);
%         hold on
%         scatter3(param_idv(fcr_id,imp(:,2),p1),param_idv(fcr_id,imp(:,2),p2),...
%             FcR_model(imp(:,2),fcr_id),200,'filled','MarkerEdgeColor','k','MarkerFaceColor',[230 159 0]/255,'linewidth',1);

%To mark patients in text
% text(param_idv(fcr_id,:,p1),param_idv(fcr_id,:,p2),...
%     FcR_model(:,fcr_id),patient_id)

%custom color map for low, med, high complex formation approx thresholds
% low_med_high_cmap = [1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;0.819607853889465,0.392156869173050,0.650980412960053;0.819607853889465,0.392156869173050,0.650980412960053;0.819607853889465,0.392156869173050,0.650980412960053;0.819607853889465,0.392156869173050,0.650980412960053;0.819607853889465,0.392156869173050,0.650980412960053;0.819607853889465,0.392156869173050,0.650980412960053;0.819607853889465,0.392156869173050,0.650980412960053;0.819607853889465,0.392156869173050,0.650980412960053;0.819607853889465,0.392156869173050,0.650980412960053;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952;0.321568638086319,0.701960802078247,0.447058826684952];
% colormap(low_med_high_cmap)

%plot surface - > currently everyother line to prevent crowding
surface(x(1:2:end),y(1:2:end),z(1:2:end,1:2:end),'LineStyle','-','FaceColor','none')%
ylabel([paramnames(p2)+" (nM)"])
xlabel([paramnames(p1)+" (nM)"])
xlim([min(x)*1 max(x)*1])
ylim([min(y)*1 max(y)*1])
yticks([1e-2 1e-1 1 1e1 1e2])
xticks([1e0 1e1 1e2 1e3])
zlabel(["Complex", "Formation (nM)"])
set(gca, 'YScale', 'log','XScale','log')%,'Zscale','log'
set(gca,'FontSize',8)
end