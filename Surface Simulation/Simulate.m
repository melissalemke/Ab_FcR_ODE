% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: March 11th, 2020

%% Run single simulation
function [yend, steadystate]  = Simulate(params, paramnames, altered,alteration) 
%% Run Normal Simulation
global tStep tArray varArray;

% Initialize Global Variables
tStep = 1; 
tArray = zeros(1e6,1); 
varArray = zeros(1e6,6);  

% Run Simulation
y0 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %Initial complex concentrations
tspan = [0 100000];
options = odeset('AbsTol',1e-50,'RelTol',1e-10);

[t,y] = ode113(@ODEs,tspan,y0,options,params);

% Post-processing 
tArray = tArray(1:min(tStep,length(y(:,1))));   % truncate unused variables
varArray = varArray(1:min(tStep,length(y(:,1))),:);

ycut = y(1:min(tStep,length(y(:,1))),:);

    % free species concentrations over time
    g1y=varArray(:,1);
    g2y=varArray(:,2);
    g3y=varArray(:,3);
    g4y=varArray(:,4);
    ey=varArray(:,5);
    fy=varArray(:,6);
    
    % Adding free species concentrations and complex summations to the
    % results matrix
    ycut(:,25) = g1y;
    ycut(:,26) = g2y;
    ycut(:,27) = g3y;
    ycut(:,28) = g4y;
    ycut(:,29) = ey;
    ycut(:,30) = fy;
    ycut(:,31) = ycut(:,15)+ycut(:,17)+ycut(:,22);
    ycut(:,32) = ycut(:,15)+ycut(:,16)+ycut(:,17)+ycut(:,18)+ycut(:,20)+ycut(:,22)+ycut(:,23);
    ycut(:,33) = sum(ycut(:,15:24),2); 
    
    % Concentrations at the last time point
    yend = ycut(end,:);
    
    % Steady state definition (subject to change with different baseline
    % parameters)
    steadystate = sum(abs(yend-ycut(end-floor(0.025*length(ycut(:,1))),:)) > 1e-2*yend);
end
