% Kirschner Lab, University of Michigan
% Kirschner Lab website: http://malthus.micro.med.umich.edu/lab/

% Last edit: March 11th, 2020 by
% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Modified for Ab-FcR model

function  ode = lhs_ode_igg_settings()

% This is an example of a new format ODE LHS settings file.

% To use this file for a new model do the following:
%     Copy this file and edit it to define the items as needed for the new
%     model.  That includes deleting these header comments and replacing them
%     with ones that describe the new model.
%
%     Create a Matlab file to contain the equations, which can be an edited
%     copy of file lhs_ode_predator_prey_ode.m. Be sure to define item
%     ode.odeModelHandle in this file accordingly.

% This is a settings file for a simple predator/prey model with 2 equations and
% 4 parameters.  This model is from page 6, section 3.1 of "A Methodology For
% Performing Global Uncertainty And Sensitivity Analysis In Systems Biology",
% Marino, et.  al., Journal of Theoretical Biology, 2008-09-07, doi:
% 10.1016/j.jtbi.2008.04.011.

% To run this model, use the following command from within a Matlab command
% window:
% lhs_ode_run_new('lhs_ode_settings_predator_prey_new')

% You can edit the the number of runs, time points to analyze, add or change
% ranges for the parameters, initial conditions and pulse items, as described
% in the comments below. 

% In this simple model all parameters are varied according to a normal
% distribution ('n' in the parameter items). Our production models will
% typically not vary all parameters and most tend to be varied according to a
% uniform distribution ('u').

% This settings file does not vary the initial conditions. You can experiment
% with varying initial conditions by specifying a distribution to use (ex. 'u'
% instead of '') and values that specify a range (uniform 'u' or log uniform
% 'lu' distribution) or that specify a mean and standard deviation (normal
% distribution, 'n').

% This settings file does not define pulse time points. Those are commented
% out.  You can experiment with pulse time points by uncommenting them and
% specify pulse time points and pulse values or ranges.

%===============================================================================
%
% Output:
%   The function return value is a structure named ode, with the following
%   members, which need to be define in this function.
%
%   NR
%       The number of times to run the ode model.
%
%   tspan
%       A vector of time points to produce model run results for.
%
%       This is an input to the ODE solver.
%ode.analysisTimePoints=
%   analysisTimePoints
%       A vector of time points to perform an analysis for.
%
%       Each element of this vector must also be an element of tspan.
%
%   odeModelHandle
%       The name of the Matlab function that contains the model equations.
%       This is an input to the ODE solver.
%       This is required.
%          
%   computeParamsInitCondHandle
%       The name of a Matlab function for computing some parameter values
%       and/or initial conditions based on other parameter values and/or
%       initial conditions.
%
%       This is not required. Leave undefined if not needed.
%
%  solverTimeLimit
%      The amount of time to allot for solving one set of parameters and initial
%      conditions. If the Matlab ODE solver takes longer than this amount the
%      solve is stopped and noted in an error variable.
%
%      In seconds. Can be fractional, ex. 1.5 or 0.5.
%
%  saveOnError
%      Whether or not to save the ODE solution results when an error occurs
%      during a call to the Matlab ODE solver.
%      1 = yes
%      0 = no.
%
%      If yes, pad the output matrix for the solution time steps that were not
%      defined, so the result data structures for all parameter sets will have
%      the same size and shape.Model_LHS.mat
%
%      If no, a failed solver call results in that run's results not being
%      saved in the output matrix.
%
%   parameters
%       A cell array for specifying the parameters.
%
%       There is one element of the cell array for each parameter.
%       The element for a parameter is also a cell array that has the
%       parameter's name, probability distribution and range.
%
%       The probability distribution is one of the following:
%           ''    No distribution, use the first value in the range.
%                 Both values in the range must be the same.
%
%           'u'    Use the 2 values in the range as the min and max for a
%                  uniform distribution.
%                  min must be less than max.
%
%           'lu'   Use the 2 values in the range as the min and max for a
%                  log uniform distribution.
%                  min must be less than max.
%
%           'n'    Use the 2 values in the range as the mean and std.dev.
%                  for a normal distribution.
%                   Std. dev. must not be 0.
%
%       The range always has 2 values, as described above in the probability
%       distribution.
%
%       If a parameter is a computed parameter (will be defined in the function
%       specified by compute_params_handle), then specify a probability
%       distribution of '' and a range of 0.0, 0.0. It doesn't really
%       matter, since any random value defined for it will be overwritten by
%       the computed value, but it would be confusing to specify a range if
%       that range won't be used.
%
%   initialConditions
%       A cell array for specifying the initial conditions.
%
%       This has the same format as the parameter cell array.
%
%       Initial conditions can be chosen from a range, using a probability
%       distribution, or be computed, just as for parameters.
%
%   outputLabels
%       A cell array with string labels for each model output.
%       Usually use the default of 'O1', 'O2', ...
%
%   pulseTimePoints
%       A vector of time points when to pulse the initial conditions.
%
%       This is not required. Leave undefined if not pulsing.
%
%       If defined then pulseValues and pulseHandle must also be defined.
%
%   pulseValues
%       A cell array for specifying pulse values.
%
%       There is one element of the cell array for each pulse value.
%       The element for a pulse value is also a cell array that has the
%       value's name, probability distribution and range.
%
%       The format is the same as for the parameters cell array.
%       see the comments below for a description of this format.


% Sample size.
ode.NR = 2000;

% Time span of the simulation.lhs_ode_settings_predator_prey_new
t_end=100000; % length of the simulations (in days)
ode.tspan=(0:10:t_end); % time points where the output is calculated
ode.analysisTimePoints=(0:10:t_end); % time points of interest for the uncertainty/sensitivity analysis

% The alpha value to use for PRCC.
% Undefine this (comment out) if no PRCC is to be performed.

% Doing PRCC is generally used for production runs.
%
% Not doing PRCC can be useful for testing and debugging the code which reads a
% settings file like this one and runs the model.


ode.alpha = 0.01;

% A function handle for the Matlab ODE model to analyze.
% For example: odeModelHandle = @ODEmodel;
%
% Be sure to spell this correctly, and that the file containing this
% function exists and is on the Matlab path. 
%
% Don't put it in quotes, ex. odeModelHandle = '@ODEmodel';
% is not correct.
%
ode.odeModelHandle = @lhs_ode_igg_ode;

% Optional function handle for computing parameters and/or initial conditions.
% Leave undefined if not needed.
%ode.computeParamsInitCondHandle = @lhs_ode_compute_params_predator_prey

% Redefine the ODE solver time limit, if you want a value different than the                                  
% default. In seconds. Can be fractional, ex. 1.5 or 0.5. 
ode.solverTimeLimit = '1000.0';

% Set saveOnError to define how to handle the case where the ODE solver
% fails. saveOnError = 1 means to pad the output matrix for the solution
% time steps that were not defined. saveOnError = 0 means dont' pad the
% oututp matrix.
ode.saveOnError = 0;

% Define the model parameters: {name, distribution, value1, value2}
%
% NOTE:
% Do not vary here any parameters that are based on other parameters and/or
% initial conditions.  Any variation done here will be overwritten when they
% are defined based on those parameters and/or initial conditions.
%
% No blank lines are allowed in between the cell arrays for the individual
% parameters.
%
% Comments lines must be preceded by "...".

ode.parameters = ...
{
    { 'f1', 'u', 0.04, 200 } ... %'u' = uniform, min, max
    { 'r1', 'u', 8e-7, 0.004  } ...
    { 'f2', 'u', 0.04, 200  } ...
    { 'r2', 'u', 8e-7, 0.004 } ...
    { 'f3', 'u', 0.04, 200  } ...
    { 'r3', 'u', 8e-7, 0.004 } ...
    { 'f4', 'u', 0.04, 200  } ...
    { 'r4', 'u', 8e-7, 0.004 } ...
    { 'f5', 'u', 0.08, 400  } ...
    { 'r5', 'u', 4e-5, 0.2 } ...
    { 'f6', 'u', 0.0028, 14  } ...
    { 'r6', 'u', 4e-5, 0.2 } ...
    { 'f7', 'u', 0.392, 1960  } ...
    { 'r7', 'u', 4e-5, 0.2 } ...
    { 'f8', 'u', 0.01, 20  } ...
    { 'r8', 'u', 4e-5, 0.2 } ...
    { 'g1tot', 'u', 3.4168e-07, 0.0017  } ...
    { 'g2tot', 'u', 5.3135e-10, 2.6567e-06 } ...
    { 'g3tot', 'u', 1.2280e-08, 6.1400e-05  } ...
    { 'g4tot', 'u', 4.0878e-10, 2.0439e-06 } ...
    { 'etot', 'u', 1e-7, 5e-4  } ...
    { 'ftot', 'u', 8e-8, 4e-4 } ...
};

% Define the initial conditions:  {name, distribution, value1, value2}
%
% NOTE:
% Do not vary here any initial conditions that are based on other parameters
% and/or initial conditions.  Any variation done here will be overwritten when
% they are defined based on those parameters and/or initial conditions.
%
% No blank lines are allowed in between the cell arrays for the individual
% parameters.
%
% Comments lines must be preceded by "...".

ode.initialConditions = ...
{
    { 'e1', '',  0.0, 0.0 } ...
    { 'e2', '',   0.0,  0.0 } ...
    { 'e3', '',  0.0, 0.0 } ...
    { 'e4', '',   0.0,  0.0 } ...
    { 'e11', '',  0.0, 0.0 } ...
    { 'e12', '',   0.0,  0.0 } ...
    { 'e13', '',  0.0, 0.0 } ...
    { 'e14', '',   0.0,  0.0 } ...
    { 'e22', '',  0.0, 0.0 } ...
    { 'e23', '',   0.0,  0.0 } ...
    { 'e24', '',  0.0, 0.0 } ...
    { 'e33', '',   0.0,  0.0 } ...
    { 'e34', '',  0.0, 0.0 } ...
    { 'e44', '',   0.0,  0.0 } ...
     { 'fe11', '',  0.0, 0.0 } ...
    { 'fe12', '',   0.0,  0.0 } ...
    { 'fe13', '',  0.0, 0.0 } ...
    { 'fe14', '',   0.0,  0.0 } ...
    { 'fe22', '',  0.0, 0.0 } ...
    { 'fe23', '',   0.0,  0.0 } ...
    { 'fe24', '',  0.0, 0.0 } ...
    { 'fe33', '',   0.0,  0.0 } ...
    { 'fe34', '',  0.0, 0.0 } ...
    { 'fe44', '',   0.0,  0.0 } ...
    { 'sum_FcR_complexes', '', 0.0, 0.0}...
};

% Labels to use for model outputs.
%
% There should be one for each model equation.
%
% The default labels are 'O1', 'O2', ..., 'Oneq' where neq is the number of
% model equations, which should be the same as the number of initial
% conditions.
icCount = length(ode.initialConditions);
ode.outputLabels = lhs_ode_default_output_labels_new(icCount);

% Pulse information - time points and values, each of which may or may not be
% varied.
%ode.pulseTimePoints = [];

% NOTE:
% Do not vary here any pulse values that are based on parameters and/or initial
% conditions.  Any variation done here will be overwritten when they are
% defined based on those parameters and/or initial conditions.
%{
ode.pulseValues = ...
{
   ... % Pulse Q with value 0.5 for each pulse time point of each run.
   % { 'Q', '', 0.5, 0.5 }, ... % For pulse time point 1.
   % { 'Q', '', 0.5, 0.5 }, ... % For pulse time point 2.
   ...
   ... % Pulse Q in the range 0.5 to 2.5 for each pulse time
   ... % point of each run.
   % { 'Q', 'u', 0.5, 2.5 }, ... % For pulse time point 1.
   % { 'Q', 'u', 0.5, 2.5 }, ... % For pulse time point 2.
   ...
   ... % Pulse Q with its initial value for each pulse
   ... % time point of each run.
   % { 'Q' }, ... % For pulse time point 1.
   % { 'Q' }, ... % For pulse time point 2.
};
%}
end % function ode = lhs_ode_settings_predator_prey_new
