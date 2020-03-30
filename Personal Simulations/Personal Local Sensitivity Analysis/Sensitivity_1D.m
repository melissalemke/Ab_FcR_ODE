% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: March 11th, 2020

function [YOutput, filettl] = Sensitivity_1D(params, dParams, paramnames,ycol)
% If not inside a PBS job, use 12 processors


% Create Output Matrices
    [YOutput] = iterate(1:length(params),dParams,params,paramnames,ycol);
    clear Outputs;
    filettl = strcat(datestr(today()),'_','1D_Sensitivity.mat');

%% Iterate Through Parameter Changes
    function [YOutput] = iterate(jloop,dParams,params,paramnames,ycol)
        YOutput = zeros(length(jloop)*length(dParams),ycol);
        lendParams = length(dParams);
        lenjloop = length(jloop);
        
        parfor r = 1:lendParams*lenjloop
            i = ceil(r/lenjloop);
            if rem(r,lenjloop)
                j=rem(r,lenjloop);
            else
                j=lenjloop;
            end
            newP = params;
            newP(jloop(j)) = dParams(i)*params(jloop(j)); %#ok<PFBNS>
            [yout steadystate] = Simulate(newP,paramnames, [],[]);
            YOutput(r,:) = yout;
            if steadystate>0
              disp([paramnames(jloop(j))+" at "+num2str(dParams(i))+"X did not run to steady state"])
            else
              disp([paramnames(jloop(j))+" at "+num2str(dParams(i))+"X ran to steady state"])  
        end
        end
    end
end