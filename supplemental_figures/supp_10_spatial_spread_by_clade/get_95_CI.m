function [CI] = get_95_CI(datapoints)
SEM = std(datapoints)/sqrt(length(datapoints));               % Standard Error
ts = tinv([0.025  0.975],length(datapoints)-1);      % T-Score
CI = mean(datapoints) + ts*SEM;                      % Confidence Intervals
end