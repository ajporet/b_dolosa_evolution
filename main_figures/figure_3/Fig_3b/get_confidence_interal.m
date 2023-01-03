function [mean_val, confidence_interval] = get_confidence_interal(datapoints, p_cutoff)

SEM = std(datapoints)/sqrt(max(size(datapoints)));              % Standard Error
ts = tinv([p_cutoff  1-p_cutoff],max(size(datapoints))-1);      % T-Score
CI = mean(datapoints) + ts*SEM;                                 % Confidence Intervals
mean_val = mean(datapoints);
confidence_interval = [CI];
end