function [mean_val, SEM] = get_SEM(datapoints)
SEM = std(datapoints)/sqrt(max(size(datapoints)));              % Standard Error                                % Confidence Intervals
mean_val = mean(datapoints);
end