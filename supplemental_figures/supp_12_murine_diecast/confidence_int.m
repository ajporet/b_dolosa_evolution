function CI = confidence_int(x)
    length_of_x = sum(~isnan(x));
    SEM = std(x,'omitnan')/sqrt(length_of_x);               % Standard Error
    ts = tinv([0.025  0.975],length_of_x-1);      % T-Score
    CI = nanmean(x) + ts*SEM;                      % Confidence Intervals
end