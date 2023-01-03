function [p_final,rejectsnull]  = binomial_proportion(sucess1,total1,sucess2,total2,alpha)

p1 = sucess1./total1;
p2 = sucess2./total2;

total_successes = sucess1 + sucess2;
total_trials = total1 + total2;
phat = total_successes/total_trials; 

z_score = (p1-p2)/(sqrt(phat*(1-phat)*(1/total1+1/total2)));
p_final = normcdf(z_score);
rejectsnull = p_final<alpha;

end
