function summed_every_other_fun =  every_other_sum(ex_array)

summed_every_other_fun = [sum(ex_array(:,1:2:end),2), sum(ex_array(:,2:2:end),2)];

end