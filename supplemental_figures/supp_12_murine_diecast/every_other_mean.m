function summed_every_other_mean =  every_other_mean(ex_array)

summed_every_other_mean = [mean(ex_array(:,1:2:end),2), mean(ex_array(:,2:2:end),2)];

end