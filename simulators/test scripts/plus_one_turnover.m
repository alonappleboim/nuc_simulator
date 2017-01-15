function [ avg_turnover_time ] = plus_one_turnover( s_hist, time, plus_one_pos, slide_len)
%plus_one_turnover a function for getting the average turnover time of the
%plus one nucleosome.
%   The function goes over the given state history and finds when the +1
%   nuc has been evicted and assembled, finally calculating how many
%   evictions were in total and the average turnover rate.

index = 1;
assem_flag = 0;
num_of_turnovers = 0;

while (index <= length(s_hist(:,1)))
    
    nuc_flag = sum(s_hist(index,plus_one_pos));
    
    if (nuc_flag > 0)
        if (assem_flag == 0)
            assem_flag = 1;
            num_of_turnovers = num_of_turnovers + 1;
        end
    end
    if (nuc_flag == 0)
        if (assem_flag == 1)
            assem_flag = 0;
        end           
    end
    
    index = index + 1;
end

avg_turnover_time = time(end) / num_of_turnovers;

end

