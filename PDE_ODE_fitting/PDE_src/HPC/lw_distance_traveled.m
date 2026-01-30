function [total_distance_traveled temp1_max] = lw_distance_traveled(la, N, R, temp1_max, total_distance_traveled)
% Reed Ogrosky
% University of North Carolina at Chapel Hill
% 10-12-12

% This function keeps track of how far the interface has traveled by
% following a local maxima

L=2*pi*la;
dx=L/N;

if temp1_max>5 && temp1_max<N-4
    old_temp1_max = temp1_max;
    [m1 temp1_max2] = min(R(temp1_max-5:temp1_max+5));
    temp1_max = temp1_max2-1+temp1_max-5;
    total_distance_traveled = total_distance_traveled+abs(temp1_max-old_temp1_max)*dx;
end
if temp1_max<=5
    old_temp1_max = temp1_max;
    index1_list=[1:temp1_max+5];
    index1_list=[index1_list N-5+temp1_max:N];
    [m1 temp1_max2] = min(R(index1_list));
    if temp1_max2<=temp1_max+5
        temp1_max = temp1_max2;
    else
        temp1_max = temp1_max2+N-11;
    end
    if temp1_max<10
        total_distance_traveled = total_distance_traveled+abs(temp1_max-old_temp1_max)*dx;
    else            
        total_distance_traveled = total_distance_traveled+abs(temp1_max-old_temp1_max-N)*dx;
    end
end
if temp1_max>=N-4
    old_temp1_max = temp1_max;
    index1_list=[temp1_max-5:N];
    index1_list=[1:temp1_max-N+5 index1_list];
    [m1 temp1_max2] = min(R(index1_list));
    if temp1_max2<=temp1_max-N+5
        temp1_max = temp1_max2;
    else
        temp1_max = temp1_max2+N-11;
    end
    if temp1_max<10
        total_distance_traveled = total_distance_traveled+abs(temp1_max-old_temp1_max+N)*dx;
    else            
        total_distance_traveled = total_distance_traveled+abs(temp1_max-old_temp1_max)*dx;
    end                
end


end

