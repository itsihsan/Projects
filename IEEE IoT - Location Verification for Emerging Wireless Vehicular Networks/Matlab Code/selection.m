function [ToA, U_Good, Y_Good, Y_Bad,V_Bad, var_t] = selection(r, standard_deviation, NLoS)

% standard_deviation for thermal noise.
% NLOS is the standard deviation for NLoS
% r is the minimum distant constraint

toa_Array_NN = [];
a=1/1000;
num_1 = exp(-a.*r); % r goes from 1->inf, num_1 goes from 1->0

for m=1:5000
    num_2 = rand(1);
    if num_1 >= num_2
        [toa_Array_NN] = lower_bound2_near(r, standard_deviation, NLoS);
    else
        [toa_Array_NN] = lower_bound2_far(standard_deviation, NLoS);
    end
    ToA_pre(m,:)=toa_Array_NN;
    toa_Array_NN = [];
    num_2 = [];
end

ToA_good = ToA_pre(:,1:9);
ToA_bad  = ToA_pre(:,10:18);
ToA = [ToA_good;ToA_bad];
U_Good = ToA_pre(:,1:4);
Y_Good = ToA_pre(:,5:8);
Y_Bad = ToA_pre(:,14:17);
V_Bad = ToA_pre(:,19:22);
var_t = unique(ToA_pre(:,23));
end