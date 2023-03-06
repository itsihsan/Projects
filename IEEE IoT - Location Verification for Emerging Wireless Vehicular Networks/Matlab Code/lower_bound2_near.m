function [toa_Array_NN] = lower_bound2_near(r, standard_deviation, NLoS)

[toa_Array, std_t, var_t] = lower_bound_near(r , standard_deviation);

%%% Below for 4 verifiers
U_Good = toa_Array(:,1:4);
% Y_Good = toa_Array(:,5:8);
Y_Good_old = toa_Array(:,5:8);
U_Bad = toa_Array(:,9:12);
V_Bad = toa_Array(:,13:16);
% Y_Bad = toa_Array(:,17:20);
Y_Bad_old = toa_Array(:,17:20);

ii = size(U_Good,1);
% % [bias,std_i,var_i] = nonLOS(ii, mu, number_of_verifiers);
[bias,~,~] = nonLOS(ii, NLoS, 4);
Y_Good = Y_Good_old + bias; % NLoS biases to Good User only.
Y_Bad  = Y_Bad_old  + bias; % NLoS biases to Bad User only.

tag_good_users = zeros(ii,1);
tag_bad_users  =  ones(ii,1);

% toa_Array_NN  = [U_Good,Y_Good,tag_good_users; U_Bad,Y_Bad,tag_bad_users];
toa_Array_NN  = [U_Good,Y_Good,tag_good_users, U_Bad,Y_Bad,tag_bad_users,V_Bad,var_t];

end