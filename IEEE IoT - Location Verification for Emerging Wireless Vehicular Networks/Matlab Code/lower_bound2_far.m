function [toa_Array_NN] = lower_bound2_far(standard_deviation, NLoS)

[toa_Array, std_t, var_t] = lower_bound_far(standard_deviation);

%%% Below for 4 verifiers
U_Good = toa_Array(:,1:4);
% Y_Good = toa_Array(:,5:8);
Y_Good_old = toa_Array(:,5:8);
U_Bad = toa_Array(:,9:12);
V_Bad = toa_Array(:,13:16);
% Y_Bad = toa_Array(:,17:20);
Y_Bad_old = toa_Array(:,17:20);

% % Below for 6 verifiers
% U_Good = toa_Array(:,1:6); % % Y_Good = toa_Array(:,1:6); % Y_Good_old = toa_Array(:,7:12);
% U_Bad = toa_Array(:,13:18); % V_Bad = toa_Array(:,19:24); % % Y_Bad = toa_Array(:,25:30); % Y_Bad_old = toa_Array(:,25:30);

ii = size(U_Good,1);
% % [bias,std_i,var_i] = nonLOS(ii, mu, number_of_verifiers);
[bias,~,~] = nonLOS(ii, NLoS, 4);
Y_Good = Y_Good_old + bias; % NLoS biases to Good User.
Y_Bad  = Y_Bad_old; % Far field approximation.

tag_good_users = zeros(ii,1);
tag_bad_users  =  ones(ii,1);

% toa_Array_NN  = [U_Good,Y_Good,tag_good_users; U_Bad,Y_Bad,tag_bad_users];
toa_Array_NN  = [U_Good,Y_Good,tag_good_users, U_Bad,Y_Bad,tag_bad_users,V_Bad,var_t];

end