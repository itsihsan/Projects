function [toa_Array, toa_Array_NN, bias, std_t, var_t, std_i, var_i, U_Good, Y_Good, U_Bad,V_Bad, Y_Bad] = lower_bound2(standard_deviation, mu)

[toa_Array, std_t, var_t] = lower_bound(standard_deviation);

% toa_Array = [UGoodUser,YGoodUser,UBadUser,VBadUser,YBadUser];

% % %%% Below for 2 verifiers
% U_Good = toa_Array(:,1:2);
% % Y_Good = toa_Array(:,3:4);
% Y_Good_old = toa_Array(:,3:4);
% U_Bad = toa_Array(:,5:6);
% V_Bad = toa_Array(:,7:8);
% % Y_Bad = toa_Array(:,9:10);
% Y_Bad_old = toa_Array(:,9:10);

% %%% Below for 3 verifiers
% U_Good = toa_Array(:,1:3);
% % Y_Good = toa_Array(:,4:6);
% Y_Good_old = toa_Array(:,4:6);
% U_Bad = toa_Array(:,7:9);
% V_Bad = toa_Array(:,10:12);
% % Y_Bad = toa_Array(:,13:15);
% Y_Bad_old = toa_Array(:,13:15);

%%% Below for 4 verifiers
U_Good = toa_Array(:,1:4);
% Y_Good = toa_Array(:,5:8);
Y_Good_old = toa_Array(:,5:8);
U_Bad = toa_Array(:,9:12);
V_Bad = toa_Array(:,13:16);
% Y_Bad = toa_Array(:,17:20);
Y_Bad_old = toa_Array(:,17:20);

% %%% Below for 5 verifiers
% U_Good = toa_Array(:,1:5);
% % Y_Good = toa_Array(:,6:10);
%  Y_Good_old = toa_Array(:,6:10);
% U_Bad = toa_Array(:,11:15);
% V_Bad = toa_Array(:,16:20);
% % Y_Bad = toa_Array(:,21:25);
% Y_Bad_old = toa_Array(:,21:25);

% % Below for 6 verifiers
% U_Good = toa_Array(:,1:6);
% % Y_Good = toa_Array(:,1:6);
% Y_Good_old = toa_Array(:,7:12);
% U_Bad = toa_Array(:,13:18);
% V_Bad = toa_Array(:,19:24);
% % Y_Bad = toa_Array(:,25:30);
% Y_Bad_old = toa_Array(:,25:30);

ii = length(U_Good);

% % [bias,std_i,var_i] = nonLOS(ii, mu, number_of_verifiers);
[bias,std_i,var_i] = nonLOS(ii, mu, 4);
Y_Good = Y_Good_old + bias; % NLoS biases to Good User only.
% Y_Good = Y_Good_old;
Y_Bad  = Y_Bad_old  + bias; % NLoS biases to Bad User only.
% Y_Bad = Y_Bad_old;

tag_good_users = zeros(ii,1);
tag_bad_users  =  ones(ii,1);

% toa_Array_NN  = [U_Good,Y_Good,tag_good_users; V_Bad,Y_Bad,tag_bad_users];
% toa_Array_LRT = [U_Good,Y_Good,tag_good_users, V_Bad,Y_Bad,tag_bad_users];

toa_Array_NN  = [U_Good,Y_Good,tag_good_users; U_Bad,Y_Bad,tag_bad_users];

end