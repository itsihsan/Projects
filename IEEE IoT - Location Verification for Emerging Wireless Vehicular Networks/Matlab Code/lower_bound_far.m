% % % % % % Modified Function
function [toa_Array, std_t, var_t] = lower_bound_far(standard_deviation)

toa_Array = []; d_ClaimedGoodUser = []; % d_TrueGoodUser = []; % VGoodUser = []; % x_value =[]; % y_value=[]; %index=[];
UGoodUser = []; YGoodUser = []; d_ClaimedBadUser = []; UBadUser = []; d_TrueBadUser=[];
WBadUser = []; T_WbadUser = []; YBadUser = [];

std_t = standard_deviation;  % Real Standard Deviation of Thermal noise in 'micro seconds'.
var_t = (std_t)^2;           % variance of Thermal noise

[good_user, ~] = locations(1,1);
bad_user = good_user;

verifier =  [0,125;0,375;1000,125;1000,375];           % 4 verifiers
SpeedOfLight = 0.3; % Speed of Light in meters per nano second

ClaimedGoodUser = good_user;
ClaimedBadUser  = bad_user;

N=4; % Number of Verifiers
sigma = noise(std_t, N);  % Thermal noise of the reciever in micro seconds

for j = 1:N
    d_ClaimedGoodUser(1,j) = norm(verifier(j,:) - ClaimedGoodUser);    % Distance of ClaimedGoodUser (Good User's claimed Location)from j-th verifier
    UGoodUser (1,j) = d_ClaimedGoodUser(1,j)/SpeedOfLight  ;           % U under Ho using good user's claimed location.
    YGoodUser (1,j) = UGoodUser (1,j) + sigma(1,j)       ;             %
    
    d_ClaimedBadUser (1,j) = norm(verifier(j,:) - ClaimedBadUser);     % Distance of ClaimedBadUser (Bad User's claimed Location)from j-th verifier.
    UBadUser (1,j)         = d_ClaimedBadUser (1,j)/SpeedOfLight ;     % U using bad user's claimed location. (Same as good user's true location in our experiment).
end

%   UBadUser
Avg_U = mean(UGoodUser,2);
VBadUser = repmat(Avg_U,1,N);
YBadUser = VBadUser + sigma;

toa_Array = [UGoodUser,YGoodUser,UBadUser,VBadUser,YBadUser];
end