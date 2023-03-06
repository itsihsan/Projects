% % % % % % Modified Function
function [toa_Array, std_t, var_t] = lower_bound_near(r , standard_deviation)

toa_Array = []; d_ClaimedGoodUser = []; % d_TrueGoodUser = []; % VGoodUser = []; % x_value =[]; % y_value=[]; %index=[];
UGoodUser = []; YGoodUser = []; d_ClaimedBadUser = []; UBadUser = []; d_TrueBadUser=[];
WBadUser = []; T_WbadUser = []; YBadUser = [];

std_t = standard_deviation;  % Real Standard Deviation of Thermal noise in 'micro seconds'.
var_t = (std_t)^2;           % variance of Thermal noise

[good_user, ~] = locations(1,1);
bad_user = good_user;
verifier =  [0,125;0,375;1000,125;1000,375];                           % 4 verifiers
SpeedOfLight = 0.3; % Speed of Light in meters per nano second

ClaimedGoodUser = good_user;
ClaimedBadUser  = bad_user;

N=4; % Number of Verifiers
sigma = noise(std_t, N);  % Thermal noise of the reciever in micro seconds

for j = 1:N
    d_ClaimedGoodUser(1,j) = norm(verifier(j,:) - ClaimedGoodUser);    % Distance of ClaimedGoodUser (Good User's claimed Location)from j-th verifier
    UGoodUser (1,j) = d_ClaimedGoodUser(1,j)/SpeedOfLight;             % U under Ho using good user's claimed location.
    YGoodUser (1,j) = UGoodUser (1,j) + sigma(1,j);                    %
    
    d_ClaimedBadUser (1,j) = norm(verifier(j,:) - ClaimedBadUser);     % Distance of ClaimedBadUser (Bad User's claimed Location)from j-th verifier.
    UBadUser (1,j)         = d_ClaimedBadUser (1,j)/SpeedOfLight;      % U using bad user's claimed location. (Same as good user's true location in our experiment).
end

%  We now make calculations for Tx_star and WbadUser
x0 =  -2500:(2500-(-2500))/5001:2500;
y0 =  -2500 :(2500-(-2500))/5001 :2500;

rmin = r;
rmax = r + 10;

N_0 = length(x0);
x0_p = x0'*ones(1,N_0);
y0_p = ones(N_0,1)*y0;
dct = sqrt((x0_p - ClaimedBadUser(1,1)).^2 + (y0_p - ClaimedBadUser(1,2)).^2);

I = find(dct>=rmin & dct<=rmax);
JJ = min(dct(I));
JJJ = find(dct==JJ);
x_coordinate = x0_p(JJJ);
y_coordinate = y0_p(JJJ);
TrueBadUser = [x_coordinate , y_coordinate];

for j = 1:N
    d2t(1,j) = norm(verifier(j,:)-TrueBadUser);
    T_WbadUser(1,j) = d2t(1,j)/SpeedOfLight;
end

Tx_star = (1/N) * sum(UBadUser - T_WbadUser);

d_TrueBadUser = [];
for jj = 1:N
    d_TrueBadUser(1,jj) = norm(verifier(jj,:) - TrueBadUser);       % Distance of TrueBadUser from jj-th verifier
    WBadUser (1,jj)     = d_TrueBadUser(1,jj)/SpeedOfLight;         % V using good user's True location. Same as his claimed Location.
end

VBadUser = Tx_star + WBadUser;
YBadUser = VBadUser + sigma;

toa_Array = [UGoodUser,YGoodUser,UBadUser,VBadUser,YBadUser];
end