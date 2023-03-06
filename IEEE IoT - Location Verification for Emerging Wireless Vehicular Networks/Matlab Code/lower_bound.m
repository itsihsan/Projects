% % % % % % Modified Function
function [toa_Array, std_t, var_t] = lower_bound(standard_deviation)
%  function TOA(number_of_verifiers,standard_deviation)

toa_Array = [];
d_ClaimedGoodUser = [];
% d_TrueGoodUser = [];
UGoodUser = [];
% VGoodUser = [];
YGoodUser = [];
d_ClaimedBadUser = [];
UBadUser = [];
d_TrueBadUser=[];
WBadUser = [];
T_WbadUser = [];
index=[];
YBadUser = [];
% x_value =[];
% y_value=[];


std_t = standard_deviation;  % Real Standard Deviation of Thermal noise in 'micro seconds'.
var_t = (std_t)^2;           % variance of Thermal noise
clf('reset');

for k =1:1:2500
    
    %   [good_user, bad_user] = locations(good, bad) % gives claimed locations of both users.
    
    %     [good_user, bad_user] = locations(1,1);
    [good_user, ~] = locations(1,1);
    bad_user = good_user;
    
%         verifier =  [0,250;1000,250];                % 2 verifiers
%     verifier =  [0,125;0,375;1000,250];                  % 3 verifiers
    verifier =  [0,125;0,375;1000,125;1000,375];           % 4 verifiers
%       verifier =  [0,125;0,250;0,375;1000,188;1000,312]; % 5 verifiers
%       verifier =  [0,125;0,250;0,375;1000,125;1000,250;1000,375]; % 6 verifiers
    %     verifier =  [0,150;0,450;1000,150;1000,450];

    v_x = verifier(:,1); v_y = verifier(:,2);
    g_x = good_user(:,1); g_y = good_user(:,2);
    b_x = bad_user(:,1); b_y = bad_user(:,2);
    scatter(v_x,v_y);
    hold on;
    plot(g_x,g_y,'g--o','MarkerSize',10)
    hold on;
    plot(b_x,b_y,'bs','MarkerSize',10)
    hold on;
    
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
    x0 =  250:(750-(250))/1001:750;
    y0 =   0 :(500-(0))/1001 :500;
    
    rmin = 200.00;
    rmax = 200.05;
    N_0 = length(x0);
    x0_p = x0'*ones(1,N_0);
    y0_p = ones(N_0,1)*y0;
    dct = sqrt((x0_p - ClaimedBadUser(1,1)).^2 + (y0_p - ClaimedBadUser(1,2)).^2);
    
    I = find(dct>=rmin & dct<=rmax);
    JJ = min(dct(I));
    JJJ = find(dct==JJ);
    % %        dct(JJJ) % just checking if the above condition is working
    x_coordinate = x0_p(JJJ);
    y_coordinate = y0_p(JJJ);
    TrueBadUser = [x_coordinate , y_coordinate];
    plot(x_coordinate,y_coordinate,'rs','MarkerSize',10)
    
    for j = 1:N
        %         d2t2(1,j)  = sqrt((verifier(j,1) - x_coordinate)^2 + (verifier(j,2) - y_coordinate)^2)
        %         T_WbadUser2 (1,j) = d2t2(1,j)/SpeedOfLight
        d2t(1,j) = norm(verifier(j,:)-TrueBadUser);
        T_WbadUser(1,j) = d2t(1,j)/SpeedOfLight;
    end
    Tx_star = (1/N) * sum(UBadUser - T_WbadUser);
    
    %          plot(x_coordinate,y_coordinate,'ro')
    % %       distance_check = norm(bad_user - TrueBadUser)
    % %       distance_check2 = norm(TrueBadUser - bad_user)
    
    %   I = find(dct>=r);
    % % %     x1 = x0_p(I);
    % % %     y1 = y0_p(I);
    % % %
    % % %     %     scatter(x1,y1,'*')
    % % %
    % % %     N_1 = length(x1);
    % % %     KL_Divergance = [];
    % % %
    % % %     for i = 1:N_1
    % % %         for j = 1:N
    % % %             d2t(1,j)  = sqrt((verifier(j,1) - x1(i))^2 + (verifier(j,2) - y1(i))^2);
    % % %             T_WbadUser (1,j) = d2t(1,j)/SpeedOfLight;
    % % %         end
    % % %         KL_Divergance(i) = (1/N) * sum(UBadUser - T_WbadUser);
    % % %     end
    % % %     %   size_KL_Divergance = size(KL_Divergance)
    % % %     %   T_WbadUser
    % % %     %   size_T_WbadUser = size(T_WbadUser)
    % % %     %     [~ , index] = min(KL_Divergance);
    % % %     [value , index] = min(KL_Divergance)
    % % %     Tx_star = value
    % % %     x_value = x1(index);
    % % %     y_value = y1(index);
    % % %     TrueBadUser = [x_value,y_value]
    % % %     plot(x_value,y_value,'rs','MarkerSize',10)
    % % %     %    distance = norm(TrueBadUser - bad_user)
    d_TrueBadUser = [];
    %     WBadUser = [];
    for jj = 1:N
        %         d_TrueBadUser(1,jj)  = sqrt((verifier(jj,1) - TrueBadUser(1,1))^2 + (verifier(jj,2) - TrueBadUser(1,2))^2);
        d_TrueBadUser(1,jj) = norm(verifier(jj,:) - TrueBadUser);       % Distance of TrueBadUser from jj-th verifier
        WBadUser (1,jj)     = d_TrueBadUser(1,jj)/SpeedOfLight;         % V using good user's True location. Same as his claimed Location.
    end
    %     WBadUser
    %  VBadUser
    VBadUser = Tx_star + WBadUser;
    %     VBadUser = WBadUser; % Tx = 0 Trail Reasons
    YBadUser = VBadUser + sigma;
    
    %     toa_Array(k,:) = [UGoodUser, VGoodUser, YGoodUser, UBadUser, VBadUser, YBadUser];
    toa_Array(k,:) = [UGoodUser,YGoodUser,UBadUser,VBadUser,YBadUser];
end
end