

%Unchanged Code



clf;
rng(0);
% Load data in the below column wise order.
% BS1 is the middle BS and considered as a reference.
% Date,Time,BS1(Long&Lat),BS2(Long&Lat),BS3(Long&Lat),User(Long&Lat),RSS(V1,V2,V3)
load 3BS-Linksys-Data-pre-export-05032019.mat

% Extracting x and y coordinates for the BS and removing the duplicates.
BS  = data(:,3:8);%  Long Lat
verifier = unique(BS,'rows');
% verifier_pre = unique(BS,'rows');
% verifier = round(verifier_pre,6);
verifier_1 = verifier(:,1:2); % Manzoor
verifier_2 = verifier(:,3:4); % Hung
verifier_3 = verifier(:,5:6); % Ihsan
RSS_ref = -32.51;
d_ref = 1;
pathloss_constant = 2.27;
% Getting the the users true locations
claimed_location = data(:,9:10); % Long Lat

%claimed_location = round(claimed_location1,6);

size_users = size(claimed_location,1);

% verifier_1 is choosen as the reference in cartesian coordinate system i.e. [0,0].
% Converting verifier_2 and verifier_3 geodatic coordinates to cartesian coordinates w.r.t verifier_1.
[xEast_2,yNorth_2,zUp_2] = geodetic2enu(verifier_2(1,2),verifier_2(1,1),0,verifier_1(1,2),verifier_1(1,1),0,wgs84Ellipsoid);
verifier_2_cartesian= [xEast_2,yNorth_2];
[xEast_3,yNorth_3,zUp_3] = geodetic2enu(verifier_3(1,2),verifier_3(1,1),0,verifier_1(1,2),verifier_1(1,1),0,wgs84Ellipsoid);
verifier_3_cartesian= [xEast_3,yNorth_3];

% Converting the user's geodatic coordinates to cartesian coordinates w.r.t BS1.
% Calculating the distance between the users to the verifiers

for n = 1:size_users
    [xEast(n),yNorth(n),zUp(n)] = geodetic2enu(claimed_location(n,2),claimed_location(n,1),0,verifier_1(1,2),verifier_1(1,1),0,wgs84Ellipsoid);
    claimed_location_cartesian(n,1)= xEast(n);
    claimed_location_cartesian(n,2)= yNorth(n);
    distance_to_verifier_1(n,1) = norm([0,0] - claimed_location_cartesian(n,:));
    distance_to_verifier_2(n,1) = norm(verifier_2_cartesian - claimed_location_cartesian(n,:));
    distance_to_verifier_3(n,1) = norm(verifier_3_cartesian - claimed_location_cartesian(n,:));
    
    RSS_theory_V1(n,1) = RSS_ref + (10* pathloss_constant * log10(d_ref/distance_to_verifier_1(n,1)));
    RSS_theory_V2(n,1) = RSS_ref + (10* pathloss_constant * log10(d_ref/distance_to_verifier_2(n,1)));
    RSS_theory_V3(n,1) = RSS_ref + (10* pathloss_constant * log10(d_ref/distance_to_verifier_3(n,1)));
end

% Plotting the verifiers.
plot(0,0,'go','MarkerSize',15);
hold on;
plot(xEast_2,yNorth_2,'go','MarkerSize',15);
hold on;
plot(xEast_3,yNorth_3,'go','MarkerSize',15);
hold on;

for m=1:size_users
    plot(claimed_location_cartesian(m,1),claimed_location_cartesian(m,2),'r*','MarkerSize',15);
    hold on;
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_x = mean(claimed_location_cartesian(:,1));
min_y = mean(claimed_location_cartesian(:,2));

% The below set used mostly.
% x_min = mean_x - 120;
% x_max = mean_x;
% y_min = min_y - 120;
% y_max = min_y + 120;


r = 40;
rmax = r+0.1;

for nn=1:size_users
    
% %     xmax1 = claimed_location_cartesian(nn,1);
% %     x_max = xmax1 + (0.5*r);
% %     % x_max = xmax1;
% %     x_min = xmax1  - r;
% %     y_max = claimed_location_cartesian(nn,2) + (0.5*r);
% %     y_min = y_max - 2*r;
    
        x_max = 400;
        x_min = -400;
        y_max = 400;
        y_min = -400;
% %         
        
% %     x_min = mean_x - 150;
% %     x_max = mean_x;
% %     y_min = min_y - 100;
% %     y_max = min_y + 100;

%     x_min = -300;
%     x_max = -25;
%     y_min = -300;
%     y_max = 40;


    x0 = x_min:(x_max-x_min)/1000:x_max;
    y0 = y_min:(y_max-y_min)/1000:y_max;
    N_0 = length(x0);
    x0_p = x0'*ones(1,N_0);
    y0_p = ones(N_0,1)*y0;
    
    dct = sqrt((x0_p - claimed_location_cartesian(nn,1)).^2 + (y0_p - claimed_location_cartesian(nn,2)).^2);
    I = find(dct==min(dct(find(dct>=r & dct<=rmax))));
    I = I(1);
    x_coordinate = x0_p(I);
    y_coordinate = y0_p(I);
    
    mal_location (nn,:)= [x_coordinate,y_coordinate];
    [mal_lat(nn,1),mal_long(nn,1),mal_h(nn,1)] = enu2geodetic(x_coordinate,y_coordinate,0,verifier_1(1,2),verifier_1(1,1),0,wgs84Ellipsoid);
    
    % Calculating distance of Mal Location to various reference points
    % Distance of malicious Location to the claimed locations.
    mal_location_distance_to_True_Location(nn,1) = norm(claimed_location_cartesian(nn,:) - mal_location (nn,:) );
    % Distance of malicious Location to the verifier 1.
    mal_location_distance_to_v1(nn,1) = norm([0,0] - mal_location (nn,:) );
    mal_RSS_v1(nn,1) = RSS_ref + (10* pathloss_constant * log10(d_ref/mal_location_distance_to_v1(nn,1)));
    % Distance of malicious Location to the verifier 2.
    mal_location_distance_to_v2(nn,1) = norm(verifier_2_cartesian - mal_location (nn,:));
    mal_RSS_v2(nn,1) = RSS_ref + (10* pathloss_constant * log10(d_ref/mal_location_distance_to_v2(nn,1)));
    % Distance of malicious Location to the verifier 3.
    mal_location_distance_to_v3(nn,1) = norm(verifier_3_cartesian - mal_location (nn,:));
    mal_RSS_v3(nn,1) = RSS_ref + (10* pathloss_constant * log10(d_ref/mal_location_distance_to_v3(nn,1)));
    
    plot(mal_location(nn,1),mal_location(nn,2),'b^','MarkerSize',15);
    hold on;
    
    
    xmax1 = [];
    x_max  = [];
    x_min  = [];
    y_max  = [];
    y_min  = [];
    x0 = [];
    y0 = [];
    N_0 = [];
    x0_p = [];
    y0_p = [];
    dct = [];
    I = [];
    x_coordinate = [];
    y_coordinate = [];
end






% Exporting in order - Date,Time,V1-GPS(Lat/Long),V1-Cartesian,V2-GPS(Lat/Long),V2-Cartesian,V3-GPS(Lat/Long),V3-Cartesian,
%                      User's-Claimed-Location-GPS(Lat/Long),User's-Claimed-Location-Cartesian,distance-of-user-to-V1,
%                      distance-of-user-to-V2,distance-of-user-to-V3,RSS(V1,V2,V3),malicious-location-lat,malicious-location-long,
%                      malicious-location-Cartesian-Coordinates,malicious-location-distance-to-v1,malicious-location-distance-to-v2,
%                      malicious-location-distance-to-v3,malicious-location-distance-to-claimed-location
%
%
T_Export = [data(:,1:2),data(:,4),data(:,3),zeros(size_users,2),data(:,6),data(:,5),repmat(verifier_2_cartesian,size_users,1),...
    data(:,8),data(:,7),repmat(verifier_3_cartesian,size_users,1),...
    data(:,11:13),claimed_location(:,2),claimed_location(:,1),claimed_location_cartesian,round(distance_to_verifier_1,2),round(distance_to_verifier_2,2),...
    round(distance_to_verifier_3,2),RSS_theory_V1,RSS_theory_V2,RSS_theory_V3,mal_lat,mal_long,mal_location,round(mal_location_distance_to_v1,2),...
    round(mal_location_distance_to_v2,2),round(mal_location_distance_to_v3,2),mal_RSS_v1,mal_RSS_v2,mal_RSS_v3,mal_location_distance_to_True_Location];







% % % % Code for optimal attack location.
% % %
% % % RSS_theory = []; % Extracting theoratical RSS based on users true locations
% % % variance_matrix = [];
% % % std_matrix= [];
% % % user_x = [];
% % % user_y = [];
% % % x_max = [];
% % % x_min = [];
% % % y_max = [];
% % % y_min = [];
% % % x0 = [];
% % % y0 = [];
% % % r = [];
% % % N_0 = [];
% % % x0_p = [];
% % % y0_p = [];
% % % dct = [];
% % % I = [];
% % % x1 = [];
% % % y1 = []
% % % N_1 = [];
% % % pathloss_constant = [];
% % % verifier_pre = [];
% % % verifier = [];
% % % d2t=[]
% % % optimal_RSS=[];
% % % KL_Divergance=[];
% % %
% % %
% % % % c = 3e8; % speed of light
% % % % f = 2.437*10^9; %2.4 Ghz wifi band - channel 6
% % % % lembda = c/f;
% % % % p = 18; %transmit power of the user in dBm.
% % % RSS_theory = RSS_t(:,9:11); % Extracting theoratical RSS based on users true locations. For malicious it will be v
% % % RSS_theory=RSS_theory';
% % % verifier_pre = unique(RSS_t(:,1:6),'rows'); % verifiers coordinates
% % % verifier = [verifier_pre(1,1:2);verifier_pre(1,3:4);verifier_pre(1,5:6)];
% % % user_x = RSS_t(:,7);
% % % user_y = RSS_t(:,8);
% % % variance_matrix = [10.4184531663998,9.77605704307982,17.3572130860042];
% % % std_matrix= sqrt(variance_matrix);
% % % RSS_ref = -32.51;
% % % d_ref = 1;
% % % pathloss_constant = [2.27,2.27,2.27];
% % %
% % % x_max = 0;
% % % x_min = -150;
% % % y_max = 45;
% % % y_min = -105;
% % % x0 = x_min:(x_max-x_min)/2000:x_max;
% % % y0 = y_min:(y_max-y_min)/2000:y_max;
% % % r = 100; %the minimum distance between true location and claimed location
% % % N_0 = length(x0);
% % % x0_p = x0'*ones(1,N_0);
% % % y0_p = ones(N_0,1)*y0;
% % %
% % % len = length(RSS_theory);
% % %
% % %
% % % N = 3; % The number of verifiers
% % %
% % % for m=1:len
% % %     dct = sqrt((x0_p - user_x(m)).^2 + (y0_p - user_y(m)).^2);
% % %     I = find(dct>=r);
% % %     x1 = x0_p(I);
% % %     y1 = y0_p(I);
% % %     N_1 = length(x1);
% % %
% % %     for i = 1:N_1
% % %         for j = 1:N
% % %             d2t = sqrt((verifier(j,1) - x1(i))^2 + (verifier(j,2) - y1(i))^2);
% % %             optimal_RSS(j,1) = RSS_ref + (10* pathloss_constant(1,j) * log10(d_ref/d2t));
% % %             d2t=[];
% % %         end
% % %         KL_Divergance(i) = (1/2) * ((RSS_theory(:,m) - optimal_RSS)'/std_matrix); % 0.5*((v-u)/rho)^2
% % %         %         KL_Divergance(i) = (1/2) * ((optimal_RSS - RSS_theory(:,m))'/std_matrix);
% % %         optimal_RSS=[];
% % %     end
% % %     [~ , index] = min(KL_Divergance);
% % %     x_value(m) = x1(index);
% % %     y_value(m)= y1(index);
% % %     distance_check(m) = sqrt((user_x(m) - x_value(m)).^2 + (user_y(m) - y_value(m)).^2);
% % %     KL_Divergance = [];
% % % end
