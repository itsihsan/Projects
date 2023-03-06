clf;
rng(0);

load 3BS-Linksys-Data-pre-export-05032019.mat
% data = data(randperm(end),:);
% data = data(1:50,:);

verifier_pre = [];
verifier = [];
user_x = [];
user_y = [];
variance_matrix = [];
std_matrix= [];
RSS_theory = []; % Extracting theoratical RSS based on users true locations
x_max = [];
x_min = [];
y_max = [];
y_min = [];
x0 = [];
y0 = [];
r = [];
N_0 = [];
x0_p = [];
y0_p = [];
dct = [];
I = [];
x1 = [];
y1 = [];
N_1 = [];
pathloss_constant = [];
d2t=[];
optimal_RSS=[];
KL_Divergance=[];

% Load data in the below column wise order.
% BS1 is the middle BS and considered as a reference.
% Date,Time,BS1(Long&Lat),BS2(Long&Lat),BS3(Long&Lat),User(Long&Lat),RSS(V1,V2,V3)


variance_matrix = [10.4184531663998,9.77605704307982,17.3572130860042];
std_matrix= sqrt(variance_matrix);
RSS_ref = -32.51;
d_ref = 1;
pathloss_constant = [2.27,2.27,2.27];

% Extracting x and y coordinates for the BS and removing the duplicates.
BS  = data(:,3:8);%  Long Lat
verifier = unique(BS,'rows');
verifier_1 = verifier(:,1:2); % Manzoor
verifier_2 = verifier(:,3:4); % Hung
verifier_3 = verifier(:,5:6); % Ihsan

% Getting the claimed locations for the users
claimed_location = data(:,9:10); % Long Lat

%claimed_location = round(claimed_location1,6);

size_users = size(claimed_location,1);

% verifier_1 is choosen as the reference in cartesian coordinate system i.e. [0,0].
% Converting verifier_2 and verifier_3 geodatic coordinates to cartesian coordinates w.r.t verifier_1.
[xEast_2,yNorth_2,zUp_2] = geodetic2enu(verifier_2(1,2),verifier_2(1,1),0,verifier_1(1,2),verifier_1(1,1),0,wgs84Ellipsoid);
verifier_2_cartesian= [xEast_2,yNorth_2];
[xEast_3,yNorth_3,zUp_3] = geodetic2enu(verifier_3(1,2),verifier_3(1,1),0,verifier_1(1,2),verifier_1(1,1),0,wgs84Ellipsoid);
verifier_3_cartesian= [xEast_3,yNorth_3];

verifier = [0,0;verifier_2_cartesian;verifier_3_cartesian]; % All verifiers in cartesian coordinates.

% Converting the user's geodatic coordinates to cartesian coordinates w.r.t BS1.
% Calculating the distance between the users to the verifiers

for n = 1:size_users
    [xEast(n),yNorth(n),zUp(n)] = geodetic2enu(claimed_location(n,2),claimed_location(n,1),0,verifier_1(1,2),verifier_1(1,1),0,wgs84Ellipsoid);
    claimed_location_cartesian(n,1)= xEast(n);
    claimed_location_cartesian(n,2)= yNorth(n);
    distance_to_verifier_1(n,1) = norm([0,0] - claimed_location_cartesian(n,:));
    distance_to_verifier_2(n,1) = norm(verifier_2_cartesian - claimed_location_cartesian(n,:));
    distance_to_verifier_3(n,1) = norm(verifier_3_cartesian - claimed_location_cartesian(n,:));
    
    % Calculating theoratical RSS using user's true distance to verifier.
    RSS_theory_V1(n,1) = RSS_ref + (10* pathloss_constant(1) * log10(d_ref/distance_to_verifier_1(n,1)));
    RSS_theory_V2(n,1) = RSS_ref + (10* pathloss_constant(1) * log10(d_ref/distance_to_verifier_2(n,1)));
    RSS_theory_V3(n,1) = RSS_ref + (10* pathloss_constant(1) * log10(d_ref/distance_to_verifier_3(n,1)));
end
user_x = xEast;  % x coordinate of the user true location.
user_y = yNorth; % y coordinate of the user true location.


% Plotting the verifiers and users true locations.
plot(0,0,'go','MarkerSize',15);
hold on;
plot(xEast_2,yNorth_2,'go','MarkerSize',15);
hold on;
plot(xEast_3,yNorth_3,'go','MarkerSize',15);
hold on;

for mm=1:size_users
    plot(claimed_location_cartesian(mm,1),claimed_location_cartesian(mm,2),'r*','MarkerSize',15);
    hold on;
end


% Code for optimal attack location.

RSS_theory = [RSS_theory_V1,RSS_theory_V2,RSS_theory_V3]; % Extracting theoratical RSS based on users true locations. For malicious it will be v
RSS_theory=RSS_theory';

% x_max = 0;
% x_min = -150;
% y_max = 45;
% y_min = -105;
% x0 = x_min:(x_max-x_min)/4000:x_max;
% y0 = y_min:(y_max-y_min)/4000:y_max;
r = 50; %the minimum distance between true location and claimed location
% N_0 = length(x0);
% x0_p = x0'*ones(1,N_0);
% y0_p = ones(N_0,1)*y0;
len = length(RSS_theory);
N = 3; % The number of verifiers

for m=1:len
    
% %         x_max = claimed_location_cartesian(m,1)+(0.5*r);
% %         %     x_max = -25;
% %         x_min = x_max  - (1.5*r);
% %         %         y_max = claimed_location_cartesian(m,2) + r;
% %         y_max = claimed_location_cartesian(m,2) + (0.5*r);
% %         y_min = y_max - (1.5*r);
% %         %     y_max = 40;
% %         %     y_min = y_max - (1.5*r);
    
    x_max = 250;
    x_min = -250;
    y_max = 250;
    y_min = -250;
    
    x0 = x_min:(x_max-x_min)/501:x_max;
    y0 = y_min:(y_max-y_min)/501:y_max;
    N_0 = length(x0);
    x0_p = x0'*ones(1,N_0);
    y0_p = ones(N_0,1)*y0;
    
    dct = sqrt((x0_p - user_x(m)).^2 + (y0_p - user_y(m)).^2);
    I = find(dct>=r);
    x1 = x0_p(I);
    y1 = y0_p(I);
    N_1 = length(x1);
    
    for i = 1:N_1
        for j = 1:N
            d2t = sqrt((verifier(j,1) - x1(i))^2 + (verifier(j,2) - y1(i))^2);
            optimal_RSS(j,1) = RSS_ref + (10* pathloss_constant(1,j) * log10(d_ref/d2t));
            d2t=[];
        end
%         KL_Divergance(i) = (1/2) * ((RSS_theory(:,m) - optimal_RSS)'/std_matrix)^2; % 0.5*((v-u)/rho)^2
%        KL_Divergance(i) = (1/2) * (Delta_V - Delta_U)' * D^(-1) * (Delta_V - Delta_U);
       KL_Divergance(i) = (1/2) * (RSS_theory(:,m) - optimal_RSS)' * (diag(variance_matrix))^(-1) * (RSS_theory(:,m) - optimal_RSS);
       optimal_RSS=[];
    end
    
    [~ , index] = min(KL_Divergance)
    x_value(m,1) = x1(index)
    y_value(m,1)= y1(index)
    index=[];
    distance_check(m) = sqrt((user_x(m) - x_value(m,1)).^2 + (user_y(m) - y_value(m,1)).^2)
    KL_Divergance = [];
    
    mal_location (m,:)= [x_value(m,1),y_value(m,1)];
    [mal_lat(m,1),mal_long(m,1),mal_h(m,1)] = enu2geodetic(x_value(m,1),y_value(m,1),0,verifier_1(1,2),verifier_1(1,1),0,wgs84Ellipsoid);
    
    % Calculating distance of Mal Location to various reference points
    % Distance of malicious Location to the claimed locations.
    mal_location_distance_to_True_Location(m,1) = norm(claimed_location_cartesian(m,:) - mal_location (m,:) );
    % Distance of malicious Location to the verifier 1.
    mal_location_distance_to_v1(m,1) = norm([0,0] - mal_location (m,:) );
    mal_RSS_v1(m,1) = RSS_ref + (10* pathloss_constant(1) * log10(d_ref/mal_location_distance_to_v1(m,1)));
    % Distance of malicious Location to the verifier 2.
    mal_location_distance_to_v2(m,1) = norm(verifier_2_cartesian - mal_location (m,:));
    mal_RSS_v2(m,1) = RSS_ref + (10* pathloss_constant(1) * log10(d_ref/mal_location_distance_to_v2(m,1)));
    
    % Distance of malicious Location to the verifier 3.
    mal_location_distance_to_v3(m,1) = norm(verifier_3_cartesian - mal_location (m,:));
    mal_RSS_v3(m,1) = RSS_ref + (10* pathloss_constant(1) * log10(d_ref/mal_location_distance_to_v3(m,1)));
    
    plot(mal_location(m,1),mal_location(m,2),'b^','MarkerSize',15);
    hold on;
    
end


% Exporting in order -
% Date,Time,V1-GPS(Lat/Long),V1-Cartesian,V2-GPS(Lat/Long),V2-Cartesian,V3-GPS(Lat/Long),V3-Cartesian,RSS(V1,V2,V3),
%                      User's-Claimed-Location-GPS(Lat/Long),User's-Claimed-Location-Cartesian,distance-of-user-to-V1,
%                      distance-of-user-to-V2,distance-of-user-to-V3,Theoratical_RSS_true_location(V1,V2,V3),malicious-location-lat,malicious-location-long,
%                      malicious-location-Cartesian-Coordinates,malicious-location-distance-to-v1,malicious-location-distance-to-v2,
%                      malicious-location-distance-to-v3,Mal_loc_RSS(V1,V2,V3),malicious-location-distance-to-claimed-location
%
%
T_Export = [data(:,1:2),data(:,4),data(:,3),zeros(len,2),data(:,6),data(:,5),repmat(verifier_2_cartesian,len,1),...
    data(:,8),data(:,7),repmat(verifier_3_cartesian,len,1),...
    data(:,11:13),claimed_location(:,2),claimed_location(:,1),claimed_location_cartesian,round(distance_to_verifier_1,2),round(distance_to_verifier_2,2),...
    round(distance_to_verifier_3,2),RSS_theory_V1,RSS_theory_V2,RSS_theory_V3,mal_lat,mal_long,mal_location,round(mal_location_distance_to_v1,2),...
    round(mal_location_distance_to_v2,2),round(mal_location_distance_to_v3,2),mal_RSS_v1,mal_RSS_v2,mal_RSS_v3,round(mal_location_distance_to_True_Location,2)];














% % % %
% % % %
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % mean_x = mean(claimed_location_cartesian(:,1));
% % % % min_y = mean(claimed_location_cartesian(:,2));
% % % %
% % % % % The below set used mostly.
% % % % % x_min = mean_x - 120;
% % % % % x_max = mean_x;
% % % % % y_min = min_y - 120;
% % % % % y_max = min_y + 120;
% % % %
% % % %
% % % % r = 50;
% % % % rmax = r+0.1;
% % % %
% % % % for nn=1:size_users
% % % %
% % % %     xmax1 = claimed_location_cartesian(nn,1);
% % % %     x_max = xmax1 + r;
% % % %     x_max = xmax1;
% % % %     x_min = xmax1  - r;
% % % %     y_max = claimed_location_cartesian(nn,2) + r;
% % % %     y_min = y_max - 2*r;
% % % %
% % % %     %     x_max = 400;
% % % %     %     x_min = -400;
% % % %     %     y_max = 400;
% % % %     %     y_min = -400;
% % % %     %     x_min = mean_x - 100;
% % % %     % x_max = mean_x;
% % % %     % y_min = min_y - 100;
% % % %     % y_max = min_y + 100;
% % % %     x0 = x_min:(x_max-x_min)/4000:x_max;
% % % %     y0 = y_min:(y_max-y_min)/4000:y_max;
% % % %     N_0 = length(x0);
% % % %     x0_p = x0'*ones(1,N_0);
% % % %     y0_p = ones(N_0,1)*y0;
% % % %
% % % %     dct = sqrt((x0_p - claimed_location_cartesian(nn,1)).^2 + (y0_p - claimed_location_cartesian(nn,2)).^2);
% % % %     I = find(dct==min(dct(find(dct>=r & dct<=rmax))));
% % % %     I = I(1);
% % % %     x_coordinate = x0_p(I);
% % % %     y_coordinate = y0_p(I);
% % % %
% % % %     mal_location (nn,:)= [x_coordinate,y_coordinate];
% % % %     [mal_lat(nn,1),mal_long(nn,1),mal_h(nn,1)] = enu2geodetic(x_coordinate,y_coordinate,0,verifier_1(1,2),verifier_1(1,1),0,wgs84Ellipsoid);
% % % %
% % % %     % Calculating distance of Mal Location to various reference points
% % % %     % Distance of malicious Location to the claimed locations.
% % % %     mal_location_distance_to_True_Location(nn,1) = norm(claimed_location_cartesian(nn,:) - mal_location (nn,:) );
% % % %     % Distance of malicious Location to the verifier 1.
% % % %     mal_location_distance_to_v1(nn,1) = norm([0,0] - mal_location (nn,:) );
% % % %     % Distance of malicious Location to the verifier 2.
% % % %     mal_location_distance_to_v2(nn,1) = norm(verifier_2_cartesian - mal_location (nn,:));
% % % %     % Distance of malicious Location to the verifier 3.
% % % %     mal_location_distance_to_v3(nn,1) = norm(verifier_3_cartesian - mal_location (nn,:));
% % % %
% % % %     plot(mal_location(nn,1),mal_location(nn,2),'b^','MarkerSize',15);
% % % %     hold on;
% % % %
% % % %
% % % %     xmax1 = [];
% % % %     x_max  = [];
% % % %     x_min  = [];
% % % %     y_max  = [];
% % % %     y_min  = [];
% % % %     x0 = [];
% % % %     y0 = [];
% % % %     N_0 = [];
% % % %     x0_p = [];
% % % %     y0_p = [];
% % % %     dct = [];
% % % %     I = [];
% % % %     x_coordinate = [];
% % % %     y_coordinate = [];
% % % % end
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % % % Exporting in order - Date,Time,V1-GPS(Lat/Long),V1-Cartesian,V2-GPS(Lat/Long),V2-Cartesian,V3-GPS(Lat/Long),V3-Cartesian,
% % % % %                      User's-Claimed-Location-GPS(Lat/Long),User's-Claimed-Location-Cartesian,distance-of-user-to-V1,
% % % % %                      distance-of-user-to-V2,distance-of-user-to-V3,RSS(V1,V2,V3),malicious-location-lat,malicious-location-long,
% % % % %                      malicious-location-Cartesian-Coordinates,malicious-location-distance-to-v1,malicious-location-distance-to-v2,
% % % % %                      malicious-location-distance-to-v3,malicious-location-distance-to-claimed-location
% % % % %
% % % % %
% % % % T_Export = [data(:,1:2),data(:,4),data(:,3),zeros(size_users,2),data(:,6),data(:,5),repmat(verifier_2_cartesian,size_users,1),...
% % % %     data(:,8),data(:,7),repmat(verifier_3_cartesian,size_users,1),...
% % % %     claimed_location(:,2),claimed_location(:,1),claimed_location_cartesian,distance_to_verifier_1,distance_to_verifier_2,...
% % % %     distance_to_verifier_3,data(:,11:13),mal_lat,mal_long,mal_location,mal_location_distance_to_v1,...
% % % %     mal_location_distance_to_v2,mal_location_distance_to_v3,mal_location_distance_to_True_Location];
% % % %
