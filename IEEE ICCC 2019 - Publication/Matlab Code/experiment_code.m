clf;
load Linksys_28Feb_1.mat
% % % R = 6371000; % Radius of the earth.
% pathloss_constant = 3;
% transmit_power = 30; % in dBm
% do = 1; % Reference Distance
% c = 3*10^8; %Speed of light.
%
% f = 2.412*10^9; %2.4 Ghz wifi band - channel 1 - iphone
% f = 2.437*10^9; %2.4 Ghz wifi band - channel 6 - linksys
% lembda = c/f; % Wavelength for 900 MHz band
% p = transmit_power + 20 * log10(lembda/(4*pi*do)) % RSS at the reference distance

% Extracting x and y coordinates for the BS and removing the duplicates.
BS  = data(:,3:4);%  Long Lat
verifier1 = unique(BS,'rows');
verifier = round(verifier1,6);
% Getting the claimed locations for the users
claimed_location1 = data(:,6:7); % Long Lat
claimed_location = round(claimed_location1,6);

% % % Converting BS coordinaes from Degrees to Radians
% % % ver_lat_rad = verifier1(1,2) * pi/180;
% % % ver_lon_rad = verifier1(1,1) * pi/180;

size_users = size(claimed_location,1);

% Calculating the distance between the users and the BS

for n = 1:size_users
    % % %         % Converting users coordinaes from Degrees to Radians
    % % %         user_lat_rad = claimed_location1(n,2) * pi/180;
    % % %         user_lon_rad = claimed_location1(n,1) * pi/180;
    % % %         % Haversine formula
    % % %         del_phi = user_lat_rad - ver_lat_rad;
    % % %         del_lambda = user_lon_rad - ver_lon_rad;
    % % %         a = (sin(del_phi/2) * sin(del_phi/2)) + (cos(user_lat_rad)*cos(ver_lat_rad)*sin(del_lambda/2)*sin(del_lambda/2));
    % % %         c = 2*atan2(sqrt(a),sqrt(1-a));
    % % %         d(n) = R*c;
    % % %
    % % %             % Distance verification check using law of cosines.
    % % %             d2(n)= acos( (sin(user_lat_rad) * sin(ver_lat_rad)) + (cos(user_lat_rad) * cos(ver_lat_rad) * cos(del_lambda))) * R;
    % % %
    % % %     % Calculating the signal strength using claimed location.
    % % % %     U(n) = p - 10 * pathloss_constant * log10(d(n)/do);
    % % % %     U(end)
    [xEast(n),yNorth(n),zUp(n)] = geodetic2enu(claimed_location(n,2),claimed_location(n,1),0,verifier(1,2),verifier(1,1),0,wgs84Ellipsoid);
    claimed_location_cartesian(n,1)= xEast(n);
    claimed_location_cartesian(n,2)= yNorth(n);
    distance_to_verifier(n,1) = norm([0,0] - claimed_location_cartesian(n,:));
    % % %     distance_to_verifier_check(n,1)= d(n);
    % % %     distance_to_verifier_check2(n,1)= d2(n);
    % % %     user_lat_rad = [];
    % % %     user_lon_rad = [];
    % % %     del_phi = [];
    % % %     del_lambda = [];
    % % %     a = [];
    % % %     c = [];
    %     plot(claimed_location_cartesian(n,1),claimed_location_cartesian(n,2))
    %     hold on;
end

plot(0,0,'g*','MarkerSize',15);
hold on;
for hh=1:size_users
    plot(claimed_location_cartesian(hh,1),claimed_location_cartesian(hh,2),'r*','MarkerSize',15);
    hold on
end

mean_x = mean(claimed_location_cartesian(:,1));
min_y = min(claimed_location_cartesian(:,2));

% x_min = mean_x - 100;
% x_max = mean_x + 100;
% y_min = min_y;
% y_max = min_y - 100;

% The below set used mostly.
% x_min = mean_x - 100;
% x_max = mean_x;
% y_min = min_y - 100;
% y_max = min_y + 100;

% Free choice of choosing the claimed location anywhere even around BS
x_min = -400;
x_max = 400;
y_min = -400;
y_max = 400;


x0 = x_min:(x_max-x_min)/4000:x_max;
y0 = y_min:(y_max-y_min)/4000:y_max;
r = 50; %the minimum distance between true location and claimed location
% rmin = r-0.002;
rmax = r+0.1;
N_0 = length(x0);
x0_p = x0'*ones(1,N_0);
y0_p = ones(N_0,1)*y0;

for nn=1:size_users
    %     for mm=1:length(x0_p)
    dct = sqrt((x0_p - claimed_location_cartesian(nn,1)).^2 + (y0_p - claimed_location_cartesian(nn,2)).^2);
    %         I = find(dct>=r & dct<=rmax)
    %         dct(I)
    %         JJ = min(dct(I))
    %         JJJ = find(dct==JJ)
    I = find(dct==min(dct(find(dct>=r & dct<=rmax))));
    I = I(1);
    x_coordinate = x0_p(I);
    %         x_coordinate = unique(x_coordinate)
    y_coordinate = y0_p(I);
    %         y_coordinate = unique(y_coordinate)
    mal_location (nn,:)= [x_coordinate,y_coordinate];
    [mal_lat(nn,1),mal_long(nn,1),mal_h(nn,1)] = enu2geodetic(x_coordinate,y_coordinate,0,verifier(1,2),verifier(1,1),0,wgs84Ellipsoid);
    mal_location_distance_to_True_Location(nn,1) = norm(claimed_location_cartesian(nn,:) - mal_location (nn,:) );
    mal_location_distance_to_BS(nn,1) = norm([0,0] - mal_location (nn,:) );
    plot(mal_location(nn,1),mal_location(nn,2),'b^','MarkerSize',15);
    hold on;
    dct = [];
    I = [];
    x_coordinate = [];
    y_coordinate = [];
end

Export = [data(:,1:2),data(:,4),data(:,3),zeros(size_users,2),claimed_location(:,2),claimed_location(:,1),claimed_location_cartesian,distance_to_verifier,mal_lat,mal_long,mal_location,mal_location_distance_to_BS,mal_location_distance_to_True_Location,data(:,5)];

% I = find(dct>=r);
% x1 = x0_p(I);
% y1 = y0_p(I);
% N_1 = length(x1);
%
