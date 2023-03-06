
%%%%%%%%%%%Constructing Covariance R matrixs Begin      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_Base_stations=4;
pathloss_constant = 3;
transmit_power = 30;

shadowing_variance = 2;  % Shadowing Variance
do = 1;                    % Reference Distance
lembda = (3*10^8)/(9*10^9); % lembda = c/f for 900 MHz band
p = transmit_power + 20 * log10(lembda/(4*pi*do));
N = N_Base_stations;
standard_deviation= 5;   
 Gt=1;
 Gr=1;                 
 PL_Ref_dB = -10*log10(Gt*Gr*lembda^2./(4*pi*do).^2);

xb = 150 * ( - 1 + 2 * rand ( 1 ,N_Base_stations ) );
yb = 150 * ( - 1 + 2 * rand ( 1 ,N_Base_stations ) );

Base_station = [10 10; 10 250; 250 10; 250 250];   % base stations

% Claimed Location for malicious or legitimate vehicle (Randomly Generated)
xc = rand(1,1)*200;  %%claimed location  (malicious/legitimate)
yc = rand(1,1)*200;

xt = rand(1,1)*100;   %Legitimate nodes
yt= rand(1,1)*100;

true_location = [xt yt];
claimed_location = [xc yc];

%%%%%  Calculating distances  of all the points in 200*200 m area from a claimed location %%%%%%%
    Q= zeros(200,200);
   
  claim_x = claimed_location(1);
  claim_y = claimed_location(2);
  
  for x=1:1: 200
        for y=1:1: 200
            Q(x,y)= sqrt((claim_x - x)^2 + (claim_y - y)^2);
        end
  end
r=100;  %range for threat model 
filtered_range = ((Q>=r) & (Q<(r+0.5)));   %%% filtering all the distance which is 100 meters away from the claimed location
[filtered_x,filtered_y] = find(filtered_range);

 Malicious_loc=[filtered_x filtered_y];
 
 r1=50;  % range is in meters
filtered_range1 = ((Q<=r1));   %%% filtering all the distance which is 100 meters away from the claimed location
[filtered_x1,filtered_y1] = find(filtered_range1);

Legitimate_loc=[filtered_x1 filtered_y1];
LegitimateX=sum(filtered_x1);



%%%%%%%%%%%%%%%%%% RSS from True Locations %%%%%%%%%%%%%%%%%%
for i = 1 : size (Base_station,1)
      d1(i) = norm (true_location - Base_station(i,:))
      PL_d_dB_tr(i) = PL_Ref_dB + 10*pathloss_constant*log10(d1(i)/do)   % no shadowing
      RSS_measured(i) = transmit_power - PL_d_dB_tr(i);              % RSS measured from true location
    
end
 RSS_measured;
  %%%calculate Ixx
    f= [10/log(10)];
    for a=1:N
      n1=pathloss_constant^2*(Base_station(a,1)-xt)^2*f;
     d11=((Base_station(a,1)+xt)^2+(Base_station(a,2)+yt)^2)^2;
     n2= 2*pathloss_constant*(Base_station(a,1)+xt)^2*(RSS_measured-pathloss_constant*log(d1));
     d22=((Base_station(a,1)+xt)^2+(Base_station(a,2)+yt)^2)
     n3=2*pathloss_constant*(RSS_measured-pathloss_constant*log(d1))
     d3=((Base_station(a,1)-xt)^2+(Base_station(a,2)-yt)^2)
      K=  sum(n1/d11+n2/d22+n3/d3)
      Ixx= (-1)* 1/shadowing_variance^2*K
    
      %%%calculate Iyy%%%%%%%%
      n12=pathloss_constant^2*(Base_station(a,2)-yt)^2*f;
      d12=((Base_station(a,1)+xt)^2+(Base_station(a,2)+yt)^2)^2;
      n22= 2*pathloss_constant*(Base_station(a,2)-yt)^2*(RSS_measured-pathloss_constant*log(d1));
      d222=((Base_station(a,1)-xt)^2+(Base_station(a,2)-yt)^2)
      n13=2*pathloss_constant*(RSS_measured-pathloss_constant*log(d1))
      d13=((Base_station(a,1)-xt)^2+(Base_station(a,2)-yt)^2)
      Iyy=(-1)*1/shadowing_variance^2*sum(n12/d12+n22/d222+n13/d13)
      
      
      %%Calculating Ixy%%%%
       n31=pathloss_constant^2*(Base_station(a,1)-xt)*(Base_station(a,2)-yt)*f;
      d31=((Base_station(a,1)+xt)^2+(Base_station(a,2)+yt)^2)^2;
      n32= 2*pathloss_constant*(Base_station(a,2)-yt)*(Base_station(a,1)-xt)*(RSS_measured-pathloss_constant*log(d1));
      d32=((Base_station(a,1)-xt)^2+(Base_station(a,2)-yt)^2)
      
      Ixy=1/shadowing_variance^2*sum(n31/d31+n32/d32)
    Ixy=Ixy
    end
  R11=[Ixx Ixy; Ixy Iyy]
 
  %RSS_measured_shadowing=add_shadowing( RSS_measured, R11);
 % Calculate DRSS at 4 different verifiers (From true locations)
for a = 1:N-1
     %   d_y_m(a) = RSS_measured_shadowing(a) - RSS_measured_shadowing(N);
         d_y_m(a) = RSS_measured(a) - RSS_measured(N);
end
 %%%%%%%%%%%  RSS from a claimed location which could be a legitmate or malicious%%%%
  for i = 1 : size (Base_station,1)
      d2(i) = norm (claimed_location - Base_station(i,:));
      PL_d_dB(i) = PL_Ref_dB + 10*pathloss_constant*log10(d2(i)/do);                 % no shadowing
      RSS_claimed(i) = transmit_power - PL_d_dB(i) ;             % RSS measured from claimed distance
   
  end
  mean_u=RSS_claimed
  %RSS_claimed_shadowing=add_shadowing(RSS_claimed, R11);
for b = 1:N-1
        %d_y(b) = RSS_claimed_shadowing(b) - RSS_claimed_shadowing(N);
        d_y(b) = RSS_claimed(b) - RSS_claimed(N);
end
   %%%%%%%%%%%%%%%%%%%%%%%%% RSS from locations r meters away (optimal attack location)from the true location%%%%%%%%%% 
  for i=1 :size(Base_station,1)
       d_1(i) = norm (Malicious_loc - Base_station(i,:))
      PL_d(i) = PL_Ref_dB + 10*pathloss_constant*log10(d_1(i)/do)                 % no shadowing
      RSS_m(i) = transmit_power - PL_d(i)     %% RSS measured from r meters away from claimed location
  end
  RSS_m;
  %RSS_Malicious_shadowing=add_shadowing(RSS_m, R11);
for c = 1:N-1
       % d_y(c) =  RSS_Malicious_shadowing(c) -  RSS_Malicious_shadowing(N);
            d_y(c) =  RSS_m(c) -  RSS_m(N);
end

 % Constructing the (N-1) dimensional covariance matrix D
    for m = 1:N-3
        for n = 1:N-3
            D(m,n) = R11(N-2,N-2) + R11(m,n) - R11(m,N-2) - R11(n,N-2);
        end
    end
    

%%%%%%%%%%%%%%%% RSS From claimed and True locations and attack locations from Base station after GPS realization faction%%%%%%
  % Obtaining GPS Realization Factor from a guassain distribution through

GPS_Realization_Factor_x = standard_deviation * randn(1,1);
GPS_Realization_Factor_y = standard_deviation * randn(1,1);

estimated_loc = [];

for hh=1:length(GPS_Realization_Factor_x)
    for kk=1:length(Base_station)
        Base_station_GPS(kk,1) = Base_station(kk,1)-GPS_Realization_Factor_x(hh);
        Base_station_GPS(kk,2) = Base_station(kk,2)-GPS_Realization_Factor_y(hh);
    end
end
%%%%%%%%%%%%%%% Calculating distances from BASE STATION_GPS to claimed %%%%%%%%%%%%%%% locations and Malicius location
 for i = 1 : size (Base_station_GPS,1)
      d3(i) = norm (claimed_location - Base_station_GPS(i,:));
      PL_d_dB(i) = PL_Ref_dB + 10*pathloss_constant*log10(d3(i)/do);                
      RSS_claimed_GPS(i) = transmit_power - PL_d_dB(i) ;
      
      d4(i) = norm (Malicious_loc - Base_station_GPS(i,:));
      PL_d_dB1(i) = PL_Ref_dB + 10*pathloss_constant*log10(d4(i)/do);                
      RSS_Malicious_GPS(i) = transmit_power - PL_d_dB1(i) ;          
 end
 
for m = 1:N-1
        d_u(m) = RSS_claimed_GPS(m) - RSS_claimed_GPS(N);
end

for n = 1:N-1
        d_v(n) = RSS_Malicious_GPS(n) - RSS_Malicious_GPS(N);
end

plot(claimed_location(:,1), claimed_location(:,2), '^', 'markersize', 10,'markerfacecolor', 'k');
hold on
plot(Base_station(:,1), Base_station(:,2), 'o', 'markersize', 10,'markerfacecolor', 'b');
plot(Legitimate_loc(3,1),Legitimate_loc(3,2),'square', 'markersize', 10, 'markerfacecolor', 'r');
 
plot(Legitimate_loc(1:1,1),Legitimate_loc(1:1,2),'diamond', 'markersize', 10, 'markerfacecolor', 'b');
plot(Legitimate_loc(3,1),Legitimate_loc(3,2),'square', 'markersize', 10, 'markerfacecolor', 'r');
plot(Legitimate_loc(4,1),Legitimate_loc(4,2),'square', 'markersize', 10, 'markerfacecolor', 'r');
plot(Legitimate_loc(7,1),Legitimate_loc(7,2),'square', 'markersize', 10, 'markerfacecolor', 'r');
plot(Legitimate_loc(10,1),Legitimate_loc(10,2),'square', 'markersize', 10, 'markerfacecolor', 'r');
plot(Legitimate_loc(20,1),Legitimate_loc(20,2),'square', 'markersize', 10, 'markerfacecolor', 'r');
plot(Legitimate_loc(30,1),Legitimate_loc(30,2),'square', 'markersize', 10, 'markerfacecolor', 'r');
plot(Legitimate_loc(15,1),Legitimate_loc(15,2),'square', 'markersize', 10, 'markerfacecolor', 'r');
plot(Legitimate_loc(28,1),Legitimate_loc(28,2),'square', 'markersize', 10, 'markerfacecolor', 'r');
plot(Legitimate_loc(25,1),Legitimate_loc(25,2),'square', 'markersize', 10, 'markerfacecolor', 'r');

% plot(Malicious_loc(1:1,1), Malicious_loc(1:1,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Legitimate_loc(1:1,1),Legitimate_loc(1:1,2),'diamond', 'markersize', 10, 'markerfacecolor', 'b');
% plot(Malicious_loc(5,1), Malicious_loc(5,2),   'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Malicious_loc(6,1), Malicious_loc(6,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Malicious_loc(7,1), Malicious_loc(7,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Malicious_loc(8,1), Malicious_loc(8,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Malicious_loc(10,1), Malicious_loc(10,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Malicious_loc(12,1), Malicious_loc(12,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Malicious_loc(11,1), Malicious_loc(11,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Malicious_loc(13,1), Malicious_loc(13,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Malicious_loc(14,1), Malicious_loc(14,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Malicious_loc(16,1), Malicious_loc(16,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Malicious_loc(17,1), Malicious_loc(17,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Malicious_loc(4,1), Malicious_loc(4,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Malicious_loc(3,1), Malicious_loc(3,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Malicious_loc(22,1), Malicious_loc(22,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');
% plot(Malicious_loc(23,1), Malicious_loc(23,2), 'square', 'markersize', 10, 'markerfacecolor', 'r');

draw_elipse(claimed_location, R11, claimed_location)


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Random Test start%%%%%%%%%%%%%%%%%%%
 
  

function pdf_value = mv_gaussian_pdf(x,u,C)
% Description
%
% Syntax: pdf_value = mv_gaussian_pdf(x,u,C)
% x: multi-variables vector 1*N
% u: mean vector 1*N
% C: Covariance Matrix N*N
% k: Dimension
x=x';
u=u';
k=length(x);
pdf_value = (1/sqrt ( (2*pi) ^ k * det(C) )) * exp (-1/2 * ( x - u)' * inv(C) * (x - u) );
    
end
%%%%% Calculate RSS after adding shadowing%%%%%%%%%%%%%%%%%%%%%%

function  [RSS]= add_shadowing(u,R11)

% Calculate the correlated shadowing and add it to the RSS mean valu
RSS_mean = u';
num = length(u);
correlated_shadowing = 5 * mvnrnd ( zeros (1,num) , R11 );
RSS = RSS_mean + correlated_shadowing;
end
function draw_elipse(mean, covar, datapoint)
    m = mean; % mean value
    n = length(m);
    covariance = covar;   %error in x variable and y variables
    
    x = datapoint;   %data point with some error 
    
    y1=x(:,1);
    y2=x(:,2);
    data = [y1 y2];

    % Calculate the eigenvectors and eigenvalues
    
    [eigenvec, eigenval ] = eig(covariance);

    % Get the index of the largest eigenvector
    [largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
    largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

    % Get the largest eigenvalue
    largest_eigenval = max(max(eigenval));

    % Get the smallest eigenvector and eigenvalue
    if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
    else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
    end

    % Calculate the angle between the x-axis and the largest eigenvector
    phi = atan2(largest_eigenvec(2), largest_eigenvec(1));

    % This angle is between -pi and pi.
    % Let’s shift it such that the angle is between 0 and 2pi
    if(phi < 0)
    phi = phi + 2*pi;
    end

    % Get the coordinates of the data mean

    % Get the 95% confidence interval error ellipse
    %chisquare_val = 2.4477;
    alpha = 0.95;
    b = chi2inv(alpha, n);

    % Get the 33% confidence interval error ellipse
    %chisquare_val = 0.8010;
   alpha1 = 0.33;
    b1 = chi2inv(alpha1, n);  %n is 2DoF here i.e., 2 unknowns
 % Get the 66% confidence interval error ellipse
    alpha2 = 0.66;
    b2 = chi2inv(alpha2, n);  %n is 2DoF here i.e., 2 unknowns

    theta = linspace(0,2*pi,100);  % one row   100 col (0 to 6.28),

    X0=m(1);
    Y0=m(2);
    a=sqrt(b*largest_eigenval);  %sigma_x  x-axes lenght
    b=sqrt(b*smallest_eigenval);  %sigm_y  y_axes lenght

    a1=sqrt(b1*largest_eigenval);  %2sigma_x  x-axes lenght
   b1=sqrt(b1*smallest_eigenval);  %2sigm_y  y_axes lenght

    a2=sqrt(b2*largest_eigenval);  %3sigma_x  x-axes lenght
    b2=sqrt(b2*smallest_eigenval);  %3sigm_y  y_axes lenght

    %Define a rotation matrix
    R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
    Q=[ a*cos( theta);b*sin( theta)]';

    Q1=[ a1*cos( theta);b1*sin( theta)]';
    Q2=[ a2*cos( theta);b2*sin( theta)]';


    %let's rotate the ellipse to some angle phi
    r_ellipse = Q * R;  %95
    r_ellipse1 = Q1 * R;  % confidence interval is 33%
    r_ellipse2 = Q2 * R;  % confidence interval is 66%


    % Draw the error ellipse
    plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'g','linewidth',1.5)
    hold on;
    plot(r_ellipse1(:,1) + X0,r_ellipse1(:,2) + Y0,'b','linewidth',1.5)
    plot(r_ellipse2(:,1) + X0,r_ellipse2(:,2) + Y0,'k','linewidth',1.5)
    
   
   % Plotting  the original data
   % plot(data(:,1), data(:,2), '+');  %plots data points randomly
    
    %g=plot(m(1),m(2),'o','markersize', 8);  % plots the mean point
    %set(g,'MarkerEdgeColor','r','MarkerFaceColor','r')
     
    hold on;

    %Plot the eigenvectors
    quiver(X0, Y0, largest_eigenvec(1)*sqrt(largest_eigenval), largest_eigenvec(2)*sqrt(largest_eigenval), '-m', 'LineWidth',2);
    quiver(X0, Y0, smallest_eigenvec(1)*sqrt(smallest_eigenval), smallest_eigenvec(2)*sqrt(smallest_eigenval), '-g', 'LineWidth',2);
    hold on;
    h= legend({'Claimed Position','Base Stations', 'Estimated Position', 'True Position'}, 'fontsize', 8 );
    rect = [0.2, 0.7, 0.2, 0.2];
    set(h, 'Position', rect)
     %Set the axis labels
    xlabel(' Xdata ','interpreter','latex','fontsize',10)
    ylabel('Ydata ','interpreter','latex','fontsize',10)
    set(gca, 'fontsize',10);


end

