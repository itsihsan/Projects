function [verifier_pre, verifier, claimed_location , optimal_attack_location, U, Delta_U, R, D] = Optimum_DRSS(standard_deviation, N_verifiers, distant_constaint)

%[verifier_pre, verifier, claimed_location , optimal_attack_location, U, Delta_U, R, D] = Optimum_DRSS(standard_deviation, N_verifiers, distant_constaint)
pathloss_constant = 3;
transmit_power = 30; % in dBm
shadowing_variance = 5;  % Shadowing Variance [Sigma]
Dc = 400; % Constant indicating degree of shadowing correlation in a
% specific enviornment. Used for evaluating covariance matrix R
do = 1; % Reference Distance

lembda = (3*10^8)/(900*10^6); % lembda = c/f for 900 MHz band
p = transmit_power + 20 * log10(lembda/(4*pi*do));
N = N_verifiers;
% Verifier True Locations (Randomly Generated)
x = 150 * ( - 1 + 2 * rand ( 1 , N_verifiers ) );
y = 150 * ( - 1 + 2 * rand ( 1 , N_verifiers ) );

verifier_pre =  [ x ; y ]';

% Claimed Location for malicious or legitimate vehicle (Randomly Generated)
xc = 150 * ( - 1 + 2 * rand ( 1 , 1 ) );
yc = 150 * ( - 1 + 2 * rand ( 1 , 1 ) ); 
claimed_location = [xc yc];

% Obtaining GPS Realization Factor from a guassain distribution through
% the standard_deviation value.
GPS_Realization_Factor_x = standard_deviation * randn(1,N);
GPS_Realization_Factor_y = standard_deviation * randn(1,N);
x_1 = x - GPS_Realization_Factor_x;
y_1 = y - GPS_Realization_Factor_y;
verifier = [x_1 ; y_1]';

optimal_attack_location = [];

 for hh=1:1
%     for kk=1:length(verifier_pre)
%         verifier(kk,1) = verifier_pre(kk,1)-GPS_Realization_Factor_x(hh);
%         verifier(kk,2) = verifier_pre(kk,2)-GPS_Realization_Factor_y(hh);
%     end
    
    N = length(verifier); % Number of verifiers
    
    for n = 1:N
        % Calculating distances of verifiers from Claimed Location
        d2t(n,:) = sqrt((claimed_location(1,1) - verifier(n,1))^2 + (claimed_location(1,2) - verifier(n,2))^2);
        % Calculating Vector U
        U(n,:) = p - 10 * pathloss_constant * log10(d2t(n,1)/do);
    end
    
    % Calculating Delta_U Mean Vector
    for m = 1:N-1
        Delta_U1(m) = U(m) - U(N);
    end
    Delta_U = Delta_U1';
    
    % Constructing the covariance matrix R
    for a = 1:N
        for b = 1:N
            d(a,b) = sqrt((verifier(a,1) - verifier(b,1))^2 + (verifier(a,2) - verifier(b,2))^2 );
            R0(a,b) = exp(-d(a,b)/Dc*log(2));
        end
    end
    
    R = shadowing_variance * R0;  % Correlated covariance matrix of shadowing noise
    
    % Constructing the (N-1) dimensional covariance matrix D
    for m = 1:N-1
        for n = 1:N-1
            D(m,n) = R(N,N) + R(m,n) - R(m,N) - R(n,N);
        end
    end
    
    %x0 = -50:(150-(-50))/2100:150;
    %y0 = -50:(100-(-50))/2100:100;
    x0 = -150:(150-(-150))/100:150;
    y0 = -150:(150-(-150))/100:150;
    r = distant_constaint; %the minimum distance between true location and claimed location
    N_0 = length(x0);
    x0_p = x0'*ones(1,N_0);
    y0_p = ones(N_0,1)*y0;
    dct = sqrt((x0_p - claimed_location(1,1)).^2 + (y0_p - claimed_location(1,2)).^2);
    I = find(dct>=r);
    x1 = x0_p(I);
    y1 = y0_p(I);
    N_1 = length(x1);
    
    for i = 1:N_1
        for j = 1:N
            d2t = sqrt((verifier(j,1) - x1(i))^2 + (verifier(j,2) - y1(i))^2);
            V(j,1) = p - 10* pathloss_constant * log10(d2t/do);
        end
        
        for m = 1:N-1
            Delta_V1(m) = V(m) - V(N);
        end
        
        Delta_V = Delta_V1';
        KL_Divergance_DRSS(i) = (1/2) * (Delta_V - Delta_U)' * D^(-1) * (Delta_V - Delta_U);
    end
    [~ , index] = min(KL_Divergance_DRSS);
    x_value(hh,1) = x1(index);
    y_value(hh,1)= y1(index);
end

% Optimum Attack Location, where detecion of the malicious vehicle will be
% minimum, after introducing GPS_Realization_Factor to the verifier location is given below

optimal_attack_location =[optimal_attack_location;x_value,y_value];

%U
%Delta_U
%r

%clf('reset');

%for hh=1:length(verifier)
    %plot(verifier(hh,1),verifier(hh,2),'r*','MarkerSize',15);
    %hold on
%end
%plot (claimed_location(1) , claimed_location(2),'--d','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',[0,0,0]);
%plot(x_value(1),y_value(1),'ko','MarkerSize',10,'MarkerEdgeColor','g','MarkerFaceColor',[0,1,0]);
%xlabel('X-Coordinate');
%ylabel('Y-Coordinate');
%title(['Graph Showing Optimal Attack Location at Std. Deviation ', num2str(standard_deviation)]);
end