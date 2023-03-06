% function [verifier, claimed_location , attack_location, U, Delta_U, R, D] = Optimum_DRSS(distant_constraint)
function [data] = nn_data(samples, distant_constraint)

pathloss_constant = 3;
transmit_power = 30; % in dBm
do = 1; % Reference Distance
lembda = (3*10^8)/(9*10^9); % lembda = c/f for 900 MHz band
p = transmit_power + 20 * log10(lembda/(4*pi*do));
% verifier =  [0,125;0,375;1000,125;1000,375];           % 4 verifiers
verifier =  [1,1;1,250;1,500; 1000,1; 1000,250; 1000,500];           % 6 verifiers

attack_location = [];

for hh=1:samples
    clear claimed_location; clear d2t; clear U; clear x0, clear y0; clear N_0;
    clear x0_p; clear y0_p; clear dct; clear I; clear x_value; clear y_value;
    clear attack_location; clear n;
    
    [claimed_location, ~] = locations(1,1);
    N = length(verifier); % Number of verifiers
    
    for n = 1:N
        % Calculating distances of verifiers from Claimed Location
        d2t(n,:) = sqrt((claimed_location(1,1) - verifier(n,1))^2 + (claimed_location(1,2) - verifier(n,2))^2);
        % Calculating Vector U
        U(n,:) = p - 10 * pathloss_constant * log10(d2t(n,1)/do);
    end
    
    x0 = 250:(750-250)/1000:750;
    y0 = 0  :(500-0)  /1000:500;
    r = distant_constraint; %the minimum distance between true location and claimed location
    N_0 = length(x0);
    x0_p = x0'*ones(1,N_0);
    y0_p = ones(N_0,1)*y0;
    dct = sqrt((x0_p - claimed_location(1,1)).^2 + (y0_p - claimed_location(1,2)).^2);
    I = find(dct>=r & dct<r+0.01); % min(a(find(a>2&a<14)))
    x_value = x0_p(I(1));
    y_value = y0_p(I(1));
    %attack_location =[attack_location;x_value,y_value];
    attack_location =[x_value,y_value];
    
    data_pre(hh,:) = [U', claimed_location, attack_location];
    
%     clf('reset');
    for hh=1:length(verifier);
        plot(verifier(hh,1),verifier(hh,2),'r*','MarkerSize',15);
        hold on
    end
    plot (claimed_location(1) , claimed_location(2),'--d','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',[0,0,0]);
    plot(x_value(1),y_value(1),'ko','MarkerSize',10,'MarkerEdgeColor','g','MarkerFaceColor',[0,1,0]);
    xlabel('X-Coordinate');
    ylabel('Y-Coordinate');
    title('Graph Showing Optimal Attack Location at Std. Deviation');
    % title(['Graph Showing Optimal Attack Location at Std. Deviation ', num2str(standard_deviation)]);
end
data = [repmat(verifier(1,:),samples,1), repmat(verifier(2,:),samples,1), repmat(verifier(3,:),samples,1), repmat(verifier(4,:),samples,1), repmat(verifier(5,:),samples,1), repmat(verifier(6,:),samples,1), data_pre];