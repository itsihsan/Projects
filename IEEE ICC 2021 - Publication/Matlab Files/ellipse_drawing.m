% Test_noise=readmatrix('5N_Test_half.csv','Range','K1:Q15000');
% test_rss = Test_noise(:,1:5);
% test_coordinates = Test_noise(:,6:7);
% clear Test_noise

noise_std = 3 ; % standard deviation of noise in dBm.
% noise_std = 3.16227766; % standard deviation of noise in mW.
gamma = 3; % free pathloss exponent.
% % % % verifier_loc =  [100,100;
% % % %     100,400;
% % % %     400,100;
% % % %     400,400;
% % % %     250,250];

% % 5 verifiers used for half/quarter trial
% verifier_loc =  [50,50;
%     50,200;
%     200,50;
%     200,200;
%     125,125];

% % 4 verifiers used for half/quarter trial
% verifier_loc =  [50,50;
%                  200,50;
%                  200,200;
%                  125,125];

% 4 verifiers - ROB task
verifier_loc =  [0,0;
    0,200;
    200,0;
    200,200];

num_examples = length(test_rss);
num_verifiers = length(verifier_loc);

% a = 10/log(exp(1)); % 10/ln(10)
a = 10/log(10); % 10/ln(10)

term_1 = (a*gamma)/(noise_std)^2;

nnInside_1 = 0;
nnInside_2 = 0;
nnInside_3 = 0;
fpInside_1 = 0;
fpInside_2 = 0;
fpInside_3 = 0;

filter = [];
for i = 1:num_examples
    for j = 1:num_verifiers
        
        d = sqrt((verifier_loc(j,1) - test_coordinates(i,1))^2 + (verifier_loc(j,2) - test_coordinates(i,2))^2);
        
        term_2_pre(1,j) = 1/(d^2);
        
        term_2_xx(1,j) = (a*gamma*((verifier_loc(j,1) - test_coordinates(i,1))^2))/(d^2);
        term_3_xx(1,j) = (test_rss(i,j) + (a*gamma*log(d))) * ( 1 - ((2/(d^2))*((verifier_loc(j,1) - test_coordinates(i,1))^2)));
        
        term_2_yy(1,j) = (a*gamma*((verifier_loc(j,2) - test_coordinates(i,2))^2))/(d^2);
        term_3_yy(1,j) = (test_rss(i,j) + (a*gamma*log(d))) * ( 1 - ((2/(d^2))*((verifier_loc(j,2) - test_coordinates(i,2))^2)));
        
        term_2_xy(1,j) = ((verifier_loc(j,1) - test_coordinates(i,1))*(verifier_loc(j,2) - test_coordinates(i,2)))/(d^4);
        term_3_xy(1,j) = ( (a * gamma) -  (2 * (test_rss(i,j)+(a*gamma*log(d)))) );
        
        clear d;
    end
    
    theta_xx = sum((term_1 * term_2_pre) .* (term_2_xx + term_3_xx));
    theta_yy = sum((term_1 * term_2_pre) .* (term_2_yy + term_3_yy));
    theta_xy = sum((term_1 * term_2_xy) .* term_3_xy);
    
    %     theta_xx = term_1 * (term_2 * term_3_xx');
    %     theta_yy = term_1 * (term_2 * term_3_yy');
    %     theta_xy = term_1 * sum(term_2_xy - term_3_xy);
    theta_yx = theta_xy;
    
    covariance_matrix = [1/theta_xx, 1/theta_xy; 1/theta_yx, 1/theta_yy];
    %     pre_covariance_matrix=[theta_xx, theta_xy; theta_yx, theta_yy];
    %     covariance_matrix = inv(pre_covariance_matrix);
    
    [eigenvec, eigenval] = eig(covariance_matrix);
    
    % Get the largest eigenvalue
    largest_eigenval = max(max(eigenval));
    
    % Get the index of the largest eigenvector
    [largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
    largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);
    
    % Get the smallest eigenvector and eigenvalue
    if(largest_eigenvec_ind_c == 1)
        smallest_eigenval = max(eigenval(:,2));
        smallest_eigenvec = eigenvec(:,2);
    else
        smallest_eigenval = max(eigenval(:,1));
        smallest_eigenvec = eigenvec(1,:);
    end
    
    if smallest_eigenval==0
        continue
    else
        filter = [filter;i];
        % Calculate the angle between the x-axis and the largest eigenvector
        phi = atan2(largest_eigenvec(2), largest_eigenvec(1));
        
        % This angle is between -pi and pi.
        % Let’s shift it such that the angle is between 0 and 2pi
        if(phi < 0)
            phi = phi + 2*pi;
        end
        
        % Get the 39% confidence interval error ellipse
        nn = length(test_coordinates(i,:));
        alpha1 = 0.39;
        b1 = chi2inv(alpha1, nn);
        semi_major_1=sqrt(b1*largest_eigenval);  % Semi Major Axis - 1 sigma ellipse
        semi_minor_1=sqrt(b1*smallest_eigenval); % Semi Minor Axis - 1 sigma ellipse
        
        %Define a rotation matrix
        R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
        theta = linspace(0,2*pi,100);  % one row   100 col (0 to 6.28),
        Q1 = [semi_major_1*cos(theta) ; semi_minor_1*sin(theta)]';
        %let's rotate the ellipse to some angle phi
        r_ellipse_1 = Q1 * R;  % 33
        
        % Get the 66% confidence interval error ellipse
        %chisquare_val = 0.8010;
        alpha2 = 0.66;
        b2 = chi2inv(alpha2, nn);  %n is 2DoF here i.e., 2 unknowns
        semi_major_2 = sqrt(b2*largest_eigenval);  % Semi Major Axis - 2 sigma ellipse
        semi_minor_2 = sqrt(b2*smallest_eigenval); % Semi Minor Axis - 2 sigma ellipse
        Q2 = [semi_major_2*cos(theta) ; semi_minor_2*sin(theta)]';
        r_ellipse_2 = Q2 * R;  % confidence interval is 66%
        
        alpha3 = 0.95;
        b3 = chi2inv(alpha3, nn);  %n is 2DoF here i.e., 2 unknowns
        semi_major_3 = sqrt(b3*largest_eigenval);  % Semi Major Axis - 3 sigma ellipse
        semi_minor_3 = sqrt(b3*smallest_eigenval); % Semi Minor Axis - 3 sigma ellipse
        Q3=[ semi_major_3*cos(theta) ; semi_minor_3*sin(theta)]';
        r_ellipse_3 = Q3 * R;  % confidence interval is 95%
        
        % Draw the error ellipses
        mean_point = test_coordinates(i,:); % m is the mean value of the ellipse
        x_mean = mean_point(1);
        y_mean = mean_point(2);
        g=plot(mean_point(1),mean_point(2),'o','markersize', 8);  % plots the mean point
        set(g,'MarkerEdgeColor','r','MarkerFaceColor','r')
        hold on;
        
        % Plotting 1st confidence Ellipse - 1 sigma
        plot(r_ellipse_1(:,1) + x_mean,r_ellipse_1(:,2) + y_mean,'g','linewidth',1.0)
        
        fp_estimated = fpEstLoc(i,:);
        nn_estimated = nnEstLoc_50nn(i,:);
        hh=plot(fp_estimated(1),fp_estimated(2),'*','markersize', 8);  % plots the fp est loc
        set(hh,'MarkerEdgeColor','b','MarkerFaceColor','b')
        ii=plot(nn_estimated(1),nn_estimated(2),'d','markersize', 8);  % plots the nn est loc
        set(ii,'MarkerEdgeColor','g','MarkerFaceColor','g')
        
        %%% Analysis for point inside/outside the 1st confidence Ellipse - 1 sigma
        % Checking if FP estimated location is within the 1st confidence ellipse.
        fp_11 = (((cosd(rad2deg(phi))*(fp_estimated(1)-mean_point(1)))+(sind(rad2deg(phi))*(fp_estimated(2)-mean_point(2))))^2/semi_major_1^2);
        fp_12 = (((sind(rad2deg(phi))*(fp_estimated(1)-mean_point(1)))-(cosd(rad2deg(phi))*(fp_estimated(2)-mean_point(2))))^2/semi_minor_1^2);
        fp_total_1 = fp_11+fp_12;
        if fp_total_1 <= 1
            fprintf('The FP estimated location is inside the 1st confidence ellipse. \n')
            fpInside_1 = fpInside_1 + 1;
        else
            fprintf('The FP estimated location is outside the 1st confidence ellipse. \n')
        end
        % Checking if NN estimated location is within the 1st confidence ellipse.
        nn_11 = (((cosd(rad2deg(phi))*(nn_estimated(1)-mean_point(1)))+(sind(rad2deg(phi))*(nn_estimated(2)-mean_point(2))))^2/semi_major_1^2);
        nn_12 = (((sind(rad2deg(phi))*(nn_estimated(1)-mean_point(1)))-(cosd(rad2deg(phi))*(nn_estimated(2)-mean_point(2))))^2/semi_minor_1^2);
        nn_total_1 = nn_11+nn_12;
        if nn_total_1 <= 1
            fprintf('The NN estimated location is inside the 1st confidence ellipse. \n')
            nnInside_1 = nnInside_1 + 1;
        else
            fprintf('The NN estimated location is outside the 1st confidence ellipse. \n')
        end
        
        % Plotting the 2nd confidence Ellipse - 2 sigma
        plot(r_ellipse_2(:,1) + x_mean,r_ellipse_2(:,2) + y_mean,'b-.','linewidth',1.0)
        %%% Analysis for point inside/outside the 2nd confidence Ellipse - 2 sigma
        % Checking if FP estimated location is within the 2nd confidence ellipse.
        fp_21 = (((cosd(rad2deg(phi))*(fp_estimated(1)-mean_point(1)))+(sind(rad2deg(phi))*(fp_estimated(2)-mean_point(2))))^2/semi_major_2^2);
        fp_22 = (((sind(rad2deg(phi))*(fp_estimated(1)-mean_point(1)))-(cosd(rad2deg(phi))*(fp_estimated(2)-mean_point(2))))^2/semi_minor_2^2);
        fp_total_2 = fp_21+fp_22;
        if fp_total_2 <= 1
            fprintf('The FP estimated location is inside the 2nd confidence ellipse. \n')
            fpInside_2 = fpInside_2 + 1;
        else
            fprintf('The FP estimated location is outside the 2nd confidence ellipse. \n')
        end
        % Checking if NN estimated location is within the 2nd confidence ellipse.
        nn_21 = (((cosd(rad2deg(phi))*(nn_estimated(1)-mean_point(1)))+(sind(rad2deg(phi))*(nn_estimated(2)-mean_point(2))))^2/semi_major_2^2);
        nn_22 = (((sind(rad2deg(phi))*(nn_estimated(1)-mean_point(1)))-(cosd(rad2deg(phi))*(nn_estimated(2)-mean_point(2))))^2/semi_minor_2^2);
        nn_total_2 = nn_21+nn_22;
        if nn_total_2 <= 1
            fprintf('The NN estimated location is inside the 2nd confidence ellipse. \n')
            nnInside_2 = nnInside_2 + 1;
        else
            fprintf('The NN estimated location is outside the 2nd confidence ellipse. \n')
        end
        
        % Plotting the 3rd confiednce Ellipse - 3 sigma
        plot(r_ellipse_3(:,1) + x_mean,r_ellipse_3(:,2) + y_mean,'k--','linewidth',1.0)
        %%% Analysis for point inside/outside the 3rd confidence Ellipse - 3 sigma
        % Checking if FP estimated location is within the 3rd confidence ellipse.
        fp_31 = (((cosd(rad2deg(phi))*(fp_estimated(1)-mean_point(1)))+(sind(rad2deg(phi))*(fp_estimated(2)-mean_point(2))))^2/semi_major_3^2);
        fp_32 = (((sind(rad2deg(phi))*(fp_estimated(1)-mean_point(1)))-(cosd(rad2deg(phi))*(fp_estimated(2)-mean_point(2))))^2/semi_minor_3^2);
        fp_total_3 = fp_31+fp_32;
        if fp_total_3 <= 1
            fprintf('The FP estimated location is inside the 3rd confidence ellipse. \n')
            %             keyboard
            fpInside_3 = fpInside_3 + 1;
        else
            fprintf('The FP estimated location is outside the 3rd confidence ellipse. \n')
        end
        % Checking if NN estimated location is within the 3rd confidence ellipse.
        nn_31 = (((cosd(rad2deg(phi))*(nn_estimated(1)-mean_point(1)))+(sind(rad2deg(phi))*(nn_estimated(2)-mean_point(2))))^2/semi_major_3^2);
        nn_32 = (((sind(rad2deg(phi))*(nn_estimated(1)-mean_point(1)))-(cosd(rad2deg(phi))*(nn_estimated(2)-mean_point(2))))^2/semi_minor_3^2);
        nn_total_3 = nn_31+nn_32;
        if nn_total_3 <= 1
            fprintf('The NN estimated location is inside the 3rd confidence ellipse. \n')
            %             keyboard
            nnInside_3 = nnInside_3 + 1;
        else
            fprintf('The NN estimated location is outside the 3rd confidence ellipse. \n')
        end
        
        clear j; clear term_2_pre; clear term_2_xx; clear term_3_xx; clear term_2_yy; clear term_3_yy; clear term_2_xy; clear term_3_xy; clear theta_xx; clear theta_xy; clear theta_yx; clear theta_yy;
        %         clear pre_covariance_matrix;
        clear covariance_matrix; clear eigenvec; clear eigenval; clear largest_eigenval; clear largest_eigenvec_ind_c; clear r; clear largest_eigenvec;
        clear smallest_eigenval; clear smallest_eigenvec; clear phi; clear alpha1; clear semi_major_1; clear b1; clear semi_minor_1; clear R; clear theta; clear Q1; clear r_ellipse_1;
        clear alpha2; clear semi_major_2; clear b2; clear semi_minor_2; clear Q2; clear r_ellipse_2;
        clear alpha3; clear semi_major_3; clear b3; clear semi_minor_3; clear Q3; clear r_ellipse_3;
        clear mean_point; clear x_mean; clear y_mean; clear g; clear hh; clear ii; clear nn;
        clear fp_estimated; clear fp_11; clear fp_12; clear fp_21; clear fp_22; clear fp_31; clear fp_32; clear fp_total_1; clear fp_total_2; clear fp_total_3;
        clear nn_estimated; clear nn_11; clear nn_12; clear nn_21; clear nn_22; clear nn_31; clear nn_32; clear nn_total_1; clear nn_total_2; clear nn_total_3;
        clf;
    end
end
