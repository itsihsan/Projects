% Test_noise=readmatrix('5N_Test_half.csv','Range','K1:Q15000');
% test_rss = Test_noise(:,1:5);
% test_coordinates = Test_noise(:,6:7);
% clear Test_noise

gamma = 3; % free pathloss exponent.
noise_std = 4; % standard deviation of noise in dBm.

% 5 verifiers used for half/quarter trial
verifier_loc =  [50,50;
    50,200;
    200,50;
    200,200;
    125,125];

% % % % % 4 verifiers - ROB task
% % % % verifier_loc =  [0,0; 0,200; 200,0; 200,200];

num_examples = length(test_rss);
num_verifiers = length(verifier_loc);

a = 10/log(10); % 10/ln(10)
term_1 = (a*gamma).^2/(noise_std)^2;

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
        term_2_xx(1,j) = (1/(d^4))*(verifier_loc(j,1) - test_coordinates(i,1))^2;
        term_2_yy(1,j) = (1/(d^4))*(verifier_loc(j,2) - test_coordinates(i,2))^2;
        term_2_xy(1,j) = (1/(d^4))*(verifier_loc(j,1) - test_coordinates(i,1))*(verifier_loc(j,2) - test_coordinates(i,2));
        clear d;
    end
    
    theta_xx = term_1* sum(term_2_xx);
    theta_yy = term_1* sum(term_2_yy);
    theta_xy = term_1* sum(term_2_xy);
    
    theta_yx = theta_xy;
    
    pre_covariance_matrix=[theta_xx, theta_xy; theta_yx, theta_yy];
    covariance_matrix = inv(pre_covariance_matrix);
    
    % New code Start
    lambda = eig(covariance_matrix);
    [S, D] = eig(covariance_matrix);
    scale = 2.447;  %for 95% confidense ellipse from sigma1
    [a, a_i] = max(lambda); %find maximum eigenvalue
    [b, b_i] = min(lambda); %find minimum eigenvalue
    a = scale*sqrt(a); %scale according to confidence interval
    b = scale*sqrt(b);
    
    if covariance_matrix(1,1) > covariance_matrix(2,2) %resolve tilt of the ellipse
        phi = atan2(S(2,a_i),S(1,a_i));
        x_axis = a;
        y_axis= b;
    else
        phi = atan2(S(2,b_i),S(1,b_i));
        x_axis = b;
        y_axis = a;
    end
    
    if(phi < 0)
        phi = phi + 2*pi;
    end

    % New code finish
    
    %Define a rotation matrix
    R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
    theta = linspace(0,2*pi,100);  % one row   100 col (0 to 6.28),
    
    Q3=[ x_axis*cos(theta) ; y_axis*sin(theta)]';
    r_ellipse_3 = Q3 * R;  % confidence interval is 95%
    
    % Draw the error ellipses
    mean_point = test_coordinates(i,:); % m is the mean value of the ellipse
    x_mean = mean_point(1);
    y_mean = mean_point(2);
    g=plot(mean_point(1),mean_point(2),'o','markersize', 8);  % plots the mean point
    set(g,'MarkerEdgeColor','r','MarkerFaceColor','r')
    hold on;
    
    fp_estimated = fpEstLoc(i,:);
    % % % % %         nn_estimated = nnEstLoc_08nn(i,:);
    hh=plot(fp_estimated(1),fp_estimated(2),'*','markersize', 8);  % plots the fp est loc
    set(hh,'MarkerEdgeColor','b','MarkerFaceColor','b')
    
    % Plotting the 3rd confiednce Ellipse - 3 sigma
    plot(r_ellipse_3(:,1) + x_mean,r_ellipse_3(:,2) + y_mean,'k--','linewidth',1.0)
    %%% Analysis for point inside/outside the 3rd confidence Ellipse - 3 sigma
    % Checking if FP estimated location is within the 3rd confidence ellipse.
    fp_31 = (((cosd(rad2deg(phi))*(fp_estimated(1)-mean_point(1)))+(sind(rad2deg(phi))*(fp_estimated(2)-mean_point(2))))^2/x_axis^2);
    fp_32 = (((sind(rad2deg(phi))*(fp_estimated(1)-mean_point(1)))-(cosd(rad2deg(phi))*(fp_estimated(2)-mean_point(2))))^2/y_axis^2);
    fp_total_3 = fp_31+fp_32;
    if fp_total_3 <= 1
        fprintf('The FP estimated location is inside the 3rd confidence ellipse. \n')
        %             keyboard
        fpInside_3 = fpInside_3 + 1;
    else
        fprintf('The FP estimated location is outside the 3rd confidence ellipse. \n')
    end
    
    clear j; clear term_2_pre; clear term_2_xx; clear term_2_yy; clear term_2_xy; clear theta_xx; clear theta_xy; clear theta_yx; clear theta_yy;
    clear pre_covariance_matrix; clear covariance_matrix;
    clear eigenvec; clear eigenval; clear largest_eigenval; clear largest_eigenvec_ind_c; clear r; clear largest_eigenvec;
    clear smallest_eigenval; clear smallest_eigenvec; clear phi; clear alpha1; clear semi_major_1; clear b1; clear semi_minor_1; clear R; clear theta; clear Q1; clear r_ellipse_1;
    clear alpha2; clear semi_major_2; clear b2; clear semi_minor_2; clear Q2; clear r_ellipse_2;
    clear alpha3; clear semi_major_3; clear b3; clear semi_minor_3; clear Q3; clear r_ellipse_3;
    clear mean_point; clear x_mean; clear y_mean; clear g; clear hh; clear ii; clear nn;
    clear fp_estimated; clear fp_11; clear fp_12; clear fp_21; clear fp_22; clear fp_31; clear fp_32; clear fp_total_1; clear fp_total_2; clear fp_total_3;
    clear nn_estimated; clear nn_11; clear nn_12; clear nn_21; clear nn_22; clear nn_31; clear nn_32; clear nn_total_1; clear nn_total_2; clear nn_total_3;
    clear lambda;clear S;clear D;clear scale;clear a;clear a_i;clear b;clear b_i;clear x_axis;clear y_axis;
    clf;
end
