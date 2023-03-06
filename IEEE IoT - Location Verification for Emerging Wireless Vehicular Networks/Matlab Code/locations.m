%=========================================================================
% Modified Version to resemble the DRSS paper. The verifiers and the true
% and claimed locations lie in the same space.
%=========================================================================
function [good_user, bad_user] = locations(good, bad)

% Generate x and y for the good user between 250 and 750 on x-axis and
% between 0 and 500 on y-axis
g_x_a = 250;      % Point a for defining range on x-axis
g_x_b = 750;      % Point b for defining range on x-axis
good_user_x = (g_x_b - g_x_a) .* rand(1,good) + g_x_a ;

g_y_a = 0;      % Point a for defining range on y-axis
g_y_b = 500;     % Point b for defining range on y-axis
good_user_y = (g_y_b - g_y_a) .* rand(1,good) + g_y_a ;

good_user = [good_user_x' , good_user_y'];

% Generate x and y coordinates for 'bad user' in a 500 X 500 area.
v_x_a = 250;      % Point a for defining range on x-axis
v_x_b = 750;      % Point b for defining range on x-axis
bad_x = (v_x_b - v_x_a) .* rand(1,bad) + v_x_a;

v_y_a = 0;      % Point a for defining range on y-axis
v_y_b = 500;    % Point b for defining range on y-axis
bad_y = (v_y_b - v_y_a) .* rand(1,bad) + v_y_a ;
bad_user = [bad_x' , bad_y'];

end


% % % % % % % % % %=========================================================================
% % % % % % % % % % Modified Version to resemble the DRSS paper. The verifiers and the true
% % % % % % % % % % and claimed locations lie in the same space.
% % % % % % % % % %=========================================================================
% % % % % % % % % function [verifier , good_user] = locations(number_of_verifiers)
% % % % % % % % % % function [good_user] = locations(number_of_verifiers)
% % % % % % % % % % function [verifier , good_user , bad_user] = locations(number_of_verifiers)
% % % % % % % % % % verifier_x = 250  * ( - 1 + 2 * rand ( 1 , number_of_verifiers ) );
% % % % % % % % % % verifier_y = 15   * ( - 1 + 2 * rand ( 1 , number_of_verifiers ) );
% % % % % % % % % 
% % % % % % % % % % number_of_verifiers_1 = 3;
% % % % % % % % % % number_of_verifiers_2 = 3;
% % % % % % % % % number_of_verifiers_1 = abs(number_of_verifiers/2);
% % % % % % % % % number_of_verifiers_2 = number_of_verifiers - number_of_verifiers_1;
% % % % % % % % % 
% % % % % % % % % % Generate x and y coordinates for 'number of verifiers' in a 100 X 30 in area - 1.
% % % % % % % % % v_x_a = 0;      % Point a for defining range on x-axis
% % % % % % % % % v_x_b = 150;    % Point b for defining range on x-axis
% % % % % % % % % verifier_x_1 = (v_x_b - v_x_a) .* rand(1, number_of_verifiers_1) + v_x_a;
% % % % % % % % % 
% % % % % % % % % v_y_a = 0;      % Point a for defining range on y-axis
% % % % % % % % % v_y_b = 30;    % Point b for defining range on y-axis
% % % % % % % % % verifier_y_1 = (v_y_b - v_y_a) .* rand(1, number_of_verifiers_1) + v_y_a ;
% % % % % % % % % verifier_1 = [verifier_x_1' , verifier_y_1'];
% % % % % % % % % 
% % % % % % % % % % Generate x and y coordinates for 'number of verifiers' in a 100 X 30 in area - 2.
% % % % % % % % % v_x_c = 450;      % Point a for defining range on x-axis
% % % % % % % % % v_x_d = 600;      % Point b for defining range on x-axis
% % % % % % % % % verifier_x_2 = (v_x_d - v_x_c) .* rand(1, number_of_verifiers_2) + v_x_c;
% % % % % % % % % 
% % % % % % % % % v_y_c = 0;        % Point a for defining range on y-axis
% % % % % % % % % v_y_d = 30;       % Point b for defining range on y-axis
% % % % % % % % % verifier_y_2 = (v_y_d - v_y_c) .* rand(1, number_of_verifiers_1) + v_y_c;
% % % % % % % % % verifier_2 = [verifier_x_2' , verifier_y_2'];
% % % % % % % % % 
% % % % % % % % % verifier = [verifier_1 ; verifier_2];
% % % % % % % % % 
% % % % % % % % % % good_user_x = 250 * ( - 1 + 2 * rand ( 1 , 1 ) );
% % % % % % % % % % good_user_y =  15 * ( - 1 + 2 * rand ( 1 , 1 ) );
% % % % % % % % % 
% % % % % % % % % % Generate x and y for the good user between 200 and 400 on x-axis and
% % % % % % % % % % between 0 and 30 on y-axis
% % % % % % % % % g_x_a = 300;      % Point a for defining range on x-axis
% % % % % % % % % g_x_b = 450;      % Point b for defining range on x-axis
% % % % % % % % % good_user_x = (g_x_b - g_x_a) .* rand(1,1) + g_x_a ;
% % % % % % % % % 
% % % % % % % % % g_y_a = 0;      % Point a for defining range on y-axis
% % % % % % % % % g_y_b = 30;     % Point b for defining range on y-axis
% % % % % % % % % good_user_y = (g_y_b - g_y_a) .* rand(1,1) + g_y_a ;
% % % % % % % % % 
% % % % % % % % % good_user = [good_user_x , good_user_y];
% % % % % % % % % 
% % % % % % % % % % % % % % % % % % % Generate x and y for the BAD user between 500 and 550 on x-axis and
% % % % % % % % % % % % % % % % % % % between 0 and 20 on y-axis
% % % % % % % % % % % % % % % % % % b_x_a = 500;    % Point a for defining range on x-axis
% % % % % % % % % % % % % % % % % % b_x_b = 550;    % Point b for defining range on x-axis
% % % % % % % % % % % % % % % % % % b_y_a = 0;      % Point a for defining range on y-axis
% % % % % % % % % % % % % % % % % % b_y_b = 20;     % Point b for defining range on y-axis
% % % % % % % % % % % % % % % % % % bad_user_x = (b_x_b - b_x_a) .* rand(1,1) + b_x_a ;
% % % % % % % % % % % % % % % % % % bad_user_y = (b_y_b - b_y_a) .* rand(1,1) + b_y_a ;
% % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % bad_user = [bad_user_x , bad_user_y];
% % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % Noise(sigma) Generation from a guassian function
% % % % % % % % % % % % % % % % % % % sigma = standard_deviation * randn(1,1);
% % % % % % % % % 
% % % % % % % % % end