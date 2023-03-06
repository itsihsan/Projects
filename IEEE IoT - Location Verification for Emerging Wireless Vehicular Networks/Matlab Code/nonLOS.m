% Below function generates NLoS using a truncated exponential
% distribution.
function [bias,std_i,var_i] = nonLOS(length_U_Good, mu, number_of_verifiers)

% mu is the standard deviation term.
% '1' is the row and the 'number_of_verifiers' are the columns.

std_i = mu;
% row_i = 1/std_i;
var_i = (std_i)^2;

for i=1:length_U_Good
    bias(i,:) = exprnd(std_i,1,number_of_verifiers);
end
end

% % Below function generates NLoS using a truncated GUASSIAN distribution.
% function [sigma2,std_i,var_i] = nonLOS(length_U_Good, mu, number_of_verifiers)
% 
% std_i = mu;
% var_i = (std_i)^2;
% for i=1:length_U_Good
%     sigma2(i,:) = abs(std_i * randn(1, number_of_verifiers));
%     %     sigma2(i,:) = std_i * randn(1, number_of_verifiers);
% end

%%% Below function generates NLoS using a UNIFORM distribution in a range.
% function [bias,std_i,var_i] = nonLOS(length_U_Good, mu, number_of_verifiers)
% 
% % mu is the standard deviation term.
% % '1' is the row and the 'number_of_verifiers' are the columns.
% 
% std_i = 0;
% var_i = (std_i)^2;
% minimum = 0;
% maximum = 0;
% for i=1:length_U_Good
%     bias(i,:) = minimum + (maximum - minimum)*rand(1,number_of_verifiers);
% end
