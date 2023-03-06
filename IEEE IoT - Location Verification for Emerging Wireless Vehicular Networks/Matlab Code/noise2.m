function [sigma] = noise2(standard_deviation, number_of_verifiers)

% '1' is the row and the 'number_of_verifiers' are the columns.

std_nlos = standard_deviation;
     sigma= exprnd(std_nlos,1,number_of_verifiers);
end
