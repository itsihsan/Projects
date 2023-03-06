function [sigma] = noise(standard_deviation, number_of_verifiers)

sigma = standard_deviation * randn(1, number_of_verifiers);

end