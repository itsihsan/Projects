function [receiver_noise] = noise(standard_deviation, n)
% standard deviation of the reciever noise is in n dBs (RSS)
receiver_noise = (standard_deviation) * randn(1, n); % Gaussian/Normal distribution- RSS error - 1 rows and n columns
end