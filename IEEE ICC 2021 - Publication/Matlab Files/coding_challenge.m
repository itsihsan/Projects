drawnValues = 20; % Sample space
N=20; % Total coins
for i = 1:factorial(N)
    y = (randsample(N,drawnValues,false))';
    absDiff = abs(y(1:drawnValues-1) - y(2:drawnValues));
    totPayment(i) = y(1)+sum(absDiff);
    y=[];
    absDiff=[];
end
digitsOld = digits(10)
meanTotPayment = vpa(mean(totPayment))
stdTotPayment = vpa(std(totPayment))