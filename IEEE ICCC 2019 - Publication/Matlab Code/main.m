function [a,b,c,xscale] = main (N,r,min,step,max,iteration)
tic
tick = 0
for t = min:step:max
    tt = round(t/step)+1;
    false_detection_record = [];
    detection_record = [];
    LRT_value_record = [];
    LRT_value_record1 = [];
    %output = [];
        for i = 1 : iteration

            [verifier_pre, verifier, claimed_location , optimal_attack_location, U, Delta_U, R, D] = Optimum_DRSS(t,N,r);

            [falsepositive, detection,LRT_value,LRT_value_1] = DRSSLV (verifier_pre,verifier,claimed_location,optimal_attack_location,R,D);

            false_detection_record = [false_detection_record , falsepositive];
            detection_record = [detection_record , detection];
            %output(i,:) = [LRT_value,detection,LRT_value_1,falsepositive];
            tick = tick + 1 % iteration indicator
        end 

    False_Positive_Rate = numel(find(false_detection_record == 1))/numel(false_detection_record);
    Detection_Rate = numel(find(detection_record == 1))/numel(detection_record);
    a (tt) = False_Positive_Rate;
    b (tt) = Detection_Rate;
    c (tt) = 0.5 * a (tt) + 0.5 * (1 - b (tt));
end
toc
xscale = min:step:max;
figure (1)
plot (xscale,a,'r') % false positive
grid on;
hold on;
saveas(gcf,'false_detection.fig','fig')

figure (2)
plot (xscale,b,'b') % detection
grid on;
hold on;
saveas(gcf,'detection.fig','fig')

figure (3)
plot (xscale,c,'g') % cost
grid on;
hold on;
saveas(gcf,'cost.fig','fig')