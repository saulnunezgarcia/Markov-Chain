%% Chernobyl
clc
clear
close all
%Defining Variables
s = log(2) / 30.15 / 365;
x = 1/500;
v = 1/20;
g = 1/5;
a = 1/30;
c = 1/15;
z = 1/12;
y = 1/5;
%Building Transition Matrix
Names = ["Stable", "Humans","Atmosphere", "Cows", "Vegtables", "Ground"];
M = [1 s s s s s; 0 1-s x y z 0; 0 0 1-s-2*x-v-g 0 0 0; 0 0 x 1-s-y c 0; 0 0 v 0 1-s-c-z a; 0 0 g 0 0 1-s-a];
mc = dtmc(M','StateNames',Names);
%% Absorbing Matrix
P = flip(flip(M',2));
Q = P(1:5,1:5);
R = P(1:5,6);
N = inv(eye(5) - Q);
M2 = (eye(5) - Q) \ R; %Absorbing Matrix
H = (N - eye(5)) / ((eye(5) .* diag(N))); %Transient Probabilities
%% Plot graph
figure(1);
graphplot(mc, 'ColorEdges', true);
title('Markov Chain Representation')
%% Plot Heatmap
figure(2);
imagesc(M);
colormap(jet);
colorbar;
axis square
xticklabels(Names);
yticklabels(Names);
title('Probabilities Heatmap')
%% Evolution
NumDays = 365*100;
Evolution = zeros(6,NumDays);
Evolution(3,1) = 1;
Humans = zeros(3,NumDays);
for t = 2:NumDays
    Humans(:,t-1) = Evolution(3:5,t-1) .* [x y z]';
    Evolution(:,t) = M * Evolution(:,t-1);
end
%% Plot Evolution
figure(3)
subplot(2,3,1)
plot(Evolution(:,1:15)');
title('15 Days')
subplot(2,3,2)
plot(Evolution(:,1:30)');
title('30 Days')
subplot(2,3,3)
plot(Evolution(:,1:200)');
title('200 Days')
subplot(2,3,4)
plot(Evolution(:,1:365)');
title('One Year')
subplot(2,3,5)
plot(Evolution(:,1:365*10)');
title('10 Years')
subplot(2,3,6)
plot(Evolution');
title('100 Years')
legend(Names);
sgtitle('Evolution')
%% Plot Human Contamination
figure(4)
plot(Humans(:,1:30)');
legend(Names(3:5));
title('Sources of Human Contamination')
%% Plot Total Human Contamination
figure(5)
bar(sum(Humans,2));
xticklabels(Names(3:5));
title('Total Human Contamination After 100 years')

