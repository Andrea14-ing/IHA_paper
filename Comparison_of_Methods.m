clear;
close all
clc
addpath(genpath('lib'));


% Comparison of the IHA and AHA calculation methods


% Load Data
load('.\DATA\\Sampledata_Gait_Hinge.mat');





%% 1 - Hinge

T1 = hinge.T1_0;
T2 = hinge.T2_0;
fs = hinge.fs;

% get GA
[~ , n] = mvarray( hinge.markers.ax1 - hinge.markers.ax2) ;
s =  (hinge.markers.ax1 + hinge.markers.ax2) ./ 2;
GA_0 = [n , s];


% Calculate the screw twist
T0_2 = Hinv(T2);
T2_1 = Tcat(T0_2,T1);
tw1 = T2twist(T1,fs,'screw');
tw2 = T2twist(T2,fs,'screw');
tw2_1_0 = tw2-tw1;
T0_1 = Hinv(T1);
S0_1 = T2S(T0_1);
tw2_1_1 = transformTwist(tw2_1_0, S0_1);

nF = size(T2_1, 3);
f = 1:1:nF;
t = f./fs;

% Plot the screw twist
figure;
plotTwist(t,tw2_1_1,1,'b');
sgtitle('Screw Twist of the hinge');


% plot the CS in the {0}
figure;
hold on; grid on; box on; grid minor;
plotH(0.1,eye(4),'k');
nF = size(T1,3); % number of samples
h1 = plotH(0.1,T1(:,:,1),'b');
h3 = plotSA(GA_0(1,:) , 0.3, 'b','LineWidth',2);
for k=1:5:200
   h2 = plotH(0.1,T2(:,:,k),'r');
end
axis equal
view([30,15]);
xlabel('x - [m]');
ylabel('y - [m]');
zlabel('z - [m]');
legend([h1, h2, h3],{'CS 1','CS 2','G.A.'});
title('Poses of hinge segments');


% Calc IHA and IHA in {0} 
[IHA , AHA] = calcIHA(tw2_1_0); 
IHA_0 = [IHA.n, IHA.s];
AHA_0 = [AHA.n, AHA.s];


% calc FHA and AHA in {0}
if 0
    [FHA, MFHA] = calcFHA(T2_1(:,:,500:1000)); % limit to the first samples for speed
    FHA_1 = [FHA.n, FHA.s];
    MFHA_1 = [MFHA.n, MFHA.s];
    save FHA FHA_1 MFHA_1
end
% load previously processed for speed
load FHA FHA_1 MFHA_1
n = transfVect(MFHA_1(:,1:3), T2RO(T1) );
s = transfPoint(MFHA_1(:,4:6), T1);
MFHA_0 = mean( [n , s] , 1 , 'omitnan');


% calc ISA and ASA == IHA with Erwin's Method in {0}
[n,s1,T] = twist2ISA(tw2_1_0);
ISA = [n,s1];
refDir = mean(GA_0(:,1:3) , 1, 'omitnan');
[n,s1,T] = twist2ASA(tw2_1_0 , refDir , 1);
ASA = [n,s1];


% Plot the axes
figure;
hold on; grid on; box on; grid minor; axis equal;
nF = size(T1,3);
RF = 1000;
h1 = plotHingeMarkers(hinge.markers,RF,'b');
h2 = plotSA2(GA_0(RF,:),0.5,'b');
h3 = plotH(0.1,T1(:,:,RF),'g');
plotH(0.1,T2(:,:,RF),'g');
h4 = plotSA2(IHA_0,0.3,'c');
h5 = plotSA2(AHA_0,0.5,'g');
h6 = plotSA2(MFHA_0,0.6,'r');
h7 = plotSA2(ASA,0.6,'y');
view([35,15]);
axis equal
xlim([-0.8743    0.3925]);
ylim([-0.5461    0.4530]);
zlim([-0.0857    0.9134]);
xlabel('x - [m]');
ylabel('y - [m]');
zlabel('z - [m]');
legend([h1,h2,h3,h4,h5,h6,h7],{'Markers','GA','Segment CSs','IHA','AHA','MFHA','ASA'});



% Dispersion analysis
out = dispersionAnalysis(IHA_0, AHA.T, 1, 1);



% Compare the axes
% GA to AHA
[ang, dist1, dist2] = compareAxes( mean(GA_0,1,'omitnan') , AHA_0 )

% GA to MFHA
[ang, dist1, dist2] = compareAxes( mean(GA_0,1,'omitnan') , mean(MFHA_0,1,'omitnan' ))

% MFHA to AHA
[ang, dist1, dist2] = compareAxes( mean(MFHA_0,1,'omitnan' ) , AHA_0)










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2 - Gait


%% Gait data
% cut in between events
if strcmp(gait.side , 'R')
    hs = gait.fluoro.EV.rhs.f;
else
    hs = gait.fluoro.EV.lhs.f;
end
frange = hs(1):hs(end);

T1_0 = gait.fluoro.T1_0(:,:,frange);
T2_0 = gait.fluoro.T2_0(:,:,frange);


% resample to gait cycle
T1_0 = FrameRange2GaitCycle_matrices(T1_0);
T2_0 = FrameRange2GaitCycle_matrices(T2_0);

nF = size(T1_0,3);
fs = gait.fluoro.fs;

GA_0 = [ squeeze( T1_0(1:3,3,:) )' , squeeze( T1_0(1:3,4,:) )' ];

% Plot the frames in the world
figure;
hold on; grid on; box on; grid minor;
plotH(0.1,eye(4),'k');
for k=1:3:nF
    h1 = plotH(0.1,T1_0(:,:,k),'b');
    h2 = plotH(0.1,T2_0(:,:,k),'r');
    h3 = plotSA(GA_0(k,:) , 0.3, 'c','LineWidth',2);
end
axis equal
view([-60,25]);
legend([h1, h2, h3],{'CS 1','CS 2','G.A.'});
title('Poses of femur and tibia');


% Calculate the screw twist of the knee
t1 = T2twist(T1_0,fs);
t2 = T2twist(T2_0,fs);
t21_0 = t2 - t1;
T0_1 = Hinv(T1_0);
S0_1 = T2S(T0_1);
t21_1 = transformTwist(t21_0,S0_1);


% Plot the screw twist of the knee
figure;
x=0:1:nF-1;
plotTwist(x,t21_1,1,'b');
sgtitle('Screw twist of tibia wrt femur seen in femur');


% change CS for GA
R0_1 = T2RO(T0_1);
GA_1 = nan(nF,6);
GA_1(:,1:3) = transfVect(GA_0(:,1:3) , R0_1);
GA_1(:,4:6) = transfPoint(GA_0(:,4:6) , T0_1);
GA_1 = nanmean(GA_1 , 1);

% calculate IHA and AHA
[IHA , AHA] = calcIHA(t21_1);



% Plot
figure;
hold on; grid on; box on; grid minor;
SA = [IHA.n,IHA.s];
plotSA(SA , 0.2, 'b');
SA = [AHA.n,AHA.s];
plotSA(SA , 0.3, 'g','LineWidth',3);
plotH(0.05,AHA.T,'g');
plotSA(GA_1, 0.3, 'r', 'LineWidth',3);
axis equal;
view([115 , 15]);
title('Gait');






