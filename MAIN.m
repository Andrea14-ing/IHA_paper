% This script demonstrates the calculation of the:
% Instantaneous Helical Axis IHA
% Finite Helical Axis
% Average Helical Axis, based on IHA and FHA
% Coordinate system attached to the AHA
% Analysis of the dispersion of the IHA or FHA

% The analysis is run on the data contained in the DATA folder and
% described in the following. 

% This script requires the functions in the 'lib' folder

% The code is divided into section. It is recommended to run the sections
% sequentially, one at a time.

% Andrea Ancillao, Ph.D.
% KU Leuven, Leuven ,BE.
% Last edit: 07/2021


% DISCLAIMER:
% This code is for demonstration purposes only.
% You running this script/function means you will not blame the author(s) 
% if this software breaks your stuff, destroys your computer or kills your cat. 
% This script/function is provided AS IT IS, without warranty of any kind. 
% Author(s) disclaim all implied warranties including, without limitation, 
% any implied warranties of merchantability or of fitness for a particular purpose. 
% The entire risk arising out of the use or performance of the sample scripts 
% and documentation remains with you. In no event shall author(s) be held liable for any damages whatsoever 
% (including, without limitation, damages for loss of business profits, business interruption, 
% loss of business information, or other pecuniary loss) arising out of the use of or inability to use 
% the script or documentation. Neither this script/function, nor any part of it may be republished 
% without author(s) express written permission. 
% Author(s) retain the right to alter this disclaimer at any time. 


clear;
close all
clc
addpath(genpath('lib'));




%% 1 Load and visualize data
load('.\DATA\\HingeSampleData.mat');
% The data consists in the motion capture recording of an artificial hinge
% joint in a free swing after an initial perturbation. For more info see:
% Ancillao, A., et al.,(2020). Estimating the instantaneous screw axis and the screw axis invariant descriptor of motion by means of inertial sensors: An experimental study with a mechanical hinge joint and comparison to the optoelectronic system. Sensors (Switzerland), 20(1). https://doi.org/10.3390/s20010049

% Load a structure named "hinge" containing the following fields:
% fs: 100                   The sampling frequency
% T1_0: [4×4×3941 double]   Pose of proximal segment
% T2_0: [4×4×3941 double]   Pose of distal segment
% GA_0: [3941×6 double]     Geometric rotation axis, direction and origin
% markers: [1×1 struct]     Coordinates of measured markers
% tw1_0: [3941×6 double]    Screw twist of proximal segment seen in {0}
% tw2_0: [3941×6 double]    Screw twist of distal segment seen in {0}
% tw21_1: [3941×6 double]   Relative Screw twist of proximal segment w.r.t the distal one, seen in {1}
%
T1 = hinge.T1_0;
T2 = hinge.T2_0;
GA_0 = hinge.GA_0;
tw21_1 = hinge.tw21_1;

% Visualize the motion in the CS{0}
figure;
hold on; grid on; box on; grid minor;
plotH(0.1,eye(4),'k');
nF = size(T1,3); % number of samples
h1 = plotH(0.1,T1(:,:,1),'b');
h3 = plotSA(GA_0(1,:) , 0.3, 'c','LineWidth',2);
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

if 0
% Animated plot of the motion (press ctr+c to interrupt)
figure;
hold on; grid on; box on; grid minor; axis equal;
nF = size(T1,3);
for k=100:5:600
    RF = k;
    cla;
    plotHingeMarkers(hinge.markers,RF,'b');
    plotSA(GA_0(RF,:) , 0.3, 'c','LineWidth',2);
    plotH(0.1,T1(:,:,RF),'b');
    plotH(0.1,T2(:,:,RF),'r');
    view([35,15]);
    xlim([-0.8743    0.3925]);
    ylim([-0.5461    0.4530]);
    zlim([-0.0857    0.9134]);
    xlabel('x - [m]');
    ylabel('y - [m]');
    zlabel('z - [m]');
    pause(0.03);
end
end


%% 2 Calculate the IHA, the AHA based on IHA, an the CA attached to the AHA
[IHA , AHA] = calcIHA(tw21_1); % calculate in {1}
IHA_1 = [IHA.n, IHA.s];
AHA_1 = [AHA.n, AHA.s];
Taha_1 = AHA.T;

% Transform to {0} for plotting
n1 = IHA_1(:,1:3);
S1 = IHA_1(:,4:6);
R1 = T2RO(T1);
n0 = transfVect(n1,R1);
S0 = transfPoint(S1,T1);
IHA_0 = [n0,S0];
n1 = AHA_1(:,1:3);
S1 = AHA_1(:,4:6);
n0 = transfVect(n1,R1);
S0 = transfPoint(S1,T1);
AHA_0 = [n0,S0];
Taha_0 = Tcat(T1,Taha_1);


% Plot the IHA and AHA in {0}
figure;
hold on; grid on; box on; grid minor; axis equal;
nF = size(T1,3);
RF = 100;
h1 = plotHingeMarkers(hinge.markers,RF,'b');
h2 = plotSA(AHA_0(RF,:), 1, 'g','LineWidth',3);
h3 = plotSA2(IHA_0 , 0.3, 'b','LineWidth',1);
h4 = plotH(0.1,T1(:,:,RF),'b');
h5 = plotH(0.1,T2(:,:,RF),'r');
h6 = plotH(0.1,Taha_0(:,:,RF),'g');
view([35,15]);
xlim([-0.8743    0.3925]);
ylim([-0.5461    0.4530]);
zlim([-0.0857    0.9134]);
xlabel('x - [m]');
ylabel('y - [m]');
zlabel('z - [m]');
legend([h1, h4, h5,h3,h2,h6],{'Markers','CS 1','CS 2','IHA','AHA','CS AHA'});


%% 3 Dispersion analysis of the IHA
% The dispersion analysis can be run in any CS. Here we run the analysis in
% the AHA CS for easier / better visualization
T1_aha = Hinv(Taha_1);
AHA_aha = [0 0 1 0 0 0];
n1 = IHA_1(:,1:3);
S1 = IHA_1(:,4:6);
R1_aha = T2RO(T1_aha);
n_aha = transfVect(n1,R1_aha);
S_aha = transfPoint(S1,T1_aha);
IHA_aha = [n_aha,S_aha];

% run the dispersion analysis 
d = 1; % plane at unitary distance from origin
dispParam = dispersionAnalysis(IHA_aha, eye(4), d, 0);

% plot in {aha}
figure;
subplot 121
hold on; grid on; box on; grid minor; axis equal;
h1 = plotSA(AHA_aha, 1.5, 'g','LineWidth',3);
h2 = plotH(0.5,eye(4),'g');
h3 = plotSA2(IHA_aha , 1.3, 'b','LineWidth',1);
n = dispParam.Tpl(1:3,3)';
O = dispParam.Tpl(1:3,4)';
h4 = plotPlane(n,O,0.5,'c',0.5);
view([35,15]);
xlabel('x - [m]');
ylabel('y - [m]');
zlabel('z - [m]');
legend([h1, h2,h3,h4],{'AHA','CS AHA','IHA','Ref. Plane'});
title('Dispersion analysis of the IHA seen in CS \{AHA\}');

subplot 122
hold on; grid on; box on; grid minor;
PlaneInters_pl = [dispParam.x , dispParam.y];
nP = size(PlaneInters_pl,1);
for k=1:nP
    plot(PlaneInters_pl(k,1), PlaneInters_pl(k,2) , 'b*');
end
plot(dispParam.o(1),dispParam.o(2) , 'r*', 'LineWidth',4);
plot_ellipse(dispParam.o,dispParam.va*dispParam.ma,dispParam.vb*dispParam.mb,1000,'r','LineWidth',2);
title('Intersection with ref. plane and confidence ellipse');
xlabel('x');
ylabel('y');
axis equal


%% 4 Parameters
Ellipse_Ax1 = dispParam.ma
Ellipse_Ax2 = dispParam.mb
Ellipse_Ratio = dispParam.ma / dispParam.mb

% Compare IHA to AHA
[ang, dist] = compareAxes( IHA_1 , AHA_1 );
figure;
subplot 211
hold on; grid on; grid minor;
plot(ang);
xlabel('Sample n.');
ylabel('[deg]');
title('Angle between IHA and AHA');

subplot 212
hold on; grid on; grid minor;
plot(dist);
xlabel('Sample n.');
ylabel('[m]');
title('DIstance between IHA and AHA');


%% 5 Now calculate the FHA


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







