% SYSTEMS MODELING AND SIMULATION
% Assignment 1 - Summer Semester 2020/2021
% Kavelidis Frantzis Dimitrios - AEM 9351 - kavelids@ece.auth.gr - ECE AUTH

%% Exercise 1 - LS
% Simulating a mass-spring-dumper system 
% Estimating unknown parameters using Least-Squares Method

%% Clearing
clear all;
close all;
clc;
%% Real values of parameters
global m b k w_force F0 F1
m = 15;         % Defining mass of system
b = 0.2;        % Defining damping
k = 2;          % Defining stiffness of the spring
w_force = 2;    % Defining the frequency of the sin / external force
F0 = 5;         % Defining the amplitude of the sin / external force
F1 = 10.5;      % Defining a constant external force
%% Simulation of the real system
tStart = 0;     % Time span
tStep = 0.1;
tEnd = 10;
tspan = tStart:tStep:tEnd; % simulation for 10 seconds
y0 = 0;   % initial values
yd0 = 0;
[t,ysol] = ode45(@(t,y) msdSyst(t,y),tspan,[y0 yd0]); % Solving DE
plot(t,ysol(:,2),'r',t,ysol(:,1))                     % Plotting simulation
title("Mass-Spring-Dumper System Simulation")
xlabel("Time [s]")
ylabel("Displacement, Speed  [m, m/s]")
legend("Speed","Displacement")
grid on
%% Estimating with Least-Squares method
% MestEr = NaN(100,1);
% KestEr = NaN(100,1);
% BestEr = NaN(100,1);
% for i = 1:1:100
% % Choosing filter poles by trial and error
% i
% p1 = i/100;
% p2 = i/100;
p1 = 0.47         % Pole Placement
p2 = 0.47
% Creating phi matrix
u = F0*sin(w_force*t) + F1;       % Recreating u force to use in Phi matrix
figure
plot(t,u);
phi1 = lsim(tf([-1 0],[1 (p1+p2) p1*p2]),ysol(:,1),t);     % zeta(1) = phi1
phi2 = lsim(tf(-1,[1 (p1+p2) p1*p2]),ysol(:,1),t);         % zeta(2) = phi2
phi3 = lsim(tf(1,[1 (p1+p2) p1*p2]),u,t);                  % zeta(3) = phi3
phi = zeros(length(t),3);                                  % Phi Matrix
phi(:,1) = phi1;
phi(:,2) = phi2;
phi(:,3) = phi3;
phiTphi = phi.'*phi;                                       % Phi squared
YTphi = ysol(:,1).'*phi;                                   % Y^T*Phi
theta0 = YTphi/phiTphi;                       % Solution of a xA = B system
mest = 1/theta0(3)                            % Estimation of m / Display
kest = mest*(theta0(2)+p1*p2)                % Estimation of k / Display
best = mest*(theta0(1)+p1+p2)                 % Estimation of b / Display
% MestEr(i) = (m-mest);
% KestEr(i) = (k-kest);
% BestEr(i) = (b-best);
% end
% Errors = [MestEr KestEr BestEr]

%% ------------------------- End of exercise 1 ---------------------------S
