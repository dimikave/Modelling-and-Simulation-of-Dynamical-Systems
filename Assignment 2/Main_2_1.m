%% SYSTEMS MODELING AND SIMULATION
% Assignment 2 - Summer Semester 2020/2021
% Kavelidis Frantzis Dimitrios - AEM 9351 - kavelids@ece.auth.gr - ECE AUTH

%% Description:
% On-line estimation of unknown parameters using:
% i) Gradient Descent Method
% ii) Lyapunov Method : 1) Parallel Configuration
%                       2) Series-Parallel Configuration
% iii) Lyapunov method / Parallel Configuration on MIMO system (2x2)
%% Part i) - Real time estimator / Gradient descent Method
% System : dx/dt = -a*x + b*u, x(0) = 0 & u = 5*sin(3t)

%% Clearing
clear all;
close all;
clc;
format longG
tic;                        % Start clock for code evaluation
%% Real values of parameters
global a b c d am gamma
a = 2;                                      % Real value of ODE parameter a
b = 1;                                      % Real value of ODE parameter b
c = 5;                                      % Magnitude of input sin
d = 3;                                      % Omega of input sin
am = 2;                                     % Filter parameter
gamma = 100;                                % Method parameter

%% Gradient Descent Method Estimator
% Time Span
tStart = 0;   
tStep = 0.1;
tEnd = 10;
tspan = tStart:tStep:tEnd;
% Initial value
initCond = [0 0 0 0 0];                                    
% Solving ODE
[tg,xg] = ode45(@(tg,xg) GraDE(tg,xg), tspan, initCond);    

% Getting results -> easy-to-understand variables
x = xg(:,1);
phi1 = xg(:,2);
phi2 = xg(:,3);
th1est = xg(:,4);
th2est = xg(:,5);

% Constructing theta and phi vectors
th_est = [th1est th2est];
phi = [phi1 phi2];

% Estimation of x
xest = th1est.*phi1 + th2est.*phi2;

% Estimation of system parameters
aest = am - th1est
best = th2est

% Error e
e = x - xest;
% Lyapunov function & derivative
V = (1/2)*(aest-a).^2 + (1/2)*(best-b).^2;
Vdot = - gamma*e.^2;
[V Vdot]

%% Plot simulation
figure
subplot(2,2,1)
plot(tg,x,'r')    
title('System Simulation - x')
xlabel('Time [s]')
ylabel('x')
grid on

subplot(2,2,2)
plot(tg,xest,'r')    
title('System Simulation - x_e_s_t')
xlabel('Time [s]')
ylabel('x_e_s_t')
grid on

subplot(2,2,3)
plot(tg,x,'r')
hold on
plot(tg,xest,'b')
hold off
title('System Simulation - x & x_e_s_t')
xlabel('Time [s]')
ylabel('x , x_e_s_t')
grid on
legend('x','x_e_s_t')

subplot(2,2,4)
plot(tg,x-xest,'r')    
title('System Simulation - Error e = x - x_e_s_t')
xlabel('Time [s]')
ylabel('e')
grid on


%% Lyapunov Plots
figure
subplot(1,2,1)
plot(tg,V,'r')    
title('Lyapunov Function')
xlabel('Time [s]')
ylabel('V')
grid on

subplot(1,2,2)
plot(tg,Vdot,'r')    
title('Lyapunov Function Derivative')
xlabel('Time [s]')
ylabel('V_d')
grid on

%% a , b estimation plots
figure
subplot(1,2,1)
plot(tg,aest,'r')
title('Estimation of a')
xlabel('Time [s]')
ylabel('a_e_s_t')
grid on

subplot(1,2,2)
plot(tg,best,'r')    
title('Estimation of b')
xlabel('Time [s]')
ylabel('b_e_s_t')
grid on

toc;                                            % Stop clock
%% ------------------------- End of Part i) -------------------------------