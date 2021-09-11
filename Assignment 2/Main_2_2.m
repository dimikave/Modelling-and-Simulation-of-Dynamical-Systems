%% SYSTEMS MODELING AND SIMULATION
% Assignment 2 - Summer Semester 2020/2021
% Kavelidis Frantzis Dimitrios - AEM 9351 - kavelids@ece.auth.gr - ECE AUTH

%% Description:
% On-line estimation of unknown parameters using:
% i) Gradient Descent Method
% ii) Lyapunov Method : 1) Parallel Configuration
%                       2) Series-Parallel Configuration
% iii) Lyapunov method / Parallel Configuration on MIMO system (2x2)
%% Part ii) - Real time estimator / Lyapunov Method
                           % 1) Parallel Configuration
                           % 2) Mixed/Series-Parallel Configuration
% System : dx/dt = -a*x + b*u, x(0) = 0 & u = 5*sin(3t)

%% Clearing
clear all;
close all;
clc;
format longG
tic;                        % Start clock for code evaluation
%% Real values of parameters and configuration parameters
global a b c d pargamma mixgamma thetam h0 f
a = 2;                     % Real value of ODE parameter a
b = 1;                     % Real value of ODE parameter b
c = 5;                     % Magnitude of input sin
d = 3;                     % Omega of input sin
pargamma = [5 1];         % Parallel Configuration parameters
mixgamma = [5 1];          % Mixed/Series-Parallel Configuration Parameters
thetam = 2;               % Method Parameter
h0 = 0.15;                 % Magnitude of noise sin
f = 20;                    % Frequency of noise sin

%% Lyapunov Method - Parallel Configuration
% Time Span
tStart = 0;   
tStep = 0.001;
tEnd = 20;
tspan = tStart:tStep:tEnd;
% Initial value
initCond = zeros(1,4);                                    
% Solving ODE for noise free and noise output
[tpar,xpar] = ode45(@(tpar,xpar) LyapPar(tpar,xpar), tspan, initCond);
[tparn,xparn] = ode45(@(tpar,xpar) LyapParNoise(tpar,xpar), tspan, initCond);

% Getting results -> easy-to-understand variables
% Noise Free
x = xpar(:,1);
th1est = xpar(:,2);
th2est = xpar(:,3);
xest = xpar(:,4);

% With Noise
xn = xparn(:,1);
th1estn = xparn(:,2);
th2estn = xparn(:,3);
xestn = xparn(:,4);

% Unknown parameters estimation
aest = th1est;
best = th2est;
aestn = th1estn;
bestn = th2estn;

[aest aestn best bestn]
% Error e
e = x - xest;
en = xn - xestn;

% Lyapunov function & Derivative - Noise free
V = (1/2)*(e.^2 + (1/pargamma(1))*(aest-a).^2 ...
    + (1/pargamma(2))*(best-b).^2);
Vdot = - a*e.^2;
[V Vdot];         % If we want to check for <= 0, >= 0 (Lyapunov stability)

% Lyapunov function & Derivative - Noise Presence
Vn = (1/2)*(en.^2 + (1/pargamma(1))*(aestn-a).^2 ...
    + (1/pargamma(2))*(bestn-b).^2);
Vdotn = - thetam*en.^2;
[Vn Vdotn];       % If we want to check for <= 0, >= 0 (Lyapunov stability)
%% Plot simulation - Parallel Configuration
figure
suptitle('Parallel Configuration')

subplot(2,2,1)
plot(tpar,x,'r')
hold on
plot(tpar,xest,'b')
hold off
title('System Simulation - x & x_e_s_t - Noise Free')
xlabel('Time [s]')
ylabel('x , x_e_s_t')
grid on
legend('x','x_e_s_t')

subplot(2,2,2)
plot(tpar,xn,'r')
hold on
plot(tpar,xestn,'b')
hold off
title('System Simulation - x & x_e_s_t - Noise Presence')
xlabel('Time [s]')
ylabel('x , x_e_s_t')
grid on
legend('x','x_e_s_t')

subplot(2,2,3)
plot(tpar,xest,'r')    
hold on
plot(tpar,xestn,'b')
hold off
title('Comparison of estimations with / without noise')
xlabel('Time [s]')
ylabel('x')
grid on
legend('x_e_s_t without noise', 'x_e_s_t with Noise')

subplot(2,2,4)
plot(tpar,e,'r')    
hold on
plot(tpar,en,'b')
hold off
title('Error e = x - x_e_s_t - Comparison with / without noise')
xlabel('Time [s]')
ylabel('e')
grid on
legend('Error without noise', 'Error with Noise')

%% Lyapunov Plots - Parallel Configuration
figure
subplot(2,2,1)
plot(tpar,V,'r')    
title('Lyapunov Function - Parallel Configuration - Noise Free')
xlabel('Time [s]')
ylabel('V')
grid on

subplot(2,2,2)
plot(tpar,Vdot,'r')    
title('Lyapunov Function Derivative - Parallel Configuration - Noise Free')
xlabel('Time [s]')
ylabel('V_d')
grid on

subplot(2,2,3)
plot(tpar,V,'r')    
title('Lyapunov Function - Parallel Configuration - Noise Presence')
xlabel('Time [s]')
ylabel('V')
grid on

subplot(2,2,4)
plot(tpar,Vdot,'r')    
title('Lyapunov Function Derivative - Parallel Configuration - Noise Presence')
xlabel('Time [s]')
ylabel('V_d')
grid on

%% Lyapunov Method - Mixed / Series-Parallel Configuration
% Initial conditions
initCondm = zeros(1,4);

% Solving ODE
[tmix,xmix] = ode45(@(tmix,xmix) LyapMix(tmix,xmix), tspan, initCond);
[tmixn,xmixn] = ode45(@(tmix,xmix) LyapMixNoise(tmix,xmix), tspan, initCond);
% Getting results -> easy-to-understand variables
% Noise Free
xm = xmix(:,1);
th1estm = xmix(:,2);
th2estm = xmix(:,3);
xestm = xmix(:,4);
% With Noise
xmn = xmixn(:,1);
th1estmn = xmixn(:,2);
th2estmn = xmixn(:,3);
xestmn = xmixn(:,4);

aestm = th1estm;
bestm = th2estm;
aestmn = th1estmn;
bestmn = th2estmn;

[aestm aestmn bestm bestmn]

% Error e
em = xm - xestm;
emn = xmn - xestmn;

% Lyapunov function & Derivative
Vm = (1/2)*(em.^2 + (1/mixgamma(1))*(aestm-a).^2 ...
    + (1/mixgamma(2))*(bestm-b).^2);
Vdotm = - thetam*em.^2;
[Vm Vdotm];       % If we want to check for <= 0, >= 0 (Lyapunov stability)

% Lyapunov function & Derivative - Noise Presence
Vmn = (1/2)*(emn.^2 + (1/pargamma(1))*(aestmn-a).^2 ...
    + (1/pargamma(2))*(bestmn-b).^2);
Vdotmn = - thetam*emn.^2;
[Vmn Vdotmn];     % If we want to check for <= 0, >= 0 (Lyapunov stability)
%% Plot simulation - Mixed Configuration
figure
suptitle('Mixed Configuration')

subplot(2,2,1)
plot(tmix,xm,'r')
hold on
plot(tmix,xestm,'b')
hold off
title('System Simulation - x & x_e_s_t - Noise Free')
xlabel('Time [s]')
ylabel('x , x_e_s_t')
grid on
legend('x','x_e_s_t')

subplot(2,2,2)
plot(tmix,xmn,'r')
hold on
plot(tmix,xestmn,'b')
hold off
title('System Simulation - x & x_e_s_t - Noise Presence')
xlabel('Time [s]')
ylabel('x , x_e_s_t')
grid on
legend('x','x_e_s_t')

subplot(2,2,3)
plot(tmix,xestm,'r')    
hold on
plot(tmix,xestmn,'b')
hold off
title('Comparison of estimations with / without noise')
xlabel('Time [s]')
ylabel('x')
grid on
legend('x_e_s_t without noise', 'x_e_s_t with Noise')

subplot(2,2,4)
plot(tmix,em,'r')    
hold on
plot(tmix,emn,'b')
hold off
title('Error e = x - x_e_s_t - Comparison with / without noise')
xlabel('Time [s]')
ylabel('e')
grid on
legend('Error without noise', 'Error with Noise')

%% Lyapunov Plots - Mixed Configuration
figure
subplot(2,2,1)
plot(tpar,Vm,'r')    
title('Lyapunov Function - Mixed Configuration - Noise Free')
xlabel('Time [s]')
ylabel('V')
grid on

subplot(2,2,2)
plot(tpar,Vdotm,'r')    
title('Lyapunov Function Derivative - Mixed Configuration - Noise Free')
xlabel('Time [s]')
ylabel('V_d')
grid on

subplot(2,2,3)
plot(tpar,Vmn,'r')    
title('Lyapunov Function - Mixed Configuration - Noise Presence')
xlabel('Time [s]')
ylabel('V')
grid on

subplot(2,2,4)
plot(tpar,Vdotmn,'r')    
title('Lyapunov Function Derivative - Mixed Configuration - Noise Presence')
xlabel('Time [s]')
ylabel('V_d')
grid on

%% Comparing methods plot
figure
subplot(1,2,1)
plot(tmix,e,'r')    
hold on
plot(tmix,em,'b')
hold off
title('Parallel / Mixed Comparison - Noise Free')
xlabel('Time [s]')
ylabel('e')
grid on
legend('Error from Parallel Conf', 'Error from Mixed Conf')

subplot(1,2,2)
plot(tmix,en,'r')    
hold on
plot(tmix,emn,'b')
hold off
title('Parallel / Mixed Comparison - Noise Presence')
xlabel('Time [s]')
ylabel('e')
grid on
legend('Error from Parallel Conf', 'Error from Mixed Conf')

%% Theta estimation plots
figure
subplot(2,2,1)
plot(tpar,th1est,'r')
hold on
plot(tpar,th1estn)
title('Theta1 - Parallel Configuration - With/Without Noise')
xlabel('Time [s]')
ylabel('theta_1_,_e_s_t')
grid on
legend('Theta1 without noise','Theta1 with noise')

subplot(2,2,2)
plot(tpar,th2est,'r')
hold on
plot(tpar,th2estn)
title('Theta2 - Parallel Configuration - With/Without Noise')
xlabel('Time [s]')
ylabel('theta_2_,_e_s_t')
grid on
legend('Theta2 without noise','Theta2 with noise')

subplot(2,2,3)
plot(tmix,th1estm,'r')
hold on
plot(tmix,th1estmn)
title('Theta1 - Mixed Configuration - With/Without Noise')
xlabel('Time [s]')
ylabel('theta_1_,_e_s_t')
grid on
legend('Theta1 without noise','Theta1 with noise')

subplot(2,2,4)
plot(tpar,th2estm,'r')
hold on
plot(tpar,th2estmn)
title('Theta2 - Mixed Configuration - With/Without Noise')
xlabel('Time [s]')
ylabel('theta_2_,_e_s_t')
grid on
legend('Theta2 without noise','Theta2 with noise')

toc;                                            % Stop clock
%% ----------------------- End of Part ii) --------------------------------