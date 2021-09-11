%% SYSTEMS MODELING AND SIMULATION
% Assignment 2 - Summer Semester 2020/2021
% Kavelidis Frantzis Dimitrios - AEM 9351 - kavelids@ece.auth.gr - ECE AUTH

%% Description:
% On-line estimation of unknown parameters using:
% i) Gradient Descent Method
% ii) Lyapunov Method : 1) Parallel Configuration
%                       2) Series-Parallel Configuration
% iii) Lyapunov method / Parallel Configuration on MIMO system (2x2)
%% Part iii) - Real time estimator / Lyapunov Method - Second order system

% System : xdot = [a11 a12; a21 a22]*x+[b1;b2]*u , x(:,:0) = [0;0]

%% Clearing
clear all;
close all;
clc;
format longG
tic;                        % Start clock for code evaluation

%% Real values of parameters and configuration parameters
global A b u1 u2 w1 w2 gamma
% A Matrix real parameters
a11 = -0.25;
a12 = 3;
a21 = -5;
a22 = -1;
A = [a11 a12; a21 a22];
% b vector real parameters
b1 = 1;
b2 = 2.2;
b = [b1;b2];
% Input parameters
u1 = 10;
w1 = 2;
u2 = 5;
w2 = 7.5;
gamma = [1 1.5];

%% Lyapunov Method / Real time estimator / Parallel Configuration
% Time Span
tStart = 0;   
tStep = 0.01;
tEnd = 100;
tspan = tStart:tStep:tEnd;
% Initial value
initCond = zeros(1,10);                                    
% Solving ODE 
[t,x] = ode45(@(tpar,xpar) LyapEst(tpar,xpar), tspan, initCond);

% Getting results -> easy to understand variables
x1 = x(:,1);
x2 = x(:,2);
x1est = x(:,3);
x2est = x(:,4);
a11est = x(:,5);
a12est = x(:,6);
a21est = x(:,7);
a22est = x(:,8);
b1est = x(:,9);
b2est = x(:,10);

Aest = [a11est a12est; a21est a22est];
best = [b1est;b2est];

% Errors
e1 = x1 - x1est;
e2 = x2 - x2est;
e = [e1;e2];

% Atilda = [a11est-A(1,1) a12est-A(1,2); a21est-A(2,1) a22est-A(2,2)];
% btilda = [b1est-b(1);b2est-b(2)];

%% Lyapunov function & Derivative
% Finding the Lyapunov equation (AP+PA'=-I) solution using A' so we instead
% solve the PA+A'P=-I that we use in theory.
P = lyap(A',[1 0;0 1]);     
V = zeros(length(t),1);
Vdot = zeros(length(t),1);
for i = 1:length(t)
    Atilda = [a11est(i)-A(1,1) a12est(i)-A(1,2); a21est(i)-A(2,1) a22est(i)-A(2,2)];
    btilda = [b1est(i)-b(1);b2est(i)-b(2)];
    V(i) = [e1(i) e2(i)]*P*[e1(i) e2(i)]'...
        + trace(Atilda'*P*Atilda)/gamma(1)...
        +trace(btilda'*P*btilda)/gamma(2);
    Vdot(i) = -[e1(i) e2(i)]*[e1(i) e2(i)]';
end
[V Vdot]
%% Simulation Plots
figure
subplot(2,2,1)
plot(t,x1,'r',t,x2,'b');
title('System Simulation - x_1 & x_2')
xlabel('Time [s]')
ylabel('x_1 ,x_2')
grid on
legend('x_1','x_2')

subplot(2,2,2)
plot(t,x1est,'r',t,x2est,'b');
title('System Simulation - x_1_,_e_s_t & x_2_,_e_s_t')
xlabel('Time [s]')
ylabel('x_1_,_e_s_t , x_2_,_e_s_t')
grid on
legend('x_1_,_e_s_t','x_2_,_e_s_t')

subplot(2,2,3)
plot(t,x1,'r',t,x1est,'b');
title('System Simulation - x_1 & x_1_,_e_s_t')
xlabel('Time [s]')
ylabel('x_1 ,x_1_,_e_s_t')
grid on
legend('x_1','x_1_,_e_s_t')

subplot(2,2,4)
plot(t,x2,'r',t,x2est,'b');
title('System Simulation - x_2 & x_2_,_e_s_t')
xlabel('Time [s]')
ylabel('x_2 , x_2_,_e_s_t')
grid on
legend('x_2','x_2_,_e_s_t')

%% Error plots
figure
subplot(1,2,1)
plot(t,e1);
title('Error of x_1 estimation (e_1)')
xlabel('Time [s]')
ylabel('e_1')
grid on

subplot(1,2,2)
plot(t,e2);
title('Error of x_2 estimation (e_2)')
xlabel('Time [s]')
ylabel('e_2')
grid on

%% Plots of estimated parameters
figure
subplot(1,2,1)
plot(t,a11est,'r',t,a12est,'g',t,a21est,'b',t,a22est,'k');
title('Real time stimation of A matrix')
xlabel('Time [s]')
ylabel('a_i_j')
grid on
legend('a_1_1','a_1_2','a_2_1','a_2_2')

subplot(1,2,2)
plot(t,b1est,'r',t,b2est,'b')
title('Real time estimation of b vector')
xlabel('Time [s]')
ylabel('b_i')
grid on
legend('b_1','b_2')

%% Lyapunov Function plot
figure
subplot(1,2,1)
plot(t,V,'r')    
title('Lyapunov Function - Parallel Configuration')
xlabel('Time [s]')
ylabel('V')
grid on

subplot(1,2,2)
plot(t,Vdot,'r')    
title('Lyapunov Function Derivative - Parallel Configuration')
xlabel('Time [s]')
ylabel('V_d')
grid on
toc;                                            % Stop clock
%% ----------------------- End of Part iii) -------------------------------