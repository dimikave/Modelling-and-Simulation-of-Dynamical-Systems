%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 1 - Summer Semester 2020/2021
% Kavelidis Frantzis Dimitrios - AEM 9351 - kavelids@ece.auth.gr - ECE AUTH

%% Exercise 2 - LS
% Estimating unknown parameters using Least-Squares Method
% on an electric circuit with two inputs/two outputs, two loops.

%% Clearing
clear all;
close all;
clc;
%% Changing Format
format longG
%% a 
% % Vector for holding the maximums for each pole placement test
%maxdif = NaN(700,1);
%% Estimating with Least-Square Method
% Choosing filter poles by trial and error
% for i = 1:1:700
% i
% p1 = i;
% p2 = i;
p1 = 370;
p2 = 370;
% Simulating for 0.1 with a time steps of 0.0001,0.00001 and 0.000001
t = 0:0.000001:1;
Vout = zeros(length(t),2);
for i = 1:length(t)
    Vout(i,:) = v(t(i));
end
VR = Vout(:,2); % Defining VR
VC = Vout(:,1); % Defining VC
% Plot of our measurements
figure()
plot(t,VR,"r")
hold on
plot(t,VC)
title("V_R and V_C measurements vs time")
xlabel("Time [s]")
ylabel("V_R,V_C  [V,V]")
legend("V_R","V_C")

% Creating inputs
global wmega amp1 amp2
wmega = 1;
amp1 = 2;
amp2 = 1;
u1 = amp1*sin(wmega*t)';
u2 = amp2*ones(length(t),1);
% Creating phi matrix
phi1 = lsim(tf([-1 0],[1 (p1+p2) p1*p2]),VC,t);     
phi2 = lsim(tf(-1,[1 (p1+p2) p1*p2]),VC,t);         
phi3 = lsim(tf([1 0],[1 (p1+p2) p1*p2]),u1,t);      
phi4 = lsim(tf(1,[1 (p1+p2) p1*p2]),u1,t);
phi5 = lsim(tf([1 0],[1 (p1+p2) p1*p2]),u2,t);
phi6 = lsim(tf(1,[1 (p1+p2) p1*p2]),u2,t);
phi = zeros(length(t),6);
phi(:,1) = phi1;
phi(:,2) = phi2;
phi(:,3) = phi3;
phi(:,4) = phi4;
phi(:,5) = phi5;
phi(:,6) = phi6;
% Solving / Estimating the parameters by finding theta*
phiTphi = phi.'*phi;                                       % Phi squared
VCTphi = VC.'*phi;                                   % Y^T*Phi
theta0 = VCTphi/phiTphi;                 % vector that is argmin_theta
theta = theta0 + [p1+p2 p1*p2 0 0 0 0]   % vector that holds the unknown parameters

% Assigning value to the estimated parameters by calculating the mean
global RC LC
RC = 1/((1/3)*(theta(1)+theta(3)+theta(5)));
LC = 1/((1/2)*(theta(2)+theta(6)));
%% Simulation using ode
[time,ysol] = ode45(@(time,y) rlcSyst(time,y),t,[VC(1) 0]);
% Defining VC,VR that came from the simulation on the model
% with the estimated parameters, and estimating their errors
VCnew = ysol(:,1);
VCdif = [VC VCnew (VC-VCnew)];
VRnew = u1+u2-VCnew;        % Using the equations of the system
VRdif = [VR VRnew (VR-VRnew)];

% Comparing plots
figure()
% VR VC
subplot(2,2,1)
plot(t,VRnew,'r',t,VCnew);
legend("V_R","V_C")
title("The measurements of voltages using our model's estimated parameters")
ylabel("Voltage [V]")
xlabel("Time [s]")

% VRnew VCnew
subplot(2,2,2)
plot(t,VR,'r',t,VC);
title("The measurements of voltages taken initially")
ylabel("Voltage [V]")
xlabel("Time [s]")
legend("V_R_n_e_w","V_C_n_e_w")

% VR VRnew
subplot(2,2,3)
plot(t,VR,'r',t,VRnew);
legend("V_R","V_R_n_e_w")
title("Comparing V_R / V_R_n_e_w")
ylabel("Voltage [V]")
xlabel("Time [s]")

% VC VCnew
subplot(2,2,4)
plot(t,VC,'r',t,VCnew);
legend("V_C","V_C_n_e_w")
title("Comparing V_R / V_C_n_e_w")
ylabel("Voltage [V]")
xlabel("Time [s]")

% Matrix holding VC and VR errors
DIF = [VCdif(:,3) VRdif(:,3)];

% Plot of VC error
figure()
subplot(2,1,1)
plot(t,DIF(:,1),'r')
title("Error of V_C")
ylabel("Voltage [V]")
xlabel("Time [s]")

% Plot of VR error
subplot(2,1,2)
plot(t,DIF(:,2))
title("Error of V_R")
ylabel("Voltage [V]")
xlabel("Time [s]")

% % MinMax: finding the index (pole value) of min of the vector holding the
% % maximums
% maxdif(i) = max(DIF(:,1));
% end
% [a,bestPole] = min(maxdif)
fprintf("Transfer Function Matrix:\n")
G1 = tf([1 0 1/LC],[1 1/RC 1/LC]);
G2 = tf([1 0 0],[1 1/RC 1/LC]);
G3 = tf([1/RC 0] ,[1 1/RC 1/LC]);
G4 = tf([1/RC 1/LC],[1 1/RC 1/LC]);
G = [G1 G2;G3 G4]
fprintf("Pause before proceeding to part b... Press any button to continue:\n")
pause;
%% b / Theta* with noise
% Copying vectors
VCb = VC;
VRb = VR;
% Adding noise on 3 random times
for i = 1:1:3
    j = randi(length(VCb));
    l = randi(length(VRb));
    VCb(j) = VCb(j) + randi([1000,20000]);
    VRb(l) = VRb(l) + randi([1000,20000]);
end

% Plot of VC with noise
figure()
subplot(2,1,1)
plot(t,VCb,'r')
title("V_C with noise")
ylabel("Voltage [V]")
xlabel("Time [s]")

% Plot of VR with noise
subplot(2,1,2)
plot(t,VRb)
title("V_R with noise")
ylabel("Voltage [V]")
xlabel("Time [s]")

% Creating Phi matrix
phi1b = lsim(tf([-1 0],[1 (p1+p2) p1*p2]),VCb,t);     
phi2b = lsim(tf(-1,[1 (p1+p2) p1*p2]),VC,t);         
phi3b = lsim(tf([1 0],[1 (p1+p2) p1*p2]),u1,t);     
phi4b = lsim(tf(1,[1 (p1+p2) p1*p2]),u1,t);
phi5b = lsim(tf([1 0],[1 (p1+p2) p1*p2]),u2,t);
phi6b = lsim(tf(1,[1 (p1+p2) p1*p2]),u2,t);
phib = zeros(length(t),6);
phib(:,1) = phi1;
phib(:,2) = phi2;
phib(:,3) = phi3;
phib(:,4) = phi4;
phib(:,5) = phi5;
phib(:,6) = phi6;
% Solving / Estimating the parameters by finding theta*
phiTphib = phib.'*phib;                                       % Phib squared
VCTphib = VCb.'*phib;                                   % Y^T*Phib
theta0b = VCTphib/phiTphib;                 % vector that is argmin_thetab
thetab = theta0b + [p1+p2 p1*p2 0 0 0 0]   % vector that holds the unknown parameters

%% ---------------------------- End of Exercise 2 -------------------------