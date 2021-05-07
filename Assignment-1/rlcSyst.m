%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 1 - Summer Semester 2020/2021
% Kavelidis Frantzis Dimitrios - AEM 9351 - kavelids@ece.auth.gr - ECE AUTH

% Function describing the Differential Equation of the RLC system
function ydot = rlcSyst(t,y)
global wmega amp1 amp2 RC LC
u1 = amp1*sin(wmega*t);
u1d = amp1*wmega*cos(wmega*t);
u2 = amp2;
u2d = 0;
ydot(1) = y(2);
ydot(2) = -(1/RC)*y(2) -(1/LC)*y(1)+(1/RC)*u1d + (1/RC)*u2d + (1/LC)*u2;
ydot = ydot';
end