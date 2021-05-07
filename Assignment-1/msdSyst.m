%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 1 - Summer Semester 2020/2021
% Kavelidis Frantzis Dimitrios - AEM 9351 - kavelids@ece.auth.gr - ECE AUTH

% Function describing the differential equation of the system
function ydot = msdSyst(t,y)
    global m b k w_force F0 F1
    u = F0*sin(w_force*t) + F1;
    ydot(1) = y(2);
    ydot(2) = u/m - y(2)*b/m - y(1)*k/m;
    ydot = ydot';
end