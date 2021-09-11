%% SYSTEMS MODELING AND SIMULATION
% Assignment 2 - Summer Semester 2020/2021
% Kavelidis Frantzis Dimitrios - AEM 9351 - kavelids@ece.auth.gr - ECE AUTH

% Lyapunov Method / Parallel Configuration Estimator
function xxdot = LyapMix(t,xx)

    xxdot = zeros(size(xx));
    global a b c d mixgamma thetam
    
    % xx description
    x = xx(1);
    theta1est = xx(2);
    theta2est = xx(3);
    xest = xx(4);
    
    % Input u
    u = c*sin(d*t);
    
    % Error e
    e = x - xest;
    
    % System description
    xxdot(1) = - a*x + b*u;                                 % xdot
    xxdot(2) = - mixgamma(1)*e*x;                           % th1^dot
    xxdot(3) = mixgamma(2)*e*u;                             % th2^dot
    xxdot(4) = -theta1est*xest + theta2est*u + thetam*e;    % x^dot
end