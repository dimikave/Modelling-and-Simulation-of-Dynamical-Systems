%% SYSTEMS MODELING AND SIMULATION
% Assignment 2 - Summer Semester 2020/2021
% Kavelidis Frantzis Dimitrios - AEM 9351 - kavelids@ece.auth.gr - ECE AUTH

% Gradient Descent Estimator
function xxdot = GraDE(t,xx)

    xxdot = zeros(size(xx));
    global a b c d am gamma
    
    % xx description
    x = xx(1);
    phi1 = xx(2);
    phi2 = xx(3);
    theta1est = xx(4);
    theta2est = xx(5);
    
    % Input u
    u = c*sin(d*t);
    
    % Error e
    e = x - (theta1est*phi1 + theta2est*phi2);
    
    % System description
    xxdot(1) = - a*x + b*u;
    xxdot(2) = - am*phi1 + x;
    xxdot(3) = - am*phi2 + u;
    xxdot(4) = gamma*e*phi1;
    xxdot(5) = gamma*e*phi2;
end