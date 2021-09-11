%% SYSTEMS MODELING AND SIMULATION
% Assignment 2 - Summer Semester 2020/2021
% Kavelidis Frantzis Dimitrios - AEM 9351 - kavelids@ece.auth.gr - ECE AUTH

% Lyapunov Method / Parallel Configuration Estimator
function xxdot = LyapEst(t,xx)

    global A b u1 u2 w1 w2 gamma
    
    % Input u
    u = u1*sin(w1*t) + u2*sin(w2*t);
    
    % xx description
    x1 = xx(1);
    x2 = xx(2);
    x1est = xx(3);
    x2est = xx(4);
    a11est = xx(5);
    a12est = xx(6);
    a21est = xx(7);
    a22est = xx(8);
    b1est = xx(9);
    b2est = xx(10);
  
    xxdot = zeros(size(xx));
    
    % Error
    e1 = x1 - x1est;
    e2 = x2 - x2est;
    
    % System description
    xxdot(1) = A(1,1)*x1 + A(1,2)*x2 + b(1)*u;          % x1dot
    xxdot(2) = A(2,1)*x1 + A(2,2)*x2 + b(2)*u;          % x2dot
    xxdot(3) = a11est*x1est + a12est*x2est + b1est*u;   % x1^dot
    xxdot(4) = a21est*x1est + a22est*x2est + b2est*u;   % x2^dot
    xxdot(5) = gamma(1)*x1est*e1;                       % a11^dot
    xxdot(6) = gamma(1)*x2est*e1;                       % a12^dot
    xxdot(7) = gamma(1)*x1est*e2;                       % a21^dot
    xxdot(8) = gamma(1)*x2est*e2;                       % a22^dot
    xxdot(9) = gamma(2)*u*e1;                           % b1^dot
    xxdot(10) = gamma(2)*u*e2;                          % b2^dot
end