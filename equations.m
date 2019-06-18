%% Falling without slipping
clear;
syms R theta2 theta theta1 ue g ii

eqn = -R*(theta2*sin(theta)+theta1^2*cos(theta)) == ue*(g+R*(theta2*cos(theta)-theta1^2*sin(theta)))*ii;

eqn2 = isolate(eqn, theta2);

pretty(eqn2);

%% Falling with slipping
clear;
syms theta2 m R Ig g R theta theta1 uc ii;

eqn = theta2 == (m*R/Ig)*(g+R*(theta2*cos(theta) - theta1^2*sin(theta)))*(uc*sin(theta)*ii-cos(theta));

eqn2 = isolate(eqn, theta2);

pretty(eqn2);

