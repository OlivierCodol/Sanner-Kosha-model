function [ddqd,ddXd] = second_derivative(t,qd,dXd,p)

% Unpack parameters
tar = p.tar;
c = p.c;
tm = p.tm;
l1 = p.l1;
l2 = p.l2;


% Get 1st & 2nd derivatives of desired cartesian coordinates
if t<=tm
    ddXd(1,1) = (tar(1)-c(1))*(2*3*10/(tm^3)*(t)-(3*4*15/(tm^4))*(t)^2+(4*5*6/(tm^5))*(t)^3);
    ddXd(2,1) = (tar(2)-c(2))*(2*3*10/(tm^3)*(t)-(3*4*15/(tm^4))*(t)^2+(4*5*6/(tm^5))*(t)^3);
else
    ddXd = zeros(2,1);
end



% Get Jacobian
J = [-l1*sin(qd(1)), -l2*sin(qd(2)); l1*cos(qd(1)), l2*cos(qd(2))];

% Get derivative of desired joint coordinates
dqd = J\dXd;

% Get derivative of the jacobian (dJ/dt   =   dJ/dQ * dQ/dt)
dJ = [-l1*dqd(1)*cos(qd(1)), -l2*dqd(2)*cos(qd(2));
    -l1*dqd(1)*sin(qd(1)), -l2*dqd(2)*sin(qd(2))];




% get second derivative of desired joint coordinates
ddqd = J\(ddXd - dJ*dqd);

end