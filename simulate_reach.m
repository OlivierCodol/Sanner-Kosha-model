function xdot = simulate_reach(t,x,p)

% Unpack required parameters
l1 = p.l1;
l2 = p.l2;
a = p.a;
tm = p.tm;
tar = p.tar;
c = p.c;
Lambda = p.Lambda;
Kd = p.Kd;
B = p.B;
w = p.w;
node = p.node;
h = p.h;

% Unpack ODE variables
q(:,1) = [x(1);x(2)];
dq(:,1) = [x(3);x(4)];
dX(:,1) = [x(11);x(12)];
Xd(:,1) = [x(13);x(14)];



%% GET DESIRED TRAJECTORY

% Get desired current joint coordinates
tmp = acos((Xd(1).^2 + Xd(2).^2 - l1.^2 - l2.^2)./(2*l1*l2));
qd(1,1) = atan2(Xd(2),Xd(1)) - asin((l2.*sin(tmp)) ./ sqrt(sum(Xd.^2)));
qd(2,1) = qd(1) + tmp;

% Get 1st & 2nd derivatives of desired cartesian coordinates
if t<=tm
    dXd(1,1) = (tar(1)-c(1))*(3*10/(tm^3)*(t)^2-(4*15/(tm^4))*(t)^3+(5*6/(tm^5))*(t)^4);
    dXd(2,1) = (tar(2)-c(2))*(3*10/(tm^3)*(t)^2-(4*15/(tm^4))*(t)^3+(5*6/(tm^5))*(t)^4);
    ddXd(1,1) = (tar(1)-c(1))*(2*3*10/(tm^3)*(t)-(3*4*15/(tm^4))*(t)^2+(4*5*6/(tm^5))*(t)^3);
    ddXd(2,1) = (tar(2)-c(2))*(2*3*10/(tm^3)*(t)-(3*4*15/(tm^4))*(t)^2+(4*5*6/(tm^5))*(t)^3);
else
    dXd = zeros(2,1);
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
    



%% GET CONTROL LAW

e = qd - q;                 % Error
de = dqd - dq;              % Derivative of error
dqr = dqd + Lambda*e;        % Filtered error
ddqr = ddqd + Lambda*de;     % Derivative of filtered error
s = (dqd - dq) + Lambda*e;

% Get dynamics matrices
V = a(3)*sin(q(2))*[-dq(2) -sum(dq); dq(1) 0];
H(1,1) = a(1) + 2*a(2)*cos(q(2)) + a(3);
H(1,2) = a(2) + a(3)*cos(q(2));
H(2,1) = H(1,2);
H(2,2) = a(2);

% Get neural network contribution
states = [ddqr;dq(1)*dqr(1);dq(1)*dqr(2);dq(2)*dqr(1);dq(2)*dqr(2);dq];
g = exp(-abs(h*q' - node)/2);
T_nn(1,1) = (w(:,:,1)*g(:,1))'*states;
T_nn(2,1) = (w(:,:,2)*g(:,2))'*states;

% Get desired torque
% Kd = zeros(2,2);
T = Kd*s + H*ddqr + V*dqr + T_nn;


    
%% RUN THE PLANT

% Get jacobian of joint to cartesian coordinates transformation
J = [-l1*sin(q(1)), -l2*sin(q(2)); l1*cos(q(1)), l2*cos(q(2))];
% Get derivative of the jacobian (dJ/dt   =   dJ/dQ * dQ/dt)
dJ = [-l1*dq(1)*cos(q(1)), -l2*dq(2)*cos(q(2));
    -l1*dq(1)*sin(q(1)), -l2*dq(2)*sin(q(2))];

% Apply force field
% B = zeros(2,2);
T = T + J'*B*J*dq;

% Apply control command to the plant
ddX = J*H^-1*(T-V*J^-1*dX) + dJ*J^-1*dX;
ddq = J\(ddX - dJ*dq);


%% PACK OUTPUT

xdot(1,1) = dq(1);
xdot(2,1) = dq(2);
xdot(3,1) = ddq(1);
xdot(4,1) = ddq(2);
xdot(5,1) = dqd(1);
xdot(6,1) = dqd(2);
xdot(7,1) = ddqd(1);
xdot(8,1) = ddqd(2);
xdot(9,1) = dX(1);
xdot(10,1) = dX(2);
xdot(11,1) = ddX(1);
xdot(12,1) = ddX(2);
xdot(13,1) = dXd(1);
xdot(14,1) = dXd(2);
xdot(15,1) = ddXd(1);
xdot(16,1) = ddXd(2);

end