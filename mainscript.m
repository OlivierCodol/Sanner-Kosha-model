clear all; close all %#ok<CLALL>

hf(1) = figure(1); clf;

%% Define parameters

% Time variables------------------------------------
tm = .65;                               % Mouvement duration (sec)
t = 0:.004:(3*tm);                      % Simulation time (sec)

% Environment variables-----------------------------
tar_dist = .1;                          % Target distance from to start (m)
c = [-.1; .375];                        % Hand start position (m)
B = .5*[-10.1, -11.2; -11.2, 11.1];        % Force-field

% Arm variables-------------------------------------
l1 = .33;                               % Length of arm (m)
l2 = .34;                               % Length of forearm (m)
m1 = 1.93;                              % Mass of limb 1 & 2
m2 = 1.52;                              % (Katayama Kawato 1991: .9 & 1.1)
i1 = 0;%.0141;                              % Inertia for limb 1
i2 = 0;%.0188;                              % Inertia for limb 2
lc1 = .165;                             % Centre of gravity for limb 1
lc2 = .19;                              % Centre of gravity for limb 2
a = [m1*lc1*lc1 + m2*l1*l1 + i1;...     % Parameters for inertia and
     m2*lc2*lc2 + i2;...                % Coriolis matrices
     m2*l1*lc2];
 
% Get initial joint coordinates---------------------
tmp = acos((c(1).^2 + c(2).^2 - l1.^2 - l2.^2)./(2*l1*l2));
q0(1,1) = atan2(c(2),c(1)) - asin((l2.*sin(tmp)) ./ sqrt(sum(c.^2)));
q0(2,1) = q0(1) + tmp;

% Control Matrices----------------------------------
Lambda = [6.5 .064; .064 6.67];         % Feedback gain matrix
Kd = [2.3, .9; .9, 2.4];                % Joint viscosity matrix

% NN parameters-------------------------------------
node = sort(reshape((-5:8)'*(-3:9),[],1));
h = 2;
gamma = [repmat(.01,6,1); .04; .04];    % Learning rate
cmax = 75;                              % Bound for update rule



%% Initialise simulation

p = ws2struct;                          % Pack parameters into a struct
ntrial = 32;                            % Number of trials per target
ntar = 8;                               % Number of targets
cc = jet(ntrial);                       % Plotting colormap
err = NaN(ntrial,2,ntar);               % Error over trials
X = NaN(numel(t),2,ntrial,ntar);        % Actual trajectory performed
dX = NaN(numel(t),2,ntrial,ntar);       % Actual trajectory derivative


%% Run simulation

for itar = 1:ntar
    
    w = zeros(8,182,2);                     % NN initial weights
    p.w = w;
    
    for k = 1:ntrial
        
        fprintf(['Target ' num2str(itar) ' Trial ' num2str(k) ' \n'])
        
        % Define target
        tar_ang = 45*(itar-1);                      % Tar. angle (degrees)
        tar(1) = c(1) + tar_dist*cosd(tar_ang);     % Tar. cartesian coord.
        tar(2) = c(2) + tar_dist*sind(tar_ang);
        p.tar = tar;
        
        % Evaluate ODE solution
        % Init = [q; dq; qd; dqd; x; dx; xd; dxd];
        init = [q0; 0;0; q0; 0;0; c; 0;0; c; 0;0];  % Initial values
        tint = [t(1) t(end)];                       % Time interval
        odearm = @(t,in,p) simulate_reach(t,in,p);  % Define fcn handle
        mvt = ode45(odearm, tint , init,[],p);      % Solve mouvement ODE
        clear thismvt
        thismvt = deval(mvt,t);                     % Evaluate ODE
        
        % Save trajectory
        q(:,1:2) = thismvt(1:2,:)';
        dq(:,1:2) = thismvt(3:4,:)';
        qd(:,1:2) = thismvt(5:6,:)';
        dqd(:,1:2) = thismvt(7:8,:)';
        X(:,:,k,itar) = thismvt(9:10,:)';
        dX(:,:,k,itar) = thismvt(11:12,:)';
        dXd(:,1:2) = thismvt(15:16,:)';
        
        % Get parameters for update law
        s = (dqd - dq)' + Lambda*(qd - q)';
        [~,ix] = max(sum(abs(s)));
        err(k,:,itar) = s(:,ix)'; % saving for plotting later
        g = exp(-abs(h*q(ix,:) - node)/2);
        [ddqd,~] = second_derivative(t(ix), qd(ix,:)', dXd(ix,:)', p);
        dqr = dqd(ix,:)' + Lambda*(qd(ix,:) - q(ix,:))';
        ddqr = ddqd + Lambda*(dqd(ix,:) - dq(ix,:))';
        states = [ddqr; dq(ix,1)*dqr(1); dq(ix,1)*dqr(2);...
            dq(ix,2)*dqr(1); dq(ix,2)*dqr(2); dq(ix,:)'];
        
        % Apply update law
        for n = 1:2
            dw = -s(n,ix)*(gamma.*states)*g(:,n)';
            them = (~(w(:,:,n)<cmax)  & ~(dw<0)) |...
                (~(w(:,:,n)>-cmax) & ~(dw>0));
            dw(them) = 0;
            w(:,:,n) = w(:,:,n) - dw;
        end
        p.w = w;
        
        % Plot hand trajectory
        if mod(k-1,floor(ntrial/15))==0 || k==ntrial
            xplot = X(:,1,k,itar);
            yplot = X(:,2,k,itar);
            plot(xplot,yplot,'.','MarkerSize',2,'MarkerEdgeColor',cc(k,:));
            hold on
            plot(p.tar(1),p.tar(2),'ro'); hold on
            axis equal
            xlim([c(1)-tar_dist*1.3,   c(1)+tar_dist*1.3])
            ylim([c(2)-tar_dist*1.3,   c(2)+tar_dist*1.3])
        end
        
    end
end

xlabel('x position')
ylabel('y position')


%% Plot error

hf(2) = figure(2); clf
set(hf(2),'position',hf(1).Position + [hf(1).Position(3) 0 0 0])
if ntar > 4
    ncol = 4;
    nrow = 2;
else
    ncol = ntar;
    nrow = 1;
end
for itar = 1:ntar
    subplot(nrow,ncol,itar)
    plot(1:ntrial,err(:,1,itar),'r-'); hold on
    plot(1:ntrial,err(:,2,itar),'b-'); hold on
end

