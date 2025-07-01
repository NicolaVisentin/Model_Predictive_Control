clear
close all
clc

addpath("functions\")
addpath("data\")

%% Nonlinear system and reference trajectory

% Load reference trajectory

load('ref_data4.mat')
t_ref=ref.t_ref;
x_ref=ref.x_ref;
u_ref=ref.u_ref;

tf_ref=t_ref(end);       % final time where reference has been computed
N_ref=length(t_ref)-1;   % discrete time steps in the reference
dt=tf_ref/N_ref;         % time resolution for the controller (30 Hz)
nx=size(x_ref,2);    
nu=size(u_ref,2);

% Initial state of the robot

e0=[-0.05 0.01 -0.1]';  % initial error wrt reference
IC=x_ref(1,:)'+e0;      % initial state

% Contraints

v_max=0.55;               % max linear velocity
v_min=0.2;                % min linear velocity
w_max=1;                  % max yaw rate
w_min=-1;                 % min yaw rate

e_ub=[0.5 0.5 0.5]';      % error upper bound
e_lb=[-0.5 -0.5 -0.5]';   % error lower bound

%% MPC controller

% Discrete linearized model for the MPC

A=@(x_ref, u_ref) [1, 0, -u_ref(1)*sin(x_ref(3)*dt); 0, 1, u_ref(1)*cos(x_ref(3)*dt); 0, 0, 1 ];
B=@(x_ref) [cos(x_ref(3))*dt 0; sin(x_ref(3))*dt 0; 0 dt];

A_ref=zeros(nx,nx,N_ref);
B_ref=zeros(nx,nu,N_ref);

for k=1:N_ref
    A_ref(:,:,k)=A(x_ref(k,:)',u_ref(k,:)');
    B_ref(:,:,k)=B(x_ref(k,:)');
end

LTV.A_ref=A_ref;
LTV.B_ref=B_ref;
LTV.C=eye(nx);

% Optimization horizon and weights

N=20;                                           % prediction horizon
Q=diag([1/e_ub(1)^2 1/e_ub(2)^2 1/e_ub(3)^2]);  % state stage weight matrix
R=diag([10/v_max^2 10/w_max^2]);                % control stage weight matrix

% Constraints definition

b_tilde=[e_ub; -e_lb];       % prepare matrices for "state" (error) constraints
b_bar=repmat(b_tilde,N,1);   % ...

u_max=[v_max w_max]';
u_min=[v_min w_min]';
ub_ub=u_max-u_ref(1,:)';                        % ub upper bound
ub_lb=u_min-u_ref(1,:)';                        % ub lower bound
Ub_ub=repmat(u_max,N_ref,1)-Umat2Uvect(u_ref);  % prepare whole sequence (N_ref) of ub upper bounds (input constraints)
Ub_lb=repmat(u_min,N_ref,1)-Umat2Uvect(u_ref);  % prepare whole sequence (N_ref) of ub lower bounds (input constraints)

% Terminal cost computation

for k=1:50:N_ref

    % Extract matrices
    A=A_ref(:,:,k);
    B=B_ref(:,:,k);

    % Solve Riccati equation
    [K,~]=dlqr(A,B,Q,R);
    K=-K;

    % Find Xf
    [Xf,Z]=FindXf(A,B,K,e_lb,e_ub,ub_lb,ub_ub,'lqr','linpr');
    Xf=Polyhedron(Xf.A,Xf.b);

    % "Update" the intersection
    if k==1
        Xf_inter=Xf;
    else
        Xf_inter = Xf & Xf_inter;
    end

    clc
    fprintf('Computing Xf: %.1f%%',k/N_ref*100)

end
clc
fprintf('Be patient...\n\n')
Xf=Xf_inter.minHRep();
clc

% Xf_proj=Xf.projection([1 2]);
% figure
% Xf_proj.plot('color','blue')
% xlabel('$e_x$','Interpreter','latex') 
% ylabel('$e_y$','Interpreter','latex')
% grid on; box on
% title('Final X_f (projection)')

H_ter=Xf.A;                  % to send to simulink
h_ter=Xf.b;                  % to send to simulink
% terminalXf.H=H_ter;          % save Xf (for simulations in ROS, to avoid downloading MPT3 also in the other OS)
% terminalXf.h=h_ter;          % ...
% save('Xf.mat','terminalXf')  % ...

% Terminal weight computation

vertices=Xf.V;
nv=size(vertices,1);
maxw=-inf;
for k=1:N_ref

    % Extract matrices
    A=A_ref(:,:,k);
    B=B_ref(:,:,k);

    % Solve Riccati equation and solve maximization
    [~,P_xi]=dlqr(A,B,Q,R);
    L=chol(P_xi,'lower');
    for ii=1:nv
        w=L*vertices(ii,:)';
        if w'*w >= maxw
            maxw=w'*w;
            P=P_xi;
        end
    end

    clc
    fprintf('Computing Vf: %.1f%%',k/N_ref*100)

end
fprintf('\n\n')

% save('Vf','P')  % save Vf (for simulations in Ubuntu)

% Prediction and cost matrices

[T_sim,S_sim]=PredMatGen(LTV,N);    
[H_sim,h_sim]=CostFuncGen(T_sim,S_sim,Q,R,P);

%% Extended Kalman filter

% Noise and covariances

d=0.1;

R_obs=diag([1 1 1]);   % noise covariance
Qw=diag([1 1]);        % disturbance covariance
Qp=10;                 % artificial noise covariance for d estimation
Q_obs=blkdiag(Qw,Qp);  % augmented disturbance covariance

% Filter initialisation

IC_obs_state=IC+[-0.5 -0.5 0.1]';  % initial guess for the filter (state estimation)
IC_obs_d=-0.1;                     % initial guess for the filter (offset estimation)
IC_obs=[IC_obs_state; IC_obs_d];   % initial guess for the filter

%% Simulink simulation

% Set up simulation

N_sim=N_ref-N;       % simulation duration in terms od time steps
tf_sim=N_sim*dt;     % simulation final time
t_sim=0:dt:tf_sim;   % time vector for the controller
time_scale=100;      % time accuracy of the simulation (dt_simulink=dt/time_scale)

% Launch simulation

out=sim("Simulink.slx");

%% Simulation results (controller)

% Extract results

t_SIM=out.tout;
u_SIM=squeeze(out.u)';
u_ref_SIM=squeeze(out.u_ref)';
x_SIM=squeeze(out.x);
x_ref_SIM=squeeze(out.x_ref)';
x_hat_SIM=squeeze(out.x_hat)';
x_meas_SIM=squeeze(out.x_meas);

e_SIM=x_SIM(1:time_scale:end,:)-x_ref_SIM(1:time_scale:end,:);

% Control input u

figure
subplot(211)
hold on
plot(t_SIM,u_SIM(:,1),'b')
plot(t_SIM,u_ref_SIM(:,1),'b--')
grid on; box on
xline(0)
yline(0)
xlabel('t [s]','Interpreter','latex')
ylabel('v [m/s]','Interpreter','latex')
title('Robot control inputs')
legend('$v$','$v_{ref}$','interpreter','latex','location','best')
hold off

subplot(212)
hold on
plot(t_SIM,u_SIM(:,2),'b')
plot(t_SIM,u_ref_SIM(:,2),'b--')
grid on; box on
xline(0)
yline(0)
xlabel('t [s]','Interpreter','latex')
ylabel('$\omega$ [rad/s]','Interpreter','latex')
legend('$\omega$','$\omega_{ref}$','interpreter','latex','location','best')
hold off

% States x

figure
subplot(311)
hold on
plot(t_SIM,x_SIM(:,1),'b')
plot(t_SIM,x_ref_SIM(:,1),'b--')
grid on; box on
xline(0)
yline(0)
xlabel('t [s]','Interpreter','latex')
ylabel('x [m]','Interpreter','latex')
title('Robot states')
legend('$x$','$x_{ref}$','interpreter','latex','location','best')
hold off

subplot(312)
hold on
plot(t_SIM,x_SIM(:,2),'b')
plot(t_SIM,x_ref_SIM(:,2),'b--')
grid on; box on
xline(0)
yline(0)
xlabel('t [s]','Interpreter','latex')
ylabel('$y$ [m]','Interpreter','latex')
legend('$y$','$y_{ref}$','interpreter','latex','location','best')
hold off

subplot(313)
hold on
plot(t_SIM,x_SIM(:,3),'b')
plot(t_SIM,x_ref_SIM(:,3),'b--')
grid on; box on
xline(0)
yline(0)
xlabel('t [s]','Interpreter','latex')
ylabel('$\theta$ [rad]','Interpreter','latex')
legend('$\theta$','$\theta_{ref}$','interpreter','latex','location','best')
hold off

% Trajectory (x,y)

figure
hold on
plot(x_SIM(1,1),x_SIM(1,2),'bo','MarkerFaceColor','b')
plot(x_SIM(:,1),x_SIM(:,2),'b')
plot(x_ref_SIM(1,1),x_ref_SIM(1,2),'bo')
plot(x_ref_SIM(:,1),x_ref_SIM(:,2),'b--')
xline(0)
yline(0)
grid on; box on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('Robot trajectory')
legend('actual starting point','traj','reference starting point','reference traj','interpreter','latex','location','best')
hold off

% "Control input" (what MPC sees as control): ub=u-u_ref

figure
subplot(211)
plot(t_SIM,u_SIM(:,1)-u_ref_SIM(:,1),'b')
grid on; box on
xline(0)
yline(0)
xlabel('t [s]','Interpreter','latex')
ylabel('$v-v_{ref}$ [m/s]','Interpreter','latex')
title('MPC "control input" (u_b)')

subplot(212)
plot(t_SIM,u_SIM(:,2)-u_ref_SIM(:,2),'b')
grid on; box on
xline(0)
yline(0)
xlabel('t [s]','Interpreter','latex')
ylabel('$\omega-\omega_{ref}$ [rad/s]','Interpreter','latex')

% Error (what MPC sees as states): e=x-x_ref

figure
subplot(311)
plot(t_SIM(1:time_scale:end),e_SIM(:,1),'b')
grid on; box on
xline(0)
yline(0)
xlabel('t [s]','Interpreter','latex')
ylabel('$e_x$ [m]','Interpreter','latex')
title('MPC "state" (error)')

subplot(312)
plot(t_SIM(1:time_scale:end),e_SIM(:,2),'b')
grid on; box on
xline(0)
yline(0)
xlabel('t [s]','Interpreter','latex')
ylabel('$e_y$ [m]','Interpreter','latex')

subplot(313)
plot(t_SIM(1:time_scale:end),e_SIM(:,3),'b')
grid on; box on
xline(0)
yline(0)
xlabel('t [s]','Interpreter','latex')
ylabel('$e_\theta$ [rad]','Interpreter','latex')

% Error trajectory

figure
hold on
plot3(e_SIM(1,1),e_SIM(1,2),e_SIM(1,3),'bo','MarkerFaceColor','b')
plot3(e_SIM(:,1),e_SIM(:,2),e_SIM(:,3),'b')
grid on; box on
plot3([0 0.05],[0 0],[0 0],'k')
plot3([0 0],[0 0.05],[0 0],'k')
plot3([0 0],[0 0],[0 0.1],'k')
xlabel('$e_x$','Interpreter','latex')
ylabel('$e_y$','Interpreter','latex')
zlabel('$e_\theta$','Interpreter','latex')
title('MPC trajectory: error trajectory')
view(3)
hold off

%% Simulation results (filter)

figure
subplot(311)
hold on
plot(t_SIM,x_SIM(:,1),'b')
plot(t_SIM,x_hat_SIM(:,1),'r')
plot(t_SIM,x_meas_SIM(:,1),'Color',[0 1 0 0.4])
grid on; box on
xline(0)
yline(0)
xlabel('t [s]','Interpreter','latex')
ylabel('$x_1$ [m]','Interpreter','latex')
title('EKF: state vs observed state vs measured state')
legend('$x$','$\hat{x}$','$x_{meas}$','Interpreter','latex','Location','best')
hold off

subplot(312)
hold on
plot(t_SIM,x_SIM(:,2),'b')
plot(t_SIM,x_hat_SIM(:,2),'r')
plot(t_SIM,x_meas_SIM(:,2),'Color',[0 1 0 0.4])
grid on; box on
xline(0)
yline(0)
xlabel('t [s]','Interpreter','latex')
ylabel('$x_2$ [m]','Interpreter','latex')
legend('$y$','$\hat{y}$','$y_{meas}$','Interpreter','latex','Location','best')
hold off

subplot(313)
hold on
plot(t_SIM,x_SIM(:,3),'b')
plot(t_SIM,x_hat_SIM(:,3),'r')
plot(t_SIM,x_meas_SIM(:,3),'Color',[0 1 0 0.4])
grid on; box on
xline(0)
yline(0)
xlabel('t [s]','Interpreter','latex')
ylabel('$x_3$ [m]','Interpreter','latex')
legend('$\theta$','$\hat{\theta}$','$\theta_{meas}$','Interpreter','latex','Location','best')
hold off

% Offset estimation

figure
hold on
plot(t_SIM,d(1)*ones(length(t_SIM),1),'b')
plot(t_SIM,x_hat_SIM(:,4),'r')
grid on; box on
xlabel('t [s]','Interpreter','latex')
ylabel('$d$','Interpreter','latex')
title('EKF: disturbance estimation')
legend('actual $d$','$\hat{d}$','Interpreter','latex','Location','best')
hold off