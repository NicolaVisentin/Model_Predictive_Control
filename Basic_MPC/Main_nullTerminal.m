clear
close all
clc

addpath("functions\")
addpath("data\")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Trajectory tracking with Xf=0 and Vf=0                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load reference trajectory
% Refer to the corresponding script "Reference_generation14.m", "Reference_generation23.m" and 
% "Reference_simulation.m" for generation and simulation of the various references.

choose=4;   % choose reference trajectory (1, 2, 3 or 4)

switch choose
    case 1
        load('ref_data1.mat')   % 1 --> trajectory n2 in the report 
    case 2
        load('ref_data2.mat')   % 2 --> trajectory n3 in the report
    case 3
        load('ref_data3.mat')   % 3 --> useless, basically
    otherwise
        load('ref_data4.mat')   % 4 --> trajectory n1 in the report
end

%% Compute LTV model

% Extract data from ref.mat

x_ref=ref.x_ref;
u_ref=ref.u_ref;
t_ref=ref.t_ref;

tf_ref=t_ref(end);       
N_ref=length(t_ref)-1;
dt=tf_ref/N_ref;
nx=size(x_ref,2);
nu=size(u_ref,2);

% Define A and B matrices (they are functions of the linearisation point, 
% so they depend on x_ref,u_ref, so they depend on the time instant)

A=@(x_ref, u_ref) [1, 0, -u_ref(1)*sin(x_ref(3)*dt); 0, 1, u_ref(1)*cos(x_ref(3)*dt); 0, 0, 1 ];
B=@(x_ref) [cos(x_ref(3))*dt 0; sin(x_ref(3))*dt 0; 0 dt];

% Compute LTV matrices (allocate them in a 3D matrix)

A_ref=zeros(nx,nx,N_ref);
B_ref=zeros(nx,nu,N_ref);

for k=1:N_ref
    A_ref(:,:,k)=A(x_ref(k,:)',u_ref(k,:)');
    B_ref(:,:,k)=B(x_ref(k,:)');
end
LTV.A_ref=A_ref;
LTV.B_ref=B_ref;
LTV.C=eye(nx);

% Define constraints

v_max=0.55;            % max linear speed
v_min=0.2;             % min linear speed
w_lim=1;               % max yaw rate (both directions)  
e_lim=[0.5 0.5 0.5]';  % max "state" (error) (both positive and negative)

%% MPC dense formulation (prediction model and cost function)
% ! Everything denoted by "ref" concerns the "virtual robot" behaviour
% computed for a certain time t_ref (containing N_ref discrete timesteps). 
% The prediction horizon is defined over N timesteps.
% The simulation horizon will be made of N_sim=N_ref-N time instants.

% Define parameters

N=20;                                                 % prediction horizon
Q=diag([10/e_lim(1)^2 10/e_lim(2)^2 10/e_lim(3)^2]);  % state stage weight matrix
R=diag([1/v_max^2 1/w_lim^2]);                        % control stage weight matrix
P=diag([0 0 0]);                                      % final weight matrix

% Compute predicted model and cost function matrices

[T_sim,S_sim]=PredMatGen(LTV,N);                
[H_sim,h_sim]=CostFuncGen(T_sim,S_sim,Q,R,P); 

%% MPC simulation
% ! Remember that the MPC is formulated on the error's dynamics, so the
% "state" for the MPC will be the error e=x-x_ref and the "control" will be
% ub=u-u_ref.

N_sim=N_ref-N;             % simulation time instants
e0=[-0.05 0.01 -0.1]';     % initial "state" (error)

% Constraints definition

u_max=[v_max w_lim]';
u_min=[v_min -w_lim]';
Ub_ub=repmat(u_max,N_ref,1)-Umat2Uvect(u_ref);  % prepare whole sequence (N_ref) of ub upper bounds (input constraints)
Ub_lb=repmat(u_min,N_ref,1)-Umat2Uvect(u_ref);  % prepare whole sequence (N_ref) of ub lower bounds (input constraints)

b_tilde=[e_lim; -(-e_lim)];  % prepare matrices for "state" (error) constraints
b_bar=repmat(b_tilde,N,1);   % ...

% MPC iterations

e_hist=zeros(N_sim+1,nx);  % initialize "state" (error) history
ub_hist=zeros(N_sim,nu);   % initialize "control" (ub) history
e_hist(1,:)=e0';           % initial condition in e_hist

optons=optimoptions('quadprog','Display','none');
for k=1:N_sim    % WITH QUADPROG

    % Extract current "state" (error)
    e=e_hist(k,:)';

    % Extract current matrices
    H=H_sim(:,:,k);
    h=h_sim(:,:,k);
    T=T_sim(:,:,k);
    S=S_sim(:,:,k);
    A_ref_pred=A_ref(:,:,k:k+N-1);
    B_ref_pred=B_ref(:,:,k:k+N-1);

    % Generate constraints matrices
    U_LB=Ub_lb((k-1)*nu+1:(k-1+N)*nu);  % estract lower bounds for ub in the current prediction horizon (input constraints)
    U_UB=Ub_ub((k-1)*nu+1:(k-1+N)*nu);  % estract upper bounds for ub in the current prediction horizon (input constraints)
    A_eq=S(end-nx+1:end,:);             % terminal constraint for the final error e(N+1)=0 in the current horizon (terminal constraint)
    b_eq=-T(end-nx+1:end,:)*e;          % ...
    A_temp=squeeze(num2cell(cat(1,A_ref_pred,-A_ref_pred),[1 2]));   % generate matrices for inequality constraints ("state" (error) constraints)
    A_bar=blkdiag(A_temp{:});                                        % ...
    B_temp=squeeze(num2cell(cat(1,B_ref_pred,-B_ref_pred),[1 2]));   % ...
    B_bar=blkdiag(B_temp{:});                                        % ...
    A_ineq=A_bar*S(1:end-nx,:)+B_bar;                                % ...
    b_ineq=b_bar-A_bar*T(1:end-nx,:)*e;                              % ...

    % quadprog solver for optimization
    Ub_opt=quadprog(H,h*e,A_ineq,b_ineq,A_eq,b_eq,U_LB,U_UB,[],optons);

    % Save first control step in ub_hist and update state with it
    ub_hist(k,:)=Ub_opt(1:nu)';

    x_curr=x_ref(k,:)'+e;
    u_curr=ub_hist(k,:)'+u_ref(k,:)';
    [~,x_next]=ode45(@(t,x) NonlinContModel(t,x,[0 dt],[u_curr u_curr]'),[0 dt],x_curr);   % apply constant input for dt to the nonlinear system
    x_next=x_next(end,:)';
    x_next(end)=wrapToPi(x_next(end));
    e_hist(k+1,1:2)= x_next(1:2)'-x_ref(k+1,1:2);
    e_hist(k+1,3)=mod(x_next(3)-x_ref(k+1,3)+pi,2*pi)-pi;  % to make sure the heading error is always between [-pi,pi and with correct sign]

    % Display progress
    clc
    fprintf('Simulating MPC: %.1f %%',k/N_sim*100)

end
fprintf('\n\n')

% Reconstruct "real" states and controls and save

x_hist=x_ref(1:N_sim+1,:)+e_hist;
u_hist=u_ref(1:N_sim,:)+ub_hist;
t_sim=t_ref(1:N_sim+1);

% res.x_hist=x_hist;
% res.u_hist=u_hist;
% res.t_sim=t_sim;
% switch choose
%     case 1
%         save('data\result1_Xf_null.mat','res');
%     case 2
%         save('data\result2_Xf_null.mat','res');
%     case 3
%         save('data\result3_Xf_null.mat','res');
%     otherwise
%         save('data\result4_Xf_null.mat','res');
% end

% % Compute region of attraction XN (only for one ref, just for comparison)
% 
% [Xn,V,~]=FindXn(A_ref(:,:,1),B_ref(:,:,1),[],N,-e_lim,e_lim,u_min-u_ref(1,:)',u_max-u_ref(1,:)','termeq','linpr');
% XN=Xn{end};
% XN=Polyhedron(XN.A,XN.b);
% XN_proj=XN.projection([1 2]);
% figure
% XN_proj.plot('color','red')
% xlabel('$e_x$','Interpreter','latex') 
% ylabel('$e_y$','Interpreter','latex')
% grid on; box on
% title('X_N projection (example for the first iteration)')
% 
% figure
% XN.plot('color','red')
% xlabel('$e_x$','Interpreter','latex') 
% ylabel('$e_y$','Interpreter','latex')
% zlabel('$e_\theta$','Interpreter','latex')
% grid on; box on
% title('X_N (example for the first iteration)')

%% Plots

% Control input u

figure
subplot(211)
hold on
stairs(t_sim(1:end-1),u_hist(:,1),'b')
plot(t_sim(1:end-1),u_ref(1:N_sim,1),'b--')
grid on; box on
xlabel('t [s]','Interpreter','latex')
ylabel('v [m/s]','Interpreter','latex')
title('Robot control input: v')
legend('$v$','$v_{ref}$','interpreter','latex','location','best')
hold off

subplot(212)
hold on
stairs(t_sim(1:end-1),u_hist(:,2))
plot(t_sim(1:end-1),u_ref(1:N_sim,2),'b--')
grid on; box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\omega$ [rad/s]','Interpreter','latex')
title('Robot control input: \omega')
legend('$\omega$','$\omega_{ref}$','interpreter','latex','location','best')
hold off

% States x

figure
subplot(211)
hold on
plot(t_sim,x_hist(:,1),'b-')
plot(t_sim,x_hist(:,2),'r-')
plot(t_sim,x_ref(1:N_sim+1,1),'b--')
plot(t_sim,x_ref(1:N_sim+1,2),'r--')
grid on; box on
xlabel('t [s]','Interpreter','latex')
ylabel('x, y [m]','Interpreter','latex')
title('Robot state: position (x,y)')
legend('$x$','$y$','$x_{ref}$','$y_{ref}$','interpreter','latex','location','best')
hold off

subplot(212)
hold on
plot(t_sim,x_hist(:,3),'b-')
plot(t_sim,x_ref(1:N_sim+1,3),'b--')
grid on; box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\theta$ [rad]','Interpreter','latex')
title('Robot state: heading \theta')
legend('$\theta$','$\theta_{ref}$','interpreter','latex','location','best')
hold off

% Trajectory (x,y)

figure
hold on
plot(x_hist(1,1),x_hist(1,2),'bo','MarkerFaceColor','b')
plot(x_hist(:,1),x_hist(:,2),'b')
plot(x_ref(1,1),x_ref(1,2),'bo')
plot(x_ref(1:N_sim+1,1),x_ref(1:N_sim+1,2),'b--')
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
stairs(t_sim(1:end-1),u_hist(:,1)-u_ref(1:N_sim,1),'b')
grid on; box on
xlabel('t [s]','Interpreter','latex')
ylabel('$u_{b,1}$ [m/s]','Interpreter','latex')
title('MPC "control input": v-v_r_e_f')

subplot(212)
stairs(t_sim(1:end-1),u_hist(:,2)-u_ref(1:N_sim,2),'b')
grid on; box on
xlabel('t [s]','Interpreter','latex')
ylabel('$u_{b,2}$ [rad/s]','Interpreter','latex')
title('MPC "control input": \omega-\omega_r_e_f')

% Error (what MPC sees as states): e=x-x_ref

figure
subplot(211)
hold on
plot(t_sim,x_hist(:,1)-x_ref(1:N_sim+1,1),'b-')
plot(t_sim,x_hist(:,2)-x_ref(1:N_sim+1,2),'r-')
grid on; box on
xlabel('t [s]','Interpreter','latex')
ylabel('$e_x \; ,\; e_y$ [m]','Interpreter','latex')
title('MPC "state": position error')
legend('$e_x$','$e_y$','interpreter','latex','location','best')
hold off

subplot(212)
plot(t_sim,x_hist(:,3)-x_ref(1:N_sim+1,3),'b-')
grid on; box on
xlabel('t [s]','Interpreter','latex')
ylabel('$e_\theta$ [rad]','Interpreter','latex')
title('MPC "state": heading error')

% Error trajectory

figure
hold on
plot3(x_hist(1,1)-x_ref(1,1),x_hist(1,2)-x_ref(1,2),x_hist(1,3)-x_ref(1,3),'bo','MarkerFaceColor','b')
plot3(x_hist(:,1)-x_ref(1:N_sim+1,1),x_hist(:,2)-x_ref(1:N_sim+1,2),x_hist(:,3)-x_ref(1:N_sim+1,3),'b')
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