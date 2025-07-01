clear
close all
clc

addpath("functions\")
addpath("data\")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Trajectory tracking with non-null Vf and Xf                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load reference trajectory
% Refer to the corresponding script "Reference_generation14.m", "Reference_generation23.m" and 
% "Reference_simulation.m" for generation and simulation of the various references.

load('ref_data4.mat')

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

% Constraints

v_max=0.55;                % max linear speed
v_min=0.2;                 % min linear speed
w_max=1;                   % max yaw rate
w_min=-1;                  % min yaw rate
e_ub=[0.5 0.5 0.5]';       % max "state" (error)
e_lb=[-0.5 -0.5 -0.5]';    % min "state" (error)

u_max=[v_max w_max]';
u_min=[v_min w_min]';
ub_ub=u_max-u_ref(1,:)';
ub_lb=u_min-u_ref(1,:)';

%% Compute terminal constraint Xf and terminal weight Vf
% The idea is to compute all Xf and then find the biggest set which is
% invariant for all the configurations.
% Concerning Vf, we want to consider che worst case scenario (higher cost).

% Horizon and optimization weights

N=6;                                               % prediction horizon
Q=diag([10/e_ub(1)^2 10/e_ub(2)^2 10/e_ub(3)^2]);  % state stage weight matrix
R=diag([1/v_max^2 1/w_max^2]);                     % control stage weight matrix

% Compute Xf

for k=1:50:N_ref

    % Extract matrices
    A=A_ref(:,:,k);
    B=B_ref(:,:,k);

    % Solve Riccati equation
    [K,~]=dlqr(A,B,Q,R);
    K=-K;

    % Find Xf
    [Xf,Z]=FindXf(A,B,K,e_lb,e_ub,ub_lb,ub_ub,'lqr','linpr');       % only computes Xf
    %[Xn,~,Z]=FindXn(A,B,K,N,e_lb,e_ub,ub_lb,ub_ub,'lqr','linpr');   % also computes all Xn
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

% % Plot 2D projection of Xf and 3D Xf
% 
% Xf_proj=Xf.projection([1 2]);
% figure
% Xf_proj.plot('color','blue')
% xlabel('$e_x$','Interpreter','latex') 
% ylabel('$e_y$','Interpreter','latex')
% grid on; box on
% title('Final X_f (projection)')
% 
% figure
% Xf.plot('color','blue')
% xlabel('$e_x$','Interpreter','latex') 
% ylabel('$e_y$','Interpreter','latex')
% zlabel('$e_\theta$','Interpreter','latex')
% grid on; box on
% title('Final X_f')

clc

% Compute Vf

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

%% MPC simulation
% ! Remember that the MPC is formulated on the error's dynamics, so the
% "state" for the MPC will be the error e=x-x_ref and the "control" will be
% ub=u-u_ref.

% Build prediction/cost matrices

N_sim=N_ref-N;          % simulation time instants
e0=[-0.05 0.01 -0.1]';  % initial "state" (error)

[T_sim,S_sim]=PredMatGen(LTV,N);    
[H_sim,h_sim]=CostFuncGen(T_sim,S_sim,Q,R,P);

% Prepare contraints matrices

Ub_ub=repmat(u_max,N_ref,1)-Umat2Uvect(u_ref);  % prepare whole sequence (N_ref) of ub upper bounds (input constraints)
Ub_lb=repmat(u_min,N_ref,1)-Umat2Uvect(u_ref);  % prepare whole sequence (N_ref) of ub lower bounds (input constraints)

b_tilde=[e_ub; -e_lb];       % prepare matrices for "state" (error) constraints
b_bar=repmat(b_tilde,N,1);   % ...

% MPC iterations

e_hist=zeros(N_sim+1,nx);  % initialize "state" (error) history
ub_hist=zeros(N_sim,nu);   % initialize "control" (ub) history
e_hist(1,:)=e0';           % initial condition in e_hist

options=optimoptions('quadprog','Display','none');
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
    A_temp=squeeze(num2cell(cat(1,A_ref_pred,-A_ref_pred),[1 2]));   % generate matrices for inequality constraints ("state" (error) constraints)
    A_bar=blkdiag(A_temp{:});                                        % ...
    B_temp=squeeze(num2cell(cat(1,B_ref_pred,-B_ref_pred),[1 2]));   % ...
    B_bar=blkdiag(B_temp{:});                                        % ...
    A_ineq1=A_bar*S(1:end-nx,:)+B_bar;                               % ...
    b_ineq1=b_bar-A_bar*T(1:end-nx,:)*e;                             % ...
    A_ineq2=Xf.A*S(end-nx+1:end,:);         % terminal constraint for the final error e(N) in the current horizon
    b_ineq2=Xf.b-Xf.A*T(end-nx+1:end,:)*e;  % ...
    A_ineq=[A_ineq1; A_ineq2];   % "total" inequality constraints
    b_ineq=[b_ineq1; b_ineq2];   % ...

    % quadprog solver for optimization
    Ub_opt=quadprog(H,h*e,A_ineq,b_ineq,[],[],U_LB,U_UB,[],options);

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
% save('data\result4_Xf_lqr.mat','res');

% % Compute region of attraction XN (only for one ref, just for comparison)
% 
% [K,~]=dlqr(A_ref(:,:,1),B_ref(:,:,1),Q,R);
% K=-K;
% 
% [Xn,V,~]=FindXn(A_ref(:,:,1),B_ref(:,:,1),K,N,e_lb,e_ub,Ub_lb(1:nu),Ub_ub(1:nu),'lqr','linpr');
% XN=Xn{end};
% Xf1=Xn{1};
% XN=Polyhedron(XN.A,XN.b);
% Xf1=Polyhedron(Xf1.A,Xf1.b);
% XN_proj=XN.projection([1 2]);
% Xf1_proj=Xf1.projection([1 2]);
% figure
% hold on
% XN_proj.plot('color','red')
% Xf1_proj.plot('color','cyan')
% Xf_proj.plot('color','blue')
% xlabel('$e_x$','Interpreter','latex') 
% ylabel('$e_y$','Interpreter','latex')
% grid on; box on
% title('X_N projection (example for the first iteration)')
% legend('$\mathcal{X}_{N,1}$','$\mathcal{X}_{f,1}$','$\mathcal{X}_f$','Interpreter','latex')
% hold off

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