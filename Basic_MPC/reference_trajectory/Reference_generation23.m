clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Generates the reference trajectory and reference inputs for a unicycle  %
% robot through analytical derivation.                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data and parameters

% Choose reference trajectory (2 or 3)

choose=2;

% Define desired trajectory

if choose==2
    % Traj 2
    a=0.15;
    b=1.5;
    traj=@(t) 3*[cos(a*t); sin(a/b*t)];                    % (x_ref(t),y_ref(t))

    trajp=@(t) 3*[-a*sin(a*t); a/b*cos(a/b*t)];            % (xp_ref(t),yp_ref(t))
    trajpp=@(t) 3*[-a^2*cos(a*t); -(a/b)^2*sin(a/b*t)];    % (xpp_ref(t),ypp_ref(t))
else
    % Traj 3
    traj=@(t) 3*[cos(0.1*t); sin(0.2*t)];                  % (x_ref(t),y_ref(t))

    trajp=@(t) 3*[-0.1*sin(0.1*t); 0.2*cos(0.2*t)];        % (xp_ref(t),yp_ref(t))
    trajpp=@(t) 3*[-0.1^2*cos(0.1*t); -0.2^2*sin(0.2*t)];  % (xpp_ref(t),ypp_ref(t))
end

% Define time parameters

if choose==2
    tf=123;     % final time (time to percur the trajectory)
    dt=1/30;    % discretisation time resolution (for traj 2)
else
    tf=10;      % final time (time to percur the trajectory)
    dt=1/30;    % discretisation time resolution (for traj 3)
end

tvect=0:dt:tf;          % vector time
N_ref=length(tvect)-1;  % referemce timesteps

%% Compute reference trajectory and corresponding theta and inputs

x_ref=zeros(N_ref+1,3);
u_ref=zeros(N_ref,2);
trajp_ref=zeros(N_ref+1,2);
trajpp_ref=zeros(N_ref+1,2);
for k=1:N_ref+1

    % Current time 
    t=tvect(k);

    % Compute everything
    x_ref(k,1:2)=traj(t)';
    trajp_ref(k,:)=trajp(t)';
    trajpp_ref(k,:)=trajpp(t)';

    theta_ref=atan2(trajp_ref(k,2),trajp_ref(k,1));
    x_ref(k,3)=theta_ref;

    if k<=N_ref
        u_ref(k,1)=sqrt(trajp_ref(k,1)^2+trajp_ref(k,2)^2);
        u_ref(k,2)=(trajp_ref(k,1)*trajpp_ref(k,2)-trajp_ref(k,2)*trajpp_ref(k,1))/(trajp_ref(k,1)^2+trajp_ref(k,2)^2);
    end

end

%% Save data

ref.x_ref=x_ref;
ref.u_ref=u_ref;
ref.trajp_ref=trajp_ref;
ref.trajpp_ref=trajpp_ref;
ref.t_ref=tvect;
if choose==2
    save('ref_data2.mat','ref')
else
    save('ref_data3.mat','ref')
end

%% Plots to visualise

% 2D trajectory

figure
plot(x_ref(:,1),x_ref(:,2))
grid on; box on
xlabel('x [m]')
ylabel('y [m]')
title('Reference trajectory')

% State evolution

figure

subplot(311)
plot(tvect,x_ref(:,1))
grid on; box on
xlabel('t [s]')
ylabel('x [m]')
legend('x(t)')
title('Reference states')

subplot(312)
plot(tvect,x_ref(:,2))
grid on; box on
xlabel('t [s]')
ylabel('y [m]')
legend('y(t)')

subplot(313)
plot(tvect,x_ref(:,3))
grid on; box on
xlabel('t [s]')
ylabel('\theta [rad]')
legend('\theta(t)')

% Input evolution

figure

subplot(211)
plot(tvect(1:end-1),u_ref(:,1))
grid on; box on
xlabel('t [s]')
ylabel('v [m/s]')
legend('v(t)')
title('Reference controls')

subplot(212)
plot(tvect(1:end-1),u_ref(:,2))
grid on; box on
xlabel('t [s]')
ylabel('\omega [rad/s]')
legend('\omega(t)')