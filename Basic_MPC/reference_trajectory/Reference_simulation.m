clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Simulates a certain open-loop reference trajectory (solving nonlinear   %
% dynamics with a certain initial condition and control)                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation

% Load and extract time vector, control and initial condition

choose=2;                       % choose reference to simulate
switch choose
    case 2
        load('ref_data2.mat')   % analytical reference ('Reference_generation23.m') --> trajectory n3 in the report        
    case 3
        load('ref_data3.mat')   % analytical reference ('Reference_generation23.m') --> useless, basically
    case 4
        load('ref_data4.mat')   % simple reference through simulations ('Reference_generation14.m') --> trajectory n1 in the report
    otherwise
        load('ref_data1.mat')   % simple reference through simulations ('Reference_generation14.m') --> trajectory n2 in the report
end

t_ref=ref.t_ref;
x_ref=ref.x_ref;
u_ref=ref.u_ref;

tf=t_ref(end);    % final time (duration)
x0=x_ref(1,:)';   % initial state

% Actual simulation

[t,x]=ode45(@(t,x) NonlinContModel(t,x,t_ref(1:end-1),u_ref),[0 tf],x0);

%% Plots

figure

subplot(211)
hold on
plot(t,x(:,1),'b')
plot(t,x(:,2),'r')
plot(t,x(:,3),'g')
grid on; box on
xlabel('t [s]','Interpreter','latex')
ylabel('states','Interpreter','latex')
title('States')
legend('$x$','$y$','$\theta$','Interpreter','latex','Location','best')

subplot(212)
plot(t_ref(1:end-1),u_ref(:,1),t_ref(1:end-1),u_ref(:,2))
grid on; box on
xlabel('t [s]','Interpreter','latex')
ylabel('controls','Interpreter','latex')
title('Controls')
legend('v','$\omega$','Interpreter','latex','Location','best')

figure
plot(x(:,1),x(:,2))
grid on; box on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('Trajectory')

%% Animation

l=0.3;   % length of the robot (adjust for visualisation)

if true

pos=x(:,1:2);
theta=x(:,3);

pause(1)
figure
hold on
tr=plot(nan,nan,'k--');
h1=plot(pos(1,1),pos(1,2),'ro');
h2=plot([pos(1,1) pos(1,1)+l*cos(theta(1))],[pos(1,2) pos(1,2)+l*cos(theta(1))],'r');
grid on; box on
xlabel('X [m]')
ylabel('Y [m]')
dynamic_title=title('Animation');
axis([1.2*min(pos(:,1)) 1.2*max(pos(:,1)) 1.2*min(pos(:,2)) 1.2*max(pos(:,2))])
axis equal

jj=2;
for ii=1:length(t)

    front_x=pos(ii,1)+l*cos(theta(ii));
    front_y=pos(ii,2)+l*sin(theta(ii));

    set(tr,'XData',pos(1:ii,1),'YData',pos(1:ii,2))
    set(h1,'XData',pos(ii,1),'YData',pos(ii,2))
    set(h2,'Xdata',[pos(ii,1) front_x],'YData',[pos(ii,2) front_y])
    set(dynamic_title,'String',sprintf('Animation (t = %.1f s)',t(ii)));

    drawnow
    pause(0.1)

end

end