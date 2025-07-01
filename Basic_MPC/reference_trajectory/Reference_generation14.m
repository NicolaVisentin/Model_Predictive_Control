clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Generates the reference trajectory and reference inputs for a unicycle  %
% robot through simulation of a given control input.                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulate a control

% Choose trajectory to save (1 or 4)

choose=4; 

% Parameters

if choose==1
    tf=10;           % final time (duration)
    x0=[0 0 0]';     % initial state
else
    tf=40;           % final time (duration)
    x0=[0 0 pi/4]';  % initial state
end

% Control definition

if choose==1
    v=@(t) 0.4+0*t;
    w=@(t) sin(t-pi/2);
else
    v=@(t) 0.4+0.1*sin(t);
    w=@(t) 0.1*sin(t);
end

u={v;w};

% Actual simulation

dt=1/30;
options=odeset(RelTol=1e-6);
[t,x]=ode45(@(t,x) NonlinContModel14(t,x,u),0:dt:tf,x0,options);

%% Save data (reference trajectory and control)

% Sample control

u=[v(t) w(t)];
u(end,:)=[];

% Save in a struct

ref.x_ref=x;
ref.u_ref=u;
ref.t_ref=t;

if choose==1
    save('ref_data1.mat','ref')
else
    save('ref_data4.mat','ref')
end

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
plot(t(1:end-1),u(:,1),t(1:end-1),u(:,2))
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

%% Animations

if true

pos=x(:,1:2);
theta=x(:,3);

pause(1)
figure
hold on
tr=plot(nan,nan,'k--');
h1=plot(pos(1,1),pos(1,2),'ro');
h2=plot([pos(1,1) pos(1,1)+cos(theta(1))],[pos(1,2) pos(1,2)+cos(theta(1))],'r');
grid on; box on
xlabel('X [m]')
ylabel('Y [m]')
dynamic_title=title('Animation');
axis([1.2*min(pos(:,1)) 1.2*max(pos(:,1)) 1.2*min(pos(:,2)) 1.2*max(pos(:,2))])
axis equal

jj=2;
for ii=1:10:length(t)

    front_x=pos(ii,1)+0.7*cos(theta(ii));
    front_y=pos(ii,2)+0.7*sin(theta(ii));

    set(tr,'XData',pos(1:ii,1),'YData',pos(1:ii,2))
    set(h1,'XData',pos(ii,1),'YData',pos(ii,2))
    set(h2,'Xdata',[pos(ii,1) front_x],'YData',[pos(ii,2) front_y])
    set(dynamic_title,'String',sprintf('Animation (t = %.1f s)',t(ii)));

    drawnow
    pause(0.1)

end

end