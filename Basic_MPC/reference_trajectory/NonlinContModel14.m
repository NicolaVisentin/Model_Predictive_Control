function xp=NonlinContModel14(t,x,u)

% Extract control
v=u{1};
v=v(t);
w=u{2};
w=w(t);
u=[v w]';

% Dynamics
xp=zeros(3,1);

xp(1)=u(1)*cos(x(3));
xp(2)=u(1)*sin(x(3));
xp(3)=u(2);

end