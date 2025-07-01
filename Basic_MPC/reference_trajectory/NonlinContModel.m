function xp=NonlinContModel(t,x,tu,u)

% Resample control at current time instant
u=interp1(tu,u,t);

% Dynamics
xp=zeros(3,1);

xp(1)=u(1)*cos(x(3));
xp(2)=u(1)*sin(x(3));
xp(3)=u(2);

end