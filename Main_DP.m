function [T,X]=Main_DP()
% script for running the double pendulum problem
clear

myode = @Double_Pendulum;

[T,X]=ode45(myode,[0,10],[pi/6,0,pi/3,0]);

for i=1:1:length(X)  %USE THIS TO COMPUTE Energy conservation
    
    [~ , Ox,Oy,Ax,Ay,props] = Double_Pendulum( T(i),X(i,:),'dummy');
    
    FA(i)=norm([Ax,Ay]);  %Pythagorus
    FO(i)=norm([Ox,Oy]);  %Pythagorus
end

figure(1)
plot(T,[X(:,1),X(:,2),X(:,3),X(:,4)])
legend('Theta1','d\Theta1/dt','Theta2','d\Theta2/dt')
title('Motion')
ylabel('Angle [radians]');
xlabel('Time [seconds]');

 PP1 = spline(T,X(:,1));
 PP2 = spline(T,X(:,3));

for t=0:0.01:T(end)
    figure(2)
    
    theta1=ppval(PP1,t);
    theta2=ppval(PP2,t);
    [x1,y1]=pol2cart(theta1-pi/2,1);
    [x2,y2]=pol2cart(theta2-pi/2,1);
    x2=x1+x2;
    y2=y1+y2;
    plot( [0,x1,x2],[0,y1,y2],'linewidth',5 )
    axis([-2.5,2.5,-2.5,2.5])
    text(1.0,1.0,sprintf('Time= %2.2f',t))
    pause(0.001)
end

[maxFo, p1] = max(FO);
[maxFa, p2] = max(FA);

figure(2)
plot(T,[FO;FA])
legend('Point O','Point A')
title('Reaction forces')
ylabel('Force [newtons]');
xlabel('Time [seconds]');
end

function [dydt, Ox,Oy,Ax,Ay, props] = Double_Pendulum( t,y, version )
%This function models a double pendulum (two rigid links)
%
% y=[theta1, theta1', theta2, theta2'];
% x=[Ox,Oy,Ax,Az,ax,ay,bx,by,theta1'',theta2''];

%Define the system parameters
L1=1; %m
L2=1; %m
m1=1; %kg
m2=1; %kg
g=9.81; %m/s^2
I0=m1*L1^2/3;
I1=m1*L1^2/12;
I2=m2*L2^2/12;

%Define theta
 theta1=y(1);
 theta2=y(3);
 C1M2=cos(theta1-theta2);
 S1M2=sin(theta1-theta2);
 S1=sin(theta1);
 S2=sin(theta2);
 C1=cos(theta1);
 C2=cos(theta2);
 
%Define constant
 c1 = m1*(L1/2)^2/2 + I1/2 + m2*L1^2/2;
 c2 = m1*(L2/2)^2/2 + I2/2;
 c3 = m2*L1*L2/2;
 c4 = g*(m1*L1/2+m2*L1);
 c5 = g*m2*L2/2;
 
%Define the matrix A
%[theta1'',         theta2'', Ax,          Ay,     Ox,      Oy,     ax,      ay;
A=[-I0,             0,    -L1*C1,     -L1*S1,       0,       0,     0,         0;...
     0, -1/12*m2*L2^2,    -L2/2*C2,   -L2/2*S2,     0,       0,     0,         0;...
     0,             0,        -1,          0,       1,       0,   -m1,         0;...
     0,             0,        0,          -1,       0,       1,     0,       -m1;...
L1/2*C1,            0,        0,           0,       0,       0,    -1,         0;...
L1/2*S1,            0,        0,           0,       0,       0,     0,        -1;...
   2*c1,      c3*C1M2,        0,           0,       0,       0,     0,         0;...
c3*C1M2,         2*c2,        0,           0,       0,       0,     0,         0];

%theta1'^2,    theta2'^2,  
B=[m1*L1*g*S1/2; 0; 0; m1*g;...
   L1*y(2)^2*S1/2; -L1*y(2)^2*C1/2; -c4*S1-c3*y(4)^2*S1M2; -c5*S2+c3*y(2)^2*S1M2];

%A
%B
%pause
X=A\B;

dydt=[y(2);X(1);y(4);X(2)];
%These lines are for later
if nargin==3
Ax=X(3);
Ay=X(4);
Ox=X(5);
Oy=X(6);

%create a data structure with the physical properties
props.L1=L1;
props.L2=L2;
props.m1=m1;
props.m2=m2;
props.I1=I1;
props.I2=I2;
props.g=g;
end

end
