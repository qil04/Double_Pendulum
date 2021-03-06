function [T,X]=Main_DP()
% script for running the double pendulum problem
clear

myode = @Double_Pendulum;

[T,X]=ode45(myode,[0,10],[pi/2,0,0,0]);

for i=1:1:length(X)  %USE THIS TO COMPUTE Energy conservation
    
    [ ~ , Ox,Oy,Ax,Ay,props ] = Double_Pendulum( T(i),X(i,:),'dummy');
    
    FO(i)=norm([Ox,Oy]);  %Pythagorus
    FA(i)=norm([Ax,Ay]);  %Pythagorus
    
    %Computation of Kinetic Energy
    
    [KE(i),PE(i),List]=Kinetic_Energy(props,X(i,:));
   
    
      KE1L(i)=List.KE1L;
      KE2L(i)=List.KE2L;
      KE1R(i)=List.KE1R;
      KE2R(i)=List.KE2R;
        P1(i)=List.P1;
        P2(i)=List.P2;
end

% figure(1)
% plot(T,[X(:,1),X(:,3)])
% legend('Theta1','Theta2')
% title('Motion')
% ylabel('Angle [radians]');
% xlabel('Time [seconds]');
% 
% 
%  PP1 = spline(T,X(:,1));
%  PP2 = spline(T,X(:,3));
% 
% for t=0:0.01:T(end)
%     figure(2)
%     
%     theta1=ppval(PP1,t);
%     theta2=ppval(PP2,t);
%     [x1,y1]=pol2cart(theta1-pi/2,props.L1);
%     [x2,y2]=pol2cart(theta2-pi/2,props.L2);
%     x2=x1+x2;
%     y2=y1+y2;
%     plot( [0,x1,x2],[0,y1,y2],'linewidth',5 )
%     axis([-2.5,2.5,-2.5,2.5])
%     text(1.5,1.5,sprintf('Time= %2.2f',t))
%     pause(0.001)
%      
% end
 

figure(2)
plot(T,[FO;FA])
legend('Point O','Point A')
title('Reaction forces')
ylabel('Force [newtons]');
xlabel('Time [seconds]');

figure(3)
plot(T,[KE;PE;KE+PE]);
legend('Kinetic','Potential','Total')
title('Conservation of Energy')
ylabel('Energy [Joules]');
xlabel('Time [seconds]');

figure(4)
plot(T,[KE1L;KE1R;KE2L;KE2R;P1;P2]);
legend('KE1L','KE1R','KE2L','KE2R','P1','P2')
title('Conservation of Energy')
ylabel('Energy [Joules]');
xlabel('Time [seconds]');

      KE1L(i)=List.KE1L;
      KE2L(i)=List.KE2L;
      KE1R(i)=List.KE1R;
      KE2R(i)=List.KE2R;
        P1(i)=List.P1;
        P2(i)=List.P2;
        
        
        end
        
        
function [ dydt, Ox,Oy,Ax,Ay, props ] = Double_Pendulum( t,y, version )
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
I0= 1/3*m1*L1^2;   I1=1/12*m1*L1^2;
I2=1/12*m2*L2^2;

%Define theta
 theta1=y(1);
 theta2=y(3);

 S1=sin(theta1);
 C1=cos(theta1);
 S2=sin(theta2);
 C2=cos(theta2);

%Define the matrix A
%[Ox,Oy,       Ax,       Ay,  ax,  ay,  bx,  by, theta1'',      theta2'';
A=[1, 0,       -1,        0, -m1,   0,   0,   0,        0,             0;...
   0, 1,        0,       -1,   0, -m1,   0,   0,        0,             0;...
   0, 0,   -L1*C1,   -L1*S1,   0,   0,   0,   0,      -I0,             0;...
   0, 0,        0,        0,  -1,   0,   0,   0,   L1/2*C1,             0;...
   0, 0,        0,        0,   0,  -1,   0,   0,   L1/2*S1,             0;...
   0, 0,        1,        0,   0,   0, -m2,   0,        0,             0;...
   0, 0,        0,        1,   0,   0,   0, -m2,        0,             0;...
   0, 0, -L2/2*C2, -L2/2*S2,   0,   0,   0,   0,        0, -1/12*m2*L2^2;...
   0, 0,        0,        0,   2,   0,  -1,   0,        0, L2/2*C2;...
   0, 0,        0,        0,   0,   2,   0,  -1,        0, L2/2*S2];
   
B=[0; m1*g; m1*L1*g*S1/2; L1*y(2)^2*S1/2; -L1*y(2)^2*C1/2; ...
   0; m2*g;            0; L2*y(4)^2*S2/2; -L2*y(4)^2*C2/2];

%A
%B
%pause
X=A\B;

dydt=[y(2);X(9);y(4);X(10)];

%These lines are for later
if nargin==3
Ox=X(1);
Oy=X(2);
Ax=X(3);
Ay=X(4);

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

