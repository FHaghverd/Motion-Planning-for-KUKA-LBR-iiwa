function [ je,Te ] = jacobin( q1,q2,q3,q4,q5,q6,q7 )
d3=0.4;
d5=0.39;
d0=0.31;
d7=0.078;
f=[pi/2;-pi/2;-pi/2;pi/2;pi/2;-pi/2;0];
a=[0;0;0;0;0;0;0];
d=[0;0;d3;0;d5;0;0];
t=[q1;q2;q3;q4;q5;q6;q7];
%Transformation matrix i about i-1
A0=[1 0 0 0;0 1 0 0;0 0 1 d0;0 0 0 1];  
A1=[cos(t(1,1)) -sin(t(1,1))*cos(f(1,1)) sin(t(1,1))*sin(f(1,1)) a(1,1)*cos(t(1,1));sin(t(1,1)) cos(t(1,1))*cos(f(1,1)) -cos(t(1 ,1))*sin(f(1,1)) a(1,1)*sin(t(1,1));0 sin(f(1,1)) cos(f(1,1)) d(1,1);0 0 0 1];
A2=[cos(t(2,1)) -sin(t(2,1))*cos(f(2,1)) sin(t(2,1))*sin(f(2,1)) a(2,1)*cos(t(2,1));sin(t(2,1)) cos(t(2,1))*cos(f(2,1)) -cos(t(2 ,1))*sin(f(2,1)) a(2,1)*sin(t(2,1));0 sin(f(2,1)) cos(f(2,1)) d(2,1);0 0 0 1];
A3=[cos(t(3,1)) -sin(t(3,1))*cos(f(3,1)) sin(t(3,1))*sin(f(3,1)) a(3,1)*cos(t(3,1));sin(t(3,1)) cos(t(3,1))*cos(f(3,1)) -cos(t(3 ,1))*sin(f(3,1)) a(3,1)*sin(t(3,1));0 sin(f(3,1)) cos(f(3,1)) d(3,1);0 0 0 1];
A4=[cos(t(4,1)) -sin(t(4,1))*cos(f(4,1)) sin(t(4,1))*sin(f(4,1)) a(4,1)*cos(t(4,1));sin(t(4,1)) cos(t(4,1))*cos(f(4,1)) -cos(t(4 ,1))*sin(f(4,1)) a(4,1)*sin(t(4,1));0 sin(f(4,1)) cos(f(4,1)) d(4,1);0 0 0 1];
A5=[cos(t(5,1)) -sin(t(5,1))*cos(f(5,1)) sin(t(5,1))*sin(f(5,1)) a(5,1)*cos(t(5,1));sin(t(5,1)) cos(t(5,1))*cos(f(5,1)) -cos(t(5 ,1))*sin(f(5,1)) a(5,1)*sin(t(5,1));0 sin(f(5,1)) cos(f(5,1)) d(5,1);0 0 0 1];
A6=[cos(t(6,1)) -sin(t(6,1))*cos(f(6,1)) sin(t(6,1))*sin(f(6,1)) a(6,1)*cos(t(6,1));sin(t(6,1)) cos(t(6,1))*cos(f(6,1)) -cos(t(6 ,1))*sin(f(6,1)) a(6,1)*sin(t(6,1));0 sin(f(6,1)) cos(f(6,1)) d(6,1);0 0 0 1];
A7=[cos(t(7,1)) -sin(t(7,1))*cos(f(7,1)) sin(t(7,1))*sin(f(7,1)) a(7,1)*cos(t(7,1));sin(t(7,1)) cos(t(7,1))*cos(f(7,1)) -cos(t(7 ,1))*sin(f(7,1)) a(7,1)*sin(t(7,1));0 sin(f(7,1)) cos(f(7,1)) d(7,1);0 0 0 1];
Ae=[1 0 0 0;0 1 0 0;0 0 1 d7;0 0 0 1];
%Transformation matrix i about 0
T1=A0*A1;
T2=A0*A1*A2;
T3=A0*A1*A2*A3;
T4=A0*A1*A2*A3*A4;
T5=A0*A1*A2*A3*A4*A5;
T6=A0*A1*A2*A3*A4*A5*A6;
Te=A0*A1*A2*A3*A4*A5*A6*A7*Ae;
% %Rotation matrix i about 0
% R0=eye(3);
% R1=T1(1:3,1:3);
% R2=T2(1:3,1:3);
% R3=T3(1:3,1:3);
% R4=T4(1:3,1:3);
% R5=T5(1:3,1:3);
% R6=T6(1:3,1:3);
% R7=Te(1:3,1:3);
%jacobian matrix
z0=[0;0;1];
z1=T1(1:3,3);
z2=T2(1:3,3);
z3=T3(1:3,3);
z4=T4(1:3,3);
z5=T5(1:3,3);
z6=T6(1:3,3);
% z7=Te(1:3,3);
p0=[0;0;0];
p1=T1(1:3,4);
p2=T2(1:3,4);
p3=T3(1:3,4);
p4=T4(1:3,4);
p5=T5(1:3,4);
p6=T6(1:3,4);
pe=Te(1:3,4);
je=[cross(z0,(pe-p0)) cross(z1,(pe-p1)) cross(z2,(pe-p2)) cross(z3,(pe-p3)) cross(z4,(pe-p4)) cross(z5,(pe-p5)) cross(z6,(pe-p6));z0 z1 z2 z3 z4 z5 z6]; %jacobian matrix of end_effector
end

