function [M,M_inv,g] = kukalwrdynamic(q,dq)
% syms q1 q2 q3 q4 q5 q6 q7
% syms dq1 dq2 dq3 dq4 dq5 dq6 dq7
% q=[q1;q2;q3;q4;q5;q6;q7];
% dq=[dq1;dq2;dq3;dq4;dq5;dq6;dq7];
q1=q(1);q2=q(2);q3=q(3);q4=q(4);q5=q(5);q6=q(6);q7=q(7);
ri=[-0.31/2;0.1;0.3;-0.39/4;0.75*0.39;0;0.05];
rm=[0;0;0;0;0;0;0];
d3=0.4;
d5=0.39;
d0=0.31;
d7=0.078;
kr1=100;kr2=100;kr3=100;kr4=100;kr5=100;kr6=100;kr7=100;
m1=0.7;m2=0.7;m3=0.7;m4=0.7;m5=0.7;m6=0.7;m7=0.7;
M1=0.7;M2=0.7;M3=0.7;M4=0.7;M5=0.7;M6=0.7;M7=0.7;
I1=10;I2=10;I3=10;I4=10;I5=10;I6=10;I7=10;
IM1=10;IM2=10;IM3=10;IM4=10;IM5=10;IM6=10;IM7=10;
g0=[0;0;-9.8];
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
%Rotation matrix i about 0
R0=eye(3);
R1=T1(1:3,1:3);
R2=T2(1:3,1:3);
R3=T3(1:3,1:3);
R4=T4(1:3,1:3);
R5=T5(1:3,1:3);
R6=T6(1:3,1:3);
Re=Te(1:3,1:3);
%position of centeral mass link i about frame i-1
rc1=[0;0;ri(1)];
rc2=[0;ri(2);0];
rc3=[0;0;ri(3)];
rc4=[0;ri(4);0];
rc5=[0;0;ri(5)];
rc6=[0;ri(6);0];
rc7=[0;0;ri(7)];
%position of motor on the link i about frame i-1
rm1=[0;0;rm(1)];
rm2=[0;rm(2);0];
rm3=[0;0;rm(3)];
rm4=[0;rm(4);0];
rm5=[0;0;rm(5)];
rm6=[0;rm(6);0];
rm7=[0;0;rm(7)];
%position of central mass of link i about 0
Rc1=R0*rc1;
Rc2=R1*rc2;
Rc3=R2*rc3;
Rc4=R3*rc4;
Rc5=R4*rc5;
Rc6=R5*rc6;
Rc7=R6*rc7;
%position of motor of link i about 0
Rm1=R0*rm1;
Rm2=R1*rm2;
Rm3=R2*rm3;
Rm4=R3*rm4;
Rm5=R4*rm5;
Rm6=R5*rm6;
Rm7=R6*rm7;
%jacobian matrix
z0=[0;0;1];
z1=T1(1:3,3);
z2=T2(1:3,3);
z3=T3(1:3,3);
z4=T4(1:3,3);
z5=T5(1:3,3);
z6=T6(1:3,3);
z7=Te(1:3,3);
p0=[0;0;0];
p1=T1(1:3,4);
p2=T2(1:3,4);
p3=T3(1:3,4);
p4=T4(1:3,4);
p5=T5(1:3,4);
p6=T6(1:3,4);
p7=Te(1:3,4);
zero=[0;0;0];
% j7=[cross(z0,(p7-p0)) cross(z1,(p7-p1)) cross(z2,(p7-p2)) cross(z3,(p7-p3)) cross(z4,(p7-p4)) cross(z5,(p7-p5)) cross(z6,(p7-p6));z0 z1 z2 z3 z4 z5 z6]; %jacobian matrix of end_effector
% j6=[cross(z0,(p6-p0)) cross(z1,(p6-p1)) cross(z2,(p6-p2)) cross(z3,(p6-p3)) cross(z4,(p6-p4)) cross(z5,(p6-p5)) zero;z0 z1 z2 z3 z4 z5 zero]; 
% j5=[cross(z0,(p5-p0)) cross(z1,(p5-p1)) cross(z2,(p5-p2)) cross(z3,(p5-p3)) cross(z4,(p5-p4)) zero zero;z0 z1 z2 z3 z4 zero zero]; 
% j4=[cross(z0,(p4-p0)) cross(z1,(p4-p1)) cross(z2,(p4-p2)) cross(z3,(p4-p3)) zero zero zero;z0 z1 z2 z3 zero zero zero]; 
% j3=[cross(z0,(p3-p0)) cross(z1,(p3-p1)) cross(z2,(p3-p2)) zero zero zero zero ;z0 z1 z2 zero zero zero zero]; 
% j2=[cross(z0,(p2-p0)) cross(z1,(p2-p1)) zero zero zero zero zero;z0 z1 zero zero zero zero zero]; 
% j1=[cross(z0,(p1-p0)) zero zero zero zero zero zero;z0 zero zero zero zero zero zero]; 
% j0=[zero zero zero zero zero zero zero;[0;0;1] zero zero zero zero zero zero]; 
%link's center of mass's position
pc1=p0+Rc1;
pc2=p1+Rc2;
pc3=p2+Rc3;
pc4=p3+Rc4;
pc5=p4+Rc5;
pc6=p5+Rc6;
pc7=p6+Rc7;
%link's motor's position
pm1=p0+Rm1;
pm2=p1+Rm2;
pm3=p2+Rm3;
pm4=p3+Rm4;
pm5=p4+Rm5;
pm6=p5+Rm6;
pm7=p6+Rm7;
% jacobian matrix of central mass of link i
jc1=[cross(z0,(pc1-p0)) zero zero zero zero zero zero;z0 zero zero zero zero zero zero];
jc2=[cross(z0,(pc2-p0)) cross(z1,(pc2-p1)) zero zero zero zero zero;z0 z1 zero zero zero zero zero];
jc3=[cross(z0,(pc3-p0)) cross(z1,(pc3-p1)) cross(z2,(pc3-p2)) zero zero zero zero;z0 z1 z2 zero zero zero zero];
jc4=[cross(z0,(pc4-p0)) cross(z1,(pc4-p1)) cross(z2,(pc4-p2)) cross(z3,(pc4-p3)) zero zero zero;z0 z1 z2 z3 zero zero zero];
jc5=[cross(z0,(pc5-p0)) cross(z1,(pc5-p1)) cross(z2,(pc5-p2)) cross(z3,(pc5-p3)) cross(z4,(pc5-p4)) zero zero;z0 z1 z2 z3 z4 zero zero];
jc6=[cross(z0,(pc6-p0)) cross(z1,(pc6-p1)) cross(z2,(pc6-p2)) cross(z3,(pc6-p3)) cross(z4,(pc6-p4)) cross(z5,(pc6-p5)) zero;z0 z1 z2 z3 z4 z5 zero];
jc7=[cross(z0,(pc7-p0)) cross(z1,(pc7-p1)) cross(z2,(pc7-p2)) cross(z3,(pc7-p3)) cross(z4,(pc7-p4)) cross(z5,(pc7-p5)) cross(z6,(pc7-p6)) ;z0 z1 z2 z3 z4 z5 z6];
%motor's jacobian matrix on center of link
jm1=[cross(z0,(pm1-p0)) zero zero zero zero zero zero;kr1*z0 zero zero zero zero zero zero];
jm2=[cross(z0,(pm2-p0)) cross(z1,(pm2-p1)) zero zero zero zero zero;z0 kr2*z1 zero zero zero zero zero];
jm3=[cross(z0,(pm3-p0)) cross(z1,(pm3-p1)) cross(z2,(pm3-p2)) zero zero zero zero;z0 z1 kr3*z2 zero zero zero zero];
jm4=[cross(z0,(pm4-p0)) cross(z1,(pm4-p1)) cross(z2,(pm4-p2)) cross(z3,(pm4-p3)) zero zero zero;z0 z1 z2 kr4*z3 zero zero zero];
jm5=[cross(z0,(pm5-p0)) cross(z1,(pm5-p1)) cross(z2,(pm5-p2)) cross(z3,(pm5-p3)) cross(z4,(pm5-p4)) zero zero;z0 z1 z2 z3 kr5*z4 zero zero];
jm6=[cross(z0,(pm6-p0)) cross(z1,(pm6-p1)) cross(z2,(pm6-p2)) cross(z3,(pm6-p3)) cross(z4,(pm6-p4)) cross(z5,(pm6-p5)) zero;z0 z1 z2 z3 z4 kr6*z5 zero];
jm7=[cross(z0,(pm7-p0)) cross(z1,(pm7-p1)) cross(z2,(pm7-p2)) cross(z3,(pm7-p3)) cross(z4,(pm7-p4)) cross(z5,(pm7-p5)) cross(z6,(pm7-p6)) ;z0 z1 z2 z3 z4 z5 z6*kr7];
% gravity matrix
g1=-((m1*g0.'*jc1(1:3,1))+(M1*g0.'*jm1(1:3,1))+(m2*g0.'*jc2(1:3,1))+(M2*g0.'*jm2(1:3,1))+(m3*g0.'*jc3(1:3,1))+(M3*g0.'*jm3(1:3,1))+(m4*g0.'*jc4(1:3,1))+(M4*g0.'*jm4(1:3,1))+(m5*g0.'*jc5(1:3,1))+(M5*g0.'*jm5(1:3,1))+(m6*g0.'*jc6(1:3,1))+(M6*g0.'*jm6(1:3,1))+(m7*g0.'*jc7(1:3,1))+(M7*g0.'*jm7(1:3,1)));
g2=-((m1*g0.'*jc1(1:3,2))+(M1*g0.'*jm1(1:3,2))+(m2*g0.'*jc2(1:3,2))+(M2*g0.'*jm2(1:3,2))+(m3*g0.'*jc3(1:3,2))+(M3*g0.'*jm3(1:3,2))+(m4*g0.'*jc4(1:3,2))+(M4*g0.'*jm4(1:3,2))+(m5*g0.'*jc5(1:3,2))+(M5*g0.'*jm5(1:3,2))+(m6*g0.'*jc6(1:3,2))+(M6*g0.'*jm6(1:3,2))+(m7*g0.'*jc7(1:3,2))+(M7*g0.'*jm7(1:3,2)));
g3=-((m1*g0.'*jc1(1:3,3))+(M1*g0.'*jm1(1:3,3))+(m2*g0.'*jc2(1:3,3))+(M2*g0.'*jm2(1:3,3))+(m3*g0.'*jc3(1:3,3))+(M3*g0.'*jm3(1:3,3))+(m4*g0.'*jc4(1:3,3))+(M4*g0.'*jm4(1:3,3))+(m5*g0.'*jc5(1:3,3))+(M5*g0.'*jm5(1:3,3))+(m6*g0.'*jc6(1:3,3))+(M6*g0.'*jm6(1:3,3))+(m7*g0.'*jc7(1:3,3))+(M7*g0.'*jm7(1:3,3)));
g4=-((m1*g0.'*jc1(1:3,4))+(M1*g0.'*jm1(1:3,4))+(m2*g0.'*jc2(1:3,4))+(M2*g0.'*jm2(1:3,4))+(m3*g0.'*jc3(1:3,4))+(M3*g0.'*jm3(1:3,4))+(m4*g0.'*jc4(1:3,4))+(M4*g0.'*jm4(1:3,4))+(m5*g0.'*jc5(1:3,4))+(M5*g0.'*jm5(1:3,4))+(m6*g0.'*jc6(1:3,4))+(M6*g0.'*jm6(1:3,4))+(m7*g0.'*jc7(1:3,4))+(M7*g0.'*jm7(1:3,4)));
g5=-((m1*g0.'*jc1(1:3,5))+(M1*g0.'*jm1(1:3,5))+(m2*g0.'*jc2(1:3,5))+(M2*g0.'*jm2(1:3,5))+(m3*g0.'*jc3(1:3,5))+(M3*g0.'*jm3(1:3,5))+(m4*g0.'*jc4(1:3,5))+(M4*g0.'*jm4(1:3,5))+(m5*g0.'*jc5(1:3,5))+(M5*g0.'*jm5(1:3,5))+(m6*g0.'*jc6(1:3,5))+(M6*g0.'*jm6(1:3,5))+(m7*g0.'*jc7(1:3,5))+(M7*g0.'*jm7(1:3,5)));
g6=-((m1*g0.'*jc1(1:3,6))+(M1*g0.'*jm1(1:3,6))+(m2*g0.'*jc2(1:3,6))+(M2*g0.'*jm2(1:3,6))+(m3*g0.'*jc3(1:3,6))+(M3*g0.'*jm3(1:3,6))+(m4*g0.'*jc4(1:3,6))+(M4*g0.'*jm4(1:3,6))+(m5*g0.'*jc5(1:3,6))+(M5*g0.'*jm5(1:3,6))+(m6*g0.'*jc6(1:3,6))+(M6*g0.'*jm6(1:3,6))+(m7*g0.'*jc7(1:3,6))+(M7*g0.'*jm7(1:3,6)));
g7=-((m1*g0.'*jc1(1:3,7))+(M1*g0.'*jm1(1:3,7))+(m2*g0.'*jc2(1:3,7))+(M2*g0.'*jm2(1:3,7))+(m3*g0.'*jc3(1:3,7))+(M3*g0.'*jm3(1:3,7))+(m4*g0.'*jc4(1:3,7))+(M4*g0.'*jm4(1:3,7))+(m5*g0.'*jc5(1:3,7))+(M5*g0.'*jm5(1:3,7))+(m6*g0.'*jc6(1:3,7))+(M6*g0.'*jm6(1:3,7))+(m7*g0.'*jc7(1:3,7))+(M7*g0.'*jm7(1:3,7)));
g=[g1;g2;g3;g4;g5;g6;g7];
%Mass matrix
b1=m1*jc1(1:3,:).'*jc1(1:3,:)+jc1(4:6,:).'*R1*I1*R1.'*jc1(4:6,:)+M1*jm1(1:3,:).'*jm1(1:3,:)+jm1(4:6,:).'*(R1)*IM1*R1.'*jm1(4:6,:);
b2=m2*jc2(1:3,:).'*jc2(1:3,:)+jc2(4:6,:).'*R2*I2*R2.'*jc2(4:6,:)+M2*jm2(1:3,:).'*jm2(1:3,:)+jm2(4:6,:).'*(R2)*IM2*R2.'*jm2(4:6,:);
b3=m3*jc3(1:3,:).'*jc3(1:3,:)+jc3(4:6,:).'*R3*I3*R3.'*jc3(4:6,:)+M3*jm3(1:3,:).'*jm3(1:3,:)+jm3(4:6,:).'*(R3)*IM3*R3.'*jm3(4:6,:);
b4=m4*jc4(1:3,:).'*jc4(1:3,:)+jc4(4:6,:).'*R4*I4*R4.'*jc4(4:6,:)+M4*jm4(1:3,:).'*jm4(1:3,:)+jm4(4:6,:).'*(R4)*IM4*R4.'*jm4(4:6,:);
b5=m5*jc5(1:3,:).'*jc5(1:3,:)+jc5(4:6,:).'*R5*I5*R5.'*jc5(4:6,:)+M5*jm5(1:3,:).'*jm5(1:3,:)+jm5(4:6,:).'*(R5)*IM5*R5.'*jm5(4:6,:);
b6=m6*jc6(1:3,:).'*jc6(1:3,:)+jc6(4:6,:).'*R6*I6*R6.'*jc6(4:6,:)+M6*jm6(1:3,:).'*jm6(1:3,:)+jm6(4:6,:).'*(R6)*IM6*R6.'*jm6(4:6,:);
b7=m7*jc7(1:3,:).'*jc7(1:3,:)+jc7(4:6,:).'*Re*I7*Re.'*jc7(4:6,:)+M7*jm7(1:3,:).'*jm7(1:3,:)+jm7(4:6,:).'*(Re)*IM7*Re.'*jm7(4:6,:);
M=b1+b2+b3+b4+b5+b6+b7;
M_inv=inv(M);
%C matrix
% c=zeros(7,7,7);
% C=zeros(7,7);
% for k=1:7
%     for i=1:7
%         for j=1:7
%             a=diff(M(i,j),q(k));
%             b=+diff(M(i,k),q(j));
%             d=diff(M(j,k),q(i));
%             h=0.5*simplify(a+b-d)
% %             c(i,j,k)=vpa(h,5);
%         end 
%     end
% end
% for i=1:7
%     for j=1:7
%         for k=1:7
%             C(i,j)=c(i,j,k)*dq(k);
%         end
%     end
% end
% % cq=C*dq;
