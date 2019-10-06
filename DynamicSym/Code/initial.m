clear all
clc
global ri rm m1 m2 m3 m4 m5 m6 m7 M1 M2 M3 M4 M5 M6 M7 ;
global I1 I2 I3 I4 I5 I6 I7 Im1 Im2 Im3 Im4 Im5 Im6 Im7 ;
global kr1 kr2 kr3 kr4 kr5 kr6 kr7 g0 ;
ri=[-0.31/2;0.1;0.3;-0.39/4;0.75*0.39;0;0.05]; %distance of the links center of mass about joint
rm=[0;0;0;0;0;0;0]; %distance of motor of each link about joint
m1=0.7;m2=0.7;m3=0.7;m4=0.7;m5=0.7;m6=0.7;m7=0.7; %link's mass
M1=0.7;M2=0.7;M3=0.7;M4=0.7;M5=0.7;M6=0.7;M7=0.7; %motor's mass
I1=10;I2=10;I3=10;I4=10;I5=10;I6=10;I7=10; %link's inersi
Im1=10;Im2=10;Im3=10;Im4=10;Im5=10;Im6=10;Im7=10; %motor's inersi
kr1=100;kr2=100;kr3=100;kr4=100;kr5=100;kr6=100;kr7=100; 
g0=[0;0;-9.8];
qi=[pi/2;0;0;pi/2;0;pi/2;pi/6];
dqi=[0;0;0;0;0;0;0];
sim simlwrkuka_dynamic_dirdyn