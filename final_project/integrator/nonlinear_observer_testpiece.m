clc;
clear all
M = 1000;
m1 = 100;
m2 = 100;
l1= 10;
l2 = 20;
g = 9.81;

A = [ 0 1 0 0 0 0 0; 0 0 -(m1*g)/M 0 -m2*g/M 0 0; 0 0 0 1 0 0 0; 0 0 -(M+m1)*g/(M*l1) 0 -m2*g/(M*l1) 0 0; 0 0 0 0 0 1 0; 0 0 -m1*g/(M*l2) 0 -(M+m2)*g/(M*l2) 0 0; -1 0 0 0 0 0 0];
B = [0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2); 0]

Kalman_C = [B A*B (A^2)*B (A^3)*B (A^4)*B (A^5)*B (A^6)*B];
Rank_K_C = rank(Kalman_C)
 
%x(t),(?1(t),?2(t)),(x(t),?2(t))or(x(t),?1(t),?2(t)).
C1= [1 0 0 0 0 0 0];
Kalman_O1 = [C1; C1*A; C1*(A^2); C1*(A^3); C1*(A^4); C1*(A^5); C1*(A^6)]; Rank_K_O1 = rank(Kalman_O1)
D1 = [0]; X01 = [0.3 0 0 0 0 0 0];

C3 = [1 0 0 0 0 0 0;0 0 1 0 0 0 0];
Kalman_O3 = [C3; C3*A; C3*(A^2); C3*(A^3); C3*(A^4); C3*(A^5); C3*(A^6)]; Rank_K_O3 = rank(Kalman_O3)
D3 =[0;0]; X03 = [0.1 0 0 0 0.0697 0 0]; 

C4=[1 0 0 0 0 0 0;0 0 1 0 0 0 0;0 0 0 0 1 0 0];
Kalman_O4 = [C4; C4*A; C4*(A^2); C4*(A^3); C4*(A^4); C4*(A^5); C4*(A^6)]; Rank_K_O4 = rank(Kalman_O4)
D4 =[0;0;0]; X04 = [0 0 0.015 0 0.0697 0 0];  

%Design of LQR
q1 =100*[1 0 0 0 0 0 0];
q2 =100*[0 1 0 0 0 0 0];
q3 =625*[0 0 1 0 0 0 0];
q4 =400*[0 0 0 1 0 0 0];
q5 =625*[0 0 0 0 1 0 0];
q6 =400*[0 0 0 0 0 1 0];
q7 =1*[0 0 0 0 0 0 1];
Q =[q1;q2;q3;q4;q5;q6;q7];% Same size as A
R = 0.01;

[K,S,e] = lqr(A,B,Q,R);

%Development of L matrix
AL = [ 0 1 0 0 0 0; 0 0 -(m1*g)/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -(M+m1)*g/(M*l1) 0 -m2*g/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -(M+m2)*g/(M*l2) 0];
BL = [0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];
CL1 = [1 0 0 0 0 0];
CL3 = [1 0 0 0 0 0 ;0 0 1 0 0 0 ];
CL4=[1 0 0 0 0 0 ;0 0 1 0 0 0 ;0 0 0 0 1 0 ];

Z = 1000*eye(6) ;
E = 0.000001*eye(6);
a = AL';
b = CL1';
F = Z*E*Z';
r = 0.0001;
[x,l,g] = care(a,b,F,r);
L1 = x*CL1'*(1/r);

%LUENBERGER OBSERVER STATE SPACE PARAMETERS FOR C1; 
AKL1 = AL - L1*CL1; 
BL1 = [BL L1]; 
DO1  = [0 0 ;0 0 ;0 0 ;0 0 ;0 0 ;0 0 ];

AKL = AKL1
BL = BL1;   
C  = C1;
CO = eye(6);
D = D1; 
DO = DO1;
X0 = X01;





