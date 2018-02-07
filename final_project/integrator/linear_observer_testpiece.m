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
CL3 = [1 0 0 0 0 0 ;0 0 0 0 1 0 ];
CL4=[1 0 0 0 0 0 ;0 0 1 0 0 0 ;0 0 0 0 1 0 ];

%Observation - Placement of Poles for A-LC doesn't affect the error in convergence :( 
p1 = 11*transpose(eig(A-B*K));
p1L = p1(1:6);
%p2 = [-5 -6 -7 -8 -9 -10];
%p3 = [-10 -11 -12 -13 -14 -15];

L1 = transpose(place(AL',CL1',p1L)); L1 = real(L1);
L3 = transpose(place(AL',CL3',p1L)); L3 = real(L3);
L4 = transpose(place(AL',CL4',p1L)); L4 = real(L4);

%LUENBERGER OBSERVER STATE SPACE PARAMETERS FOR C1; 
AKL1 = AL - L1*CL1; 
BL1 = [BL L1]; 
DO1  = [0 0 ;0 0 ;0 0 ;0 0 ;0 0 ;0 0 ];

%LUENBERGER OBSERVER STATE SPACE PARAMETERS FOR C3
AKL3 = AL - L3*CL3; 
BL3 = [BL L3];
DO3 = [0 0 0; 0 0 0; 0 0 0;0 0 0; 0 0 0; 0 0 0];  

%LUENBERGER OBSERVER STATE SPACE PARAMETERS FOR C4
AKL4 = AL - L4*CL4; 
BL4 = [BL L4]; 
DO4 = [0 0 0 0; 0 0 0 0; 0 0 0 0;0 0 0 0; 0 0 0 0; 0 0 0 0];


AKL = AKL4;
BL = BL4;   
C  = C4;
CO = eye(6);
D = D4; 
DO = DO4;
X0 = X04;





