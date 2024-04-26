
clc
clear all
%% A Membership Function-Dependent Lâˆ? Performance Based Output Interval Estimation Method for
%% Discrete-Time Takagiâ€“Sugeno Fuzzy Systems Applied to Fault Detection
% We are considering the system 
% A1=[-6.56621004566210	-31.7459613049935;
% 0	-1.13904731051587];
% A2=[-6.56621004566210	-228.294967643244;
%  0	-1.99993195827741];
% B=[1	0	2.28310502283105;0	1	0];
% D=[1;0];
%
% author:  Yi Li, Jiuxiang Dong
% contact: 1329031564@qq.com


%% Continuous model
A1=[-6.56621004566210	-31.7459613049935;
0	-1.13904731051587];
A2=[-6.56621004566210	-228.294967643244;
 0	-1.99993195827741];
B=[1	0	2.28310502283105;0	1	0];
C=[1 0;0 1];
D=[1;0];
%% Discrete model
T1=0.1;
A1=eye(2)+T1*A1;
B1=T1*B;
W1_1=T1*D;
A2=eye(2)+T1*A2;
B2=T1*B;
W1_2=T1*D;
W2=[0;0];
%% Control input
u=[350;1;1100];

G = zeros(2,2);
Z = eye(2);


%% Weight parameter
rou_1 = 1;
rou_2 = 1; 

%% Noisy signal
% noise=0.5-rand(1,10000);

%% Disturbance distribution matrix
F_w1=[0.07;0.8];
F_w2=[0.1;1];
F_v2=[0.3;1.2];
F_v1=[0.2;0.6];

%% Fault distribution matrix
F11=[1;0.2];
F12=[1;0.3];
F2=[0.1;0.1];


%% Variable definition of linear matrix inequalities
miu = sdpvar(1,1);
P1 = sdpvar(2,2);
P2 = sdpvar(2,2);
fi1=zeros(2,2);
fi21 = sdpvar(2,2,'full');
fi22 = sdpvar(2,2,'full');
Y = sdpvar(2,2,'full');

lan=0.5;
gama =1.238;


L1=[miu>=0];

L21 = [P1>=0 ];
L22 = [P2>=0 ];

% L31=[C'*C-rou_1*gama*lan*P1 C'*F_v1;...
%     (C'*F_v1)' F_v1'*F_v1-rou_1*gama*eye(1)*rou_1*gama*eye(1)+rou_1*gama*eye(1)*miu]<=0;
% L32=[C'*C-rou_2*gama*lan*P2 C'*F_v2;...
%     (C'*F_v2)' F_v2'*F_v2-rou_2*gama*eye(1)*rou_2*gama*eye(1)+rou_2*gama*eye(1)*miu]<=0;


L31=[-lan*P1 zeros(2,1) C';
    zeros(1,2) -(rou_1*gama-miu) F_v1';
    C F_v1 -rou_1*gama*eye(2)]<=0;

L32=[-lan*P2 zeros(2,1) C';
    zeros(1,2) -(rou_2*gama-miu) F_v2';
    C F_v2 -rou_2*gama*eye(2)]<=0;


L411=[(lan-1)*P1 zeros(2,1) (Y*A1-fi21*C)';...
    zeros(1,2) -miu*eye(1) (Y*F_w1-fi21*F_v1)';...
    Y*A1-fi21*C Y*F_w1-fi21*F_v1 (P1-Y-Y')]<=0;

L412=[(lan-1)*P1 zeros(2,1) (Y*A2-fi22*C)';...
    zeros(1,2) -miu*eye(1) (Y*F_w2-fi22*F_v2)';...
    Y*A2-fi22*C Y*F_w2-fi22*F_v2 (P2-Y-Y')]<=0;


L421=[(lan-1)*P2 zeros(2,1) (Y*A1-fi21*C)';...
    zeros(1,2) -miu*eye(1) (Y*F_w1-fi21*F_v1)';...
    Y*A1-fi21*C Y*F_w1-fi21*F_v1 (P1-Y-Y')]<=0;

L422=[(lan-1)*P2 zeros(2,1) (Y*A2-fi22*C)';...
    zeros(1,2) -miu*eye(1) (Y*F_w2-fi22*F_v2)';...
    Y*A2-fi22*C Y*F_w2-fi22*F_v2 (P2-Y-Y')]<=0;




LL =L1+L21+L22+L31+L32+L411+L412+L421+L422;
%% solve LMI


solvesdp(LL)


gama= double(gama)
P1 = double(P1);
P2 = double(P2);
Y = double(Y);
fi1 = double(fi1);
fi21 = double(fi21);
fi22 = double(fi22);
lan = double(lan);
gama = double(gama);
miu = double(miu);
%% Calculate observer gain
L1_tra=inv(Y)*fi1;
L21_tra=inv(Y)*fi21;
L22_tra=inv(Y)*fi22;

%% eigenvalue
eigenvalues1 = eig(P1);
eigenvalues2 = eig(P2);
maxEigenvalue1 = max(eigenvalues1);
maxEigenvalue2 = max(eigenvalues2);
maxEigenvalue_tra = max([maxEigenvalue1,maxEigenvalue2]);
%% Save result
save L1_tra L1_tra;
save L21_tra L21_tra;
save L22_tra L22_tra;
save maxEigenvalue_tra maxEigenvalue_tra;
