clc
clear all
%% A Membership Function-Dependent L∞ Performance Based Output Interval Estimation Method for
%% Discrete-Time Takagi–Sugeno Fuzzy Systems Applied to Fault Detection
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
rou_1 = 0.75;
rou_2 = 1; 
gama=1.238;
lan=0.5;
%% Noisy signal
% noise=0.5-rand(1,10000);
load noise.mat;
%% Disturbance distribution matrix
F_w1=[0.07;0.8];
F_w2=[0.1;1];
F_v2=[0.3;1.2];
F_v1=[0.2;0.6];

%% Fault distribution matrix
F11=[1;0];
F12=[2;0];
F2=[1;0];

%% Observer gain
load L1_th2.mat
load L1_tra.mat
load L21_th2.mat
load L21_tra.mat
load L22_th2.mat
load L22_tra.mat
load maxEigenvalue_th2.mat
load maxEigenvalue_tra.mat

%% initialize
x=[347.74;0.6];
xh_tra=  [347.74;0.6];
xh_th2=  [347.74;0.6];


y=C*x;

r=[0,0];
u=[350;1;880];



 e_th2_0 = x - xh_th2;
 V_th2_0 = e_th2_0' * maxEigenvalue_th2* e_th2_0;

 e_tra_0 = x - xh_tra;
 V_tra_0 = e_tra_0' * maxEigenvalue_tra* e_tra_0;

n=1;
d(:,1)=0;
w11=1;
w22=0;

for i=1:1:200
w1q=w11;
w2q=w22;
%% fault
   if i<=100
       f=0;
   else
        f=0;
   end
 %% control input
u=[350;1;860];
   %% 隶属度函数更新
      if (y(1,1)<324.4)
           w1(i)=1;
   else if y(1,1)>350
           w1(i)=0;
   else
           w1(i)=cos((y(1,1)-324.4)*pi/(25.6*2));
   end
   end

   w2(i)=1-w1(i);


w11=w1(i);
w22=w2(i);

   %% System matrix
   A=w1(i)*A1+w2(i)*A2;
   B=w1(i)*B1+w2(i)*B2;
   W1=w1(i)*W1_1+w2(i)*W1_2;
   F_w=w1(i)*F_w1+w2(i)*F_w2;
   F_v=w1(i)*F_v1+w2(i)*F_v2;
    F1=w1(i)*F11+w2(i)*F12;
 %% System update
   xq=x;
   yq=y;
   x=A*x+B*u+F_w*noise(i)+F1*f;
   y=C*xq+1*F_v*noise(i)+F2*f;

  %%  observer based on Theorem 2
      
      
      L_th2 = L21_th2 * w1(i) +L22_th2 * w2(i) ;
      xhq_th2=xh_th2;
      xh_th2 = A *xh_th2 + B * u + L_th2 *( y - C*xh_th2 )  ; 
      yh_th2=C*xhq_th2;

   %% observer based on traditional H∞ performance

      L_tra = L21_tra * w1(i) +L22_tra * w2(i) ;
      xhq_tra=xh_tra;
      xh_tra = A *xh_tra + B * u + L_tra *( y - C*xh_tra )  ; 
      yh_tra=C*xhq_tra;






   %% 记录状态和状态估计值
   x_m(:,i)=x;
  x_mh_th2(:,i)=xh_th2;
  x_mh_tra(:,i)=xh_tra;

 %% L∞
omega=0.5;
gama_final_th2(i) = (rou_1 * w1q + rou_2 * w2q ) * gama;
eee_th2 (i)= sqrt(gama_final_th2(i) * ( lan * power((1-lan),i) * V_th2_0 + 1 * gama_final_th2(i) * power(omega,2) ));  


eee_tra (i)= sqrt(gama* ( lan * power((1-lan),i) * V_tra_0 + 1 * gama * power(omega,2) ));  




x_mq(:,i)=xq;
y_m(:,i)=y;


x_mhq_th2(:,i)=xhq_th2;
x_mhq_tra(:,i)=xhq_tra;
y_m_gj_th2=C * x_mhq_th2;
y_m_gj_tra=C * x_mhq_tra;

 y_m_gj1_th2(1,i) = yh_th2(1,1);
 y_m_gj2_th2(1,i) = yh_th2(2,1);

 y_m_gj1_tra(1,i) = yh_tra(1,1);
 y_m_gj2_tra(1,i) = yh_tra(2,1);



 y_m_zj = y_m;
 
 y_m_1(1,i) = y(1,1);
 y_m_2(1,i) = y(2,1);





 e_y_x_tra(i,:) = -eee_tra(i);  
 e_y_s_tra(i,:) = +eee_tra(i); 
 
 y_x_1_tra(i,:) = y_m_gj1_tra(:,i) + e_y_x_tra(i); 
 y_s_1_tra(i,:) = y_m_gj1_tra(:,i) + e_y_s_tra(i); 

 y_x_2_tra(i,:) = y_m_gj2_tra(:,i) + e_y_x_tra(i); 
 y_s_2_tra(i,:) = y_m_gj2_tra(:,i) + e_y_s_tra(i); 

 e_y_x_th2(i,:) = -eee_th2(i);  
 e_y_s_th2(i,:) = +eee_th2(i); 
 
 y_x_1_th2(i,:) = y_m_gj1_th2(:,i) + e_y_x_th2(i); 
 y_s_1_th2(i,:) = y_m_gj1_th2(:,i) + e_y_s_th2(i); 

 y_x_2_th2(i,:) = y_m_gj2_th2(:,i) + e_y_x_th2(i); 
 y_s_2_th2(i,:) = y_m_gj2_th2(:,i) + e_y_s_th2(i); 

%% figure 5
if i>47&&i<49


uuu=-pi:0.0001:pi;
a_tra=y_m_gj1_tra(i);
b_tra=y_m_gj2_tra(i);
a_th2= y_m_gj1_th2(i);
b_th2= y_m_gj2_th2(i);
rrr_tra=eee_tra(i);
rrr_th2=eee_th2(i);
xxx_tra=rrr_tra*sin(uuu)+a_tra;
yyy_tra=rrr_tra*cos(uuu)+b_tra;
xxx_th2=rrr_th2*sin(uuu)+a_th2;
yyy_th2=rrr_th2*cos(uuu)+b_th2;
 lll=length(xxx_tra);
 a1=i*ones(lll,1)';
fill3(a1,xxx_th2,yyy_th2,[0,39,45]./255,'EdgeColor',[0,39,45]./255,'LineWidth',0.5);
alpha(0.1)
hold on
plot3(i, y(1), y(2),'xr')
hold on





fill3(a1,xxx_tra,yyy_tra,[180,39,45]./255,'EdgeColor',[180,39,45]./255,'LineWidth',0.5);
alpha(0.1)












    end



      end

b2 = xlabel('Time(k)');
b3 = ylabel('$y_1$');
b4 = zlabel('$y_2$');
b1=legend('Theorem 2 Feasible domain','$y$','[34] Feasible domain');
set([b1],'Interpreter','latex','fontsize',12)
set([b2 b3 b4],'Interpreter','latex','fontsize',14)
hold on

