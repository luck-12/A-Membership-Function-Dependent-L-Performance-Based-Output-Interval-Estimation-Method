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
        f=0*(i-100);
   end
      ff(i)=f;
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





      end

%% 取出各个状态
x_m1=x_m(1,:);
x_m2=x_m(2,:);

xh_m1_th2=x_mh_th2(1,:);
xh_m2_th2=x_mh_th2(2,:);

xh_m1_tra=x_mh_tra(1,:);
xh_m2_tra=x_mh_tra(2,:);



%% 区间观测

% 
% xxx_1 = x_x_1';% x1下限 
% xsx_1 = x_s_1';% x1上限
% xxx_2 = x_x_2';% x2下限 
% xsx_2 = x_s_2';% x2上限 


y_1_1 = y_m_zj(1,:);
y_1_2 = y_m_zj(2,:);


ysx_1_tra =  y_s_1_tra';
yxx_1_tra =  y_x_1_tra';

ysx_2_tra =  y_s_2_tra';
yxx_2_tra =  y_x_2_tra';

ysx_1_th2 =  y_s_1_th2';
yxx_1_th2 =  y_x_1_th2';

ysx_2_th2=  y_s_2_th2';
yxx_2_th2 =  y_x_2_th2';

y_gj1 = y_m_gj_th2(1,:);
y_gj2 = y_m_gj_th2(2,:);


y_gj1 = y_m_gj_tra(1,:);
y_gj2 = y_m_gj_tra(2,:);
% % % % %% y1 区间观测


%% fig 1
figure

Cblue = [0,70,222]/255;
Cgreen = [69,189,155]/255;
Cred = [255,32,39]/255;
Cbblue = [0.49,0.18,0.56];



subplot(2,1,1)
plot(gama_final_th2,'-','Color',Cblue,'LineWidth',1.2)
hold on
plot(gama*ones(1,200),'Color',Cred,'LineWidth',1.2)
b3 = ylabel('$\bar \gamma$');
% b4 = zlabel('$r_2$');
b1=legend('Theorem 2','[34]');
set([b1],'Interpreter','latex','fontsize',12)
set([b3],'Interpreter','latex','fontsize',14)
axis([0 200 0.5 2]);
subplot(2,1,2)
plot(w1,'-','Color',Cblue,'LineWidth',1.2)
hold on
 plot(w2,'-','Color',Cred,'LineWidth',1.2)
% hold on

b3 = ylabel('$\sigma$');
% b4 = zlabel('$r_2$');
b1=legend('$\sigma_1$','$\sigma_2$');
set([b1],'Interpreter','latex','fontsize',12)
set([b3],'Interpreter','latex','fontsize',14)
axis([0 200 0 1.5]);




%% fig 2

figure
plot(ysx_1_th2)
hold on
plot(yxx_1_th2)
hold on
plot(y_1_1)
hold on
plot(ysx_1_tra)
hold on
plot(yxx_1_tra)
hold on


b2 = xlabel('Time(k)');
b3 = ylabel('$y_2$');
% b4 = zlabel('$r_2$');
b1=legend('Theorem 2 $\overline{y}_{1k}^{L}$','Theorem 2 $\underline{y}_{1k}^{L}$','$y_1$','[34] $\overline{y}_{1k}^{L}$','[34] $\underline{y}_{1k}^{L}$');
set([b1],'Interpreter','latex','fontsize',12)
set([b2 b3],'Interpreter','latex','fontsize',14)
axis([0 200 320 380]);


%% fig 3
figure
plot(ysx_2_th2)
hold on
plot(yxx_2_th2)
hold on
plot(y_1_2)
hold on
plot(ysx_2_tra)
hold on
plot(yxx_2_tra)
hold on


b2 = xlabel('Time(k)');
b3 = ylabel('$y_2$');
% b4 = zlabel('$r_2$');
b1=legend('Theorem 2 $\overline{y}_{2k}^{L}$','Theorem 2 $\underline{y}_{2k}^{L}$','$y_2$','[34] $\overline{y}_{2k}^{L}$','[34] $\underline{y}_{2k}^{L}$');
set([b1],'Interpreter','latex','fontsize',12)
set([b2 b3],'Interpreter','latex','fontsize',14)

axis([0 200 -1 5]);

