%Part 1 - Finding State Space System for the given Non-Linear Model A,B,C,D
clear
clc

%Operating Conditions 
m=.120;
g=9.80;
y1o=2;
yco=12;
y2o=-2;
y12o=yco+y2o-y1o;

a=1.65;
b=6.2;
c=2.69;
d=4.2;

%Initial condition 
xi1 = [-.01 0];
xi2 = [.01 0];

%Simulation Time
ts=20;t=[0:0.1:ts];
%Let x1=y1 x2=y1' x3=y2 x4=y2' be the states of the system

%Linearizing system about operating point we get
x1o=y1o;x3o=y2o;x2o=0;x4o=0;
% x1o'=0;x2o'=0;x3o'=0;x4o'=0;
u1o=((c/((yco+x3o-x1o+d)^4))+m*g)*(a*(x1o+b)^4);
u2o=((-c/((yco+x3o-x1o+d)^4))+m*g)*(a*(-x3o+b)^4);

A_l = [0 1 0 0;((-4*u1o/(a*m*((x1o+b)^5)))-(4*c/(m*((yco+x3o-x1o+d)^5)))) 0 (4*c/(m*(yco+x3o-x1o+d)^5)) 0;0 0 0 1;(4*c/(m*(yco+x3o-x1o+d)^5)) 0 ((4*u2o/(a*m*(-x3o+b)^5))-(4*c/(m*(yco+x3o-x1o+d)^5))) 0]
B_l = [0 0;1/(a*m*(x1o+b)^4) 0;0 0;0 1/(a*m*(-x3o+b)^4)]
C_l = [1 0 0 0;0 0 1 0]
D_l = [0 0;0 0]
sys = ss(A_l,B_l,C_l,D_l)
tf(sys)

%Decoupling MIMO sytem into two SISO systems 
A1=A_l(1:2,1:2);
B1=B_l(1:2,1);
C1=C_l(1,1:2);
D1=0;
sys1 = ss(A1,B1,C1,D1);
TF1=tf(sys1);

A2=A_l(3:4,3:4);
B2=B_l(3:4,2);
C2=C_l(2,3:4);
D2=0;
sys2 = ss(A2,B2,C2,D2);
TF2=tf(sys2);

%canonical forms for the systems 
%system 1
c_obs1 = canon(sys1,'companion')
c_jordan1 = canon(sys1,'model')
c_cont1 = transpose(c_obs1)
%system 2 
c_obs2 = canon(sys2,'companion')
c_jordan2 = canon(sys2,'model')
c_cont2 = transpose(c_obs2)

%Response for system1 and system 2
figure(2)
opt=stepDataOptions('StepAmplitude',1);
subplot(2,2,1)
y1=step(sys1,opt);
plot(y1);
title('Step Response-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,2,2)
y2= step(sys2,opt);
plot(y2);
title('Step Response-System 2')
xlabel('Time')
ylabel('{\delta}y2')
subplot(2,2,3)
y1=impulse(sys1);
plot(y1);
title('Impulse Response-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,2,4)
y2=impulse(sys2);
plot(y2);
title('Impulse Response-System 2')
xlabel('Time')
ylabel('{\delta}y2')

%Root Locus
figure(3)
subplot(2,1,1)
rlocus(sys1);title('Root Locus-System 1');
subplot(2,1,2)
rlocus(sys2);title('Root Locus-System 2');

%Bode Plot
figure(4)
subplot(2,1,1)
title('Bode Plot-System 1')
bode(sys1);
subplot(2,1,2)
title('Bode Plot-System 2')
bode(sys2);

%PID control
%K1=pidtune(sys1,'PID');
P1 =pid(9300,7700,2700);
PID1 =feedback(P1*TF1,1);
P2=pid(8400,2100,2700);
PID2 =feedback(P2*TF2,1);

%step input response for closed loop systems with PID control
figure(5)
subplot(2,3,1)
y1 =step(PID1,opt,t);
plot(t,y1);
title('Step Response-System 1')
xlabel('Time')
ylabel('{\delta}y1')
y2 = step(PID2,opt,t);
subplot(2,3,4)
plot(t,y2);
title('Step Response-System 2' )
xlabel('Time')
ylabel('{\delta}y2')
r = square(t);
subplot(2,3,2)
y1 =lsim(PID1,r,t);
plot(t,y1);
title('Square Input Response-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,3,5)
y2 = lsim(PID2,r,t);
plot(t,y2);
title('Square Input Response-System 2')
xlabel('Time')
ylabel('{\delta}y2')
r = sin(t);
subplot(2,3,3)
y1 =lsim(PID1,r,t);
plot(t,y1);
title('Sinosoidal Input Response-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,3,6)
y2 = lsim(PID2, r, t);
plot(t,y2);
title('Sinosoidal Input Response-System 2')
xlabel('Time')
ylabel('{\delta}y2')

%Control Input Signal for PID compensated closed loop
figure(6)
suptitle('Control Input signal')
r = ones(size(t));
subplot(2,3,1)
u1 = feedback(sys1,P1);
u =lsim(u1,r,t);
plot(t,u);
title('Step Reference Input-System 1')
xlabel('Time')
ylabel('{\delta}u1')
subplot(2,3,4)
u2 = feedback(sys2,P2);
u =lsim(u2,r,t);
plot(t,u);
title('Step Reference Input-System 2')
xlabel('Time')
ylabel('{\delta}u2')
r = sin(t);
subplot(2,3,2)
u1 = feedback(sys1,P1);
u =lsim(u1,r,t);
plot(t,u);
title('Sinsoidal Reference Input-System 1')
xlabel('Time')
ylabel('{\delta}u1')
subplot(2,3,5)
u2 = feedback(sys2,P2);
u =lsim(u2,r,t);
plot(t,u);
title('Sinsoidal Reference Input-System 2')
xlabel('Time')
ylabel('{\delta}u2')
r = square(t);
subplot(2,3,3)
u1 = feedback(sys1,P1);
u =lsim(u1,r,t);
plot(t,u);
title('Square Reference Input-System 1')
xlabel('Time')
ylabel('{\delta}u1')
subplot(2,3,6)
u2 = feedback(sys2,P2);
u =lsim(u2,r,t);
plot(t,u);
title('Square Reference Input-System 2')
xlabel('Time')
ylabel('{\delta}u2')

%Pole placement for system 1 and system 2 with settling time of 5 seconds
p1 = -1 ;
p2 = -2 ;

%feedback Gains for closed loop system
K1 = place(A1,B1,[p1 p2]);
K2 = place(A2,B2,[p1 p2]);
sys_c1 = ss(A1-B1*K1,B1,C1,D1);
sys_c2 = ss(A2-B2*K2,B2,C2,D2);

%Nbar for sys1 and sys 2
%sys1
[a,b,c,d] = ssdata(sys1);
s =size(a,1);
Z = [zeros([1,s]) 1];
N = inv([a,b;c,d])*Z';
Nx = N(1:s);
Nu = N(1+s);
Nbar1=Nu + K1*Nx;
%sys2
[a,b,c,d] = ssdata(sys2);
s =size(a,1);
Z = [zeros([1,s]) 1];
N = inv([a,b;c,d])*Z';
Nx = N(1:s);
Nu = N(1+s);
Nbar2=Nu + K2*Nx;

%Full State Feedback System Response
figure(8)
r=zeros(size(t));
subplot(2,4,1)
suptitle('Response for Full state feedback system')
[y,t,x]=lsim(sys_c1,Nbar1*r,t,xi1);
plot(t,y);
title('Zero Input Response-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,4,5)
[y,t,x]=lsim(sys_c2,Nbar2*r,t,xi2);
plot(t,y);
title('Zero Input Response-System 2')
xlabel('Time')
ylabel('{\delta}y2')
r=ones(size(t));
subplot(2,4,2)
[y,t,x]=lsim(sys_c1,Nbar1*r,t,xi1);
plot(t,y);
title('Step Input Response-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,4,6)
[y,t,x]=lsim(sys_c2,Nbar2*r,t,xi2);
plot(t,y);
title('Step Input Response-System 2')
xlabel('Time')
ylabel('{\delta}y2')
r=sin(t);
subplot(2,4,3)
[y,t,x]=lsim(sys_c1,Nbar1*r,t,xi1);
plot(t,y);
title('Sinosoidal Input Response-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,4,7)
[y,t,x]=lsim(sys_c2,Nbar2*r,t,xi2);
plot(t,y);
title('Sinosoidal Input Response-System 2')
xlabel('Time')
ylabel('{\delta}y2')
r=square(t);
subplot(2,4,4)
[y,t,x]=lsim(sys_c1,Nbar1*r,t,xi1);
plot(t,y);
title('Square Input Response-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,4,8)
[y,t,x]=lsim(sys_c2,Nbar2*r,t,xi2);
plot(t,y);
title('Square Input Response-System 2')
xlabel('Time')
ylabel('{\delta}y2')

%Control Input Signal for Full State Feedback
figure(9)
r=zeros(size(t));
subplot(2,4,1)
suptitle('Control Input Signal for Full state feedback')
[y,t,x]=lsim(sys_c1,Nbar1*r,t,xi1);
u=Nbar1*r-x*K1';
plot(t,u);
title('Zero Input Reference-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,4,5)
[y,t,x]=lsim(sys_c2,Nbar2*r,t,xi2);
u=Nbar2*r-x*K2';
plot(t,u);
title('Zero Input Reference-System 2')
xlabel('Time')
ylabel('{\delta}y2')
r=ones(size(t));
subplot(2,4,2)
[y,t,x]=lsim(sys_c1,Nbar1*r,t,xi1);
u=Nbar1*r-x*K1';
plot(t,u);
title('Step Input Reference-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,4,6)
[y,t,x]=lsim(sys_c2,Nbar2*r,t,xi2);
u=Nbar2*r-x*K2';
plot(t,u);
title('Step Input Reference-System 2')
xlabel('Time')
ylabel('{\delta}y2')
r=sin(t);
subplot(2,4,3)
[y,t,x]=lsim(sys_c1,Nbar1*r,t,xi1);
u=Nbar1*r-x*K1';
plot(t,u)
title('Sinosoidal Input Reference-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,4,7)
[y,t,x]=lsim(sys_c2,Nbar2*r,t,xi2);
u=Nbar2*r-x*K2';
plot(t,u)
title('Sinosoidal Input Reference-System 2')
xlabel('Time')
ylabel('{\delta}y2')
r=square(t);
subplot(2,4,4)
[y,t,x]=lsim(sys_c1,Nbar1*r,t,xi1);
u=Nbar1*r-x*K1';
plot(t,u)
title('Sqaure Input Reference-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,4,8)
[y,t,x]=lsim(sys_c2,Nbar2*r,t,xi2);
u=Nbar2*r-x*K2';
plot(t,u)
title('Sqaure Input Reference-System 2')
xlabel('Time')
ylabel('{\delta}y2')

%Full Observer poles
o1 = 4*p1;
o2 = 4*p2;

L1 = place(A1',C1',[o1 o2])';
L2 = place(A2',C2',[o1 o2])';

%Full Observer 1 (Regulation Case)
Ae1 = A1-L1*C1;
Be1 = L1;
Ce1 = C1;
De1 = D1;
syse1 = ss(Ae1,Be1,Ce1,De1);

%Full Observer 2 (Regulation Case)
Ae2 = A2-L2*C2;
Be2 = L2;
Ce2 = C2;
De2 = D2;
syse2 = ss(Ae2,Be2,Ce2,De2);

%Response of Full State Oberserver to input (Regulation Case)
figure(10)
subplot(2,3,1)
title('Full order observer Response')
y = lsim(syse1,ones(size(t)),t);
plot(t,y);
title('Step Response-Observer 1')
xlabel('Time')
ylabel('$$\hat{y1}$$','Interpreter','Latex')
subplot(2,3,4)
title('Full order observer Response')
y = lsim(syse2,ones(size(t)),t);
plot(t,y);
title('Step Response-Observer 2')
xlabel('Time')
ylabel('$$\hat{y2}$$','Interpreter','Latex')
subplot(2,3,2)
title('Full order observer Response')
y = lsim(syse1,sin(t),t);
plot(t,y);
title('Sinosoidal Response-Observer 1')
xlabel('Time')
ylabel('$$\hat{y1}$$','Interpreter','Latex')
subplot(2,3,5)
title('Full order observer Response')
y = lsim(syse2,sin(t),t);
plot(t,y);
title('Sinosoidal Response-Observer 2')
xlabel('Time')
ylabel('$$\hat{y2}$$','Interpreter','Latex')
subplot(2,3,3)
title('Full order observer Response')
y = lsim(syse1,square(t),t);
plot(t,y);
title('Square Response-Observer 1')
xlabel('Time')
ylabel('$$\hat{y1}$$','Interpreter','Latex')
subplot(2,3,6)
title('Full order observer Response')
y = lsim(syse2,square(t),t);
plot(t,y);
title('Square Response-Observer 2')
xlabel('Time')
ylabel('$$\hat{y1}$$','Interpreter','Latex')

%Observer Controllers TF (Regulation)
%System 1
Aec1 = A1-B1*K1-L1*C1;
Bec1 = L1;
Cec1 = -K1;
Dec1 = D1;

sysec1= ss(Aec1,Bec1,Cec1,Dec1);
TFec1 = tf(sysec1);

%System 2
Aec2 = A2-B2*K2-L2*C2;
Bec2 = L2;
Cec2 = -K2;
Dec2 = D2;

sysec2= ss(Aec2,Bec2,Cec2,Dec2);
TFec2 = tf(sysec2)

%Reduced Observer 
%System 1 
a11_1=A1(1:1,1:1);A1e_1=A1(1:1,2:2);Ae1_1=A1(2:2,1:1);Aee_1=A1(2:2,2:2);
Be_1 = B1(2:2);b1_1 = B1(1:1);

Lr1= place(Aee_1',A1e_1',[o1]);

Ar1=(Aee_1-Lr1*A1e_1);Br1=(Ae1_1-Lr1*a11_1+Aee_1*Lr1-Lr1*A1e_1*Lr1);Cr1=1;Dr1=0;
sysre1=ss(Ar1,Br1,Cr1,Dr1);
%System 2
a11_2=A2(1:1,1:1);A1e_2=A2(1:1,2:2);Ae1_2=A2(2:2,1:1);Aee_2=A2(2:2,2:2);
Be_2 = B2(2:2);b1_2 = B2(1:1);
Lr2= place(Aee_2',A1e_2',[o1]);

Ar2=(Aee_2-Lr2*A1e_2);Br2=(Ae1_2-Lr2*a11_2+Aee_2*Lr2-Lr2*A1e_2*Lr2);
Cr2=1;Dr2=0;
sysre2=ss(Ar2,Br2,Cr2,Dr2);

figure(13)
subplot(2,3,1)
y = lsim(sysre1,ones(size(t)),t);
plot(t,y);
title('Step Response-Reduced Observer 1')
xlabel('Time')
ylabel('$$\hat{x2}$$','Interpreter','Latex')
subplot(2,3,4)
title('Full order observer Response')
y = lsim(sysre2,ones(size(t)),t);
plot(t,y);
title('Step Response-Reduced Observer 2')
xlabel('Time')
ylabel('$$\hat{x4}$$','Interpreter','Latex')
subplot(2,3,2)
title('Full order observer Response')
y = lsim(sysre1,sin(t),t);
plot(t,y);
title('Sinosoidal Response-Reduced Observer 1')
xlabel('Time')
ylabel('$$\hat{x2}$$','Interpreter','Latex')
subplot(2,3,5)
title('Full order observer Response')
y = lsim(sysre2,sin(t),t);
plot(t,y);
title('Sinosoidal Response-Reduced Observer 2')
xlabel('Time')
ylabel('$$\hat{x4}$$','Interpreter','Latex')
subplot(2,3,3)
title('Full order observer Response')
y = lsim(sysre1,square(t),t);
plot(t,y);
title('Square Response-Reducecd Observer 1')
xlabel('Time')
ylabel('$$\hat{x2}$$','Interpreter','Latex')
subplot(2,3,6)
title('Full order observer Response')
y = lsim(sysre2,square(t),t);
plot(t,y);
title('Square Response-Reduced Observer 2')
xlabel('Time')
ylabel('$$\hat{x4}$$','Interpreter','Latex')

%Closed Loop System 1 with Observer Feedback Control
Ao1 = [ A1-B1*K1            B1*K1
       zeros(size(A1))    A1-L1*C1 ];
Bo1 = [    B1*Nbar1
       zeros(size(B1)) ];
Co1 = [ C1    zeros(size(C1)) ];

sysob1 = ss(Ao1,Bo1,Co1,0);

%Closed Loop Observer 2 with Observer Feedback Control
Ao2 = [ A2-B2*K2            B2*K2
       zeros(size(A2))    A2-L2*C2 ];
Bo2 = [    B2*Nbar2
       zeros(size(B2)) ];
Co2 = [ C2    zeros(size(C2)) ];

sysob2 = ss(Ao2,Bo2,Co2,0);

% Response of closed loop system with oberserver state feedback control
figure(11)
subplot(2,4,1)
suptitle('Closed Loop Response with observer controller')
lsim(sysob1,zeros(size(t)),t,[xi1 xi1]);
title('Zero Reference Input-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,4,5)
lsim(sysob2,zeros(size(t)),t,[xi2 xi2]);
title('Zero Reference Input-System 2')
xlabel('Time')
ylabel('{\delta}y2')
subplot(2,4,2)
y = lsim(sysob1,ones(size(t)),t,[xi1 xi1]);
plot(t,y)
title('Step Reference Input-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,4,6)
y = lsim(sysob2,ones(size(t)),t,[xi2 xi2]);
plot(t,y)
title('Step Reference Input-System 2')
xlabel('Time')
ylabel('{\delta}y2')
subplot(2,4,3)
y = lsim(sysob1,sin(t),t,[xi1 xi1]);
plot(t,y)
title('Sinosoidal Reference Input-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,4,7)
y = lsim(sysob2,sin(t),t,[xi2 xi2]);
plot(t,y)
title('Sinosoidal Reference Input-System 2')
xlabel('Time')
ylabel('{\delta}y2')
subplot(2,4,4)
y = lsim(sysob1,square(t),t,[xi1 xi1]);
plot(t,y)
title('Square Reference Input-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,4,8)
y = lsim(sysob2,square(t),t,[xi2 xi2]);
plot(t,y)
title('Square Reference Input-System 2')
xlabel('Time')
ylabel('{\delta}y2')

%Response Comparison of PID with Estimated Feedback Control
figure(12)
r=ones(size(t));
subplot(2,3,1)
hold on 
y1 =step(PID1,opt,t);
plot(t,y1);
[y,t,x]=lsim(sysob1,r,t,[xi1 xi1]);
plot(t,y)
title('Step Response-System 1')
xlabel('Time')
ylabel('{\delta}y1')
y2 = step(PID2,opt,t);
subplot(2,3,4)
hold on
plot(t,y2);
[y,t,x]=lsim(sysob2,r,t,[xi2 xi2]);
plot(t,y)
title('Step Response-System 2' )
xlabel('Time')
ylabel('{\delta}y2')
r = square(t);
subplot(2,3,3)
hold on
y1 =lsim(PID1,r,t);
plot(t,y1);
[y,t,x]=lsim(sysob1,r,t,[xi1 xi1]);
plot(t,y)
title('Square Input Response-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,3,6)
hold on
y2 = lsim(PID2,r,t);
plot(t,y2);
[y,t,x]=lsim(sysob2,r,t,[xi2 xi2]);
plot(t,y)
title('Square Input Response-System 2')
xlabel('Time')
ylabel('{\delta}y2')
r = sin(t);
subplot(2,3,2)
hold on
y1 =lsim(PID1,r,t);
plot(t,y1);
[y,t,x]=lsim(sysob1,r,t,[xi1 xi1]);
plot(t,y)
title('Sinosoidal Input Response-System 1')
xlabel('Time')
ylabel('{\delta}y1')
subplot(2,3,5)
hold on
y2 = lsim(PID2, r, t);
plot(t,y2);
[y,t,x]=lsim(sysob2,r,t,[xi2 xi2]);
plot(t,y)
title('Sinosoidal Input Response-System 2')
xlabel('Time')
ylabel('{\delta}y2')