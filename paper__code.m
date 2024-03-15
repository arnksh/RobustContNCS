%% Adaptive Kalman Filtering with parameter estimation in NCS with Random Sensor Delays,Packet Dropouts and Missing Measurements.
% System Structure and Dimensions
% x(x+1)=A(th)x(k)+Bw(k)+gu(k)  ; Where A(th) is the continous function of unknown constant parameter th.
% y(k+1)=Cx(k+1)+v(k+1); A(2,2),B(2,1),g(2,1),C(1,2) and w(k),v(k),u(k) are scalar input.
clear all
%% System model matrices.
% a=[.8 0;.9 .2];
% a=[0 1;1e-5 .99998];
% a=[0 0.0590;-0.0008 0.9672];
a=[1.7240 -0.74383;1 0];
dA=a;
B=[1;.5];
dB=0;
g=[1;1];
dg=0;
C=[1 1];
dC=0;
Pw=0.5; Pv=.25;
bm=[0.6 1.4];
% l=2.5;%for randon variable theta
% l=1.4;
l=0.4; %for const. theta with contr.
[n,m]=size(B);
r=[.8 .02 .09 .09]; %probability for uncertain cases 1,2,3,4
u=0;
du=0;
% tht=ranumb(bm(1),bm(2),1,100);
% k=0:400;
% tht=1+0.4*sin(k*pi/200);
tht(1:6000)=0.7;
tht(6001:12000)=1.3;

N=length(tht);
%% Controller Design
[K,SHI]=dcontr122(a,bm,g);
% [K1,K2]=dcontr123(a,bm,g)
%% Initializations
x=[10;8]; xs=zeros(n,1); th=1;
xa=zeros(n,100); xsa=zeros(n,100); tha=zeros(1,100); tr_e=zeros(1,100);
xa(:,1)=x;
Xs=[xs;xs];
y=0;
tha(1)=th;
F=[eye(2*n) zeros(2*n,1)];
T=[eye(n) zeros(n)];
p0=20*eye(n);
P=[p0 p0;p0 p0];%initial error covariance of FX-Xs
PW=[Pw zeros(m,1) zeros(m,1);zeros(1,m) Pv zeros(1);zeros(1,m) zeros(1) Pv];
w=ranum(.8191,0,sqrt(Pw),m,40000);  %process noise
v=ranum(.2047,0,sqrt(Pv),1,40000);  %measurement noise
vk_=v(:,1);
%% Gradient initialization dxs=Gradient of xs(1/0)
% for k=0 to 1
dxs=zeros(n,1);  % diff. of xs(1/0) for all components of the theta.
A=a*th;
dAs=[dA zeros(n);zeros(n) zeros(n)];
dP=zeros(2*n,2*n);
%% Stchatic parameter for packet dropout
%Representation of random packet dropout:0=packet droped and 1=packet arrived
load packet
pas=pas(2,:);
%% Coeff. matrices of the model with no uncertainty,one step sensor delay,missing mearsurements.
B1=[B zeros(n,1) zeros(n,1);zeros(n,m) zeros(n,1) zeros(n,1);zeros(1,m) eye(1) zeros(1)];
B2=[B zeros(n,1) zeros(n,1);zeros(n,m) zeros(n,1) zeros(n,1);zeros(1,m) zeros(1) eye(1)];
B3=[B zeros(n,1) zeros(n,1);zeros(n,m) zeros(n,1) zeros(n,1);zeros(1,m) eye(1) zeros(1)];
C1=[C zeros(1,n) zeros(1)];
C2=[zeros(1,n) C zeros(1)];
C3=[zeros(1,n) zeros(1,n) zeros(1)];
I=[zeros(2*1,m) eye(2*1)];
H1=[eye(1) zeros(1)]*I;
H2=[zeros(1) eye(1)]*I;
H3=[eye(1) zeros(1)]*I;
piq=zeros(1,3);
for q=1:3
    piq(q)=r(q)/sum(r(1:3));
end
CS=piq(1)*C1+piq(2)*C2+piq(3)*C3;
HS=piq(1)*H1+piq(2)*H2+piq(3)*H3;
Cs=CS*F';
s1=piq(1)*F*B1*PW*B1'*F'+piq(2)*F*B2*PW*B2'*F'+piq(3)*F*B3*PW*B3'*F';
s2=piq(1)*H1*PW*H1'+piq(2)*H2*PW*H2'+piq(3)*H3*PW*H3';
%% Filter algorithm for k>0
for k=1:N
    As=[A zeros(n);eye(n) zeros(n)];
    A0=a*tht(k);
    xk=A0*x+B*w(:,k)+g*u;
    X=[xk' x' y']';
    W=[w(:,k)' v(:,k)' vk_']';
    
    xsk=A*xs+g*u; %One step prediction of state.
    dxsk=dA*xs+A*dxs+g*du;
    Xs=[xsk' xs']';
    dXs=[dxsk' dxs']';
    
    % Predicted covariance
    dP=dAs*P*As'+As*dP*As'+As*P*dAs';
    P=As*P*As'+s1;
    
    if pas(k)==0
        yk=y;
        Ls=zeros(n,1);
        dLs=zeros(n,1);
    else
        yk=CS*X+HS*W;
        Ls=T*P*Cs'/(Cs*P*Cs'+s2); % Gain calculation
        dLs=T*dP*Cs'/(Cs*P*Cs'+s2)+(T*P*Cs'/(Cs*P*Cs'+s2))*Cs*dP*Cs'/(Cs*P*Cs'+s2);
    end
    
    yn=yk-Cs*Xs;
    dyn=-Cs*dXs;
    xs=xsk+Ls*yn;   % Updates state
    dxs=dxsk+dLs*yn+Ls*dyn;
    
    %Update of cov. and gradient of cov.
    Rs=[Ls' zeros(1,n)]';
    dRs=[dLs' zeros(1,n)]';
    s3=piq(1)*Rs*H1*PW*H1'*Rs'+piq(2)*Rs*H2*PW*H2'*Rs'+piq(3)*Rs*H3*PW*H3'*Rs';
    dP=-dRs*Cs*P*(eye(2*n)-Rs*Cs)'-(eye(2*n)-Rs*Cs)*P*(dRs*Cs)'+(eye(2*n)-Rs*Cs)*dP*(eye(2*n)-Rs*Cs)'+dRs*s2*Rs'+Rs*s2*dRs';
    P=(eye(2*n)-Rs*Cs)*P*(eye(2*n)-Rs*Cs)'+s3;

    % Parameter update
    thp=th-l*dyn'*yn/k;
    if thp<bm(1)
       th=bm(1);
    else if thp>bm(2)
            th=bm(2);
        else
            th=thp;
        end
    end
    A=a*th; %Update of the system matrix.
%     
    u=K*xs;
    du=K*dxs;

    %Storing and updating the loop variables.
    ek=T*P*T';
    tr_e(k)=trace(ek);
    th_e(k)=tht(k)-th;
    
    y=yk;
    xa(:,k+1)=xk;
    x=xk;
    xsa(:,k+1)=xs;
    vk_=v(:,k);
    tha(k+1)=th;
end
%% Plotings
figure
for i=1:n
    subplot(1,2,i)
    hold on
    box on
    plot(xa(i,:),'b','linewidth',1.5)
    plot(xsa(i,:),'r--','linewidth',1.5)
    grid
    title(['x',int2str(i),' WITHOUT CONTR."'],'FontSize',12)
    xlabel('time k','FontSize',16)
    ylabel(['x',int2str(i)],'FontSize',16)
    legend('Actual state','Estimated state')
    hold off
end
figure
box on
% subplot(1,2,1)
plot(tr_e,'linewidth',1.5)
grid
title('\theta=Constant','FontSize',16)
xlabel('time k','FontSize',16)
ylabel('Tr. Error Cov.','FontSize',16)
% subplot(1,2,2)
% plot(tr_e1,'linewidth',1.5)
% grid
% title('\theta=Randomly changing','FontSize',16)
% xlabel('time k','FontSize',16)
% ylabel('Tr. Error Cov.','FontSize',16)
    
figure
box on
hold on
plot(tha,'r--','linewidth',1.5)
plot(tht,'b')
legend('Estimated \theta','Actual \theta')
grid
title('\theta WITH CONTR.','FontSize',12)
xlabel('time k','FontSize',16)
ylabel('\theta(k)','FontSize',16)

figure
box on
plot(th_e,'r','linewidth',1.5)
grid
title('Error between actual and estimated \theta with contr.','FontSize',12)
xlabel('time k','FontSize',16)
ylabel('\theta_e(k)','FontSize',16)
%% end
