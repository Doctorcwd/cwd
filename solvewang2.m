close all;    
clear;  
clc;                 
x0=[0.1;0.15;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];     %初始化%
tspan=0:0.01:20;
[t,x]=ode15s('wang2',tspan,x0);                             
w=1-0.35*cos(t);
%  w=-2-0.5*cos(t); w=-1-0.5*cos(t);
u=zeros(length(t),1);
event1=zeros(length(t),1);
z1=zeros(length(t),1);
M1=0;M1last=0;
%%
% figure;  
% plot(t,x(:,1),'k',t,yr,'--r','linewidth',1.5);                            
% xlabel('t(s)','fontsize',16,'FontWeight','bold','Fontname', 'Times New Roman');   
%  str1='$$x_1$$';
%   str2='$$y_r$$';
% legend({str1,str2},'interpreter','latex','fontsize',16)
% ylabel('Tracking performance','FontWeight','bold','Fontname', 'Times New Roman','fontsize',16);   
%% 
figure;  
plot(t,x(:,2),'g',t,x(:,3),'b','linewidth',1.5);                            
xlabel('t(s)','fontsize',16,'FontWeight','bold','Fontname', 'Times New Roman');  
str2='$$x_2$$';
str3='$$x_3$$';
legend({str2,str3},'interpreter','latex','fontsize',16)
ylabel('States x_2 and x_3','FontWeight','bold','Fontname', 'Times New Roman','fontsize',16)
%   axis([0 20 -15 25])
%   set(gca,'YTick',[-15:10:25])

% %%
% figure;            
% plot(t,x(:,1)-yr,'k','linewidth',1.5);                                            
% xlabel('t(s)','fontsize',16,'FontWeight','bold','Fontname', 'Times New Roman');                                                    
% ylabel('Tracking error','FontWeight','bold','Fontname', 'Times New Roman','fontsize',16); 
% u=zeros(length(t),1);
%% 
figure;  
plot(t,x(:,1),'g',t,x(:,4),'r','linewidth',1.5);                            
xlabel('t(s)','fontsize',16); 
str1='$$x_1$$';
 str2='$$\hat x_1$$';
legend({str1,str2},'interpreter','latex','fontsize',16)
ylabel('State x_1 ,\hat State x_1','FontWeight','bold','Fontname', 'Times New Roman','fontsize',16)
 axis([0 20 -4 4])
  set(gca,'YTick',[-6:2:20])
%%
figure;  
plot(t,x(:,1),'b',t,x(:,1)/w,'g','linewidth',1.5);                            
xlabel('t(s)','fontsize',16); 
str1='$$x_1$$';
 str2='$$\check x_1$$';
legend({str1,str2},'interpreter','latex','fontsize',16)
ylabel('State','FontWeight','bold','Fontname', 'Times New Roman','fontsize',16)
axis([0 20 -10 10])
% figure;            
% plot(t,x(:,1)-yr,'k','linewidth',1.5);                                            
% xlabel('t(s)','fontsize',16,'FontWeight','bold','Fontname', 'Times New Roman');                                                    
% ylabel('Tracking error','FontWeight','bold','Fontname', 'Times New Roman','fontsize',16); 
% u=zeros(length(t),1);
%%
 for kk=1:length(t);     
% k1=2;k2=0.8;l1=2;l2=2;K1=1;K2=1;c1=0.7;c2=0.7;d=0.01;
%  k1=0.5;k2=0.5;l1=2;l2=2;K1=0.5;K2=1;c1=0.57;c2=0.56;d=0.01;
%k1=0.5;k2=0.5;l1=3;l2=2;K1=0.1;K2=1;c1=0.57;c2=0.57;d=0.01;
w=-3-0.05*cos(t(kk));
k1=0.5;k2=0.5;k3=0.5;K1=0.1;K2=0.1;K3=0.1;c1=0.57;c2=0.57;c3=0.57;d=0.01;s11=0.01;s12=0.01;
% o=0.001;k1=90;k2=7;;d=0.0
% k1=0.1;k2=0.5;K1=0.1;K2=1;c1=0.57;c2=0.57;d=0.01;s11=0.01;s12=0.01;
lambda=1/(1+w);
%phi2=10;
delta1=0.01;varepsilon1=7;m1=2;barm1=11;
z1=x(kk,4);
alpha1=-1/2*z1*x(kk,7)-3/2*z1-k1*z1^(2*c1-1);
z2=x(kk,5)-alpha1;
alpha2=-k2*z2^(2*c2-1)-2*z2+1/2*z2*x(kk,8);
z3=x(kk,6)-alpha2;
alpha3=-k3*z3^(2*c3-1)-2*z3+1/2*z3*x(kk,9);
%k1=0.8;k2=0.8;l1=3;l2=2;K1=0.5;K2=1.5;c1=0.57;c2=0.56;d=0.01;
% c1=0.57;c2=0.56;
w1(kk)=-(1+delta1)*(alpha3*tanh(z3*alpha3/varepsilon1)+barm1*tanh(z3*barm1/varepsilon1));
if kk==1
   event1(kk)=0;
end
if kk>1 && abs(w1(kk)-u(kk-1))<delta1*abs(u(kk-1))+m1
      u(kk)=u(kk-1);
%       num1=num1+1;%触发次数计数
  else u(kk)=w1(kk);
%         rt11=[rt11;count];%计算触发时间间隔 并赋值给数组rt1
%         rtt11=[rtt11;count-tk11];
%         tk11=count;%把当前触发时刻重新记录tk
        M1=0.01*kk-M1last;
        event1(kk)=M1;
        M1last=0.01*kk;
end

% k1=2;k2=0.8;l1=1;l2=1;K1=1;K2=1;c1=0.6;c2=0.6;d=0.01;

% k1=0.8;k2=0.8;l1=2;l2=2;K1=1;K2=1;c1=0.57;c2=0.56;d=0.01;x还行

 end  
 %%
figure;
plot(t,u,'m','linewidth',1.5);
xlabel('t(s)','fontsize',16,'FontWeight','bold','Fontname', 'Times New Roman');                                                    
ylabel('Control input','FontWeight','bold','Fontname', 'Times New Roman','fontsize',16); 
%%
figure;
stem(t,event1,'b','linewidth',1.5)
 axis([0 20 0.001 6])
%  set(gca,'YTick',[0.5:0.5:2.5])
% legend('Agent1','interpreter','latex');
xlabel('Time(sec)','FontWeight','bold','Fontname', 'Times New Roman','fontsize',14);
ylabel('t_{k+1}-t_k','FontWeight','bold','Fontname', 'Times New Roman','fontsize',14);
%% 
figure;  
plot(t,x(:,1)-x(:,4),'r','linewidth',1.5);                          
xlabel('t(s)','fontsize',16);   
 str1='$$\delta_1$$';
legend({str1},'interpreter','latex','fontsize',16)
ylabel('Estimation error','FontWeight','bold','Fontname', 'Times New Roman','fontsize',16)
%% 
figure;  
plot(t,x(:,1),'k','linewidth',1.5);                            
xlabel('t(s)','fontsize',16);   
 str1='$$x_1$$';
legend({str1},'interpreter','latex','fontsize',16)
ylabel('State','FontWeight','bold','Fontname', 'Times New Roman','fontsize',16)
axis([0 20 -1 1])
%% 
figure;  
plot(t,x(:,7),'k',t,x(:,8),'r',t,x(:,9),'g','linewidth',1.5);                            
xlabel('t(s)','fontsize',16);   
 str1='$$\theta_1$$';
  str2='$$\theta_2$$';
   str3='$$\theta_3$$';
legend({str1,str2,str3},'interpreter','latex','fontsize',16)
ylabel('State','FontWeight','bold','Fontname', 'Times New Roman','fontsize',16)
%    %%
% figure;            
% plot(t,x(:,2),'k','linewidth',1.5);                                            
% xlabel('t(s)','fontsize',16);                                                    
% ylabel('Tracking error','FontWeight','bold','Fontname', 'Times New Roman','fontsize',16);          
% %% 
% figure;                                                                
% plot(t,v,'k','linewidth',1.5);                  
% xlabel('t(s)','fontsize',16);                 
% ylabel('Control input','FontWeight','bold','Fontname', 'Times New Roman','fontsize',16);  
%%    
 
%  figure;                                                                
% plot(t,x(:,2),'linewidth',1.5);                  
% xlabel('t(s)');                 
% ylabel('Control input');     
% % 
%  figure;                                                                
% plot(t,N1,'linewidth',1.5);                  
% xlabel('t(s)');                 
% ylabel('Control input'); 
% %%
%  figure;                                                                
% plot(t,N2,'linewidth',1.5);                  
% xlabel('t(s)');                 
% ylabel('Control input');     


                                                                                                                                                                                   
  
  
  
  
  