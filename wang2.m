function dx=wang2(t,x);
global u
k1=0.5;k2=0.1;k3=0.5;K1=0.5;K2=0.5;K3=0.1;c1=0.57;c2=0.57;c3=0.57;d=0.1;s11=0.1;s12=0.1;s13=0.1;
% o=0.001;k1=90;k2=7;;d=0.0
% k1=0.1;k2=0.5;K1=0.1;K2=1;c1=0.57;c2=0.57;d=0.01;s11=0.01;s12=0.01;
w=1-0.35*cos(t);
lambda=1/(1+w);
%phi2=10;
delta1=0.01;varepsilon1=7;m1=2;barm1=11;
z1=x(4);
alpha1=-1/2*z1*x(7)-3/2*z1-k1*z1^(2*c1-1);
z2=x(5)-alpha1;
alpha2=-k2*z2^(2*c2-1)-2*z2+1/2*z2*x(8);
z3=x(6)-alpha2;
alpha3=-k3*z3^(2*c3-1)-2*z3+1/2*z3*x(9);
w1=-(1+delta1)*(alpha3*tanh(z3*alpha3/varepsilon1)+barm1*tanh(z3*barm1/varepsilon1));
if t==0;
   u=w1;
end
 if t>0 && abs(w1-u)<delta1*abs(u)+m1;
      u=u
  else u=w1;
 end
%%  %% ********************** 模糊 1*****************
   sum1=0; %fuzzy logic system for first subsystems                       
   for i=1:5
   uf1(i)=exp(-0.5*(x(13)-3+i)^2);  %basis function vector uf1
   sum1=sum1+uf1(i); %phi1隶属度求和
   end     
   for i=1:5
    phi11(i)=uf1(i)/sum1;%行向量
   end 
     sum2=0; %fuzzy logic system for first subsystems                       
   for i=1:5
   uf2(i)=exp(-0.5*(x(13)-3+i)^2)*exp(-0.5*(x(14)-3+i)^2);  %basis function vector uf1
   sum2=sum2+uf2(i); %phi1隶属度求和
   end     
   for i=1:5
    phi12(i)=uf2(i)/sum2;
   end 
   sum3=0; %fuzzy logic system for first subsystems                       
   for i=1:5
   uf3(i)=exp(-0.5*(x(13)-3+i)^2)*exp(-0.5*(x(14)-3+i)^2);  %basis function vector uf1
   sum2=sum2+uf2(i); %phi1隶属度求和
   end     
   for i=1:5
    phi13(i)=uf3(i)/sum2;
   end 

 %%
%  for iii=1:125;
%       s1(iii)=exp(-((x(1)-1+0.08*iii)^2)/4);   %RBF设计
% end           
% for iii=1:125;
%       s2(iii)=exp(-((x(1)-1+0.08*iii)^2+(x(2)-1+0.08*iii)^2)/4);   %RBF设计
% end

% if v>=1  
%     u=(1-0.2*sin(v))*(v-1);
% else  if v<=-0.5  
%  u=(0.8-0.1*cos(v))*(v+0.5);     
%     else  
%      u=0;   
%     end 
% end  
 
  dx=[x(2)+x(1)^2;
      x(3)+x(1)^2*x(2);
      u+x(1)^2*x(2)* x(3);
      x(5)+K1*(x(1)/lambda-x(4))+[x(10) x(11) x(12) x(13) x(14)]*phi11';
      x(6)+K2*(x(1)/lambda-x(4))+[x(15) x(16) x(17) x(18) x(19)]*phi12';
      u+K3*(x(1)/lambda-x(3))+[x(20) x(21) x(22) x(23) x(24)]*phi13';
      1/2*z1^2-d*x(7);
      1/2*z2^2-d*x(8);
      1/2*z2^2-d*x(9);
        -z1*phi11(1)-s11*x(10);
        -z1*phi11(2)-s11*x(11);
         -z1*phi11(3)-s11*x(12);
         -z1*phi11(4)-s11*x(13);
        -z1*phi11(5)-s11*x(14);
         -z2*phi12(1)-s12*x(15);
         -z2*phi12(2)-s12*x(16);
         -z2*phi12(3)-s12*x(17);
        -z2*phi12(4)-s12*x(18);
         -z2*phi12(5)-s12*x(19);
         -z3*phi13(1)-s13*x(20);
         -z3*phi13(2)-s13*x(21);
         -z3*phi13(3)-s13*x(22);
        -z3*phi13(4)-s13*x(23);
         -z3*phi13(5)-s13*x(24);
      
%       z1t*dfai1;
%       z2t*dfai2;
%       z2t*dfai3;
%      1/2*z2t^2*phi2^2-o*x(6);
%      1/2*z2t^2*alpha2^2-o*x(7)
%      1/2*z2t^2*(exp(m)*sin(pi*x(3))*k1*lambdam1^2*x(2)/lambda)^2-o*x(8);
%      1/2*z2t^2*alpha1^2-o*x(9);
  
%   -x(2)+x(1)*exp(-0.5*x(2))+0.1*x(1)^2+sin(x(1)^2)*x(3);                                                                                                          
%        u+0.2*exp(-x(2))+x(1)*sin(x(2))+x(1)^2*x(3);     
%        -x(3)+x(1)^2+0.5;
%         p1*a1^2*e1^2*s1*s1'-gamma1*x(4);
%         p2*a2^2*e2^2*s2*s2'-gamma2*x(5);
%        e1^2+a1^2*e1^2*x(4)*s1*s1';
%         e2^2+a2^2*e2^2*x(5)*s2*s2'
];   



