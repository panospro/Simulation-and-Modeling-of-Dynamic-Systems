% clear;
% clc;
% close all;
% %% Topic 1
% %values of m,b,k
% m=10;
% b=0.3;
% k=1.5;
% 
% w_real=[(b/m)-5; (k/m)-3; 1/m];
% 
% %initial conditions
% x0(1)=0;
% x0(2)=0;
% 
% %time between [0,10] with step 0.1
% tspan=0:0.1:10;
% 
% %input of system 
% u=@(t)10*sin(3*t)+5;
% 
% %simulation
% [t, state]=ode15s(@(t,state)odefun(t,state,m,b,k,u),tspan,x0);
% 
% %plot
% figure(1);
% plot(t,state(:,1));
% xlabel('The position in time')
% ylabel('sec')
% 
% figure(2);
% plot(t,state(:,2));
% xlabel('The velocity in time')
% ylabel('sec')
% 
% figure(3);
% plot(t,state(:,1), t,state(:,2));
% 
% figure(4);
% plot(state(:,1), state(:,2));
% xlabel('position')
% ylabel('velocity')
% title('phaze diagrams')
% 
% 
% 
% f(1)=tf([-1,0],[1,5,3])*state(2);
% 
% f(2)=tf([-1],[1,5,3])*state(2);
% 
% f(3)=tf([1],[1,5,3])*u(0.2);
% 
% 
% f=lsim(tf([-1,0],[1,5,3]), state(:,1), tspan);
% 
% f(:,2)=lsim(tf([0,-1],[1,5,3]), state(:,1), tspan);
% 
% f(:,3)= lsim(tf([0,1],[1,5,3]), (u(tspan))', tspan);
% 
% w_estimated=(state(:,1)')*f*inv(f'*f);
% 
% %estimated values of m,b,k
% m_estimated = 1/w_estimated(1,3);
% disp('m_estimated=');
% disp(m_estimated);
% 
% b_estimated = (w_estimated(1,1)+5)*m_estimated;
% disp('b_estimated=');
% disp(b_estimated);
% 
% k_estimated = (w_estimated(1,2)+3)*m_estimated;
% disp('k_estimated=');
% disp(k_estimated);
% 
% 
% function dx = odefun(t,x,m,b,k,u)
% 
% dx(1) = x(2);
% dx(2) = 1/m * (-k*x(1) - b*(x(2)) + u(t));
% 
% dx=dx';
% 
% end
% 
% 
% 
% 

% 
% pcode v.p
% 
% result=v .p
%% Topic 2

clc;
clear;
format short g;

lamda1=3;
lamda2=5;

%input of system 
u1=@(t)3*sin(2*t);
u2=@(t)2;
u1dot=@(t)6*cos(2*t);
u2dot=@(t)0;

%time between [0,3] with step 0.00001
tspan=0:0.00001:3;

%values of VR, VC
i=1;
for t=0:0.00001:3
    [VR VC]=v(t);
    VR(i,1)=Vout(2);
    VC(i,1)=Vout(1);
    u2matrix(i,1)=u2(t);
    i=i+1;
end

%{
VC(50,1)=500000;
VC(10000,1)=800000;
VC(20000,1)=100000;

VR(50,1)=500000;
VR(10000,1)=800000;
VR(20000,1)=100000;
%}
    
%matrix F
f(:,1)=lsim(tf([-1,0],[1,lamda1,lamda2]), VC(:,1), tspan);

f(:,2)=lsim(tf([0,-1],[1,lamda1,lamda2]), VC(:,1), tspan);

f(:,3)= lsim(tf([1,0],[1,lamda1,lamda2]), (u1(tspan))', tspan);
    
f(:,4)= lsim(tf([0,1],[1,lamda1,lamda2]), (u1(tspan))', tspan);
    
f(:,5)= lsim(tf([1,0],[1,lamda1,lamda2]), (u2matrix)', tspan);
    
f(:,6)= lsim(tf([1,0],[1,lamda1,lamda2]), (u2matrix)', tspan);

    
w_estimated=(VC(:,1)')*f*inv(f'*f);
disp(w_estimated);

%initial conditions
vo(1)=VC(1);
vo(2)=0;
        
[t, state]=ode15s(@(t,state)odefun2(t,state,w_estimated(5),w_estimated(6),u1dot),tspan,vo);

%real values of VC
figure(1);
plot(t,VC(:,1));
    
%estimated values of VC
figure(2);
plot(t,state(:,1));
    
%real values of VR
figure(3);
plot(t,VR(:,1));
    
%estimated values of VR
figure(4);
plot(t,u1(t)+u2(t)-state(:,1));



function [dV] = odefun2(t,V,w3,w6,u1dot)

dV(1)=V(2);

dV(2)=w3*u1dot(t) - w3*V(2) + w6-w6*V(1);

dV=dV';

end
% function [theta0] = thetaestimate(y, phi)
% % From the theory we know that the measurements of N time points should
% % to be at least equal to d, the dimension number of ö. We consider
% % that something like this is true (otherwise we would just put an if to check it)
% % and we consider that the inverse of (sum(phi*phi^T)/N) exists so that
% % get the ready formula.
%     N = length(y);
%     for i = 1:N
% 
%     theta0 = ((sum(phi*phi'))/length(y))' * ((sum(phi*y))/N);
% 
%     end
% end

