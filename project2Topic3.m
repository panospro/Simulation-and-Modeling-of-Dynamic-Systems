clear;
clc;
close all;


%% Topic 3
%values of aij,bi,u
a11= -0.5;
a12= -3;
a21=  4;
a22= -2;
b1=1;
b2=1.4;
u=@(t)(7.5*cos(3*t)+10*cos(2*t));

%initial conditions
x0=[0;0];
theta0=[0 0];

%time between [0,10] with step 0.1
tspan=0:0.1:10;


%matrices A,B  
AA=[a11 a12;a21 a22];
BB=[b1; b2];

%values of a,b
a=-AA;
b=BB;

%theta values
theta1=a;
theta2=b;
thetaM=5;

%solving the differential
[t, x]=ode15s(@(t,x)odefun(t,x,a,b,u),tspan,x0);

%transfer function simulated response
% x = repmat(x, 50,1);
f1=lsim(tf([1],[1,thetaM]), x(:,1), tspan);     
f2=lsim(tf([1],[1,thetaM]), (u(tspan))', tspan);

%in order to be suitable for a simulation
A=ts2func(f1,'Times',tspan);
B=ts2func(f2,'Times',tspan);
C=ts2func(x,'Times',tspan);

%solving the differentials
[t,thetaHat1]= ode15s(@(t,thetaHat)dth(t,thetaHat,thetaM-theta1,A,C),tspan,theta0);
[t,thetaHat2]= ode15s(@(t,thetaHat)dth(t,thetaHat,theta2,B,u),tspan,theta0);



%chaning the value of theta0 to see what happens
theta0=[5;5];  
[t,thetaHat3]= ode15s(@(t,thetaHat3)dth(t,thetaHat3,thetaM-theta1,A),tspan,theta0);
[t,thetaHat4]= ode15s(@(t,thetaHat3)dth(t,thetaHat3,theta2,B,u),tspan,theta0);

E=ts2func(thetaHat2,'Times',tspan);
F=ts2func(thetaM-thetaHat3,'Times',tspan);
[t,xHat]= ode15s(@(t,xHat)odefun(t,xHat,F,E,u,thetaM,C),tspan,x0);


%% plots
figure(1);
plot(t,x,t,xHat);
title('Values of x');
xlabel('Time');
ylabel('x');
legend('Real x','Estimated x');

figure(2);
plot(t,thetaM-thetaHat1(:,1),t,thetaM-thetaHat3(:,1));
title('Estimated values of a');
xlabel('Time');
ylabel('a');
legend('with и0=0','with и0/=0');

figure(3);
plot(t,thetaM-thetaHat1(:,2),t,thetaM-thetaHat3(:,2));
title('Estimated values of a');
xlabel('Time');
ylabel('a');
legend('with и0=0','with и0/=0');

figure(4);
plot(t,thetaHat2(:,1),t,thetaHat4(:,1));
title('Estimated values of b1');
xlabel('Time');
ylabel('b');
legend('with и0=0','with и0/=0');

figure(5);
plot(t,thetaHat2(:,2),t,thetaHat4(:,2));
title('Estimated values of b2');
xlabel('Time');
ylabel('b');
legend('with и0=0','with и0/=0');

figure(6);
plot(t,x(:,1),t,xHat(:,1));
title('Values of x with ип/=0');
xlabel('Time');
ylabel('x');
legend('Real x1','Estimated x1');

figure(7);
plot(t,x(:,2),t,xHat(:,2));
title('Values of x with ип/=0');
xlabel('Time');
ylabel('x');
legend('Real x2','Estimated x2');

figure(8);
plot(t,x);
title('Values of x');
xlabel('Time');
ylabel('x');
legend('Real x');


%% functions
function dx = odefun(t,x,a,b,u,thetaM,xReal)
    switch nargin
        %normal odefun
        case 5
            dx = -a*x+b*u(t);
            dx(:)=dx';
            
        %odefun with different structures
        case 7
            dx=-a(t)*x(1,:)+b(t)*u(t);
            dx(:)=dx';
    end
end

function dth = dth(t,thetaHat,thetaStar,f,x)
    switch nargin
        %normal dth
        case 4
            dth=(thetaStar - thetaHat)*f(t)*thetaHat*f(t);
            dth(:)=dth';
            
        %dth with x   
        case 5
            dth=(thetaStar - thetaHat)*f(t)*x(t);
            dth(:)=dth';
    end     
end