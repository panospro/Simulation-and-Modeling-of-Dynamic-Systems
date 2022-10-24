clear;
clc;
close all;


%% Topic 1

%% a
%values of a,b,u
a=1.5;
b=2;
u=@(t)3;


%initial conditions
x0=0;

%time between [0,10] with step 0.1
tspan=0:0.1:10;

%simulation
[t, x]=ode15s(@(t,x)odefun(t,x,a,b,u),tspan,x0);

%theta values
theta1=a;
theta2=b;
thetaM=10;

for i=1:length(tspan)
    U(i)=u(t);
end

%transfer function simulated response
f1=lsim(tf([1],[1,thetaM]), x, tspan);
f2=lsim(tf([1],[1,thetaM]), U', tspan);

%in order to be suitable for a simulation
A=ts2func(f1,'Times',tspan);
B=ts2func(f2,'Times',tspan);

theta0=0;
[t,theta1Hat]= ode15s(@(t,thetaHat)dth(t,thetaHat,thetaM-theta1,A),tspan,theta0);
[t,theta2Hat]= ode15s(@(t,thetaHat)dth(t,thetaHat,theta2,B),tspan,theta0);



%% plots1
figure(1);
plot(t,x);
title('Real values of x');
xlabel('Time');
ylabel('x');
hold on;

figure(2);
plot(t,thetaM-theta1Hat);
title('Estimated values of a');
xlabel('Time');
ylabel('a');
hold on;

figure(3);
plot(t,theta2Hat);
title('Estimated values of b');
xlabel('Time');
ylabel('b');
hold on;

thetaHat=[theta1Hat theta2Hat];
f=[f1 f2];

for i=1:length(tspan)
    xHat(i)=thetaHat(i,:)*f(i,:)';
end

figure(4);
plot(t,xHat);
title('Estimated values of x');
xlabel('Time');
ylabel('x');
hold on;



%% b

u=@(t)3*cos(2*t);

%simulation
[t, x]=ode15s(@(t,x)odefun(t,x,a,b,u),tspan,x0);

theta1=a;
theta2=b;
thetaM=5;

f1=lsim(tf([1],[1,thetaM]), x, tspan);
f2=lsim(tf([1],[1,thetaM]), (u(tspan))', tspan);

%in order to be suitable for a simulation
A=ts2func(f1,'Times',tspan);
B=ts2func(f2,'Times',tspan);

theta0=0;
[t,theta1Hat]= ode15s(@(t,thetaHat)dth(t,thetaHat,thetaM-theta1,A),tspan,theta0);
[t,theta2Hat]= ode15s(@(t,thetaHat)dth(t,thetaHat,theta2,B),tspan,theta0);


%% plots2
figure(1);
plot(t,x);
title('Real values of x');
xlabel('Time');
ylabel('x');
legend('u=3','u=3cos3t')

figure(2);
plot(t,thetaM-theta1Hat);
title('Estimated values of a');
xlabel('Time');
ylabel('a');
legend('u=3','u=3cos3t')

figure(3);
plot(t,theta2Hat);
title('Estimated values of b');
xlabel('Time');
ylabel('b');
legend('u=3','u=3cos3t')

thetaHat=[theta1Hat theta2Hat];
f=[f1 f2];

for i=1:length(tspan)
    xHat(i)=thetaHat(i,:)*f(i,:)';
end

figure(4);
plot(t,xHat);
title('Estimated values of x');
xlabel('Time');
ylabel('x');
legend('u=3','u=3cos3t')


%% functions
function dx = odefun(t,x,a,b,u)
    dx = -a*x+b*u(t);
    dx=dx';
end

function dth = dth(t,thetaHat,thetaStar,f)
    gamma=100;
    dth=-(thetaHat - thetaStar)*f(t)*f(t)*gamma;
    dth=dth';
end