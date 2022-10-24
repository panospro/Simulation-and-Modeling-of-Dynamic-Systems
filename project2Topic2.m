clear;
clc;
close all;


%% Topic 2
%values of a,b,u
a=1.5;
b=2;
u=@(t)3*cos(2*t);

%initial conditions
x0=0;
theta0=0;

%time between [0,10] with step 0.1
tspan=0:0.1:10;

%theta values
theta1=a;
theta2=b;
thetaM=10;              

%noise values
n0=25;
f=30;
n=@(t)n0*sin(2*pi*f*t);

%solving the differential
[t, x]=ode15s(@(t,x)odefun(t,x,a,b,u),tspan,x0);

for i=1:length(tspan)
    for j=tspan
        x(i)=x(i)+n(j);
        break;
    end
end

%transfer function simulated response
f1=lsim(tf([1],[1,thetaM]), x, tspan);
f2=lsim(tf([1],[1,thetaM]), (u(tspan))', tspan);

%in order to be suitable for a simulation
A=ts2func(f1,'Times',tspan);
B=ts2func(f2,'Times',tspan);
C=ts2func(x,'Times',tspan);

%solving the differentials
[t,thetaHat1]= ode15s(@(t,thetaHat)dth(t,thetaHat,thetaM-theta1,A,C),tspan,theta0);
[t,thetaHat2]= ode15s(@(t,thetaHat)dth(t,thetaHat,theta2,B,u),tspan,theta0);

D=ts2func(thetaM-thetaHat1,'Times',tspan);
E=ts2func(thetaHat2,'Times',tspan);
[t,xHat]= ode15s(@(t,xHat)odefun(t,xHat,D,E,u,thetaM,C,'m'),tspan,x0);

%chaning the value of theta0 to see what happens
theta0=5;  
[t,thetaHat3]= ode15s(@(t,thetaHat3)dth(t,thetaHat3,thetaM-theta1,A),tspan,theta0);
[t,thetaHat4]= ode15s(@(t,thetaHat3)dth(t,thetaHat3,theta2,B,u),tspan,theta0);

F=ts2func(thetaM-thetaHat3,'Times',tspan);
[t,xHat2]= ode15s(@(t,xHat)odefun(t,xHat,F,E,u,thetaM,C,'p'),tspan,x0);


%% plots
figure(1);
plot(t,x,t,xHat);
title('Values of x');
xlabel('Time');
ylabel('x');
legend('Real x','Estimated x');

figure(2);
plot(t,thetaM-thetaHat1,t,thetaM-thetaHat3);
title('Estimated values of a');
xlabel('Time');
ylabel('a');
legend('with и0=0','with и0/=0');

figure(3);
plot(t,thetaHat2,t,thetaHat4);
title('Estimated values of b');
xlabel('Time');
ylabel('b');
legend('with и0=0','with и0/=0');

figure(4);
plot(t,x,t,xHat2);
title('Values of x ме ип/=0');
xlabel('Time');
ylabel('x');
legend('Real x','Estimated x');

figure(5);
plot(t,x,t,xHat,t,xHat2);
title('Values of x');
xlabel('Time');
ylabel('x');
legend('Real x','Estimated x ме и0=0','Estimated x ме и0/=0');

%% functions
function dx = odefun(t,x,a,b,u,thetaM,xReal,structure)
    switch nargin
        %normal odefun
        case 5
            dx = -a*x+b*u(t);
            dx=dx';
            
        %odefun with different structures
        case 8
            switch structure
                case "p"
                    dx=-a(t)*x+b(t)*u(t);
                    dx=dx';
                case "m"
                    dx=-a(t)*x+b(t)*u(t)+thetaM*(xReal(t)-x);
                    dx=dx';
            end
    end
end

function dth = dth(t,thetaHat,thetaStar,f,x)
    switch nargin
        %normal dth
        case 4
            dth=(thetaStar - thetaHat)*f(t)*thetaHat*f(t);
            dth=dth';
            
        %dth with x   
        case 5
            dth=(thetaStar - thetaHat)*f(t)*x(t);
            dth=dth';
    end     
end