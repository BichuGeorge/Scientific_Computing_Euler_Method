% Question 1
% Part a
clear all
h=0.01;
x=0.1:h:1;
y(1)=1;
for i = 1:length(x)-1,
    dydx= x(i).^2 - y(i);
    y(i+1) = y(i)+dydx*h;
    fprintf('="Y"\n\t   %0.02f',y(i));
end
%%fprintf('="Y"\n\t   %0.01f',y);
plot(x,y);
grid on;
hold on
clear all
analytical_matrix = [];

for i=0.1:0.01:1
 y_new = (-exp((-i)) + i.^2 - 2*i + 2);
 analytical_matrix = [analytical_matrix y_new];
end

plot(0.1:0.01:1 ,  analytical_matrix,'--')
hold on

clear all
h=0.05;
x=0.1:h:1;
y(1)=1;
for i = 1:length(x)-1,
    dydx= x(i).^2 - y(i);
    y(i+1) = y(i) + dydx*h;
    fprintf('="Y"\n\t   %0.02f',y(i));
end
%%fprintf('="Y"\n\t   %0.01f',y);
plot(x,y,'*');
grid on;

% Part b

clear all
h=0.01;
x=0.1:h:1;
y(1)=1;
for i = 1:length(x)-1,
    dydx= 3*y(i) + 3*x(i);
    y(i+1) = y(i)+dydx*h;
    fprintf('="Y"\n\t   %0.02f',y(i));
end
%%fprintf('="Y"\n\t   %0.01f',y);
plot(x,y);
grid on;
hold on
clear all
analytical_matrix = [];

for i=0.1:0.01:1
 y_new = ((4/3)*exp(3*i) - i - (1/3));
 analytical_matrix = [analytical_matrix y_new];
end

plot(0.1:0.01:1 ,  analytical_matrix,'--')
hold on

clear all
h=0.05;
x=0.1:h:1;
y(1)=1;
for i = 1:length(x)-1,
    dydx= 3*y(i) + 3*x(i);
    y(i+1) = y(i) + dydx*h;
    fprintf('="Y"\n\t   %0.02f',y(i));
end
%%fprintf('="Y"\n\t   %0.01f',y);
plot(x,y,'*');
grid on;

% % Part c

clear all
h=0.01;
x=0.1:h:1;
y(1)=1;
for i = 1:length(x)-1,
    dydx= -x(i)*y(i);
    y(i+1) = y(i)+dydx*h;
    fprintf('="Y"\n\t   %0.02f',y(i));
end
%%fprintf('="Y"\n\t   %0.01f',y);
plot(x,y);
grid on;
hold on
clear all
analytical_matrix = [];

for i=0.1:0.01:1
 y_new = (exp(-( (i^2)/2) ));
 analytical_matrix = [analytical_matrix y_new];
end

plot(0.1:0.01:1 ,  analytical_matrix,'--')
hold on

clear all
h=0.05;
x=0.1:h:1;
y(1)=1;
for i = 1:length(x)-1,
    dydx= -x(i)*y(i);
    y(i+1) = y(i) + dydx*h;
    fprintf('="Y"\n\t   %0.02f',y(i));
end
%%fprintf('="Y"\n\t   %0.01f',y);
plot(x,y,'*');
grid on;

% % Part d
% 
clear all
h=0.01;
x=0.1:h:1;
y(1)=1;
for i = 1:length(x)-1,
    dydx= 2*x(i)*y(i)*y(i);
    y(i+1) = y(i)+dydx*h;
    fprintf('="Y"\n\t   %0.02f',y(i));
end
%%fprintf('="Y"\n\t   %0.01f',y);
plot(x,y);
grid on;
hold on
clear all
analytical_matrix = [];

for i=0.1:0.01:1
 y_new = (1 / (1-i^2) );
 analytical_matrix = [analytical_matrix y_new];
end

plot(0.1:0.01:1 ,  analytical_matrix,'--')
hold on

clear all
h=0.05;
x=0.1:h:1;
y(1)=1;
for i = 1:length(x)-1,
    dydx= 2*x(i)*y(i)*y(i);
    y(i+1) = y(i) + dydx*h;
    fprintf('="Y"\n\t   %0.02f',y(i));
end
%%fprintf('="Y"\n\t   %0.01f',y);
plot(x,y,'*');
grid on;

% Part e

clear all
h=0.05;
x=0.1:h:1;
y(1)=1/10;
for i = 1:length(x)-1
    dydx= exp(-(2*x(i))) - 2*y(i);
    y(i+1) = y(i) + dydx*h;
    fprintf('="Y"\n\t   %0.02f',y(i));
end
%%fprintf('="Y"\n\t   %0.01f',y);
plot(x,y);
grid on;
hold on
clear all
analytical_matrix = [];

for i=0.1:0.05:1
 y_new = (1/10)*exp((-2)*i) + i*exp((-2)*i);
 analytical_matrix = [analytical_matrix y_new];
end

plot(0.1:0.05:1 ,  analytical_matrix,'--')
hold on

clear all
h=0.01;
x=0.1:h:1;
y(1)=1/10;
for i = 1:length(x)-1,
    dydx= exp(-(2*x(i))) - 2*y(i);
    y(i+1) = y(i) + dydx*h;
    fprintf('="Y"\n\t   %0.02f',y(i));
end
%%fprintf('="Y"\n\t   %0.01f',y);
plot(x,y,'*');
grid on;

% Question 2

% mx'' - mg + cx'^2 = 0
% x' = v , x'' = v'

% x' = v
% mv' +cv^2 = mg

h = 0.0001;
t = 0:h:4;
m = 3.4 * 10^-5;
g = 9.81;
c= 0.5;
x0=0;
v0=0;
x=zeros(length(t),1);
v=zeros(length(t),1);
x(1)=x0;
v(1)=v0;

for i=1:length(t)-1
    x(i+1) = x(i) + h*v(i);
    v(i+1) = v(i) + h*( (m*g - c*(v(i)^2))/m ) ;
    fprintf('x(%0.2f)=%0.6f\n',t(i+1), x(i+1))
    fprintf('v(%0.2f)=%0.6f\n',t(i+1), v(i+1))
end
plot(t,x)
hold on
plot(t,v,'--')

% Question 3
% y'' + (g/L) y = 0
% y(0) = pi/2
%  y'(0) = 0
% g/L = 1

% y' = v
% v' = -y

h = 0.01;
t = 0:h:500;
y0=pi/2;
v0=0;
y=zeros(length(t),1);
v=zeros(length(t),1);
y(1)=y0;
v(1)=v0;

for i=1:length(t)-1
    y(i+1) = y(i) + h*v(i);
    v(i+1) = v(i) + h*(-y(i));
    fprintf('x(%0.2f)=%0.6f\n',t(i+1), y(i+1))
    fprintf('v(%0.2f)=%0.6f\n',t(i+1), v(i+1))
end
plot(t,y)
hold on
plot(t,v,'--')
