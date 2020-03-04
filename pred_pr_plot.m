clear all
close all

h = 0.01; %step size
t = 0:h:50; %time
n = length(t);
x1 = zeros(1, n); %prey
x2 = zeros(1, n); %predator
a1 = zeros(1, n); %food  
x_w_st1 = 10; %x with star
%x1 -> x_w_st1
%T1 = 1;
%T2 = 1;
x1(1) = 1;
x2(1) = 1;
a1(1) = 1; %a1

A = load('param3.txt');
a2 = A(:,1);
b1 = A(:,2);
b2 = A(:,3);
T1 = A(:,4);
T2 = A(:,5);

for i = 1:n - h
    f1 = a1(i) .* x1(i) - b1(i) * x1(i) * x2(i);
    f2 = -a2(i) .* x2(i) + b2(i) * x1(i) * x2(i);
    f3 = 0;
    phi = -((x1(i) - x_w_st1) / (T2(i) * x1(i))) + b1(i) * x2(i); %phi(x1, x2)
    dphi = - ((x_w_st1) / (T2(i) * x1(i)^2)) * f1 + b1(i) * f2; % dphi(x1,x2)/dt
    psi = a1(i) - phi; %psi^(I)
    U(i) = -(psi/T1(i)) + dphi; 
    x1(i+1) = x1(i) + h*f1;
    x2(i+1) = x2(i) + h*f2;
    a1(i+1) = a1(i) + h*(f3 + U(i)); 
end %for euler


plot(t,x1, t,x2,  t, a1); hold on
hold off
legend('predator', 'pray', 'food')
xlabel('time');
ylabel('population');

%smth wrong
%my brain i think :\
f2 = figure;
figure(f2);
plot(x_w_st1, x1);
hold on;
hold off;
