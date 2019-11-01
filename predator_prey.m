clear all


%const
a1 = 1.0; %alpha 1
a2 = 1.0; %alpha 2
b1 = 0.5; %beta 1
b2 = 0.5; %beta 2
h = 0.1; %step size
t = 0:h:50; %time
n = length(t);
x1 = zeros(1, n); %prey
x2 = zeros(1, n); %predator
x3 = zeros(1, n); %food
x_w_st1 = zeros(1,n); %x with star
T1 = 1;
T2 = 1;
x1(1) = 1;
x2(1) = 1;
x3(1) = 1;
x_w_st1(1) = 1;


%euler's method
for i = 1:length(t) - h
    f1 = a1 * x1(i) - b1 * x1(i) * x2(i);
    f2 = -a2 * x2(i) + b2 * x1(i) * x2(i);
    f3 = 0;
    phi = -((x1(i) - x_w_st1(i)) / (T2 * x1(i))) + b1 * x2(i); %phi(x1, x2)
    dphi = - ((x_w_st1(i)) / (T2 * x1(i)^2)) * f1 + b1 * f2; % dphi(x1,x2)/dt
    psi = a1 - phi; %psi^(I)
    U(i) = -(psi/T1) + dphi; 
    x1(i+1) = x1(i) + h*f1;
    x2(i+1) = x2(i) + h*f2;
    x3(i+1) = x3(i) + h*(f3 + U(i));
end


%plot
plot(t,x1,'b.', t,x2, 'r.', t,x3, 'm.'); hold on
hold off
legend('predator', 'pray', 'food')
xlabel('time');
ylabel('population');


