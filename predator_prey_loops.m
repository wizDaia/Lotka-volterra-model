clear all
close all

fid = fopen('param3.txt','w');
h = 0.1; %step size
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
eps = 0.1 * x_w_st1;

for a2 = 0:0.1:1.5
    for b1 = 0:0.1:1.5
        for b2 = 0:0.1:1.5
            for T1 = 0:0.1:1
                for T2 = 0:0.1:1
                    %euler's method
                    for i = 1:n - h
                        f1 = a1(i) * x1(i) - b1 * x1(i) * x2(i);
                        f2 = -a2 * x2(i) + b2 * x1(i) * x2(i);
                        f3 = 0;
                        phi = -((x1(i) - x_w_st1) / (T2 * x1(i))) + b1 * x2(i); %phi(x1, x2)
                        dphi = - ((x_w_st1) / (T2 * x1(i)^2)) * f1 + b1 * f2; % dphi(x1,x2)/dt
                        psi = a1(i) - phi; %psi^(I)
                        U(i) = -(psi/T1) + dphi; 
                        x1(i+1) = x1(i) + h*f1;
                        x2(i+1) = x2(i) + h*f2;
                        a1(i+1) = a1(i) + h*(f3 + U(i));    
                    end %for euler
                    if ((x1(end) - x1(1) <= 10^3)) 
                        if((mean(x1(end - 9:end))) < eps)                   
							fprintf(fid, '%2.5f\t', a2);
							fprintf(fid, '%2.5f\t', b1);
                            fprintf(fid, '%2.5f\t', b2);
							fprintf(fid, '%2.5f\t', T1);
							fprintf(fid, '%2.5f\n', T2);
                        end   
                    end %if
                end %T2
            end %T1
        end %b2
    end %b1
end %a2


fclose(fid);
run('pred_pr_plot.m');


            