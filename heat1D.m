clear all; clc;

%%% matetial
rho = 7800;
c = 500;
lambda = 46;
a = lambda / (rho * c);  % thermal conductivity

%%% parameters
L = 0.25; % length
t = 12600;

%%% bordering
T0 = 20; % initial temperature
h = 0.001; % height
tau = h^2 / (4 * a);

Nl = round(L / h +1);
Nt = round(t / tau);

Tl = 1300; % Left border temperature
Tr = 1300; % Right border temperature

alpha = 1900;

gu = input('Please, Enter the kind of boundary conditions: ');

T = ones(Nl, Nt);

for j = 1:Nl
       T(j, 1) = T0; 
end

if gu == 1
    q1 = 0;
    Tlborder = 1200;
    T(1,1) = Tlborder;
    Trborder = 1200;
    T(Nl, 1) = Trborder;
    for i = 1:Nt
        T(1, i) = Tlborder;
        T(Nl, i) = Trborder;
    end
    alpha = 0;
    
else if gu == 2
        q1 = 5000;
        Tlborder = 0;
        Trborder = 0;
        alpha = 0;
    else if gu == 3
            q1 = 0;
            Tlborder = 0;
            Trborder = 0;
        end
    end
end

 for i = 1:Nt-1
    for  j = 1:Nl
        if gu == 1
            T(1, i+1) = Tlborder;
            T(Nl, i+1) = Trborder;
        if gu == 2 || gu == 3
            if j == 1
                T(j, i+1) = T(j, i) + ((Tlborder + q1 + alpha * (Tl - T(j, i))) +((lambda * (T(j+1, i) - T(j , i)) )/ h)) * ...
               (2 * tau / (c * rho * h));
            end
            if j == Nl
                T(j, i+1) = T(j, i) + ((Trborder + q1 + alpha * (Tr - T(j, i))) +((lambda * (T(j-1, i) - T(j , i))) / h)) * ...
               (2 * tau / (c * rho * h));                                     
            end
        T(j, i+1) = T(j, i) + tau * a * (T(j-1, i) - 2 * T(j , i) + T(j+1, i)) / h^2;
        end
    end
 end

 Tend = ones(Nl, 1);
 len = 1:1:Nl;
 for i = 1:1:Nl
    Tend(i) = T(i, Nt); 
 end

plot(len, Tend)
grid on
hold on 
xlabel('x')
ylabel('Temperature')

s = Tend(120)
s2 = Tend(130)

 
