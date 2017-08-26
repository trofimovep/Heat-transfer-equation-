clear all; clc;

%%% matetial
rho = 7800;
c = 500;
lambda = 46;

a = lambda / (rho * c);  %температуропроводность


%%% parameters
lentgh = 0.25;
%L = lentgh / 2;
L = 0.25;
% n=fix(N/2);
% T=[repmat(400,1,n),repmat(300,1,N-n)];% вектор начальных температур
% %Tr=400*(1+0.1*(cos(w*(0:tau:t_end))));% вектор температур правой границы
% Tr=1900*((T - 20)/h);% вектор температур правой границы
% x=0:h:L;% координаты узлов сетки
% % вспомогательные коэффициенты:
% k=ro*c/tau;
% A=lamda/h^2;
% B=2*A+k;
% C=A;
% f=(0.1+2*x(2:end))*ro*c;
% %-----------------------
% alfa=zeros(1,N);% подготовка массива для коэффициентов alfa
% bet=repmat(Tl,1,N);% подготовка массива для коэффициентов beta
% for I=2:length(Tr)
%     % Вычисление коэффициентов alfa и beta:
%     for J=2:N-1
%         s=B-C*alfa(J-1);
%         alfa(J)=A/s;
%         bet(J)=(T(J)*k+C*bet(J-1)+f(J))/s;
%     end
%     % Определение неизвестного поля температуры:
%     T(N)=Tr(I);
%     for J=N-1:-1:1
%         T(J)=alfa(J)*T(J+1)+bet(J);
%     end
% end
% plot(x,T),grid on


t = 12600;

% %% bordering
T0 = 20;

h = 0.001;
tau = h^2 / (4 * a);

ll = round( lentgh / h + 1);

Nl = round(L / h +1);
Nt = round(t / tau);

Tl = 1300;
Tr = 1300;

alpha = 1900;

gu = input('Введите род граничных условий: ');




T = ones(Nl, Nt);
q = ones(Nl, Nt);

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
         if j == Nl %%%%  os
                           
             
           T(j, i+1) = T(j, i) + ((Trborder + q1 + alpha * (Tr - T(j, i))) +((lambda * (T(j-1, i) - T(j , i))) / h)) * ...
               (2 * tau / (c * rho * h));                                     
            end                                             
          
                   
        T(j, i+1) = T(j, i) + tau * a * (T(j-1, i) - 2 * T(j , i) + T(j+1, i)) / h^2;
                 
            end
        end
    end
 
















%  for i = 1:1:Nt-1
%      
%     for  j = 1:1:Nl 
%      
%         if j == 1
%             
%            T(j, i+1) = T(j, i) + (alpha * (Tl - T(j, i)) +( lambda * (T(j+1, i) - T(j , i)) / h)) * ...
%                (2 * tau / (c * rho * h)); 
%            
%         else if j == Nl %%%%  os
%                                 
%                      if mod(ll,2) == 0
%                                                      
%                                                              
%                 T(j, i+1) =  T(j, i) + tau * a * (T(j-1, i) - 2 * T(j , i) + T(j, i)) / h^2;
%                      
%                      else 
%               T(j, i+1) =  T(j, i) + tau * a * (T(j-1, i) - 2 * T(j , i) + T(j-1, i)) / h^2;
%                          
%                      end 
%             else
%         
%         T(j, i+1) = T(j, i) + tau * a * (T(j-1, i) - 2 * T(j , i) + T(j+1, i)) / h^2;
%                  
%         
%                 
%             end
%         end
%     end
%  end

%  v = ones(1, Nt);
%  k = 1:1:Nt;
% for i =1:Nt
%  v(i) = T(1,i);
% end
%  plot(v,k)
 

 Tend = ones(Nl, 1);
 
 len = 1:1:Nl;
 
 for i = 1:1:Nl
    Tend(i) = T(i, Nt); 
 end
 
 
 
plot(len, Tend)
grid on
hold on 
xlabel('x')
ylabel('Температура')



s = Tend(120)
s2 = Tend(130)

 
