clc;
%%%%%%%%%%%%%%%%       2D      %%%%%%%%%%%%%%%%%%%%%%
% Matetial
rho = 7800;
c = 500;
lambda = 46;
a = lambda / (rho * c);       % thermal diffusitivy
% Parameters 
A = 0.25;                     %  height
B = 1.5;                      %  width
t = 600;                      %  time

%  Bordering
T0 = 20;
hx = 0.01;
hy = 0.01;
tau = 1 / (a * 12 * (1 / (hx * hx) + 1 / (hy * hy))) ;

Nx = round(B / (hx));
Ny = round(A / (hy));
Nt = round(t / tau);

condition = input('Please, Enter the knid of border conditions: ');

tic
for i = 1 : Nx
   for j = 1 : Ny
       T(i, j, 1) = T0; 
   end
end

%%%%%%%%%%  FIRST CONDITIONS 

if(condition == 1)

Tl = 1300; 
Tr = 1300;
Tu = 1300;
Td = 1300;

% Tl = input('Enter the LEFT temperature: ');
% Tr = input('Enter the RIGHT temperature: ');
% Tu = input('Enter the UP temperature: ');
% Td = input('Enter the DOWN temperature: ');
    
for k = 1 : Nt - 1
    for i = 1 : Nx
        for j = 1 : Ny
            
                  if i == 1 && j == 1
                T(i, j, k + 1) = T(i, j, k) + 2 * a * tau * (((T(i + 1, j, k) + T(i, j + 1, k) - 2 * T(i, j, k)) / hx^2) + ...
                      (Tl + Td ) / (lambda * hx));
                               
             else if i == 1 && j == Ny
                T(i, j, k + 1) = T(i, j, k) + 2 * a * tau * (((T(i + 1, j, k) + T(i, j - 1, k) - 2 * T(i, j, k)) / hx^2) + ...
                      (Tl + Tu) / (lambda * hx));
                 
               
             else if i == Nx && j == 1
                 T(i, j, k + 1) = T(i, j, k) + 2 * a * tau * (((T(i - 1, j, k) + T(i, j + 1, k) - 2 * T(i, j, k)) / hx^2) + ...
                      (Tr + Td) / (lambda * hx));
                  
                  
             else if i == Nx && j == Ny
                 T(i, j, k + 1) = T(i, j, k) + 2 * a * tau * (((T(i - 1, j, k) + T(i, j - 1, k) - 2 * T(i, j, k)) / hx^2) + ...
                      (Tr + Tu) / (lambda * hx));
                       
                     
             else if i == 1 
                % Left border 
                T(i, j, k + 1) = T(i, j, k) + a * tau * (T(i, j + 1, k)- 2 * T(i, j, k) + T(i, j - 1, k)) + ...
               (Tl + (lambda * (T(i + 1, j, k) - T(i, j , k)) / hx)) * (2 * tau / (c * rho * hx)); 
            
            
            else if i == Nx
                % Right border
                T(i, j, k + 1) = T(i, j, k) + a * tau * (T(i, j + 1, k)- 2 * T(i, j, k) + T(i, j - 1, k)) + ...
               (Tr + (lambda * (T(i - 1, j, k) - T(i, j , k)) / hx)) * (2 * tau / (c * rho * hx));
               
        
            else if j == 1
                        % Lower border
                T(i, j, k + 1) = T(i, j, k) + a * tau * (T(i + 1, j, k)- 2 * T(i, j, k) + T(i - 1, j, k)) + ...
               (Td + (lambda * (T(i, j + 1, k) - T(i, j , k)) / hx)) * (2 * tau / (c * rho * hx));
                  
           
            else if j == Ny
                            % Top border
                T(i, j, k + 1) = T(i, j, k) + a * tau * (T(i + 1, j, k)- 2 * T(i, j, k) + T(i - 1, j, k)) + ...
               (Tu  + (lambda * (T(i, j - 1, k) - T(i, j , k)) / hx)) * (2 * tau / (c * rho * hx)); 
                
             else 
                            % internal points
                     T(i , j, k + 1) = T(i, j, k) + a * tau * ((T(i + 1, j, k) - 2 * T(i, j, k) + T(i - 1, j, k)) / hx^2 + ...
                    (T(i, j + 1, k) - 2 * T(i, j, k) + T(i, j - 1, k)) / hy^2);
                        
                end                          
                end
                end
                end 
                end
                end
                end
                end
        end
    end
end                      
toc
    
%%%%%%%%%%  SECOND CONDITIONS 

else if(condition == 2)
ql = 1300; 
qr = 1300;
qu = 1300;
qd = 1300;

% ql = input('Enter the LEFT heat flow: ');
% qr = input('Enter the RIGHT flow: ');
% qu = input('Enter the UP flow: ');
% qd = input('Enter the DOWN flow: ');
    
  for k = 1 : Nt - 1
    for i = 1 : Nx
        for j = 1 : Ny
            
             if i == 1 && j == 1
                T(i, j, k + 1) = T(i, j, k) + 2 * a * tau * (((T(i + 1, j, k) + T(i, j + 1, k) - 2 * T(i, j, k)) / hx^2) + ...
                      (ql + qd ) / (lambda * hx));
                               
             else if i == 1 && j == Ny
                T(i, j, k + 1) = T(i, j, k) + 2 * a * tau * (((T(i + 1, j, k) + T(i, j - 1, k) - 2 * T(i, j, k)) / hx^2) + ...
                      (ql + qu) / (lambda * hx));
                 
               
                 else if i == Nx && j == 1
                 T(i, j, k + 1) = T(i, j, k) + 2 * a * tau * (((T(i - 1, j, k) + T(i, j + 1, k) - 2 * T(i, j, k)) / hx^2) + ...
                      (qr + qd) / (lambda * hx));
                  
                  
                     else if i == Nx && j == Ny
                 T(i, j, k + 1) = T(i, j, k) + 2 * a * tau * (((T(i - 1, j, k) + T(i, j - 1, k) - 2 * T(i, j, k)) / hx^2) + ...
                      (qr + qu) / (lambda * hx));
                       
                     
                    else if i == 1 
                % Left border 
                T(i, j, k + 1) = T(i, j, k) + a * tau * (T(i, j + 1, k)- 2 * T(i, j, k) + T(i, j - 1, k)) + ...
               (ql + (lambda * (T(i + 1, j, k) - T(i, j , k)) / hx)) * (2 * tau / (c * rho * hx)); 
            
            
            else if i == Nx
                % Right border
                T(i, j, k + 1) = T(i, j, k) + a * tau * (T(i, j + 1, k)- 2 * T(i, j, k) + T(i, j - 1, k)) + ...
               (qr + (lambda * (T(i - 1, j, k) - T(i, j , k)) / hx)) * (2 * tau / (c * rho * hx));
               
        
                else if j == 1
                        % Lower border
                T(i, j, k + 1) = T(i, j, k) + a * tau * (T(i + 1, j, k)- 2 * T(i, j, k) + T(i - 1, j, k)) + ...
               (qd + (lambda * (T(i, j + 1, k) - T(i, j , k)) / hx)) * (2 * tau / (c * rho * hx));
                  
           
                    else if j == Ny
                            % Top border
                T(i, j, k + 1) = T(i, j, k) + a * tau * (T(i + 1, j, k)- 2 * T(i, j, k) + T(i - 1, j, k)) + ...
               (qu  + (lambda * (T(i, j - 1, k) - T(i, j , k)) / hx)) * (2 * tau / (c * rho * hx)); 
                     
                        else 
                            % internal points
                     T(i , j, k + 1) = T(i, j, k) + a * tau * ((T(i + 1, j, k) - 2 * T(i, j, k) + T(i - 1, j, k)) / hx^2 + ...
                    (T(i, j + 1, k) - 2 * T(i, j, k) + T(i, j - 1, k)) / hy^2);
                        
                                          end
                                      end
                                  end
                              end 
                          end
                      end
                  end
              end
          end
     end
end                      
toc
    

%%%%%%%%%%  THIRD CONDITIONS 

else if(condition == 3)
            
Tl = 1300; 
Tr = 1300;
Tu = 1300;
Td = 1300;

alpha_l = 2000;
alpha_r = 2000;
alpha_u = 2000;
alpha_d = 2000;

% Tl = input('Enter the LEFT temperature: ');
% Tr = input('Enter the RIGHT temperature: ');
% Tu = input('Enter the UP temperature: ');
% Td = input('Enter the DOWN temperature: '); 

% alpha_l = input('Enter the LEFT alpha');
% alpha_r = input('Enter the RIGHT alpha');
% alpha_u = input('Enter the UP alpha');
% alpha_d = input('Enter the DOWN alpha');
            
    
for k = 1 : Nt - 1
    for i = 1 : Nx
        for j = 1 : Ny
            
             if i == 1 && j == 1
                T(i, j, k + 1) = T(i, j, k) + 2 * a * tau * (((T(i + 1, j, k) + T(i, j + 1, k) - 2 * T(i, j, k)) / hx^2) + ...
                      (alpha_l * (Tl - T(i, j, k)) + (alpha_d * (Td - T(i, j, k)))) / (lambda * hx));
                               
             else if i == 1 && j == Ny
                T(i, j, k + 1) = T(i, j, k) + 2 * a * tau * (((T(i + 1, j, k) + T(i, j - 1, k) - 2 * T(i, j, k)) / hx^2) + ...
                      (alpha_l * (Tl - T(i, j, k)) + (alpha_u * (Tu - T(i, j, k)))) / (lambda * hx));
                 
               
                 else if i == Nx && j == 1
                 T(i, j, k + 1) = T(i, j, k) + 2 * a * tau * (((T(i - 1, j, k) + T(i, j + 1, k) - 2 * T(i, j, k)) / hx^2) + ...
                      (alpha_r * (Tr - T(i, j, k)) + (alpha_d * (Td - T(i, j, k)))) / (lambda * hx));
                  
                  
                     else if i == Nx && j == Ny
                 T(i, j, k + 1) = T(i, j, k) + 2 * a * tau * (((T(i - 1, j, k) + T(i, j - 1, k) - 2 * T(i, j, k)) / hx^2) + ...
                      (alpha_r * (Tr - T(i, j, k)) + (alpha_u * (Tu - T(i, j, k)))) / (lambda * hx));
                       
                     
                    else if i == 1 
                % Left border 
                T(i, j, k + 1) = T(i, j, k) + a * tau * (T(i, j + 1, k)- 2 * T(i, j, k) + T(i, j - 1, k)) + ...
               (alpha_l * (Tl - T(i, j, k)) + (lambda * (T(i + 1, j, k) - T(i, j , k)) / hx)) * (2 * tau / (c * rho * hx)); 
            
            
            else if i == Nx
                % Right border
                T(i, j, k + 1) = T(i, j, k) + a * tau * (T(i, j + 1, k)- 2 * T(i, j, k) + T(i, j - 1, k)) + ...
               (alpha_r * (Tr - T(i, j, k)) + (lambda * (T(i - 1, j, k) - T(i, j , k)) / hx)) * (2 * tau / (c * rho * hx));
               
        
                else if j == 1
                        % Lower border
                T(i, j, k + 1) = T(i, j, k) + a * tau * (T(i + 1, j, k)- 2 * T(i, j, k) + T(i - 1, j, k)) + ...
               (alpha_d * (Td - T(i, j, k)) + (lambda * (T(i, j + 1, k) - T(i, j , k)) / hx)) * (2 * tau / (c * rho * hx));
                  
           
                    else if j == Ny
                            % Top border
                T(i, j, k + 1) = T(i, j, k) + a * tau * (T(i + 1, j, k)- 2 * T(i, j, k) + T(i - 1, j, k)) + ...
               (alpha_u * (Tu - T(i, j, k)) + (lambda * (T(i, j - 1, k) - T(i, j , k)) / hx)) * (2 * tau / (c * rho * hx)); 
                     
                        else 
                            % internal points
                     T(i , j, k + 1) = T(i, j, k) + a * tau * ((T(i + 1, j, k) - 2 * T(i, j, k) + T(i - 1, j, k)) / hx^2 + ...
                    (T(i, j + 1, k) - 2 * T(i, j, k) + T(i, j - 1, k)) / hy^2);
                        
                                          end
                                      end
                                  end
                              end 
                          end
                      end
                  end
              end
          end
     end
end
toc

else
disp('Wrong condition !')
        end
    end
end

Tend = zeros(Nx, Ny);
for i = 1 : Nx
    for j = 1 : Ny
Tend(i, j) = T(i, j, Nt);
    end
end
pcolor(Tend)
colorbar



