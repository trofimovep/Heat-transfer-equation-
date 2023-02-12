clc; clear;
% Matetial
rho = 7800;
c = 500;
lambda = 46;

% rho = input('Enter the density:');
% c = input('Enter the heat capacity: ');
% lambda = input('Enter the coefficient of thermal conductivity: ');  

a = lambda / (rho * c);       % thermal diffusitivy

% Parameters 
A = 0.2;           % i           %  height
B = 0.2;           % j           %  width  
C = 0.02 ;          % k           %  length
t = 600;            % p

% A = input('Enter the height: ');
% B = input('Enter the width: ');
% C = input('Enter the length: ');
% t = input('Enter the time: ');

%  Bordering
T0 = 0;
% T0 = input('Enter the initial Temperature T0: ');

dx = 0.01; 
dy = dx;
dz = dx;
dt = dx / 2;

Nx = round(A / (dx));
Ny = round(B / (dy));
Nz = round(C / (dz));
Nt = round(t / dt);

%Tempearture of environment
T_left = 1500; % left
T_right = 500; % right
T_up = 500; % up
T_down = 5000; % down
T_front = 0; % front
T_back = 10; % back

% flow
alpha_up = 1900; 
alpha_left = 850;
alpha_right = 850;
alpha_down = 1900/5;
alpha_front = 850;
alpha_back = 850;

parfor i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
             T(i, j, k, 1) = T0;
        end
    end
end
size(T)

tic
for k = 1 : Nt-1
    for x = 1 : Nx
        for y = 1 : Ny
            for z = 1 : Nz 
                    
%throw the corners

if x == 1 && y == 1 && z == 1  
    T(x, y, z, k+1) = T(x, y, z, k)...
        + 2 * a * dt * (T(x+1, y, z, k) + T(x, y + 1, z, k) + T(x, y, z+1, k) * 3 * T(x, y, z, k))/dx^2 ...
+ (alpha_down * (T_down - T(x, y, z, k)) + alpha_back * (T_back - T(x, y, z, k)) + alpha_left * (T_left - T(x, y, z, k))) ...
* (2 * dt / c / rho / dx);

elseif x == Nx && y == 1 && z == 1
        T(x, y, z, k+1) = T(x, y, z, k)...
        + 2 * a * dt * (T(x - 1, y, z, k) + T(x, y + 1, z, k) + T(x, y, z+1, k) * 3 * T(x, y, z, k))/dx^2 ...
+ (alpha_down * (T_down - T(x, y, z, k)) + alpha_front * (T_front - T(x, y, z, k)) + alpha_left * (T_left - T(x, y, z, k))) ...
* (2 * dt / c / rho / dx);

elseif x == 1 && y == Ny && z == 1
        T(x, y, z, k+1) = T(x, y, z, k)...
        + 2 * a * dt * (T(x +1, y, z, k) + T(x, y-1, z, k) + T(x, y, z+1, k) * 3 * T(x, y, z, k))/dx^2 ...
+ (alpha_down * (T_down - T(x, y, z, k)) + alpha_back * (T_back - T(x, y, z, k)) + alpha_right * (T_right - T(x, y, z, k))) ...
* (2 * dt / c / rho / dx);

 elseif x == 1 && y == 1  && z == Nz
        T(x, y, z, k+1) = T(x, y, z, k)...
        + 2 * a * dt * (T(x +1, y, z, k) + T(x, y+1, z, k) + T(x, y, z-1, k) * 3 * T(x, y, z, k))/dx^2 ...
+ (alpha_up * (T_up - T(x, y, z, k)) + alpha_back * (T_back - T(x, y, z, k)) + alpha_left * (T_left - T(x, y, z, k))) ...
* (2 * dt / c / rho / dx);

 elseif x == 1 && y == Ny && z == Nz
        T(x, y, z, k+1) = T(x, y, z, k)...
        + 2 * a * dt * (T(x +1, y, z, k) + T(x, y-1, z, k) + T(x, y, z-1, k) * 3 * T(x, y, z, k))/dx^2 ...
+ (alpha_up * (T_up - T(x, y, z, k)) + alpha_back * (T_back - T(x, y, z, k)) + alpha_right * (T_right - T(x, y, z, k))) ...
* (2 * dt / c / rho / dx);

 elseif x == Nx && y == 1 && z == Nz  
       T(x, y, z, k+1) = T(x, y, z, k)...
        + 2 * a * dt * (T(x-1, y, z, k) + T(x, y+1, z, k) + T(x, y, z-1, k) * 3 * T(x, y, z, k))/dx^2 ...
+ (alpha_up * (T_up - T(x, y, z, k)) + alpha_front * (T_front - T(x, y, z, k)) + alpha_left * (T_left - T(x, y, z, k))) ...
* (2 * dt / c / rho / dx);      
 
elseif x == Nx && y == Ny && z == 1  
        T(x, y, z, k+1) = T(x, y, z, k)...
        + 2 * a * dt * (T(x-1, y, z, k) + T(x, y-1, z, k) + T(x, y, z+1, k) * 3 * T(x, y, z, k))/dx^2 ...
+ (alpha_down * (T_down -  - T(x, y, z, k)) + alpha_front * (T_front - T(x, y, z, k)) + alpha_right * (T_right - T(x, y, z, k))) ...
* (2 * dt / c / rho / dx);

    elseif x == Nx && y == Ny && z == Nz
        T(x, y, z, k+1) = T(x, y, z, k)...
        + 2 * a * dt * (T(x-1, y, z, k) + T(x, y-1, z, k) + T(x, y, z-1, k) * 3 * T(x, y, z, k))/dx^2 ...
+ (alpha_up * (T_up - T(x, y, z, k)) + alpha_front * (T_front - T(x, y, z, k)) + alpha_right * (T_right - T(x, y, z, k))) ...
* (2 * dt / c / rho / dx);
             

%%% throw edges

        elseif y == 1 && z == 1
       T(x, y, z, k+1) = T(x, y, z, k) +...
   a * ((2 * (T(x, y+1, z, k) - 2 * T(x, y, z, k) + T(x, y, z+1))/dx^2) + ((T(x-1,y, z) - 2*T(x, y, z) + T(x+1, y, z))/dx^2))+...
   (alpha_down * (T_down - T(x, y, z, k)) + alpha_left * (T_left - T(x, y, z, k))) * (2 * dt / c /rho / dx);
        
            elseif y == 1 && z == Nz        
       T(x, y, z, k+1) = T(x, y, z, k) + ...
   a * ((2 * (T(x, y+1, z, k) - 2 * T(x, y, z, k) + T(x, y, z-1))/dx^2) + ((T(x-1,y, z) - 2*T(x, y, z) + T(x+1, y, z))/dx^2))+...
   (alpha_up * (T_up - T(x, y, z, k)) + alpha_left * (T_left - T(x, y, z, k))) * (2 * dt / c /rho / dx);


           elseif y == Ny && z == 1        
       T(x, y, z, k+1) = T(x, y, z, k) + ...
   a * ((2 * (T(x, y-1, z, k) - 2 * T(x, y, z, k) + T(x, y, z+1))/dx^2) + ((T(x-1,y, z) - 2*T(x, y, z) + T(x+1, y, z))/dx^2))+...
   (alpha_down * (T_down - T(x, y, z, k)) + alpha_right * (T_right - T(x, y, z, k))) * (2 * dt / c /rho / dx);
   
          elseif y == Ny && z == Nz        
       T(x, y, z, k+1) = T(x, y, z, k) + ...
   a * ((2 * (T(x, y-1, z, k) - 2 * T(x, y, z, k) + T(x, y, z-1))/dx^2) + ((T(x-1,y, z) - 2*T(x, y, z) + T(x+1, y, z))/dx^2))+...
   (alpha_up * (T_up - T(x, y, z, k)) + alpha_right * (T_right - T(x, y, z, k))) * (2 * dt / c /rho / dx);


elseif x == 1 && z == 1        
       T(x, y, z, k+1) = T(x, y, z, k) + ...
   a * ((2 * (T(x+1, y, z, k) - 2 * T(x, y, z, k) + T(x, y, z+1))/dx^2) + ((T(x,y-1, z) - 2*T(x, y, z) + T(x, y+1, z))/dx^2))+...
   (alpha_down * (T_down - T(x, y, z, k)) + alpha_back * (T_back - T(x, y, z, k))) * (2 * dt / c /rho / dx);

elseif x == 1 && z == Nz        
       T(x, y, z, k+1) = T(x, y, z, k) + ...
   a * ((2 * (T(x+1, y, z, k) - 2 * T(x, y, z, k) + T(x, y, z-1))/dx^2) + ((T(x,y-1, z) - 2*T(x, y, z) + T(x, y+1, z))/dx^2))+...
   (alpha_up * (T_up - T(x, y, z, k)) + alpha_back * (T_back - T(x, y, z, k))) * (2 * dt / c /rho / dx);

elseif x == Nx && z == 1        
       T(x, y, z, k+1) = T(x, y, z, k) + ...
   a * ((2 * (T(x-1, y, z, k) - 2 * T(x, y, z, k) + T(x, y, z+1))/dx^2) + ((T(x,y-1, z) - 2*T(x, y, z) + T(x, y+1, z))/dx^2))+...
   (alpha_down * (T_down - T(x, y, z, k)) + alpha_front * (T_front - T(x, y, z, k))) * (2 * dt / c /rho / dx);

elseif x == Nx && z == Nz        
       T(x, y, z, k+1) = T(x, y, z, k) + ...
   a * ((2 * (T(x-1, y, z, k) - 2 * T(x, y, z, k) + T(x, y, z-1))/dx^2) + ((T(x,y-1, z) - 2*T(x, y, z) + T(x, y+1, z))/dx^2))+...
   (alpha_up * (T_up - T(x, y, z, k)) + alpha_front * (T_front - T(x, y, z, k))) * (2 * dt / c /rho / dx);

elseif x == Nx && z == Nz        
       T(x, y, z, k+1) = T(x, y, z, k) + ...
   a * ((2 * (T(x-1, y, z, k) - 2 * T(x, y, z, k) + T(x, y, z-1))/dx^2) + ((T(x,y-1, z) - 2*T(x, y, z) + T(x, y+1, z))/dx^2))+...
   (alpha_up * (T_up - T(x, y, z, k)) + alpha_front * (T_front - T(x, y, z, k))) * (2 * dt / c /rho / dx);



elseif x == 1 && y == 1        
       T(x, y, z, k+1) = T(x, y, z, k) + ...
   a * ((2 * (T(x+1, y, z, k) - 2 * T(x, y+1, z, k) + T(x, y, z))/dx^2) + ((T(x,y, z-1) - 2*T(x, y, z) + T(x, y+1, z+1))/dx^2))+...
   (alpha_left * (T_left - T(x, y, z, k)) + alpha_back * (T_back - T(x, y, z, k))) * (2 * dt / c /rho / dx);

elseif x == 1 && y == Ny        
       T(x, y, z, k+1) = T(x, y, z, k) + ...
   a * ((2 * (T(x+1, y, z, k) - 2 * T(x, y-1, z, k) + T(x, y, z))/dx^2) + ((T(x,y, z-1) - 2*T(x, y, z) + T(x, y+1, z+1))/dx^2))+...
   (alpha_right * (T_right - T(x, y, z, k)) + alpha_back * (T_back - T(x, y, z, k))) * (2 * dt / c /rho / dx);

elseif x == Nx && y == 1        
       T(x, y, z, k+1) = T(x, y, z, k) + ...
   a * ((2 * (T(x-1, y, z, k) - 2 * T(x, y+1, z, k) + T(x, y, z))/dx^2) + ((T(x,y, z-1) - 2*T(x, y, z) + T(x, y+1, z+1))/dx^2))+...
   (alpha_left * (T_left - T(x, y, z, k)) + alpha_front * (T_front - T(x, y, z, k))) * (2 * dt / c /rho / dx);

elseif x == Nx && y == Ny        
       T(x, y, z, k+1) = T(x, y, z, k) + ...
   a * ((2 * (T(x-1, y, z, k) - 2 * T(x, y-1, z, k) + T(x, y, z))/dx^2) + ((T(x,y, z-1) - 2*T(x, y, z) + T(x, y+1, z+1))/dx^2))+...
   (alpha_right * (T_right - T(x, y, z, k)) + alpha_front * (T_front - T(x, y, z, k))) * (2 * dt / c /rho / dx);

%%%         throw facet


    elseif z == 1
            
T(x, y, z, k+1) = T(x, y, z, k) + ...
    a * dt * (((T(x-1, y, z, k)- 2 * T(x, y, z, k) + T(x+1, y, z, k))/dx^2 + (T(x, y-1, z, k) - 2 * T(x, y, z, k) + T(x, y+1, z, k))/dx^2)...
    + (lambda * (T(x , y, z + 1, k) - T(x, y, z, k)))/dx + alpha_down * (T_down - T(x, y, z, k)))...
    *(2 * dt /c / rho / dx);

elseif x == 1
                
                T(x, y, z, k+1) = T(x, y, z, k) + ...
    a * dt * (((T(x, y-1, z, k)- 2 * T(x, y, z, k) + T(x, y+1, z, k))/dx^2 + (T(x, y, z-1, k) - 2 * T(x, y, z, k) + T(x, y, z+1, k))/dx^2)...
    + (lambda * (T(x + 1 , y, z, k) - T(x, y, z, k)))/dx + alpha_back * (T_back - T(x, y, z, k)))...
    *(2 * dt /c / rho / dx);

elseif y == 1
                
                T(x, y, z, k+1) = T(x, y, z, k) + ...
    a * dt * (((T(x-1, y, z, k)- 2 * T(x, y, z, k) + T(x+1, y, z, k))/dx^2 + (T(x, y, z-1, k) - 2 * T(x, y, z, k) + T(x, y, z+1, k))/dx^2)...
    + (lambda * (T(x , y + 1, z, k) - T(x, y, z, k)))/dx + alpha_left * (T_left - T(x, y, z, k)))...
    *(2 * dt /c / rho / dx);



        elseif z == Nz
                
                T(x, y, z, k+1) = T(x, y, z, k) + ...
    a * dt * (((T(x-1, y, z, k)- 2 * T(x, y, z, k) + T(x+1, y, z, k))/dx^2 + (T(x, y-1, z, k) - 2 * T(x, y, z, k) + T(x, y+1, z, k))/dx^2)...
    + (lambda * (T(x , y, z - 1, k) - T(x, y, z, k)))/dx + alpha_up * (T_up - T(x, y, z, k)))...
    *(2 * dt /c / rho / dx);


elseif y == Nz
                
                T(x, y, z, k+1) = T(x, y, z, k) + ...
    a * dt * (((T(x-1, y, z, k)- 2 * T(x, y, z, k) + T(x+1, y, z, k))/dx^2 + (T(x, y, z-1, k) - 2 * T(x, y, z, k) + T(x, y, z+1, k))/dx^2)...
    + (lambda * (T(x , y - 1, z, k) - T(x, y, z, k)))/dx + alpha_right * (T_right - T(x, y, z, k)))...
    *(2 * dt /c / rho / dx);



elseif x == Nx
                
                T(x, y, z, k+1) = T(x, y, z, k) + ...
    a * dt * (((T(x, y-1, z, k)- 2 * T(x, y, z, k) + T(x, y+1, z, k))/dx^2 + (T(x, y, z-1, k) - 2 * T(x, y, z, k) + T(x, y, z+1, k))/dx^2)...
    + (lambda * (T(x - 1 , y, z, k) - T(x, y, z, k)))/dx + alpha_front * (T_front - T(x, y, z, k)))...
    *(2 * dt /c / rho / dx);

    else
        
        T(x, y, z, k+1) = T(x, y, z, k) + ...
            a * dt * (((T(x-1, y, z) - 2 * T(x, y, z, k) + T(x+1, y, z, k))/dx^2 + ...
            ((T(x, y -1, z, k) - 2 * T(x, y, z, k) + T(x, y+1, z, k))/dx^2) + ...
            ((T(x, y, z+1, k) - 2 * T(x, y, z, k) + T(x, y, z-1, k))/dx^2)));
        
end
            end
        end
    end
        
        clearvars T(x, y, z, k - 2);
end
toc
           

 %% in XY flat 
 Txy = zeros(Nx, Ny);
  zz = Nz;
  for i = 1 : Nx
    for j = 1 : Ny
          Txy(i, j) = T(i, j, 1, Nt);
     end
  end
  pcolor(Txy)
  
%   %% in YZ flat
%   Tyz = zeros(Ny, Nz);
% %   x = round(B / 2); 
%   xx = round(Nx / 2);
% for i = 1:Ny
%     for j = 1:Nz
%           Tyz(i, j) = T(xx, i, j, Nt);
%      end
% end 
% pcolor(Tyz)
% 
%   %% in XZ flat
%   Txz = zeros(Nx, Nz);
%   yy = round(Ny / 2); 
% for i = 1:Nx
%     for j = 1:Nz
%           Txz(i, j) = T(i, yy, j, Nt);
%      end
% end 
% pcolor(Txz)
                    


