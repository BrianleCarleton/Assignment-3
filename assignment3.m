clearvars; close all; clc
set(0,'DefaultFigureWindowStyle','docked')
T = 300;
q = 1.602e-19;
n_elec = 1*10e13;
kb = 1.38064852*10^-23;
kt = 1.38064852*10^-23*T;
m = (9.10938356*10^-31)*0.26;
vth = sqrt(2*kt/m);
width = 200e-9;
height = 100e-9;
n = 1000;
xp = zeros(n,1);
x = (width).*rand(n,1) ;
y = (height).*rand(n,1);
yp = zeros(n,1);
Nt = 7;
tmn = 0.2e-12;
dt = (1/200)*(width/vth) ;
sample = randi(n,n,1);
F = (0.1/200e-9)*(1.602e-19);
a = F/m;

%maxwell-boltzmann distribution
vx = randn(n,1) *vth/sqrt(2);
vy = randn(n,1) *vth/sqrt(2);
vx_thermalized = randn(n,1) *vth/sqrt(2);
vy_thermalized = randn(n,1) *vth/sqrt(2);


hold on
axis([0 width 0 height])


for t =1:Nt  

        
    %position calculations 
    xp = x;
    yp = y;
    
    vx = vx + a*dt;
    
    x = x + vx*dt+(1/2)*(a)*(dt^2) ;
    y = y + (vy)*dt ;
    

    %scattering probabilty
    Pscat = 1-exp(-dt/tmn);
    scat = Pscat > rand(n,1);
    vx(scat) = vx_thermalized(scat);
    vy(scat) = vy_thermalized(scat);


     %periodic boundary condition
     xp(x>width) = xp(x>width)-width;
     x(x>width) = x(x>width)-width;
     xp(x<0) = xp(x<0)+width;
     x(x<0) = x(x<0)+width;
     

   %particle hitting top and bottom
    vy(y>height) = -vy(y>height);
    vy(y<0) = -vy(y<0);
    
    vx_ave(t) = sum(vx)/n;
                            % 3D: C elec/m^3 m/s = A/m^2
                            %     I = J*Area
    Jx = q*n_elec*vx_ave;   % 2D: C elec/m^2 m/s = A/m
    I = Jx*width;           %     I = J*Width
    % n_particle = num of particles/L*W
    % simple example 3D:
    % n_elec = 100 elec/unit vol
    % n_particles = 10 particles/unit vol
    % elec/particle = 10


    
    %ploting 
    for i = 1:n
       plot([xp(i),x(i)],[yp(i),y(i)],'Seriesindex',i)
        hold on 
    end

       pause(0.1)
end
 

Nbinx = 40;
Nbiny = 30;
binx = ceil(x*Nbinx/width);   % 100nm -> bin 5
% ->     [ 1, 4, 7, 7 .. bin_of_last_part ]
biny = ceil(y*Nbiny/height);

for  i=1:Nbinx
    for j = 1:Nbiny
    match = binx==i & biny==j;
    ebox(i,j) = sum(match);
    end
end

for  i=1:Nbinx
    for j = 1:Nbiny
    match = binx==i & biny==j;
    temperature(i,j) = sum(((sqrt(((vx(match)).^2)+((vy(match)).^2))*m)/2)/kb);
    end
end

z = 1:1:Nt;
figure(2)
surf(ebox)
title('Density')
figure(3)
plot(z,I)
title('Time vs Current')
figure(4)
surf(temperature)
title('Temperature')
        

