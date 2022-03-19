clearvars; close all; clc
set(0,'DefaultFigureWindowStyle','docked')
nx = 100;
ny = 50;
f = zeros(1,nx*ny);
G = sparse(ny*nx);
M = zeros(nx,ny);
l = 100e-9;
w = 200e-9;
sigma = ones(nx,ny);
xx = linspace(0,l,nx);
yy = linspace(0,w,ny);
T = 300;
q = 1.602e-19;
n_elec = 1*10e13;
kb = 1.38064852*10^-23;
kt = 1.38064852*10^-23*T;
m = (9.10938356*10^-31)*0.26;
vth = sqrt(2*kt/m);
width = 200e-9;
height = 100e-9;
numberofelectrons = 1000;
xp = zeros(numberofelectrons,1);
x = (width).*rand(numberofelectrons,1);
y = (0.55e-7-0.45e-7).*rand(numberofelectrons,1)+0.45e-7;
yp = zeros(numberofelectrons,1);
Nt = 1000;
tmn = 0.2e-12;
dt = (1/200)*(width/vth) ;
sample = randi(numberofelectrons,numberofelectrons,1);
  


for i = 1:nx
   for j = 1:ny
    if xx(i)> l*0.3 & xx(i)< l*0.65 & (yy(j)< w*0.3 | yy(j)> w*0.6)
        sigma(i,j) = 10e-2;
    end
   end
end

    

for i = 1:nx 
        for j = 1:ny
        n = j +(i-1)*ny;
        nxm = j + (i-2)*ny;
        nxp = j + (i)*ny;
        nyp = (j+1) + (i-1)*ny; 
        nym = (j-1) +(i-1)*ny;
        
            if i == 1
            G(n,n) = sigma(i,j);
            f(n) = 1;
            
            elseif i == nx
            G(n,n) = sigma(i,j);   
            
            elseif j == 1 %bottom   
            G(n,nxp) = sigma(i+1,j);
            G(n,nxm) = sigma(i-1,j);
            G(n,nyp) = sigma(i,j+1);
            G(n,n) = -(G(n,nxp)+G(n,nxm)+G(n,nyp));
            
            elseif j == ny %top
            G(n,nxp) = sigma(i+1,j);
            G(n,nxm) = sigma(i-1,j);
            G(n,nym) = sigma(i,j-1);
            G(n,n) = -(G(n,nxp)+G(n,nxm)+G(n,nym));
          
            else 
            G(n,nxp) = sigma(i+1,j);
            G(n,nxm) = sigma(i-1,j);
            G(n,nym) = sigma(i,j-1);
            G(n,nyp) = sigma(i,j+1);
            G(n,n) = -(G(n,nxp)+G(n,nxm)+G(n,nym)+G(n,nyp));
            end           
        end
    end    
V = G\f';

    for i = 1:nx
        for j = 1:ny
            n = j + (i-1)*ny;
            M(i,j) = V(n);
        end
       % surf(M,'linestyle','none')  
    end
    
[Ex,Ey] = gradient(M);

Fx = (Ex)*(1.602e-19);
Fy = (Ey)*(1.602e-19);
ax = Fx/m; 
ay= Fy/m; 

%maxwell-boltzmann distribution
vx = randn(numberofelectrons,1) *vth/sqrt(2);
vy = randn(numberofelectrons,1) *vth/sqrt(2);
vx_thermalized = randn(numberofelectrons,1) *vth/sqrt(2);
vy_thermalized = randn(numberofelectrons,1) *vth/sqrt(2);

boxes = {};
boxes{1}.X = [0.8 0.4]*1e-7;
boxes{1}.Y = [0.6 0.4]*1e-7;
boxes{2}.X = [0.8 0.4]*1e-7;
boxes{2}.Y = [0 0.4]*1e-7;

%retangle(xmin ymin howwide howlong)
for z = 1:2
hold on
rectangle('position',[boxes{z}.X(1), boxes{z}.Y(1),boxes{z}.X(2),boxes{z}.Y(2)])
axis([0 width 0 height])
end

for t =1:Nt  
    %position calculations 
    xp = x;
    yp = y;
   
    xbin = discretize(x,nx);
    ybin = discretize(y,ny);
    
    axp = zeros(numberofelectrons,1);
    ayp = zeros(numberofelectrons,1);
 
    for i = 1:numberofelectrons
        axp(i) = ax(xbin(i),ybin(i));
        ayp(i) = ay(xbin(i),ybin(i));
    end

    x = x + vx*dt+(1/2)*(axp)*(dt^2) ;
    y = y + vy*dt+(1/2)*(ayp)*(dt^2)  ;
    
    vx = vx + axp*dt;
    vy = vy + ayp*dt;

    %scattering probabilty
    Pscat = 1-exp(-dt/tmn);
    scat = Pscat > rand(numberofelectrons,1);
    vx(scat) = vx_thermalized(scat);
    vy(scat) = vy_thermalized(scat);
      
    %Boundary conditions 
    %inside top and bottom box
    InTopBox = x < 1.2e-7 & x > 0.8e-7 & y > 0.6e-7;
    InBotBox = x < 1.2e-7 & x > 0.8e-7 & y < 0.4e-7;
    
    %right and left of the boxes
    reflectxtopright = xp > 1.2e-7 & yp > 0.6e-7;
    reflextxtopleft = xp < 0.8e-7 & yp > 0.6e-7;
    reflectxbotright = xp > 1.2e-7 & yp < 0.4e-7;
    reflextxbotleft = xp < 0.8e-7 & yp < 0.4e-7;
    reflectmiddle = xp < 1.2e-7 & xp > 0.8e-7;

    %reverse xspeed for topbox 
    vx(reflectxtopright & InTopBox) = -vx(reflectxtopright & InTopBox);
    vx(reflextxtopleft & InTopBox) = -vx(reflextxtopleft & InTopBox);
   
    %reverse xspeed for bottom box
    vx(reflectxbotright & InBotBox) = -vx(reflectxbotright & InBotBox);
    vx(reflextxbotleft & InBotBox) = -vx(reflextxbotleft & InBotBox);
    
    %reverse x for bottom box
    x(reflectxbotright & InBotBox) = xp(reflectxbotright & InBotBox);
    x(reflextxbotleft & InBotBox) = xp(reflextxbotleft & InBotBox);
    
    %reverse x for top box
    x(reflectxtopright & InTopBox) = xp(reflectxtopright & InTopBox);
    x(reflextxtopleft & InTopBox) = xp(reflextxtopleft & InTopBox);
    
    %reverse y for top and bottom box
    y(reflectxbotright & InBotBox) = yp(reflectxbotright & InBotBox);
    y(reflextxbotleft & InBotBox) = yp(reflextxbotleft & InBotBox);
    y(reflectxtopright & InTopBox) = yp(reflectxtopright & InTopBox);
    y(reflextxtopleft & InTopBox) = yp(reflextxtopleft & InTopBox);
    
    %reverse yspeed and y for middle 
    vy(reflectmiddle & (InTopBox|InBotBox)) = -vy(reflectmiddle & (InTopBox|InBotBox));
    y(reflectmiddle & (InTopBox|InBotBox)) = yp(reflectmiddle & (InTopBox|InBotBox));

    %periodic boundary condition
    xp(x>width) = xp(x>width)-width;
    x(x>width) = x(x>width)-width;
    xp(x<0) = xp(x<0)+width;
    x(x<0) = x(x<0)+width;
     

   %particle hitting top and bottom
   vy(y>height) = -vy(y>height);
   vy(y<0) = -vy(y<0);
    
   
    v_ave(t) = sum(sqrt((vx.^2)+(vy.^2)))/n;
                            % 3D: C elec/m^3 m/s = A/m^2
                            %     I = J*Area
    Jx = q*n_elec*v_ave;   % 2D: C elec/m^2 m/s = A/m
    I = Jx*width;           %     I = J*Width
    % n_particle = num of particles/L*W
    % simple example 3D:
    % n_elec = 100 elec/unit vol
    % n_particles = 10 particles/unit vol
    % elec/particle = 10
   
 
    %ploting 
    for i = 1:numberofelectrons
       plot([xp(i),x(i)],[yp(i),y(i)],'Seriesindex',i)
        hold on 
    end
    
     pause(0.1)
end
 


binx = ceil(x*nx/width);   % 100nm -> bin 5
% ->     [ 1, 4, 7, 7 .. bin_of_last_part ]
biny = ceil(y*ny/height);

for  i=1:nx
    for j = 1:ny
    match = binx==i & biny==j;
    ebox(i,j) = sum(match);
    end
end




for  i=1:nx
    for j = 1:ny
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
quiver(Ex,Ey)
figure(5)
surf(temperature)
title('Temperature')
figure(6)
surf(M,'linestyle','none')  
figure(7)
surf(sigma)  
