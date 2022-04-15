clearvars; close all; clc
set(0,'DefaultFigureWindowStyle','docked')
nx =  20;
ny =  10;
Vo = 0.8;
f = zeros(1,nx*ny);
G = sparse(ny*nx);
M = zeros(nx,ny);
l = 200e-9;
w = 100e-9;
sigma = ones(nx,ny);
T = 300;
q = 1.602e-19;
n_elec = 1e19; %m^-2 = 1e15 cm^-2   (per cm squared)
kb = 1.38064852*10^-23;
kt = 1.38064852*10^-23*T;
m = (9.10938356*10^-31)*0.26;
vth = sqrt(2*kt/m);
numberofelectrons = 50;
xp = zeros(numberofelectrons,1);
X = rand(numberofelectrons,1)*l;
Y = rand(numberofelectrons,1)*w;
yp = zeros(numberofelectrons,1);
Nt = 300;
tmn = 0.2e-12;
dt = (1/200)*(l/vth);
x = linspace(0,l,nx);
y = linspace(0,w,ny);
dx = l/(nx-1);
dy = w/(ny-1);


boxes = {};
boxes{1}.X = [0.8 0.4]*1e-7;
boxes{1}.Y = [0.6 0.4]*1e-7;
boxes{2}.X = [0.8 0.4]*1e-7;
boxes{2}.Y = [0 0.4]*1e-7;



it = (X < 1.2e-7) & (X > 0.8e-7) & (Y > 0.6e-7);
ib = (X < 1.2e-7) & (X > 0.8e-7) & (Y < 0.4e-7);
iXY = it | ib;
while(sum(iXY) > 0)
    X(iXY) = rand(1,sum(iXY))*l;
    Y(iXY) = rand(1,sum(iXY))*w;
    it = X < 1.2e-7 & X > 0.8e-7 & Y > 0.6e-7;
    ib = X < 1.2e-7 & X > 0.8e-7 & Y < 0.4e-7;
    iXY = it | ib;
end



for i = 1:nx
   for j = 1:ny
    if x(i)> l*0.3 & x(i)< l*0.65 & (y(j)< w*0.3 | y(j)> w*0.6)
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
            f(n) = Vo;
            
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
       %surf(M) %,'linestyle','none')  
    end
    
    
[Ey,Ex] = gradient(M);
Ex = -Ex/dx;
Ey = -Ey/dy;
Fx = (Ex)*(1.602e-19);
Fy = (Ey)*(1.602e-19);
ax = Fx/m; 
ay= Fy/m; 

%maxwell-boltzmann distribution
vx = randn(numberofelectrons,1) *vth/sqrt(2);
vy = randn(numberofelectrons,1) *vth/sqrt(2);
vx_thermalized = randn(numberofelectrons,1) *vth/sqrt(2);
vy_thermalized = randn(numberofelectrons,1) *vth/sqrt(2);



%retangle(xmin ymin howwide howlong)
for z = 1:2
hold on
rectangle('position',[boxes{z}.X(1), boxes{z}.Y(1),boxes{z}.X(2),boxes{z}.Y(2)])
axis([0 l 0 w])
end

for t =1:Nt  
    %position calculations 
    xp = X;
    yp = Y;
   
    xbin = discretize(X,nx);
    ybin = discretize(Y,ny);
    
    axp = zeros(numberofelectrons,1);
    ayp = zeros(numberofelectrons,1);
 
    for i = 1:numberofelectrons
        axp(i) = ax(xbin(i),ybin(i));
        ayp(i) = ay(xbin(i),ybin(i));
    end

    X = X + vx*dt+(1/2)*(axp)*(dt^2);
    Y = Y + vy*dt+(1/2)*(ayp)*(dt^2);
    
    vx = vx + axp*dt;
    vy = vy + ayp*dt;

    %scattering probabilty
    Pscat = 1-exp(-dt/tmn);
    scat = Pscat > rand(numberofelectrons,1);
    vx(scat) = vx_thermalized(scat);
    vy(scat) = vy_thermalized(scat);
      
    %Boundary conditions 
    %inside top and bottom box
    InTopBox = X < 1.2e-7 & X > 0.8e-7 & Y > 0.6e-7;
    InBotBox = X < 1.2e-7 & X > 0.8e-7 & Y < 0.4e-7;
    
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
    X(reflectxbotright & InBotBox) = xp(reflectxbotright & InBotBox);
    X(reflextxbotleft & InBotBox) = xp(reflextxbotleft & InBotBox);
    
    %reverse x for top box
    X(reflectxtopright & InTopBox) = xp(reflectxtopright & InTopBox);
    X(reflextxtopleft & InTopBox) = xp(reflextxtopleft & InTopBox);
    
    %reverse y for top and bottom box
    Y(reflectxbotright & InBotBox) = yp(reflectxbotright & InBotBox);
    Y(reflextxbotleft & InBotBox) = yp(reflextxbotleft & InBotBox);
    Y(reflectxtopright & InTopBox) = yp(reflectxtopright & InTopBox);
    Y(reflextxtopleft & InTopBox) = yp(reflextxtopleft & InTopBox);
    
    %reverse yspeed and y for middle 
    vy(reflectmiddle & (InTopBox|InBotBox)) = -vy(reflectmiddle & (InTopBox|InBotBox));
    Y(reflectmiddle & (InTopBox|InBotBox)) = yp(reflectmiddle & (InTopBox|InBotBox));

    %periodic boundary condition
    xp(X>l) = xp(X>l)-l;
    X(X>l) = X(X>l)-l;
    xp(X<0) = xp(X<0)+l;
    X(X<0) = X(X<0)+l;
     

   %particle hitting top and bottom
   vy(Y>w) = -vy(Y>w);
   vy(Y<0) = -vy(Y<0);
    
   
    v_x(t) = mean(vx);
                            % 3D: C elec/m^3 m/s = A/m^2
                            %     I = J*Area
    Jx = q*n_elec*v_x;   % 2D: C elec/m^2 m/s = A/m
    I = Jx*w;           %     I = J*height
    % n_particle = num of particles/L*W
    % simple example 3D:
    % n_elec = 100 elec/unit vol
    % n_particles = 10 particles/unit vol
    % elec/particle = 10
   
 
    %ploting 
    for i = 1:numberofelectrons
         plot([xp(i),X(i)],[yp(i),Y(i)],'Seriesindex',i)
          hold on 
    end
      pause(0.1)
end
 


binx = ceil(X*nx/l);   % 100nm -> bin 5
% ->     [ 1, 4, 7, 7 .. bin_of_last_part ]
biny = ceil(Y*ny/w);

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
view(60,30)
title('Density')
figure(3)
plot(z,I)
title('Time vs Current')
ylabel('Current in Amperes')
xlabel('Time Steps')
figure(4)
quiver(Ex',Ey')
figure(5)
surf(temperature)
view(60,30)
title('Temperature')
figure(6)
surf(M)
view(90,5)
figure(7)
surf(sigma)  
view(60,30)
