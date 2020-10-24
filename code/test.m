clear;
clc;
% This program is built to solve the lid driven flow in a cavity by applaying lattice Boltzmann method.

% condition constant
Lx=100; %\SI{m}
Ly=100; %\SI{m}

%number of cells in x and y direction
dL=1; %\SI{m}
dT=1; %\SI{s}
NX=round(Lx/dL);
NY=round(Ly/dL);

% D2Q9 lattice constant
global w;
w =[ 4/9,  1/9,  1/9,  1/9,  1/9, 1/36, 1/36, 1/36, 1/36];
global c;
c =[0, 0; 1, 0; 0, 1;-1, 0; 0,-1; 1, 1;-1, 1;-1,-1; 1,-1];

[x,y]=meshgrid(1:NX+1,1:NY+1);
density = zeros(NY+1,NX+1);
u       = zeros(NY+1,NX+1,2);
f       = zeros(NY+1,NX+1,9);
f_eq    = zeros(NY+1,NX+1,9);

%initialize
c_s =sqrt(1/3)*dL/dT;
density_0=1;  %\SI{kg/m^3}
nu   =0.0667; %\SI{m^2/s}
Re=1000;
U=nu*Re/Lx;
tau=nu/(c_s^2)+0.5*dT;

density(:,:)=density_0;
u(NY+1,:,1)=U;
for k=1:9
    f(:,:,k)=w(k)*density(:,:).*(1+(u(:,:,1)*c(k,1)+u(:,:,2)*c(k,2))/(c_s^2)+((u(:,:,1)*c(k,1)+u(:,:,2)*c(k,2)).^2)/(2*c_s^4)-(u(:,:,1).*u(:,:,1)+u(:,:,2).*u(:,:,2))/(2*c_s^2));
end


% lattice Boltzmann starts
% video recording starts
objvideo=VideoWriter('example.avi');
objvideo.FrameRate=round(120);
open(objvideo)
H=525*3;
W=700*3;

MaxIteration=1800;%max ieration number
for ite=1:MaxIteration
    
    % collision
    for k=1:9
        f_eq=w(k)*density(:,:).*(1+(u(:,:,1)*c(k,1)+u(:,:,2)*c(k,2))/(c_s^2)+((u(:,:,1)*c(k,1)+u(:,:,2)*c(k,2)).^2)/(2*c_s^4)-(u(:,:,1).*u(:,:,1)+u(:,:,2).*u(:,:,2))/(2*c_s^2));
        f(:,:,k)=f(:,:,k)*(1-dT/tau)+f_eq*dT/tau;
    end
    
    % streaming
    f(:,2:NX+1,2) = f(:,  1:NX,2);% left to right
    f(:,  1:NX,4) = f(:,2:NX+1,4);% right to left
    f(2:NY+1,:,3) = f(1:NY  ,:,3);% bottom to top
    f(1:NY  ,:,5) = f(2:NY+1,:,5);% top to bottom
    f(2:NY+1,2:NX+1,6) = f(  1:NY,  1:NX,6); % leftbottom to righttop
    f(2:NY+1,  1:NX,7) = f(  1:NY,2:NX+1,7); % rightbottom to lefttop
    f(1:NY  ,1:NX  ,8) = f(2:NY+1,2:NX+1,8); % righttop to leftbottom
    f(1:NY  ,2:NX+1,9) = f(2:NY+1,  1:NX,9); % lefttop to rightbottom
    
    % boundary condition
    % left bounce back
    f(2:NY+1,1,2)=f(2:NY+1,1,4);
    f(2:NY+1,1,6)=f(2:NY+1,1,8);
    f(2:NY+1,1,9)=f(2:NY+1,1,7);
    
    % right bounce back
    f(2:NY+1,NX+1,4)=f(2:NY+1,NX+1,2);
    f(2:NY+1,NX+1,7)=f(2:NY+1,NX+1,9);
    f(2:NY+1,NX+1,8)=f(2:NY+1,NX+1,6);
    
    % bottom bounce back
    f(1,:,3)=f(1,:,5);
    f(1,:,6)=f(1,:,8);
    f(1,:,7)=f(1,:,9);
    
    % moving lid
    density_temp=f(NY+1,2:NX,1)+f(NY+1,2:NX,2)+f(NY+1,2:NX,4)+2*(f(NY+1,2:NX,3)+f(NY+1,2:NX,6)+f(NY+1,2:NX,7));
    f(NY+1,2:NX,5)=f(NY+1,2:NX,3);
    f(NY+1,2:NX,9)=f(NY+1,2:NX,7)-0.5*(f(NY+1,2:NX,2)-f(NY+1,2:NX,4))+0.5*density_temp*U;
    f(NY+1,2:NX,8)=f(NY+1,2:NX,6)+0.5*(f(NY+1,2:NX,2)-f(NY+1,2:NX,4))-0.5*density_temp*U;
    
    % computate density
    density=sum(f,3);
    
    % compute velocity
    u=zeros(NY+1,NX+1,2);
    for k=1:9
        u(:,:,1)=u(:,:,1)+f(:,:,k)*c(k,1);
        u(:,:,2)=u(:,:,2)+f(:,:,k)*c(k,2);
    end
    u=u./density;
    
    % initialize the velocity of lid
    u(NY+1,2:NX,1)=U;
    u(NY+1,2:NX,2)=0;
    
    %video recording starts
    u_norm=sqrt(u(:,:,1).^2+u(:,:,2).^2);
    frame=getframe(gcf);
    frame.cdata=imresize(frame.cdata,[H W]);
    pcolor(x,y,u_norm);
    shading interp;
    hcb=colorbar;
    title(hcb,'speed (m/s)');
    axis equal;
    view([0,0,1]);
    writeVideo(objvideo,frame);
end
close(objvideo);