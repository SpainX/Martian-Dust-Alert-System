%% DUST CLOUD FORMATION AND TRANSPORT EQUATION
% Based on the interpolated velocity, pressure and temperature data,
% this script computes the movement and formation on dust storm, based
% on the discretization of the base surroundings and the saltation
% conditions

clear
clc, close all
r1=100; %far station ring radius(km)
r2=50; %close station ring radius (km)
n=8; %number of stations per ring
t=0;
tsim=500; %simulation time
dt=0.1;

rng= 150; %range in km
res= 10; %mesh resolution in km
%%MESH DEFINITION
[Xm, Ym]= meshgrid(-rng:res:rng,-rng:res:rng);
for i=1:n %station definition
    x1(i)=r1*cos(2*pi*i/n);
    x2(i)=r2*cos(2*pi*i/n+2*pi/2/n);
    y1(i)=r1*sin(2*pi*i/n);
    y2(i)=r2*sin(2*pi*i/n+2*pi/2/n);
end
i=1;
%%DATA INPUT
[Vx,Vy]=vfield([x1,x2],[y1,y2]);
P=800.*rand(1,2*n);
T=278.*rand(1,2*n);
X=[x1' y1';x2' y2'];
% clear x1 y1 x2 y2
stations= [X Vx' Vy' P' T'];
H=zeros(size(Xm,1));
G=zeros(size(Xm,1));

%% MESHED VARIABLES
[Vxm]=griddata(X(:,1),X(:,2),Vx',Xm,Ym,'v4');
[Vym]=griddata(X(:,1),X(:,2),Vy',Xm,Ym,'v4');
[Pm]=griddata(X(:,1),X(:,2),P',Xm,Ym,'v4');
[Tm]=griddata(X(:,1),X(:,2),T',Xm,Ym,'v4');

%% Dust volume boundary conditions
G(1,:)=0;
G(:,1)=0;
G(end,:)=0;
G(:,end)=0;
G0=zeros(size(Xm,1));

%sim loop
while t<=tsim
    [Vx,Vy]=vfield([x1,x2],[y1,y2]);
    [Vxm]=griddata(X(:,1),X(:,2),Vx',Xm,Ym,'v4');
    [Vym]=griddata(X(:,1),X(:,2),Vy',Xm,Ym,'v4');
    [Pm]=griddata(X(:,1),X(:,2),P',Xm,Ym,'v4');
    [Tm]=griddata(X(:,1),X(:,2),T',Xm,Ym,'v4');
    [Hm]=dustflux(norm([Vxm Vym]),Pm,Tm,150e-6);
    for i=2:size(Xm,1)-1
        for j=2:size(Xm,1)-1
            
            G(i,j)=G0(i,j)+Hm(i,j)*res*dt+G0(i+1,j)*dot([Vxm(i+1,j),Vym(i+1,j)],[0,1])./norm([Vxm(i+1,j),Vym(i+1,j)])+G0(i-1,j)*dot([Vxm(i-1,j),Vym(i-1,j)],[0,-1])./norm([Vxm(i-1,j),Vym(i-1,j)])+G0(i,j+1)*dot([Vxm(i,j+1),Vym(i,j+1)],[-1,0])./norm([Vxm(i,j+1),Vym(i,j+1)])+G0(i,j-1)*dot([Vxm(i,j-1),Vym(i,j-1)],[1,0])./norm([Vxm(i,j-1),Vym(i,j-1)]);
        end
    end
    G0=G; %initialization of cariable G
    t=t+dt; %time advance 
    % post processing
    pcolor(Xm,Ym,G);
    colormap(flipud(gray)) , hold on
    scatter(X(:,1),X(:,2),'c')
    plot(0,0,'rX')
    xlabel('x [km]')
    ylabel('y [km]')
    legend('storms','stations','base')
    drawnow
    clf
end

hold on
%mesh(Xm,Ym,zeros(31,31)-0.1)
% quiver(stations(:,1),stations(:,2),Vx',Vy')
% quiver(Xm,Ym,Vxm,Vym)
%mesh(Xm,Ym,zeros(31,31))