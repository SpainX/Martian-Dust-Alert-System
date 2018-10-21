clear all, close all, clc
%VELOCITY FIELD INTERPOLATION USING SUBSTATION SENSORS
%this script is used to prove that a typical, not too complex wind velocity
%can be interpolated using the data from the substations without commiting
%a high amount of error
global t
r1=100;
r2=50;
n=8;
t=0;
while t<100;
for i=1:n
    x1(i)=r1*cos(2*pi*i/n);
    x2(i)=r2*cos(2*pi*i/n+2*pi/2/n);
    y1(i)=r1*sin(2*pi*i/n);
    y2(i)=r2*sin(2*pi*i/n+2*pi/2/n);
end
plot(x1,y1,'ob'), hold on
plot(x2,y2,'ob')
axis([-300 300 -300 300])
grid minor
P=[x1,x2;y1,y2];
[x,y]=meshgrid([-300:30:300],[-300:30:300]); %mesh of positions
[V,W]=vfield(x,y); %real wind velocity field, unknown to the substations
quiver(x,y,V,W,'c')
view(60,60)
[V1,W1]=vfield(x1,y1);
quiver(x1,y1,V1,W1,'r','LineWidth',2)
[V2,W2]=vfield(x2,y2);
quiver(x2,y2,V2,W2,'r','LineWidth',2);
xlabel('x  [km]'); ylabel('y [km]');
[xp,yp]=meshgrid([x1,x2],[y1,y2]);
[Vp,Wp]=vfield(xp,yp);
[x,y]=meshgrid([-150:30:150],[-150:30:150]);
%WIND INTERPOLATION based on what the substations see
[Vint]=griddata(xp,yp,Vp,x,y,'v4'); 
[Wint]=griddata(xp,yp,Wp,x,y,'v4');
quiver(x,y,Vint,Wint);
legend('Sensor Location - 1st Row','Sensor Location - 2nd Row','Real Wind','Data from Sensors','-','Interpolated Data')
% alpha 0.5
drawnow
t=t+1;
clf
end

% As it can be seen, the interpolated data of the velocity field is similar
% to the real wind fields. based on this, we can use this interpolation to
% compute the dust transport and to calculate and predict the path of the
% storms