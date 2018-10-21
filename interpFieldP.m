r1=100;
r2=50;
n=8;
t=0;
for i=1:n
    x1(i)=r1*cos(2*pi*i/n);
    x2(i)=r2*cos(2*pi*i/n+2*pi/2/n);
    y1(i)=r1*sin(2*pi*i/n);
    y2(i)=r2*sin(2*pi*i/n+2*pi/2/n);
end
[xp,yp]=meshgrid([x1,x2],[y1,y2]);
[Vp,Wp]=vfield(xp,yp);
[Vint]=griddata(xp,yp,Vp,x,y);
[Wint]=griddata(xp,yp,Wp,x,y);
quiver(x,y,Vint,Wint);

