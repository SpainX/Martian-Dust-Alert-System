function [V,W]=vfield(x,y)
%DUMMY T-DEPENDANT VELOCITY FIELDS FOR INTERPOLATION
global t
% V=1+x./sqrt(x.^2+y.^2)*sin(2*pi*t/100); W=y./sqrt(x.^2+y.^2)*cos(2*pi*t/100);
V=t/100-y./sqrt(x.^2+y.^2)*sin(2*pi*t/100)+(t/100); W=x./sqrt(x.^2+y.^2)*cos(2*pi*t/100);
% V=sin(2*pi*t)*100*ones(size(x,1)); W=0;
% V=100.*x./x;W=0.*x./x;
end