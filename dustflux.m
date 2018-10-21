function H = dustflux(U_mod,P,T,Dp)

%%ENVIRONMENTAL VARIABLES
g= 3.711; %average gravity acceleration
R=8.314472*1000/44; %Ideal gas constant considering a 100% CO2 atmosphere

%%DUST VARIABLES
gamma=3e-4; %empirical factor
An=0.0123; %empirical factor
rho_dust=2500; %particle density
%Dp= 180e-6; %particle diameter
z1=2; %sensor height
z0=0.01; %roughness height

%%FORMULATION
rho=P/(R*T); %Assuming ideal gas

U_drag= 0.4.*U_mod./log(z1/z0); %drag velocity using Mulholland expression

Ut_drag=sqrt(An.*(rho_dust./rho.*g.*Dp+gamma./(rho.*Dp))); %minimum drag velocity for saltation

U_ratio=Ut_drag./U_drag;

H = 2.61.*rho./g.*U_mod.^3.*(1-U_ratio).*(1+U_ratio).^2; %calculated dust flux

H(H<0)=0; %real dust flux





