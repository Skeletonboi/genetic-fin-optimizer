load('velocityprofile.mat','yeet')
yeet = [VelocityProfile.Mach];
OUT = NaN(1,5800);

for i = 1:5800
    OUT(i) = drag_ac(yeet(i));
end

plot(yeet,OUT)

function Cd = drag_ac(M)
X = [0.1,0.1,0.1,0.1];
% Constants (all dimensions in meters):
L_n = 0.548;    % Nose Length
F_r = 0.406;    % Fin Root Chord
F_t = 0.127;    % Fin Tip Chord
F_s = 0.152;    % Fin Semi-Span
S = 0.100;      % Sweep Distance
L_r = 5.844;    % Rocket Length
L_z = 0.0762;   % Nozzle Length
N = 3;          % # of fins
Cn_n = 0.5;     % Nose Cone Coefficient
t = 0.00635;    % Max Fin Root Thickness
X_tc = 0.00635; % Distance from Fin Leading Edge to Max Thickness
L_red = 0.05995;% Length Reduction @ back
D_noz = 0.08204;% Nozzle Diameter
D_nos = 0.140;  % Nose Base Diameter
D_end = 0.1077; % End Diameter
F_w = 0.00381;  % Fin Thickness
F_fl = 0.21895; % Fin Front Length
Cg = 3.372354308; %Centre of Gravity 

% Geometric Dimensions
Afe = (1/2)*(X(1)+X(2))*X(4); % Fin Exposed Area
Afp = (1/2)*D_nos*X(1) + Afe; % Fin Planform Area
mid =  X(4)/cos(atan((X(3)+(X(2)/2)-(X(1)/2))/X(4))); % Midchord

%Reynolds Number
rho = 1.225; %sea level, 15degC (ISA) [kg/m^3]
Ma = 340.3; %sea level, 15degC [m/s]

V = M*Ma;
mu = 1.802*(10^-5); % 15degC [kg/m*s]
nu = 1.470*(10^-5);
Re = (rho*V*mid)/nu; % 1,072,936.44 > Re_c = 5x10^5
Re_c = 5*(10^5);

B = Re_c*((0.074/(Re^(1/5)))-(1.328/sqrt(Re)));
if Re > Re_c
    Cff = (0.074/(Re^(1/5)))-(B/Re); % Viscous Friction Force Coefficient
elseif Re <= Re_c
    Cff = 1.328/(sqrt(Re));
end
% Fin Drag (scales with 1st deg wrt Cff)
Cdf = 2*Cff*(1+2*(F_w/mid))*((4*N*Afp)/(pi*D_nos^2))

% Fin-Body Interference Drag (scales with 1st deg wrt Cff)
Cdi = 2*Cff*(1+2*(F_w/mid))*((4*N*(Afp-Afe))/(pi*D_nos^2))



% AOA Correction Drag Coefficient
alpha = 2; % Angle of Attack [deg]
delta = 0.825; % Constant interpolated from Fig 4. [AerodynamicCoefficients.pdf]
nu = 0.6125; % Constant interpolated from Fig 4. [AerodynamicCoefficients.pdf]
L_ts = 2*X(4) + D_nos; % Total Span
Rs = L_ts/D_nos; % Fin Section Ratio

kfb = 0.8065*Rs^2 + 1.1553*Rs % Fin-body intereference coefficient
kbf = 0.1935*Rs^2 + 0.8174*Rs + 1 % Body-fin interference coefficient

% calculation of AOA correction coefficient:
Cdf_alpha =  (alpha^2)*(((1.2*(Afp*4))/(pi*D_nos^2))+3.12*(kfb+kbf-1)*((Afe*4)/(pi*D_nos^2)))

Cd_incomp = Cdf + Cdi + Cdf_alpha;

% Compressibility Correction (extending to compressible flow regime)
% Prandtl-Glauert compressibility correction
if M < 0.8
    Cd = Cd_incomp/(sqrt(1-M^2));
elseif M > 1.1
    Cd = Cd_incomp/(sqrt(M^2-1));
elseif 0.8 <= M <= 1.1
    Cd = Cd_incomp/(sqrt(1-0.8^2)); %Ketchledge [1993] correction to avoid approaching infinity as M -> 1
end



    

end

%test = linspace(0,5,1000)
%plot(test,Cdf_alpha(test))