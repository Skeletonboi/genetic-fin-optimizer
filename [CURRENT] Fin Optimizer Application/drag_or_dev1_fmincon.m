load('velocityprofile.mat','yeet')
%TEST = [0.1,0.1,0.1,0.1];
Cd_out = NaN(1,5800);
Cd_incomp_out = NaN(1,5800);

for i = 1:5800
    [Cd_out(i),Cd_incomp_out(i)] = drag(yeet(i));
end

hold on
plot(yeet,Cd_out)
plot(yeet,Cd_incomp_out)
legend("Comp","Incomp")
hold off

area = num_int(Cd_out)

function a = num_int(vec)
h = 1; %dx size
a = 0;
n = length(vec);
for i = 1:2:n-2
    da = (1/3)*h*(vec(i)+4*vec(i+1)+vec(i+2));
    a = a + da;
end

end


function [Cd,Cd_incomp] = drag(M)
X = [0.1,0.1,0.1,0.1];
% Interference Drag and Fin-Tip Vortices are ignored as they are considered
% "relatively small"
% Also assumed a fully turbulent boundary layer
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
Afe = (1/2)*(X(1)+X(2))*X(4);
Afp = (1/2)*D_nos*X(1) + Afe; % Fin Planform Area (+ Body)
mid =  X(4)/cos(atan((X(3)+(X(2)/2)-(X(1)/2))/X(4))); % Midchord
tr = X(2)/X(1);
MAC = X(1)*(2/3)*((1+tr+tr^2)/(1+tr));
%Abody = 
Aref = N*(F_w*X(4));

%Reynolds Number
rho = 1.225 ;%sea level, 15degC (ISA) [kg/m^3]

Ma = 340.3 ;%sea level, 15degC [m/s]
V = M*Ma;
mu = 1.802*(10^-5); % 15degC [kg/m*s]
nu = 1.470*(10^-5); %15degC
Rc = 51*((5*10^-6)/L_r)^-1.039;
Re = (V*mid)/mu; % REYNOLDS NUMBER (MAC = characteristic length, no density)
%Skin Friction Coefficients for Varying Regimes of Flow
if Re < 10^4
    Cf = 1.48*(10^-2);
elseif Re < Rc
    Cf = 1/((1.5*log(Re)-5.6)^2);
elseif Re > Rc
    Cf = 0.032*((5*10^-6)/L_r)^0.2;
end
%Compressibility Corrections
if V < Ma
    Cfc = Cf*(1-0.1*M^2);
elseif V > Ma
    turb = Cf/((1+0.15*M^2)^0.58);
    rough = Cf/(1+0.18*M^2);
    if Re < Rc
        Cfc = turb;
    elseif Re > Rc
        if rough < turb
            Cfc = turb;
        elseif rough > turb
            Cfc = rough;
        end
    end
end

%Skin Friction Drag Coefficient
Cdf = (Cfc*((1+(2*F_w)/MAC)*6*Afe))/Aref;

% Fin Pressure Drag Coefficient
% firstly, the leading edge pressure drag (assuming rounded leading edge):
if M < 0.9
    Cd_le = ((1-M^2)^-0.417) - 1;
elseif M < 1
    Cd_le = 1- 1.785*(M-0.9);
elseif M >= 1
    Cd_le = 1.214 - (0.502/M^2) + (0.1095/M^4);
end
% now correcting it for sweep:
Cd_le_sweep = Cd_le*(cos(atan(X(3)/X(4))^2));
% Cd_te (or trailing edge) is 0 due to tapered geometry. The resulting Cdp
% (fin pressure drag) is thus a scaled value of Cd_le_sweep.

Cdp = (Cd_le_sweep*6*Afe)/Aref;

Cd_incomp = Cdp+Cdf;

% Compressibility Correction (extending to compressible flow regime)
% Prandtl-Glauert compressibility correction
if M < 0.8
    Cd = Cd_incomp/(sqrt(1-M^2));
elseif M > 1.1
    Cd = Cd_incomp/(sqrt(M^2-1));
elseif 0.8 <= M <= 1.1
    Cd = Cd_incomp/(sqrt(1-0.8^2)); %Ketchledge [1993] correction to avoid approaching infinity as M -> 1
end


% NOTE: Transonic aerodynamic properties are accounted for "by using some
% suitable interpolation function" [OR,pg.18] Not sure where to find this. OR
% also does not take into account Hypersonic (Ma > 5)

end

