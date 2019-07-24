
X0 = [0.1;0.1;0.1;0.1]; %Initial Conditions for [a,b,m,s]
A = [-0.75,1,0,0;
    0,-0.9,1,0;
    -1,0,0,1];
B = [0;0;0];
Aeq = [];
Beq = [];
LB = [0;0;0;0];
UB = [0.1524;Inf;Inf;Inf];

X = fmincon(@(X) FUN(X),X0,A,B,Aeq,Beq,LB,UB,@(X) NONLCON(X))
% FUN = Fin Drag Equation (Subsonic)

M = CP_FINAL(X)
Cdf_NET = DRAG(X)


function F = FUN(X)
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

%Reynolds Number
rho = 1.225 %sea level, 15degC (ISA) [kg/m^3]
Ma = 340.3 %sea level, 15degC [m/s]
V = 0.9*Ma
mu = 1.802*(10^-5); % 15degC [kg/m*s]
Re = (rho*V*mid)/mu; % 1,072,936.44 > Re_c = 5x10^5

B = (5*(10^5))*((0.074/(Re^(1/5)))-(1.328/sqrt(Re)));
Cff = (0.074/((Re)^(1/5)))-(B/Re);
Cdf = 2*Cff*(1+2*(F_w/mid))*(4*3*Afp/(pi*(D_nos^2))) 

Cdi = 2*Cff*(1+2*(F_w/mid))*(4*3*(Afp-Afe)/(pi*(D_nos^2)))

F = Cdf + Cdi;
end


function [margin_final] = CP_FINAL(X)
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
F_w = 0.00381;  % Fin Width (Thickness)
F_fl = 0.21895; % Fin Front Length
Cg = 3.372354308; %Centre of Gravity 

% X (distance from nose tip to component's Cp):
% NOTE: Normal Body and Shoulder Forces are neglected 
X_n = L_n/2;
X_b = L_r - L_red + (L_red/3)*(1+1/(1+D_nos/D_end));

% Component Coefficient of Normal Force (Cn)
Cn_n = 2;
Cn_b = 2*(((D_end/D_nos)^2)-1);

%Optimization Variables = X(1) = a, X(2) = b, X(3) = m, X(4) = s

% Using default values for all other variables,
% % we can cast final Cp position as cp(a,b,m,s):
mid =  X(4)/cos(atan((X(3)+(X(2)/2)-(X(1)/2))/X(4)));
l =  (L_r-X(1))-L_red;
int = 1+ (D_nos/2)/(X(4)+(D_nos/2));
c =  int*(4*N*(X(4)/D_nos)^2)/(1+sqrt(1+(2*mid/(X(1)+X(2)))^2));
x =  l+ (X(3)*(X(1)+2*X(2)))/(3*(X(1)+X(2))) + (1/6)*(X(1)+X(2)-(X(1)*X(2))/(X(1)+X(2)));
cp =  ((Cn_n*X_n)+(Cn_b*X_b)+(c*x))/(Cn_n + Cn_b + c);

margin_final = (cp-Cg)/D_nos;

end

function [C,Ceq] = NONLCON(X)
C = [];
% Constants (all dimensions in meters):
L_n = 0.548;    % Nose Length
F_r = 0.406;    % Fin Root Chord
F_t = 0.127;    % Fin Tip Chord
F_s = 0.152;  % Fin Semi-Span
S = 0.100;      % Sweep Distance
L_r = 5.844;  % Rocket Length
L_z = 0.0762;   % Nozzle Length
N = 3;          % # of fins
Cn_n = 0.5;     % Nose Cone Coefficient
t = 0.00635;    % Max Fin Root Thickness
X_tc = 0.00635; % Distance from Fin Leading Edge to Max Thickness
L_red = 0.05995;% Length Reduction @ back
D_noz = 0.08204;% Nozzle Diameter
D_nos = 0.140;  % Nose Base Diameter
D_end = 0.1077; % End Diameter
F_w = 0.00381;    % Fin Width (Thickness)
F_fl = 0.21895; % Fin Front Length
Cg = 3.372354308; %Centre of Gravity 

% X (distance from nose tip to component's Cp):
% NOTE: Normal Body and Shoulder Forces are neglected 
X_n = L_n/2;
X_b = L_r - L_red + (L_red/3)*(1+1/(1+D_nos/D_end));

% Component Coefficient of Normal Force (Cn)
Cn_n = 2;
Cn_b = 2*(((D_end/D_nos)^2)-1);

%Optimization Variables = X(1) = a, X(2) = b, X(3) = m, X(4) = s

% Using default values for all other variables,
% % we can cast final Cp position as cp(a,b,m,s):
mid =  X(4)/cos(atan((X(3)+(X(2)/2)-(X(1)/2))/X(4)));
l =  (L_r-X(1))-L_red;
int = 1+ (D_nos/2)/(X(4)+(D_nos/2));
c =  int*(4*N*(X(4)/D_nos)^2)/(1+sqrt(1+(2*mid/(X(1)+X(2)))^2));
x =  l+ (X(3)*(X(1)+2*X(2)))/(3*(X(1)+X(2))) + (1/6)*(X(1)+X(2)-(X(1)*X(2))/(X(1)+X(2)));
cp =  ((Cn_n*X_n)+(Cn_b*X_b)+(c*x))/(Cn_n + Cn_b + c);

margin = (cp-Cg)/D_nos;
Ceq = [margin-2];

end


function drag = DRAG(X)
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
Afe = (1/2)*(X(1)+X(2))*X(4)
Afp = (1/2)*D_nos*X(1) + Afe; % Fin Planform Area (+ Body)
mid =  X(4)/cos(atan((X(3)+(X(2)/2)-(X(1)/2))/X(4))); % Midchord
tr = X(2)/X(1);
MAC = X(1)*(2/3)*((1+tr+tr^2)/(1+tr));
Abody = 
Aref = 

%Reynolds Number
rho = 1.225 ;%sea level, 15degC (ISA) [kg/m^3]
M = 0.3;
Ma = 340.3 ;%sea level, 15degC [m/s]
V = M*Ma;
mu = 1.802*(10^-5); % 15degC [kg/m*s]
Rc = 51*((5*10^-6)/mid)^-1.039;
Re = (rho*V*mid)/mu; % 1,072,936.44 > Re_c = 5x10^5

if Re < 10^4:
    Cf = 1.48*10^-2;
elseif Re < Rc:
    Cf = 1/(1.5*log(Re)-5.6)^2;
elseif Re > Rc:
    Cf = 0.032*((5*10^-6)/mid)^0.2
end

if V < Ma:
    Cfc = Cf*(1-0.1*M^2)
elseif V > Ma:
    turb = Cf/((1+0.15*M^2)^0.58);
    rough = Cf/(1+0.18*M^2);
    if Re < Rc:
        Cfc = turb;
    elseif Re > Rc:
        if rough < turb:
            Cfc = turb;
        elseif rough > turb:
            Cfc = rough;
        end
    end
end

CDf = (Cfc*((1+(1/(2*(L_r/D_nos))))*Abody + (1+(2*F_w)/MAC)*6*Afe))/Aref



end
