 %{
X0 = [0.1;0.1;0.1;0.1]; %Initial Conditions for [a,b,m,s]
            A = [-0.75,1,0,0;
            0,-0.9,1,0;
            -1,0,0,1];
            B = [0;0;0];
            Aeq = []; 
            Beq = [];
            LB = [0;0;0;0];
            UB = [0.1524;Inf;Inf;Inf];

            X = fmincon(@(X) DRAG_AC(X),X0,A,B,Aeq,Beq,LB,UB,@(X) NONLCON(X))
            % FUN = Fin Drag Equation (Subsonic)

            % M = CP_FINAL(app,X);
            DRAG_AC(X)
            X
            
            
            graph_vec = [0,0;X(3),X(4);(X(3)+X(2)),X(4);X(1),0];
%}
X = [0.1194,0.0326,0.0294,0.1194]
DRAG_AC(X)


function F = DRAG_AC(X)
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
ma_val = 0.8;
V = ma_val*Ma;
mu = 1.802*(10^-5); % 15degC [kg/m*s]
Re = (rho*V*mid)/mu; % 1,072,936.44 > Re_c = 5x10^5
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

% Angle of Attack Correction Drag Coefficient
alpha = 5; % Angle of Attack [deg]
delta = 0.825; % Constant interpolated from Fig 4. [AerodynamicCoefficients.pdf]
nu = 0.6125; % Constant interpolated from Fig 4. [AerodynamicCoefficients.pdf]
L_ts = 2*X(4) + D_nos; % Total Span
Rs = L_ts/D_nos; % Fin Section Ratio

kfb = 0.8065*Rs^2 + 1.1553*Rs % Fin-body intereference coefficient
kbf = 0.1935*Rs^2 + 0.8174*Rs + 1 % Body-fin interference coefficient

%CORRECTION DRAG FACTOR FOR ANGLE OF ATTACK
Cdf_alpha = (alpha^2)*(((1.2*(Afp*4))/(pi*D_nos^2))+3.12*(kfb+kbf-1)*((Afe*4)/(pi*D_nos^2)))

% TOTAL DRAG COEFFICIENT W/O COMPRESSIBLE FLOW CORRECTION
C_uncompress = Cdf+Cdi+Cdf_alpha; 
%{
% Compressible Flow Correction
if ma_val < 0.8
    C_compress = C_uncompress/sqrt(1-(ma_val)^2);
elseif ma_val > 1.1
    
    C_compress = C_uncompress/sqrt(((ma_val)^2)-1);
else
    C_compress = C_uncompress/sqrt(1-0.8^2);
end

F = C_compress;
%}
F = C_uncompress; % TO BE DELETED AFTER!!! ---------------
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
            