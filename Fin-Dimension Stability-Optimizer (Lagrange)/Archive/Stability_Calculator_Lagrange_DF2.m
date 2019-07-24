%Stability Calculator (by Victor Yip)

% Constants (all dimensions in meters):
L_n = 0.548;    % Nose Length
F_r = 0.150;    % Fin Root Chord
F_t = 0.033;    % Fin Tip Chord
F_s = 0.12192;  % Fin Semi-Span
S = 0.084;      % Sweep Distance
L_r = 2.65575;  % Rocket Length
L_z = 0.0762;   % Nozzle Length
N = 3;          % # of fins
Cn_n = 0.5;     % Nose Cone Coefficient
t = 0.00635;    % Max Fin Root Thickness
X_tc = 0.00635; % Distance from Fin Leading Edge to Max Thickness
L_red = 0.05995;% Length Reduction @ back
D_noz = 0.08204;% Nozzle Diameter
D_nos = 0.140;  % Nose Base Diameter
D_end = 0.1077; % End Diameter
F_w = 0.00381;    % Fin Width
F_fl = 0.21895; % Fin Front Length

% Calculated Dimensions
F_m = F_s/cos(atan((S+(F_t/2)-(F_r/2))/F_s)); % Fin Mid-Chord Line
L_nf = L_r-F_r-L_red; % Nose Tip to Fin Chord Leading Edge

% X (distance from nose tip to component's Cp):
% NOTE: Normal Body and Shoulder Forces are neglected 
X_n = L_n/2;
X_b = L_r - L_red + (L_red/3)*(1+1/(1+D_nos/D_end));
X_fb = L_nf+ (S*(F_r+2*F_t))/(3*(F_r+F_t)) + (1/6)*(F_r+F_t-(F_r*F_t)/(F_r+F_t));

% Component Coefficient of Normal Force (Cn)
Cn_n = 2;
Cn_b = 2*(((D_end/D_nos)^2)-1);
Cn_f = (4*N*(F_s/D_nos)^2)/(1+sqrt(1+(2*F_m/(F_r+F_t))^2));
int = 1+ (D_nos/2)/(F_s+(D_nos/2));
Cn_fb = int*Cn_f;

% Net Coefficient of Normal Force
Cn_net = Cn_n + Cn_b + Cn_fb;
% Net Cp Distance from Nose Tip 
X_net = ((Cn_n*X_n)+(Cn_b*X_b)+(Cn_fb*X_fb))/Cn_net;

% LAGRANGE MULTIPLIER OPTIMIZATION
% Variables:
syms a b m s
% Using default values for all other variables,
% we can cast final Cp position as cp(a,b,m,s):
mid = @(a,b,m,s) s/cos(atan((m+(b/2)-(a/2))/s));
l = @(a) (L_r-a)-L_red;
int = @(s) 1+ (D_nos/2)/(s+(D_nos/2));
c = @(a,b,m,s) int(s)*(4*N*(s/D_nos)^2)/(1+sqrt(1+(2*mid(a,b,m,s)/(a+b))^2));
x = @(a,b,m,s) l(a)+ (m*(a+2*b))/(3*(a+b)) + (1/6)*(a+b-(a*b)/(a+b));
cp = @(a,b,m,s) ((Cn_n*X_n)+(Cn_b*X_b)+(c(a,b,m,s)*x(a,b,m,s)))/(Cn_n + Cn_b + c(a,b,m,s));
% To verify correctness, the following paramaters should = 1.7196:
cp(0.15,0.033,0.084,0.122) 
%We now take the gradient of the Cp function (which we want 
% to optimize) by first finding its partial derivatives w.r.t. to 
% each independent variable:
cp_a = diff(cp,a)
cp_b = diff(cp,b)
cp_m = diff(cp,m)
cp_s = diff(cp,s)
% We now define the constraint equation:
g = @(a,b) b-(a/2) %Tip chord < 1/2*Root Chord
g_a = diff(g,a)
g_b = diff(g,b)
