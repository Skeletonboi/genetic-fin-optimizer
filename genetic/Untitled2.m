vp = load('altmachcg.mat') ;
atm = load('atmostable.mat');

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

rho = @(h) atm.atmosalt(round(h/100)+1,4);%1. Density
Ma = @(h) atm.atmosalt(round(h/100)+1,5); %2. Speed of Sound
mu = @(h) atm.atmosalt(round(h/100)+1,6); %3. Dynamic Viscosity of Air (~1.8e-5)
nu = @(h) (atm.atmosalt(round(h/100)+1,6))/(atm.atmosalt(round(h/100)+1,4)); % 4. Kinematic Viscosity (Dyn.Visc./Density)

value_rocket = {L_n,L_r,L_z,N,Cn_n,t,X_tc,L_red,D_noz,D_nos,D_end,F_w,F_fl}; 
value_atmo = {rho,Ma,mu,nu};

struct_rocket = struct("r",value_rocket);
struct_atmo = struct("a",value_atmo);

X = [0.3048,0.1524,0.1016,0.1016]

CP_barrow(X,vp,struct_rocket,struct_atmo)


function [cop] = CP_barrow(X,vp,str_r,str_a)
% CURRENTLY IMPLEMENTED AS BARROWMAN METHOD
	L_n = str_r(1).r;
	L_red = str_r(8).r;
	D_nos = str_r(10).r;
	D_end = str_r(11).r;
	L_r = str_r(2).r;
	N = str_r(4).r;
            
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
	x_fin =  l+ (X(3)*(X(1)+2*X(2)))/(3*(X(1)+X(2))) + (1/6)*(X(1)+X(2)-(X(1)*X(2))/(X(1)+X(2)));
	cop =  ((Cn_n*X_n)+(Cn_b*X_b)+(c*x_fin))/(Cn_n + Cn_b + c);       
end