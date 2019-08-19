% GENETIC ALGORITHM OPTIMIZER FOR FIN-DIMENSIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an ordinary genetic algorithm optimizer that will attempt to find
% the most optimal fin-dimensions in terms of least drag while satisfying
% the stability caliber constraint (across the entire given velocity
% profile). 
%
% The purpose of this optimizer is to provide secondary 
%
% The optimizer must be given 1) rocket constants, and optionally,
% 2) manually-imposed inter-fin-dimensional constraints, and/or 2) a target
% stability caliber constraint.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rocket/Atmospheric Constants Initialization
vp = load('altmachcg.mat') ;
atm = load('atmostable.mat');

L_n = 1.1;    %1. Nose Length
L_r = 8.0;    %2. Rocket Length
L_z = 0.30;   %3. Nozzle Length
N = 4;          %4. # of fins
Cn_n = 0.5;     %5. Nose Cone Coefficient
t = 0.00635;    %6. Max Fin Root Thickness
X_tc = 0.00635; %7. Distance from Fin Leading Edge to Max Thickness
L_red = 0.08;%8. Length Reduction @ back
D_noz = 0.10;%9. Nozzle Diameter
D_nos = 0.30;  %10. Nose Base Diameter
D_end = 0.24; %11. End Diameter
F_w = 0.006;  %12. Fin Width (Thickness)
F_fl = 0.21895; %13. Fin Front Length


rho = @(h) atm.atmosalt(round(h/100)+1,4);%1. Density
Ma = @(h) atm.atmosalt(round(h/100)+1,5); %2. Speed of Sound
mu = @(h) atm.atmosalt(round(h/100)+1,6); %3. Dynamic Viscosity of Air (~1.8e-5)
nu = @(h) (atm.atmosalt(round(h/100)+1,6))/(atm.atmosalt(round(h/100)+1,4)); % 4. Kinematic Viscosity (Dyn.Visc./Density)

value_rocket = {L_n,L_r,L_z,N,Cn_n,t,X_tc,L_red,D_noz,D_nos,D_end,F_w,F_fl}; % <- TBD: MOVING CG WRT ALTITUDE/MACH
value_atmo = {rho,Ma,mu,nu};

struct_rocket = struct("r",value_rocket);
struct_atmo = struct("a",value_atmo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


