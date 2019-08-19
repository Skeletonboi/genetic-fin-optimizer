%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENETIC ALGORITHM OPTIMIZER FOR FIN-DIMENSIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an ordinary genetic algorithm optimizer that will attempt to find
% the most optimal fin-dimensions in terms of least drag while satisfying
% the stability caliber constraint (across the entire given velocity
% profile). 
%
% The purpose of this optimizer is to assist in the validation of the
% results attained from MATLAB's built-in "fmincon" optimizer. The current
% iteration of "fmincon" uses an interior-point optimization method which
% cannot reliably produce consistent solutions (ie. does not converge on a
% globally optimal solution) in non-convex functions. Due to the
% complicated nature of the objective (drag) function, its convexity has
% not been determined at this time (Aug 19th 2019).
% 
% As such, this genetic algorithm aims to provide a certain degree of
% randomization in the search area, which - when combined with multiple
% variable initializations - will improve the confidence interval of a
% mutually agreed optimal solution.
%
% The optimizer must be given 1) rocket constants, 2) a target
% stability caliber constraint, and optionally, 3) a manually-imposed 
% inter-fin-dimensional constraint.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

value_rocket = {L_n,L_r,L_z,N,Cn_n,t,X_tc,L_red,D_noz,D_nos,D_end,F_w,F_fl}; 
value_atmo = {rho,Ma,mu,nu};

struct_rocket = struct("r",value_rocket);
struct_atmo = struct("a",value_atmo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Population Initialization
num_pop = 100;
init_set = rand(num_pop,5); % 4 independent variables
init_score = NaN(1,num_pop);
% Evaluating fitness scores for initialization set
for i = 1:num_pop
    init_score(i) = fitness([init_set(i,1),init_set(i,2),init_set(i,3),init_set(i,4)]);
end


% FITNESS FUNCTION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The goal of the fitness function is to reward fin dimension combinations
% that yield low drag, and a stability caliber that is greater than but
% close to the input target caliber.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [score] = fitness(X)
    

end
% Compute Fitness
% LOOP
% selection
% crossover
% mutation
% compute fitness
% UNTIL POP HAS CONVERGED
% STOP


