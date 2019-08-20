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
% Initial Population Initialization
num_pop = 100;
ub = 0.75;
lb = 0;
init_set = (ub-lb).*rand(num_pop,5)+lb; % 4 independent variables
init_score = NaN(1,num_pop);
max1 = {0,0,[]};
max2 = {0,0,[]};
k_max = struct("km1",max1,"km2",max2);
% Evaluating fitness scores for initialization set
for i = 1:num_pop
    X = [init_set(i,1),init_set(i,2),init_set(i,3),init_set(i,4)];
    init_score(i) = fitness(X,vp,struct_rocket,struct_atmo);
    % Keeping two highest scoring sets
    if init_score(i) > k_max(2).km2
        if init_score(i) < k_max(2).km1
            k_max(2).km2 = init_score(i);
            k_max(1).km2 = i;
            k_max(3).km2 = X;
        else
            temp2 = k_max(2).km1;
            temp1 = k_max(1).km1;
            temp3 = k_max(3).km1;
            k_max(2).km2 = temp2;
            k_max(1).km2 = temp1;
            k_max(3).km2 = temp3;
            k_max(2).km1 = init_score(i);
            k_max(1).km1 = i;
            k_max(3).km1 = X;
        end
    end
end   
% Crossover the two best performing fin-dimension sets
parent_X = [mean([k_max(3).km1;k_max(3).km2])];
% Randomized mutation of parent dimensions (invoking exploration of domain
% space versus pure exploitation of good guesses.
%
% This step will require some fiddling (of hyperparameters such as method
% of mutation) to acquire proper balance of exploitation-to-exploration.
%
% POSSIBLE IMPROVEMENT AVENUES: Mutation("learning") momentumt (ie. in SGD
% ), variable mutation rate, etc; These methods can help move the
% convergence point out of local minima by initially overshooting with
% large learning/mutation rates. 






% Compute Fitness
% LOOP
% selection
% crossover
% mutation
% compute fitness
% UNTIL POP HAS CONVERGED
% STOP


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
% FITNESS FUNCTION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The goal of the fitness function is to reward fin dimension combinations
% that yield low drag, and a stability caliber that is greater than but
% close to the input target caliber.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRAINT-BASED STABILITY CALIBER FITNESS FUNC: 
% {    0; 0<x<2;
% {    1; 2<x;
% PREFFERED STABILITY CALIBER FITNESS FUNC:
% {    exp(50x-100); 0<x<=2;
% {    (-1/(1+exp(-0.8x+5)))+1+(1/(1+exp(-1.6+5))); 2<x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAG FITNESS FUNC(will need to tune for optimal fitness function):
% {    exp(0.0001x); 0<x<50'000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [score] = fitness(X,vp,str_r,str_a)
    D_nos = str_r(10).r;
    n = 6001;
    f_drag = exp(-0.00001*drag(X,vp,str_r,str_a));
    cop = CP_barrow(X,vp,str_r,str_a);
    for i = 1:n
        margin = (cop*100-vp.altmachcg(i,2))/(D_nos*100);
        if margin < 2
            f_stab = 0;
            break;
        elseif margin >= 2
            f_stab = 1;
        end
    end
    score = f_drag*f_stab;
end
function F = drag(X,vp,str_r,str_a)        
   % n = length(vp.altmachcg);
    n = 6001;
    Cd_out = NaN(1,n);
    Cd_incomp_out = NaN(1,n);
    for i = 1:n
        [Cd_out(i),Cd_incomp_out(i)] = drag_OR(X,vp.altmachcg(i,3),vp.altmachcg(i,1),str_r,str_a);
    end         
    F = num_int(Cd_out);
end
function a = num_int(vec)
    h = 1; %dx size
    a = 0;
    n = length(vec);
for i = 1:2:n-2
    da = (1/3)*h*(vec(i)+4*vec(i+1)+vec(i+2));
    a = a + da;
end

end
function [Cd,Cd_incomp] = drag_OR(X,M,alt,str_r,str_a)
    % Interference Drag and Fin-Tip Vortices are ignored as they are considered
    % "relatively small"
    % Also assumed a fully turbulent boundary layer

    % RETREIVE CONSTANTS FROM STRUCTS
    D_nos = str_r(10).r;
    F_w = str_r(12).r;
    N = str_r(4).r;
    L_r = str_r(2).r;
    Ma = str_a(2).a(alt);

    % Geometric Dimensions
    Afe = (1/2)*(X(1)+X(2))*X(4);
    Afp = (1/2)*D_nos*X(1) + Afe; % Fin Planform Area (+ Body)
    mid =  X(4)/cos(atan((X(3)+(X(2)/2)-(X(1)/2))/X(4))); % Midchord
    tr = X(2)/X(1);
    MAC = X(1)*(2/3)*((1+tr+tr^2)/(1+tr)); % Is actually the MGC approximated to be MAC


    %Abody = 
    Aref = N*(F_w*X(4));

    %REYNOLDS NUMBER CALCULATION
    V = M*Ma; %Instantaneous speed of rocket
    mu = str_a(3).a(alt);
    Rc = 51*((5*10^-6)/L_r)^-1.039;

    Re = (V*mid)/mu; % REYNOLDS NUMBER (MAC = characteristic length, no density)

    %Skin Friction Coefficients for Varying Regimes of Flow
    if Re < 10^4
        Cf = 1.48*(10^-2);
    elseif Re < Rc
        Cf = 1/((1.5*log(Re)-5.6)^2);
    elseif Re > Rcm 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

