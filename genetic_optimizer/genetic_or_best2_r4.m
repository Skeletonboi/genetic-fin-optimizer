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
% repeatedly found optimal solution.
%
% The optimizer must be given 1) rocket constants, 2) a target
% stability caliber constraint, and optionally, 3) a manually-imposed 
% fin-dimension constraints.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rocket/Atmospheric Constants Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vp = load('altmachcg.mat') ;
atm = load('atmostable.mat');

L_n = 0.5588;    %1. Nose Length
L_r = 5.95122;    %2. Rocket Length
L_z = 0.07366;   %3. Nozzle Length
N = 4;          %4. # of fins
Cn_n = 0.5;     %5. Nose Cone Coefficient
t = 0.00635;    %6. Max Fin Root Thickness
X_tc = 0.00635; %7. Distance from Fin Leading Edge to Max Thickness
L_red = 0.08;%8. Length Reduction @ back
D_noz = 0.062865*2;%9. Nozzle Diameter
D_nos = 0.06985*2;  %10. Nose Base Diameter
D_end = 0.062865*2; %11. End Diameter
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


for n = 1:1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Population Initialization & Evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_pop = 100;
aub = toMeters(12);
bub = toMeters(12);
mub = toMeters(12);
sub = toMeters(6);
alb = toMeters(0);
blb = toMeters(0);
mlb = toMeters(0);
slb = toMeters(0);
% Initialize pop. space of 4 fin dims.
init_set = NaN(num_pop,4); 
% Initializing each specimen's dimensions randomly (uniform dist.)
% according to relative upper & lower bounds
init_set(:,1) = (aub-alb)*rand(num_pop,1)+alb;
init_set(:,2) = (bub-blb)*rand(num_pop,1)+blb;
init_set(:,3) = (mub-mlb)*rand(num_pop,1)+mlb;
init_set(:,4) = (sub-slb)*rand(num_pop,1)+slb;
%init_set = (aub-alb).*rand(num_pop,4)+alb;
% Initialize space for population scores
init_score = NaN(1,num_pop); 
% Initialize space for best 2 specimens
km1_i = 0;
km1_s = 0;
km1_dim = [0,0,0,0];
km2_i = 0;
km2_s = 0;
km2_dim = [0,0,0,0];
% Evaluating fitness scores for the population
for i = 1:num_pop
    X = [init_set(i,1),init_set(i,2),init_set(i,3),init_set(i,4)];
    init_score(i) = fitness(X,vp,struct_rocket,struct_atmo);
    % Keeping two highest scoring specimens
    if init_score(i) > km2_s
        if init_score(i) <= km1_s
            km2_s = init_score(i);
            km2_i = i;
            km2_dim = X;
        else
            tempi = km1_i;
            temps = km1_s;
            tempdim = km1_dim;
            km2_i = tempi;
            km2_s = temps;
            km2_dim = tempdim;
            km1_s = init_score(i);
            km1_i = i;
            km1_dim = X;
        end
    end
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [INPUT] Selection process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose a method of selection:
% (1) Best mean-of-two
%parent_X = [mean([k_max(3).km1;k_max(3).km2])];
% (2) Probabilistic best (selection prob. directly proportional set fitness
% score)
%ind = p_select(init_score);
%parent_X = [init_set(ind,1),init_set(ind,2),init_set(ind,3),init_set(ind,4)];
% (3) Probabilistic best mean-of-two (selection prob. directly proportional
% to set fitness score)

% ind1 = km1_i;
% ind2 = km2_i;
% mom = [init_set(ind1,1),init_set(ind1,2),init_set(ind1,3),init_set(ind1,4)];
% dad = [init_set(ind2,1),init_set(ind2,2),init_set(ind2,3),init_set(ind2,4)];
mom = km1_dim;
dad = km2_dim;
par = [mean([mom;dad])];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [INPUT] Evolution Process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose number of evolutions to perform
num_evo = 100;

next = true;
p_par = zeros(1,4);
cei = 0;
% Looping over # of evolutions
%for e = 1:num_evo
while next
    cei = cei + 1;
    fprintf('Current Evolution: %d \n',cei);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mutation Process
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create randomized mutation population of parent dimensions (invoking 
    % exploration of domain space versus pure exploitation of good guesses)
    % 
    % This step will require some fiddling (of hyperparameters such as method
    % of mutation) to acquire proper balance of exploitation-to-exploration.
    % 
    % POSSIBLE IMPROVEMENT AVENUES: Mutation("learning") momentum (ie. in SGD
    % ), variable mutation rate, etc; These methods can help move the
    % convergence point out of local minima by initially overshooting them with
    % large learning/mutation rates to hopefully reach the global minimum.
    %
    % The chosen mutation algorithm chooses two random independent variables 
    %(out of the 4) to mutate. The mutation is accomplished by adding or
    % subtracting the existing value by a random variable whose value is 
    % bounded by the existing parent value and the upper/lower bound. The 
    % probability distribution of the random variable is an exponential
    % distribution; this distribution (hyperparameter) was chosen to be 
    % exponential as it results in a greater probability of the mutation to 
    % result with a value near the existing parent value. 
    % Again, this is decision is a hyperparameter, and can be tuned.

    % i) Initialize space for mutation population
    mut_pop = zeros(num_pop,4);
    mut_score = zeros(num_pop,1);

    km1_i = 0;
    km1_s = 0;
    km1_dim = [0,0,0,0];
    km2_i = 0;
    km2_s = 0;
    km2_dim = [0,0,0,0];
    % mut_pop = gpuArray(mut_pop_cpu);
    % mut_score = gpuArray(mut_score_cpu);

    % ii) Compute difference between existing parent values and UB/LB
    dau = aub-par(1);
    dal = par(1)-alb;
    dbu = bub-par(2);
    dbl = par(2)-blb;
    dmu = mub-par(3);
    dml = par(3)-mlb;
    dsu = sub-par(4);
    dsl = par(4)-slb;
    dmat = [dau,dbu,dmu,dsu;dal,dbl,dml,dsl];

    % iii) Compute the mean value of probability distribution for each diff.
    mmat = -dmat./log(0.0001);

    parfor p = 1:num_pop
        % 1. Choose whether mutation is additive or subtractive (50/50 split)
        if rand > 0.5
            dir1 = 1;
            sign1 = 1;
        else
            dir1 = 2;
            sign1 = -1;
        end
        if rand > 0.5
            dir2 = 1;
            sign2 = 1;
        else
            dir2 = 2;
            sign2 = -1;
        end
        if rand > 0.5
            dir3 = 1;
            sign3 = 1;
        else
            dir3 = 2;
            sign3 = -1;
        end
        if rand > 0.5
            dir4 = 1;
            sign4 = 1;
        else
            dir4 = 2;
            sign4 = -1;
        end

        % 2. Choose which two of the design variables to mutate
    %      d_dist = [0.25,0.25,0.25,0.25];
    %      i = p_select(d_dist);
    %      j = p_select(d_dist);

        % 3. Generate random pertubation according to bounded exponential dist.
        pert_a = exprnd(mmat(dir1,1));
        pert_b = exprnd(mmat(dir2,2));
        pert_m = exprnd(mmat(dir3,3));
        pert_s = exprnd(mmat(dir4,4));
        % 3. Generate random pertubation according to bounded uniform dist.
%         decay = 1/log(cei);
%         pert_a = rand*dmat(dir1,1)*decay;
%         pert_b = rand*dmat(dir2,2)*decay;
%         pert_m = rand*dmat(dir3,3)*decay;
%         pert_s = rand*dmat(dir4,4)*decay;


        % 4. Perform addition/subtraction to mutate
        par_new = par;
        par_new(1) = par(1) + sign1*pert_a;
        par_new(2) = par(2) + sign2*pert_b;
        par_new(3) = par(3) + sign3*pert_m;
        par_new(4) = par(4) + sign3*pert_s;
        % 5. Add new mutation to mutation population
        for q = 1:4
            mut_pop(p,q) = par_new(q);
        end

        % 6. Compute the fitness score of the new mutation
        mut_score(p) = fitness(par_new,vp,struct_rocket,struct_atmo);

    end

    % iv) SELECTION (Post-mutation population) 
    % Finding best performing pair
    for i = 1:num_pop
        if mut_score(i) > km2_s
            if mut_score(i) < km1_s
                km2_s = mut_score(i);
                km2_i = i;
                km2_dim = [mut_pop(i,1),mut_pop(i,2),mut_pop(i,3),mut_pop(i,4)];
            else
                tempi = km1_i;
                temps = km1_s;
                tempdim = km1_dim;
                km2_i = tempi;
                km2_s = temps;
                km2_dim = tempdim;
                km1_s =  mut_score(i);
                km1_i = i;
                km1_dim = [mut_pop(i,1),mut_pop(i,2),mut_pop(i,3),mut_pop(i,4)];
            end
        end
    end

    disp(mut_score);
    % m_ind1 = km1_i;
    % m_ind2 = km2_i;
    % m_mom = [mut_pop(m_ind1,1),mut_pop(m_ind1,2),mut_pop(m_ind1,3),mut_pop(m_ind1,4)];
    % m_dad = [mut_pop(m_ind2,1),mut_pop(m_ind2,2),mut_pop(m_ind2,3),mut_pop(m_ind2,4)];
    % Update parent dimensions IFF new parent is better than previous parent
    %p_par = par;
    m_par = mean([km1_dim;km2_dim]); % Temporary new mutated parent (for checking)

    if cei > 5
        if fitness(m_par,vp,struct_rocket,struct_atmo) > fitness(p_par,vp,struct_rocket,struct_atmo)
            p_par = par;
            par = m_par; % Establish as new mutated parent
        end
    else
        p_par = par;
        par = m_par;
    end

    fprintf('Evolution Fitness: %d \n',fitness(par,vp,struct_rocket,struct_atmo));

    % Check for convergence
    next = true;
    dims_converged = 0;
    %for i = 1:4
    %    if abs(p_par(i)-par(i)) < 0.0001
    %        dims_converged = dims_converged + 1;
    %   end
    %    if dims_converged == 4
    %        next = false;
    %        fprintf('Converged, stopping evolution');
    %    end
    %end
    hold on
    x_bar = [0,par(3),par(3)+par(2),par(1)];
    y_bar = [0,par(4),par(4),0];
    xlim([0,2]);
    ylim([0,2]);
    plot(x_bar,y_bar);
    legend('a');
end

end
% Compute Fitness
% LOOP
% selection
% crossover
% mutation
% compute fitness
% UNTIL POPULATION HAS CONVERGED
% STOP



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

function [score] = fitness(X,vp,str_r,str_a)
    % [INPUT] Set scoring function parameter by the 10th percentile;
    % ie. What drag should a score of 10% (out of 100%) be given to?
   % tenth = 100000; % This hyperparameter is completely arbitrary for now
   % expnt = (log(0.1))/tenth;
   expnt = -0.001;
    
    % Calculate fitness score of fin-set dimensions
    D_nos = str_r(10).r;
    n = 6001;
    
    f_drag = exp(expnt*(drag(X,vp,str_r,str_a)-73000));
    cop = CP_barrow(X,vp,str_r,str_a);
    f_stab = 1;
    for i = 1:n
        margin = (cop*100-vp.altmachcg(i,2))/(D_nos*100);
        if margin < 2
            f_stab = 0;
            break;
        end
    end
    
    f_const = 1;
    %for j = 1:4
    %    if X(j) > 0.3
    %        f_const = 0;
    %    end
    %end
    if  X(1) < 0.01 || X(4) < 0.01 %X(2)<0.75*X(1)|| X(1)-(X(3)+X(2)) > 5 || X(1)-(X(3)+X(2)) < 3
        f_const = 0;
    end
    
    score = f_drag*f_stab*f_const;
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
function m = toMeters(in)
    m = in*0.0254;
end
function [Cd,Cd_incomp] = drag_OR(X,M,alt,str_r,str_a)
    % Interference Drag and Fin-Tip Vortices are ignored as they are 
    % considered "relatively small"
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
    elseif Re <= Rc
        Cf = 1/((1.5*log(Re)-5.6)^2);
    else %Re > Rc
        Cf = 0.032*((5*10^-6)/L_r)^0.2;
    end
    %Compressibility Corrections
    if V < Ma
        Cfc = Cf*(1-0.1*M^2);
    else%if V >= Ma
        turb = Cf/((1+0.15*M^2)^0.58);
        rough = Cf/(1+0.18*M^2);
        if Re < Rc
            Cfc = turb;
        else%if Re > Rc
            if rough < turb
                Cfc = turb;
            else%if rough > turb
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probabilistic selection from given discrete p-dist:
function F = p_select(P)
    x = cumsum([0 P(:).'/sum(P(:))]);
    x(end) = 1e3*eps + x(end);
    [a a] = histc(rand,x);
    F = a;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
