%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


par = [0.0948475132378872,0.0586723604026952,0.0552529770761651,0.298768459000582]

hold on
x_bar = [0,par(3),par(3)+par(2),par(1)];
y_bar = [0,par(4),par(4),0];
xlim([0,0.8]);
ylim([0,0.8]);
plot(x_bar,y_bar);
legend('a');

disp(fitness(par,vp,struct_rocket,struct_atmo));

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
%         else%if margin >= 2
%             f_stab = 1;
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