vp = load('altmachcg.mat') ;
atm = load('atmostable.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant Declaration (all dimensions in meters):
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
X0 = [0.15;0.15;0.15;0.15]; %Initial Conditions for [a,b,m,s]
A = [-0.75,1,0,0;
    0,-0.9,1,0;
    -1,0,0,1];
B = [0;0;0];
Aeq = [];
Beq = [];
LB = [0;0;-0.3;0];
UB = [0.3;0.3;0.3;0.3];
%OPTIONS = optimoptions('fmincon','Algorithm','interior-point');
X = fmincon(@(X) drag(X,vp,struct_rocket,struct_atmo),X0,A,B,Aeq,Beq,LB,UB,@(X) NONLCON(X,vp,struct_rocket,struct_atmo))
% % FUN = Fin Drag Equation (Subsonic)
            
% M = CP_FINAL(app,X);
% for j = 1:4
%     FinDragCoefficient = drag([X(1)-j*0.02,X(2)-j*0.02,X(3)+j*0.02,X(4)-j*0.02],vp,struct_rocket,struct_atmo)
%     margin_vec = margin_confirm([X(1)-j*0.02,X(2)-j*0.02,X(3)+j*0.02,X(4)-j*0.02],vp,struct_rocket,struct_atmo);
%     mean(margin_vec)
%     min(margin_vec)
%     hold on
%     a = X(1)-j*0.02;
%     b = X(2)-j*0.02;
%     m = X(3)+j*0.02;
%     s = X(4)-j*0.02;
% 
%     x_bar = [0,m,m+b,a];
%     y_bar = [0,s,s,0];
%             
%     plot(x_bar,y_bar);
%     xlim([0,0.3]);
%     ylim([0,0.3]);
%     hold off
% end    


    
RootChord = X(1);
TipChord = X(2);
SweepDistance = X(3);
SemiSpan = X(4);

x_bar = [0,X(3),X(3)+X(2),X(1)];
y_bar = [0,X(4),X(4),0];
            
plot(x_bar,y_bar);
xlim([0,0.3]);
ylim([0,0.3]);

margin_vec = margin_confirm(X,vp,struct_rocket,struct_atmo);
%disp(margin_vec)

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

function [Cd,Cd_alpha] = drag_AC(X,M)
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
            V = M*Ma;
            mu = 1.802*(10^-5); % 15degC [kg/m*s]
            Re = (rho*V*mid)/mu; % 1,072,936.44 > Re_c = 5x10^5
            Re_c = 5*(10^5);
            
            B = Re_c*((0.074/(Re^(1/5)))-(1.328/sqrt(Re)));
            if Re > Re_c
                Cff = (0.074/((Re)^(1/5)))-(B/Re); % Viscous Friction Force Coefficient
            elseif Re <= Re_c
                Cff = 1.328/(sqrt(Re));
            end
            % Fin Drag (scales with 1st deg wrt Cff)
            Cdf = 2*Cff*(1+2*(F_w/mid))*(4*N*Afp/(pi*(D_nos^2)));
            
            % Fin-Body Interference Drag (scales with 1st deg wrt Cff)
            Cdi = 2*Cff*(1+2*(F_w/mid))*(4*N*(Afp-Afe)/(pi*(D_nos^2)));
            
            % Angle of Attack Correction Drag Coefficient
            alpha = 5; % Angle of Attack [deg]
            delta = 0.825; % Constant interpolated from Fig 4. [AerodynamicCoefficients.pdf]
            nu = 0.6125; % Constant interpolated from Fig 4. [AerodynamicCoefficients.pdf]
            L_ts = 2*X(4) + D_nos; % Total Span
            Rs = L_ts/D_nos; % Fin Section Ratio
            
            kfb = 0.8065*Rs^2 + 1.1553*Rs; % Fin-body intereference coefficient
            kbf = 0.1935*Rs^2 + 0.8174*Rs + 1; % Body-fin interference coefficient
            
            %CORRECTION DRAG FACTOR FOR ANGLE OF ATTACK
            Cdf_alpha = (alpha^2)*(((1.2*(Afp*4))/(pi*D_nos^2))+3.12*(kfb+kbf-1)*((Afe*4)/(pi*D_nos^2)));
            
            % TOTAL DRAG COEFFICIENT W/O COMPRESSIBLE FLOW CORRECTION
            C_uncompress = Cdf+Cdi;
            C_uncompress_alpha = Cdf + Cdi + Cdf_alpha;
            % Compressible Flow Correction
            if M < 0.8
                C_compress = C_uncompress/sqrt(1-(M)^2);
                C_compress_alpha = C_uncompress_alpha/sqrt(1-(M)^2);
            elseif M > 1.1
                
                C_compress = C_uncompress/sqrt(((M)^2)-1);
                C_compress_alpha = C_uncompress_alpha/sqrt(((M)^2)-1);
            else
                C_compress = C_uncompress/sqrt(1-0.8^2);
                C_compress_alpha = C_uncompress_alpha/sqrt(1-0.8^2);
            end
            
            Cd = C_compress;
            Cd_alpha = C_compress_alpha;
            
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
        
function [c,ceq] = NONLCON(X,vp,str_r,str_a)
	D_nos = str_r(10).r;
            
	c = [0,0];
	ceq = [];
	cop = CP_barrow(X,vp,str_r,str_a);
	%n = length(vp.altmachcg);
	n = 6001;
    for i = 1:n
    	margin = (cop*100-vp.altmachcg(i,2))/(D_nos*100);
        c(1) = c(1) + (-margin+2); % Based on c(x) <= 0;
        c(2) = c(2) + (margin-3);
        %marg_resid = max([margin,2])-min([margin,2]);
    end          
end

function [out] = margin_confirm(X,vp,str_r,str_a)
	D_nos = str_r(10).r;
            
            
	cop = CP_barrow(X,vp,str_r,str_a);
	out = NaN(1,6001);
	%n = length(vp.altmachcg);
	n = 6001;
    for i = 1:n

    	margin = (cop*100-vp.altmachcg(i,2))/(D_nos*100);
    	out(i) = margin;
    end          
end
  
