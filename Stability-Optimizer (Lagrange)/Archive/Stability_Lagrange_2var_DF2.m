%Stability Optimizer(by Victor Yip)

% LAGRANGE MULTIPLIER OPTIMIZATION
n = 10;
x0 = [0.1,0.1,2];
TEST = zeros(n,4);
TEST([1,n+1,(2*n)+1]) = x0;
for i = 1:n-1
    TEST([1+i,n+1+i,(2*n)+1+i]) = TEST([i,n+i,(2*n)+i])+lagrange(TEST([i,n+i,(2*n)+i]))
end


function [X] = lagrange(P)
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

% X (distance from nose tip to component's Cp):
% NOTE: Normal Body and Shoulder Forces are neglected 
X_n = L_n/2;
X_b = L_r - L_red + (L_red/3)*(1+1/(1+D_nos/D_end));

% Component Coefficient of Normal Force (Cn)
Cn_n = 2;
Cn_b = 2*(((D_end/D_nos)^2)-1);
int = 1+ (D_nos/2)/(F_s+(D_nos/2));

% LAGRANGE MULTIPLIER OPTIMIZATION
% Variables:
syms a m lambda
% Using default values for all other variables,
% we can cast final Cp position as cp(a,b,m,s):
mid = F_s/cos(atan((m+(F_t/2)-(a/2))/F_s));
l = (L_r-a)-L_red;
%int = @(s) 1+ (D_nos/2)/(F_s+(D_nos/2));
c =  int*(4*N*(F_s/D_nos)^2)/(1+sqrt(1+(2*mid/(a+F_t))^2));
x =  l+ (m*(a+2*F_t))/(3*(a+F_t)) + (1/6)*(a+F_t-(a*F_t)/(a+F_t));
cp =  ((Cn_n*X_n)+(Cn_b*X_b)+(c*x))/(Cn_n + Cn_b + c);
% To verify correctness, the following paramaters should = 1.7196:

%We now take the gradient of the Cp function (which we want 
% to optimize) by first finding its partial derivatives w.r.t. to 
% each independent variable:
cp_a = diff(cp,a);
%cp_b = diff(cp,b);
cp_m = diff(cp,m);
%cp_s = diff(cp,s);
% We now define the constraint equation, along with the multiplier:
g = F_t-(3*a/4); %Tip chord < 1/2*Root Chord
g_a =  lambda*diff(g,a);
%g_b = lambda*diff(g,b);
g_m =  lambda*diff(g,m); %ZERO RESULT
%g_s = lambda*diff(g,s); %ZERO RESULT

%We must now solve the following system of equations
y1 =  cp_a+g_a;
%y2 = cp_b+g_b;
y3 = cp_m+g_m;
%y4 = cp_s-g_s;
y5 = g;



%[a,m,lambda] = solve(y1,y3,y5,[a,m,lambda]) %THIS METHOD TAKES
% EVEN WITH ONLY 2 VARIABLES (3 UNK, 3 EQNS)

% NEWTONS METHOD - NUMERICAL APPROXIMATION OF NON-LINEAR SYSTEMS
% We first define the initial x0 "guess" vector:
% [a,b,lambda]

% We then define the function vector and Jacobian matrix:
F = [y1;y3;y5];
J = [diff(y1,a),diff(y1,m),diff(y1,lambda);
    diff(y3,a),diff(y3,m),diff(y3,lambda);
    diff(y5,a),diff(y5,m),diff(y5,lambda)];

% We now compute F and J values at x0:
F = eval(subs(F,{a,m,lambda},P));
J = eval(subs(J,{a,m,lambda},P));
syms k1 k2 k3 
% We now perform sudo-Gaussian Elimination to find k1,k2,k3:
ge1 = k1*J(1)+k2*J(2)+k3*J(3) + F(1);
ge2 = k1*J(4)+k2*J(5)+k3*J(6) + F(2);
ge3 = k1*J(7)+k2*J(8)+k3*J(9)+ F(3);


[k1,k2,k3]= solve(ge1,ge2,ge3,k1,k2,k3);
X = [k1,k2,k3];
end