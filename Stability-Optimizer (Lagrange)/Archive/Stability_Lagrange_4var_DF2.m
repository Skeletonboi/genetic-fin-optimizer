%Stability Optimizer(by Victor Yip)

% LAGRANGE MULTIPLIER OPTIMIZATION
n = 10;
x0 = [0.2,0.1,0.1,0.2,2,2];
TEST = zeros(n,6);
TEST([1,n+1,(2*n)+1,(3*n)+1,(4*n)+1,(5*n)+1]) = x0;
for i = 1:n-1
    TEST([1+i,n+1+i,(2*n)+1+i,(3*n)+1+i,(4*n)+1+i,(5*n)+1+i]) = TEST([i,n+i,(2*n)+i,(3*n)+i,(4*n)+i,(5*n)+i])+lagrange(TEST([i,n+i,(2*n)+i,(3*n)+i,(4*n)+i,(5*n)+i]))
end



%-------------------------------------------------------------------



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

% LAGRANGE MULTIPLIER OPTIMIZATION
% Variables:
syms a b m s lambda1 lambda2
% Using default values for all other variables,
% we can cast final Cp position as cp(a,b,m,s):
mid =  s/cos(atan((m+(b/2)-(a/2))/s));
l =  (L_r-a)-L_red;
int = 1+ (D_nos/2)/(s+(D_nos/2));
c =  int*(4*N*(s/D_nos)^2)/(1+sqrt(1+(2*mid/(a+b))^2));
x =  l+ (m*(a+2*b))/(3*(a+b)) + (1/6)*(a+b-(a*b)/(a+b));
cp =  ((Cn_n*X_n)+(Cn_b*X_b)+(c*x))/(Cn_n + Cn_b + c);
% To verify correctness, the following paramaters should = 1.7196:

%We now take the gradient of the Cp function (which we want 
% to optimize) by first finding its partial derivatives w.r.t. to 
% each independent variable:
cp_a = diff(cp,a);
cp_b = diff(cp,b);
cp_m = diff(cp,m);
cp_s = diff(cp,s);
% We now define the constraint equation, along with the multiplier:
g1 = b-(3*a/4);%Tip chord < 1/2*Root Chord
g2 = s-D_nos;

g1_a =  lambda1*diff(g1,a);
g1_b = lambda1*diff(g1,b);
g1_m =  lambda1*diff(g1,m); %ZERO RESULT
g1_s = lambda1*diff(g1,s); %ZERO RESULT

g2_a =  lambda2*diff(g2,a);
g2_b = lambda2*diff(g2,b);
g2_m =  lambda2*diff(g2,m); 
g2_s = lambda2*diff(g2,s); 

%We must now solve the following system of equations
z1 = cp_a+g1_a+g2_a;
z2 = cp_b+g1_b+g2_b;
z3 = cp_m+g1_m+g2_m;
z4 = cp_s+g1_s+g2_s;
z5 = g1;
z6 = g2;

%[a,b,m,s,lambda] = solve(y1,y2,y3,y4,y5,[a,b,m,s,lambda]) %THIS METHOD
% TAKES TOO LONG EVEN WITH ONLY 2 VARIABLES (5 UNK, 5 EQNS)

% NEWTONS METHOD - NUMERICAL APPROXIMATION OF NON-LINEAR SYSTEMS
% We define the function vector and Jacobian matrix:
F = [z1;z2;z3;z4;z5;z6];
J = [diff(z1,a),diff(z1,b),diff(z1,m),diff(z1,s),diff(z1,lambda1),diff(z1,lambda2);
    diff(z2,a),diff(z2,b),diff(z2,m),diff(z2,s),diff(z2,lambda1),diff(z2,lambda2);
    diff(z3,a),diff(z3,b),diff(z3,m),diff(z3,s),diff(z3,lambda1),diff(z3,lambda2);
    diff(z4,a),diff(z4,b),diff(z4,m),diff(z4,s),diff(z4,lambda1),diff(z4,lambda2);
    diff(z5,a),diff(z5,b),diff(z5,m),diff(z5,s),diff(z5,lambda1),diff(z5,lambda2);
    diff(z6,a),diff(z6,b),diff(z6,m),diff(z6,s),diff(z6,lambda1),diff(z6,lambda2)];

% We now compute F and J values at x0:
F = eval(subs(F,{a,b,m,s,lambda1,lambda2},P));
J = eval(subs(J,{a,b,m,s,lambda1,lambda2},P));
syms k1 k2 k3 k4 k5 k6
% We now perform "sudo-Gaussian Elimination" to find k1,k2,k3,k4:
ge1 = k1*J(1)+k2*J(7)+k3*J(13)+k4*J(19)+k5*J(25)+k6*J(31) + F(1);
ge2 = k1*J(2)+k2*J(8)+k3*J(14)+k4*J(20)+k5*J(26)+k6*J(32) + F(2);
ge3 = k1*J(3)+k2*J(9)+k3*J(15)+k4*J(21)+k5*J(27)+k6*J(33) + F(3);
ge4 = k1*J(4)+k2*J(10)+k3*J(16)+k4*J(22)+k5*J(28)+k6*J(34) + F(4);
ge5 = k1*J(5)+k2*J(11)+k3*J(17)+k4*J(23)+k5*J(29)+k6*J(35) + F(5);
ge6 = k1*J(6)+k2*J(13)+k3*J(18)+k4*J(24)+k5*J(30)+k6*J(36) + F(6);

[k1,k2,k3,k4,k5,k6]= solve([ge1,ge2,ge3,ge4,ge5,ge6],[k1,k2,k3,k4,k5,k6]);
X = [k1,k2,k3,k4,k5,k6];
end