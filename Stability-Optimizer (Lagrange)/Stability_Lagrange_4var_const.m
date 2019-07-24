%Fin Stability Optimization - 4 variable, 3 constraints
%By: Victor Yip

% MAIN FUNC. FOR LAGRANGE OPT.
n = 20;
x0 = [0.2,0.1,0.1,0.1,2,2,2,2,1,1,1,1];
TEST = zeros(n,12);
TEST([1,n+1,(2*n)+1,(3*n)+1,(4*n)+1,(5*n)+1,(6*n)+1,(7*n)+1,(8*n)+1,(9*n)+1,(10*n)+1,(11*n)+1]) = x0;
for i = 1:n-1
    TEST([1+i,n+1+i,(2*n)+1+i,(3*n)+1+i,(4*n)+1+i,(5*n)+1+i,(6*n)+1+i,(7*n)+1+i,(8*n)+1+i,(9*n)+1+i,(10*n)+1+i,(11*n)+1+i]) = TEST([i,n+i,(2*n)+i,(3*n)+i,(4*n)+i,(5*n)+i,(6*n)+i,(7*n)+i,(8*n)+i,(9*n)+i,(10*n)+i,(11*n)+i])+lagrange(TEST([i,n+i,(2*n)+i,(3*n)+i,(4*n)+i,(5*n)+i,(6*n)+i,(7*n)+i,(8*n)+i,(9*n)+i,(10*n)+i,(11*n)+i]))
end
%----------------------------------Final Centre of Pressure Position-------
% Constants (all dimensions in meters):
L_n = 0.548;    % Nose Length
F_r = 0.406;    % Fin Root Chord
F_t = 0.127;    % Fin Tip Chord
F_s = 0.152;  % Fin Semi-Span
S = 0.100;      % Sweep Distance
L_r = 5.844;  % Rocket Length
L_z = 0.0762;   % Nozzle Length
N = 3;          % # of fins
Cn_n = 0.5;     % Nose Cone Coefficient
t = 0.00635;    % Max Fin Root Thickness
X_tc = 0.00635; % Distance from Fin Leading Edge to Max Thickness
L_red = 0.05995;% Length Reduction @ back
D_noz = 0.08204;% Nozzle Diameter
D_nos = 0.140;  % Nose Base Diameter
D_end = 0.1077; % End Diameter
F_w = 0.00381;    % Fin Width (Thickness)
F_fl = 0.21895; % Fin Front Length
Cg = 3.372354308; %Centre of Gravity 

% X (distance from nose tip to component's Cp):
% NOTE: Normal Body and Shoulder Forces are neglected 
X_n = L_n/2;
X_b = L_r - L_red + (L_red/3)*(1+1/(1+D_nos/D_end));

% Component Coefficient of Normal Force (Cn)
Cn_n = 2;
Cn_b = 2*(((D_end/D_nos)^2)-1);

% LAGRANGE MULTIPLIER OPTIMIZATION
% Variables:
syms a b m s; %a = root chord, b = tip chord, m = sweep, s = semispan
% Using default values for all other variables,
% we can cast final Cp position as cp(a,b,m,s):
mid =  s/cos(atan((m+(b/2)-(a/2))/s));
l =  (L_r-a)-L_red;
int = 1+ (D_nos/2)/(s+(D_nos/2));
c =  int*(4*N*(s/D_nos)^2)/(1+sqrt(1+(2*mid/(a+b))^2));
x =  l+ (m*(a+2*b))/(3*(a+b)) + (1/6)*(a+b-(a*b)/(a+b));
cp =  ((Cn_n*X_n)+(Cn_b*X_b)+(c*x))/(Cn_n + Cn_b + c);

CP_FINAL = eval(subs(cp,{a,b,m,s},TEST([n,2*n,3*n,4*n])))
MARGIN = (CP_FINAL-Cg)/D_nos

x_axis = [0,4.0511811,8.5511811,6,0]
y_axis = [0,5.5,5.5,0,0]
plot(x_axis,y_axis)
xlim([0 9])
ylim([0 9])



%------------------------Lagrange Optimization Function-------------------------------------------
function [X] = lagrange(P)
% Constants (all dimensions in meters):
L_n = 0.548;    % Nose Length
F_r = 0.406;    % Fin Root Chord
F_t = 0.127;    % Fin Tip Chord
F_s = 0.152;  % Fin Semi-Span
S = 0.100;      % Sweep Distance
L_r = 5.844;  % Rocket Length
L_z = 0.0762;   % Nozzle Length
N = 3;          % # of fins
Cn_n = 0.5;     % Nose Cone Coefficient
t = 0.00635;    % Max Fin Root Thickness
X_tc = 0.00635; % Distance from Fin Leading Edge to Max Thickness
L_red = 0.05995;% Length Reduction @ back
D_noz = 0.08204;% Nozzle Diameter
D_nos = 0.140;  % Nose Base Diameter
D_end = 0.1077; % End Diameter
F_w = 0.00381;    % Fin Width (Thickness)
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
syms a b m s lambda1 lambda2 lambda3 lambda4 s1 s2 s3 s4
% Using default values for all other variables,
% we can cast final Cp position as cp(a,b,m,s):
mid =  s/cos(atan((m+(b/2)-(a/2))/s));
l =  (L_r-a)-L_red;
int = 1+ (D_nos/2)/(s+(D_nos/2));
c =  int*(4*N*(s/D_nos)^2)/(1+sqrt(1+(2*mid/(a+b))^2));
x =  l+ (m*(a+2*b))/(3*(a+b)) + (1/6)*(a+b-(a*b)/(a+b));
cp =  ((Cn_n*X_n)+(Cn_b*X_b)+(c*x))/(Cn_n + Cn_b + c);


%We now take the gradient of the Cp function (which we want 
% to optimize) by first finding its partial derivatives w.r.t. to 
% each independent variable:
cp_a = diff(cp,a);
cp_b = diff(cp,b);
cp_m = diff(cp,m);
cp_s = diff(cp,s);
% We now define the constraint equation, along with the multiplier:

g1 = (3*a/4)-b-(s1)^2; 
g2 = a-s-(s2)^2;
g3 = 0.9*b-m-(s3)^2;
g4 = 0.1524-a-(s4)^2;


g1_a =  lambda1*diff(g1,a);
g1_b = lambda1*diff(g1,b);
g1_m =  lambda1*diff(g1,m); 
g1_s = lambda1*diff(g1,s); 

g2_a =  lambda2*diff(g2,a);
g2_b = lambda2*diff(g2,b);
g2_m =  lambda2*diff(g2,m); 
g2_s = lambda2*diff(g2,s); 

g3_a =  lambda3*diff(g3,a);
g3_b = lambda3*diff(g3,b);
g3_m =  lambda3*diff(g3,m); 
g3_s = lambda3*diff(g3,s); 

g4_a =  lambda4*diff(g4,a);
g4_b = lambda4*diff(g4,b);
g4_m =  lambda4*diff(g4,m); 
g4_s = lambda4*diff(g4,s); 
%We must now solve the following system of equations
z1 = cp_a+g1_a+g2_a+g3_a+g4_a;
z2 = cp_b+g1_b+g2_b+g3_b+g4_b;
z3 = cp_m+g1_m+g2_m+g3_m+g4_m;
z4 = cp_s+g1_s+g2_s+g3_s+g4_s;
z5 = g1;
z6 = g2;
z7 = g3;
z8 = g4;
z9 = lambda1*diff(g1,s1);
z10 = lambda2*diff(g2,s2);
z11 = lambda3*diff(g3,s3);
z12 = lambda4*diff(g4,s4);

%[a,b,m,s,lambda] = solve(y1,y2,y3,y4,y5,[a,b,m,s,lambda]) %THIS METHOD
% TAKES TOO LONG EVEN WITH ONLY 2 VARIABLES (5 UNK, 5 EQNS)

% NEWTONS METHOD - NUMERICAL APPROXIMATION OF NON-LINEAR SYSTEMS
% We define the function vector and Jacobian matrix:
F = [z1;z2;z3;z4;z5;z6;z7;z8;z9;z10;z11;z12];
J = [diff(z1,a),diff(z1,b),diff(z1,m),diff(z1,s),diff(z1,lambda1),diff(z1,lambda2),diff(z1,lambda3),diff(z1,lambda4),diff(z1,s1),diff(z1,s2),diff(z1,s3),diff(z1,s4);
    diff(z2,a),diff(z2,b),diff(z2,m),diff(z2,s),diff(z2,lambda1),diff(z2,lambda2),diff(z2,lambda3),diff(z2,lambda4),diff(z2,s1),diff(z2,s2),diff(z2,s3),diff(z2,s4);
    diff(z3,a),diff(z3,b),diff(z3,m),diff(z3,s),diff(z3,lambda1),diff(z3,lambda2),diff(z3,lambda3),diff(z3,lambda4),diff(z3,s1),diff(z3,s2),diff(z3,s3),diff(z3,s4);
    diff(z4,a),diff(z4,b),diff(z4,m),diff(z4,s),diff(z4,lambda1),diff(z4,lambda2),diff(z4,lambda3),diff(z4,lambda4),diff(z4,s1),diff(z4,s2),diff(z4,s3),diff(z4,s4);
    diff(z5,a),diff(z5,b),diff(z5,m),diff(z5,s),diff(z5,lambda1),diff(z5,lambda2),diff(z5,lambda3),diff(z5,lambda4),diff(z5,s1),diff(z5,s2),diff(z5,s3),diff(z5,s4);
    diff(z6,a),diff(z6,b),diff(z6,m),diff(z6,s),diff(z6,lambda1),diff(z6,lambda2),diff(z6,lambda3),diff(z6,lambda4),diff(z6,s1),diff(z6,s2),diff(z6,s3),diff(z6,s4);
    diff(z7,a),diff(z7,b),diff(z7,m),diff(z7,s),diff(z7,lambda1),diff(z7,lambda2),diff(z7,lambda3),diff(z7,lambda4),diff(z7,s1),diff(z7,s2),diff(z7,s3),diff(z7,s4);
    diff(z8,a),diff(z8,b),diff(z8,m),diff(z8,s),diff(z8,lambda1),diff(z8,lambda2),diff(z8,lambda3),diff(z8,lambda4),diff(z8,s1),diff(z8,s2),diff(z8,s3),diff(z8,s4);
    diff(z9,a),diff(z9,b),diff(z9,m),diff(z9,s),diff(z9,lambda1),diff(z9,lambda2),diff(z9,lambda3),diff(z9,lambda4),diff(z9,s1),diff(z9,s2),diff(z9,s3),diff(z9,s4);
    diff(z10,a),diff(z10,b),diff(z10,m),diff(z10,s),diff(z10,lambda1),diff(z10,lambda2),diff(z10,lambda3),diff(z10,lambda4),diff(z10,s1),diff(z10,s2),diff(z10,s3),diff(z10,s4);
    diff(z11,a),diff(z11,b),diff(z11,m),diff(z11,s),diff(z11,lambda1),diff(z11,lambda2),diff(z11,lambda3),diff(z11,lambda4),diff(z11,s1),diff(z11,s2),diff(z11,s3),diff(z11,s4);
    diff(z12,a),diff(z12,b),diff(z12,m),diff(z12,s),diff(z12,lambda1),diff(z12,lambda2),diff(z12,lambda3),diff(z12,lambda4),diff(z12,s1),diff(z12,s2),diff(z12,s3),diff(z12,s4)];

% We now compute F and J values at x0:
F = eval(subs(F,{a,b,m,s,lambda1,lambda2,lambda3,lambda4,s1,s2,s3,s4},P));
J = eval(subs(J,{a,b,m,s,lambda1,lambda2,lambda3,lambda4,s1,s2,s3,s4},P));
syms k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12
% We now perform "sudo-Gaussian Elimination" to find k1,k2,k3,k4:
ge1 = k1*J(1)+k2*J(13)+k3*J(25)+k4*J(37)+k5*J(49)+k6*J(61)+k7*J(73)+k8*J(85)+k9*J(97)+k10*J(109)+k11*J(121)+k12*J(133) + F(1);
ge2 = k1*J(2)+k2*J(14)+k3*J(26)+k4*J(38)+k5*J(50)+k6*J(62)+k7*J(74)+k8*J(86)+k9*J(98)+k10*J(110)+k11*J(122)+k12*J(134) + F(2);
ge3 = k1*J(3)+k2*J(15)+k3*J(27)+k4*J(39)+k5*J(51)+k6*J(63)+k7*J(75)+k8*J(87)+k9*J(99)+k10*J(111)+k11*J(123)+k12*J(135) + F(3);
ge4 = k1*J(4)+k2*J(16)+k3*J(28)+k4*J(40)+k5*J(52)+k6*J(64)+k7*J(76)+k8*J(88)+k9*J(100)+k10*J(112)+k11*J(124)+k12*J(136) + F(4);
ge5 = k1*J(5)+k2*J(17)+k3*J(29)+k4*J(41)+k5*J(53)+k6*J(65)+k7*J(77)+k8*J(89)+k9*J(101)+k10*J(113)+k11*J(125)+k12*J(137) + F(5);
ge6 = k1*J(6)+k2*J(18)+k3*J(30)+k4*J(42)+k5*J(54)+k6*J(66)+k7*J(78)+k8*J(90)+k9*J(102)+k10*J(114)+k11*J(126)+k12*J(138) + F(6);
ge7 = k1*J(7)+k2*J(19)+k3*J(31)+k4*J(43)+k5*J(55)+k6*J(67)+k7*J(79)+k8*J(91)+k9*J(103)+k10*J(115)+k11*J(127)+k12*J(139) + F(7);
ge8 = k1*J(8)+k2*J(20)+k3*J(32)+k4*J(44)+k5*J(56)+k6*J(68)+k7*J(80)+k8*J(92)+k9*J(104)+k10*J(116)+k11*J(128)+k12*J(140) + F(8);
ge9 = k1*J(9)+k2*J(21)+k3*J(33)+k4*J(45)+k5*J(57)+k6*J(69)+k7*J(81)+k8*J(93)+k9*J(105)+k10*J(117)+k11*J(129)+k12*J(141) + F(9);
ge10 = k1*J(10)+k2*J(22)+k3*J(34)+k4*J(46)+k5*J(58)+k6*J(70)+k7*J(82)+k8*J(94)+k9*J(106)+k10*J(118)+k11*J(130)+k12*J(142) + F(10);
ge11 = k1*J(11)+k2*J(23)+k3*J(35)+k4*J(47)+k5*J(59)+k6*J(71)+k7*J(83)+k8*J(95)+k9*J(107)+k10*J(119)+k11*J(131)+k12*J(143) + F(11);
ge12 = k1*J(12)+k2*J(24)+k3*J(36)+k4*J(48)+k5*J(60)+k6*J(72)+k7*J(84)+k8*J(96)+k9*J(108)+k10*J(120)+k11*J(132)+k12*J(144) + F(12);

[k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12]= solve([ge1,ge2,ge3,ge4,ge5,ge6,ge7,ge8,ge9,ge10,ge11,ge12],[k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12]);
X = [k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12];
end