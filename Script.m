%% Coded for the fulfilment of Master's Degree at Politecnico Di Milano
% Author:: Eshwar Bhargav Bhupanam
% Course:: Modeling and Simulation of Aerospace Systems
% Topic:: Mass, Spring and Damper System
% Year:: 2020-2021

%% Initializing....
clear vars; clc; close all

%% Symbolic representation (EOM)
syms x(t) m k b v(t) omega_n zeta
v = diff(x);
EOM_undamp = children(diff(v)+(k*x/m) == 0);
m = 1;
k = 100;
EOM = diff(v,t)+(k*x/m) == 0;
% 1 - Transforming EOM in two first order equations
V = odeToVectorField(EOM);

% 2 - Save to matlab File
MSode = matlabFunction(V,'vars', {'t', 'Y'},'file','MSode.m');

% 3 - Integrating the system in time domain
tspan = [0 10];
x0 = 0.2;
v0 = 4;
Y = [x0; v0];
options = odeset;
[t,Y] = ode45(@MSode,tspan,Y,options); % where Y = [x-position v-speed]
figure()
plot(t,Y)
legend({'Position', 'Speed'});
ylabel('Position [m], Speed [m/s]')
xlabel('Time [s]')
title('Mass-Spring System');

% 4 - EOM with damping term
syms k m
EOM_damp = EOM_undamp(1)+((b/m)*v) == 0;
EOM_damp = isolate(EOM_damp,diff(x,2));
EOM_damp_new = subs(EOM_damp,[k b],...
    [m*omega_n^2 2*m*zeta*omega_n]);   % Transforming coefficients from k,m,b to omega_n, zeta
EOM_damp_new = simplify(EOM_damp_new,6);
[V_damp,Ydamp] = odeToVectorField(EOM_damp_new);
MSDode = matlabFunction(V_damp,'file','MSDode.m');

% 5 - Equilibrium point of the system
omega_n = 10;
zeta = {0.1 'Underdamped'; 0.9999 'Criticallydamped'; 1.5 'Overdamped'};
for i = 1:size(zeta,1)
    [equilibrium_point.(zeta{i,2}),~,~,~,JACOB.(zeta{i,2})]=fsolve(@(Y) MSDode(Y,omega_n,zeta{i,1}),[2,40]);
end

% 6 - EOM solution using 'dsolve' and system response w.r.t eigenvalues
% v = diff(x);
condition = [x(0) == 0.2, v(0) == 4];
Solution = dsolve(EOM_damp_new,condition);

%% Plotting
figure()
for i = 1:size(zeta,1)
    eigen.(zeta{i,2}) = eig(JACOB.(zeta{i,2}));
    hold on; grid on;
    plot(real(eigen.(zeta{i,2})),imag(eigen.(zeta{i,2})),'^')
end
xlabel ( 'Re\(({\lambda_i})\)', 'Interpreter', 'latex' )
ylabel ( 'Imag\(({\lambda_i})\)', 'Interpreter', 'latex')
title('System Equilibrium Analysis');

figure()
for i = 1:size(zeta,1)
    Sol.(zeta{i,2}) = subs(Solution, {'omega_n','zeta'}, {omega_n,[zeta{i,1}]});
    hold on; grid on;
    fplot(Sol.(zeta{i,2}),[0 10],'LineWidth',1.5)
end
ylabel('System response [x(t)]')
xlabel('Time [t]')
title('Mass-Spring-Damper System');
