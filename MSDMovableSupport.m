%% Coded for the fulfilment of Master's Degree at Politecnico Di Milano
% Author:: Eshwar Bhargav Bhupanam
% Course:: Modeling and Simulation of Aerospace Systems
% Topic:: Mass, Spring and Damper System with Movable Support
% Year:: 2020-2021

%% Initializing....
clear vars; clc; close all

%% Mass Spring System with movable support - Physical model inputs
par.m = 3; % Mass in 'kg'
par.k = 100; % Spring resistance in 'N/m'
par.x0 = 0.2; % Initial position of support in 'm'
par.omega_f = 5; % Frequency of support in 's-1'
par.omega_n = sqrt(par.k/par.m); % Natural frequency of the system in 's-1'
par.y0 = [0 0]; % Considering the system to be in static equilibrium
par.zeta_guess = 0.08;  % Damper coefficient of zeta (Underdamped condition)
                         % Initial guess for 'b'
% Mathematical model
[t_two,t_exp,y,ddot_y,ddot_yexp,b_final,t_iguess,y_iguess] = mechanicalsystem(par);

% Plotting
figure(1)

plot(t_iguess,y_iguess); grid on;
legend({'Position @ Initial guess', 'Speed @ Initial guess'}, 'Interpreter', 'latex');
ylabel('Position [m], Speed [m/s]', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')

figure(2)

plot(t_two,y); grid on;
legend({'Position @ Final', 'Speed @ Final'}, 'Interpreter', 'latex');
ylabel('Position [m], Speed [m/s]', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')

figure(3)

plot(t_two,ddot_y); hold on; grid on;
plot(t_exp',ddot_yexp);
legend({'Accleration\(_{num}\)', 'Accleration\(_{exp}\)'}, 'Interpreter', 'latex');
ylabel('Accleration [\(m/s^2\)]', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')

%% Functions

function [t,t_exp,y,ddot_y,ddot_yexp,b_final,t_iguess,y_iguess] = mechanicalsystem(par)
syms x y(t) b t
EOM = diff(y,t,2)+((b*diff(y,t))/par.m)+((par.omega_n^2)*(y-x)) == 0;
EOM_main = isolate(EOM, diff(y,2));

time.ti = 0;
time.tf = 10;
time.tspan = time.ti:0.1:time.tf;

iter = 1;
tol = 1e-2;
error(iter) = tol+100;

% Importing data from samples.txt
ddot_yexp = importdata('samples.txt');
t_exp=0:0.1:10; hold on;

while error(iter)>tol
    EOM = EOM_main;
    iter = iter + 1;
    if iter==2
        par.zeta = par.zeta_guess; 
    else
        increment = 1e-3;
        par.zeta = par.zeta+increment;
    end
    par.b(iter) = 2*par.zeta*par.omega_n*par.m;
    EOM = subs(EOM,'b',par.b(iter));
    
    [t,y,~] = damperspringmass(t,EOM,par.x0,time.tspan,par.omega_f,par.y0); % Calculates Position and Velocity
    if iter == 2
        t_iguess = t;
        y_iguess = y;
    end
    [~,~,ddot_y] = damperspringmass(t,EOM,par.x0,time.tspan,par.omega_f,par.y0,y); % Calculates Accleration
    
    error(iter) = max(abs(ddot_y-ddot_yexp));
    
    if error(iter)>error(iter-1)
        b_final = par.b(iter-1);
        break;
    end
    clear b t; syms b t;
end
figure(4)
plot(error(:,2:end));
legend({'Accleration Error'}, 'Interpreter', 'latex');
ylabel('max\((abs(\ddot{y}_{num}-\ddot{y}_{exp}))\) [\(m/s^2\)]', 'Interpreter', 'latex')
xlabel('Iterations', 'Interpreter', 'latex')

end

function [t,y,ddot_y] = damperspringmass(teval,EOM,x0,tspan,omega_f,y0,y)
u = exp(-10.*teval);
x = x0.*u.*cos(omega_f.*teval);  % dynamics of the support
EOM = subs(EOM,x);
EOM = simplify(EOM);
ddot_y = zeros(length(EOM),1);
for i = 1:length(EOM)
    V = odeToVectorField(EOM(i));
    if nargin == 6
        MSDode = matlabFunction(V,'vars', {'t', 'Y'},'file','MSDode.m');
        options = odeset;
        [t,y] = ode45(@MSDode,tspan,y0,options);
        ddot_y = [];
        delete MSDode.m;
    elseif nargin>6
        MSDode2 = matlabFunction(V,'vars', {'Y'},'file','MSDode2.m');
        Y = y(i,:);
        V = MSDode2(Y);
        ddot_y(i) = V(2);
        t = tspan;
        delete MSDode2.m;
    end
end
end
