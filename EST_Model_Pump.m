clear all, close all, clc
%% Constants
g = 9.81; %Gravitational constant [m/s^2]
rho_w = 1025;
visc_w = 1.11e-3; %Dynamic viscosity seawater [Pa*s]
r = 0.045*10^-3; %Relative roughness
n_pump = 0.6; %Efficiency

%% Test Parameters
D_pipe = 2; %Diameter of pipe/hydrogenerator [m]
    A_pipe = pi*(D_pipe^2)/4; %Cross-sectional area pipe [m^2]
N_pipe = 500; %Number of pipes/generators
L_pipe = 10; %Length of the pipe [m]
L_reservoir = 3000; %Length/width reservoir [m]
    A_reservoir = L_reservoir^2; %Cross-sectional area reservoir [m^2]
%H_reservoir = 15; %Height of reservoir above sea level [m]
    %V_reservoir = A_reservoir * H_reservoir; %Maximum volume capacity reservoir [m^3]
t_end = 13274; %End time 'experiment'
H_sea = 25; %Sea level height [m]

P_surplus = 1*10^9; %Power Surplus [W]
P_per_pump = P_surplus/N_pipe;

%% Setup
n=1:t_end;
Time = n.';
H_water = zeros(length(n),1); %Height of water in the reservoir above sea level [m]
    H_water(1) = H_sea-0.25;
V_water = zeros(length(n),1); %Volume of water in the reservoir (above sea level) [m^3]
    V_water(1) = A_reservoir*H_water(1);
Delta_P = zeros(length(n),1);
    Delta_P(1) = rho_w*g*(H_sea-H_water(1));
v = zeros(length(n),1);
    v(1) = 0;
Q = zeros(length(n),1);
    Q(1) = N_pipe*v(1)*A_pipe;
f_D = zeros(length(n),1);
    f_D(1) = 0;
Re = zeros(length(n),1);    
    Re(1) = 0;

%% For-loop
for t=2:t_end
    V_water(t) = V_water(t-1) - N_pipe*Q(t-1);
    H_water(t) = V_water(t)/A_reservoir;
    Delta_P(t) = rho_w * (g * (H_sea-H_water(t)) - (f_D(t).*((v(t)).^2) * L_pipe/(2*D_pipe)));
    Q(t) = P_per_pump/(n_pump*Delta_P(t));
    v(t) = Q(t)/A_pipe;
    Re(t) = rho_w*v(t).*L_pipe/visc_w;
        syms F
        eqn_fD = 1/sqrt(F) + 2*log10(r/3.7+(2.51/(Re(t)*sqrt(F)))) == 0;
    f_D(t) = real(vpasolve(eqn_fD,F));
end

%% Plots
T = tiledlayout(2,2);
nexttile
plot(v)
xlim([2 60])
title('Water speed over time')
xlabel('Time [s]')
ylabel('Velocity [m/s]')

nexttile
plot(Q)
xlim([2 60])
title('Flow rate over time')
xlabel('Time [s]')
ylabel('Flow rate [m^3/s]')

nexttile
plot(H_water)
xlim([2 60])
title('Water height over time')
xlabel('Time [s]')
ylabel('Height [m]')

nexttile
plot(V_water)
xlim([2 60])
title('Volume in reservoir over time')
xlabel('Time [s]')
ylabel('Volume [m^3]')


savefig('Water_Height')