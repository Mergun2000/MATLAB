clear all, close all, clc
%% Constants
rho_w = 1025; %Density seawater [kg/m^3]
g = 9.81; %Gravitational constant [m/s^2]
visc_w = 1.11e-3; %Dynamic viscosity seawater [Pa*s]
r = 0.045*10^-3; %Relative roughness
n_generator = 0.6; %Efficiency

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

%% Setup
n=1:t_end;
Time = n.';
H_water = zeros(length(n),1); %Height of water in the reservoir above sea level [m]
    H_water(1) = 0;
V_water = zeros(length(n),1); %Volume of water in the reservoir (above sea level) [m^3]
    V_water(1) = L_reservoir^2*H_water(1);
v = zeros(length(n),1); %Flow speed of water through the generator [m/s]
    v(1) = sqrt(2*g*(H_sea - H_water(1))); 
Q = zeros(length(n),1); %Flow rate through single generator [m^3/s]
    Q(1) = N_pipe*v(1)*A_pipe;
Re = zeros(length(n),1); %Reynolds number [Pas]
    Re(1) = (rho_w*v(1)*D_pipe)/visc_w;
f_D = zeros(length(n),1); %Darcy friction coefficient [-]
        syms f
        eqn_fD = 1/sqrt(f) + 2*log10(r/3.7+(2.51/(Re(1)*sqrt(f)))) == 0;
    f_D(1) = vpasolve(eqn_fD,f);
delta_P = zeros(length(n),1); %Pressure difference before/after pump [Pa]
    delta_P(1) = rho_w*g*(H_sea-H_water(1));
P_output = zeros(length(n),1); %Total power output generators [W]
    P_output(1) = n_generator*delta_P(1)*Q(1);
E_out = zeros(length(n),1);
    E_out(1) = P_output(1)/3600000000000;

%% Loop
for t=2:t_end
    V_water(t) = real(V_water(t-1) + Q(t-1)); %Volume of liquid in reservoir
    H_water(t) = real(V_water(t)/ A_reservoir); %Water height in tank
    v(t) = real(sqrt((2*g*(H_sea - H_water(t)))/(1+(f_D(t-1).*L_pipe/D_pipe)))); %Outflow velocity
    Q(t) = real(N_pipe*v(t)*A_pipe); %Total flow rate
    Re(t) = real(rho_w*v(t).*D_pipe/visc_w); %Reynolds number
        syms F
        eqn_fD = 1/sqrt(F) + 2*log10(r/3.7+(2.51/(Re(t)*sqrt(F)))) == 0;
    f_D(t) = real(vpasolve(eqn_fD,F)); %Darcy friction factor
    delta_P(t) = rho_w*(g*(H_sea-H_water(t)) - f_D(t)*v(t).^2*L_pipe/(2*D_pipe)); %Pressure difference before/after pump 
    P_output(t) = n_generator*delta_P(t)*Q(t);%Total power output generators
    E_out(t) = E_out(t-1) + (P_output(t))/3600000000000;
end 
%% Plots
T = tiledlayout(2,3);

nexttile
plot(v)
xlim([0 t_end])
title('Water speed over time')
xlabel('Time [s]')
ylabel('Velocity [m/s]')

nexttile
plot(Q)
xlim([0 t_end])
title('Flow rate over time')
xlabel('Time [s]')
ylabel('Flow rate [m^3/s]')

nexttile
plot(P_output)
xlim([0 t_end])
title('Power output over time')
xlabel('Time [s]')
ylabel('Power output [W]')

nexttile
plot(H_water)
xlim([0 t_end])
title('Water height over time')
xlabel('Time [s]')
ylabel('Height [m]')

nexttile
plot(V_water)
xlim([0 t_end])
title('Volume in reservoir over time')
xlabel('Time [s]')
ylabel('Volume [m^3]')

nexttile
plot(E_out)
xlim([0 t_end])
title('Energy over storage in time')
xlabel('Time [s]')
ylabel('Energy [Gwh]')

%savefig('Water_Height.png')