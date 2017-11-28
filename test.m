%%

% %%

I_I0 = [0.0958];
mu_Au_30 = 2.752e1; % [cm2/g]
% mu_Au_30 = 2.349e1; % mu_en
mu_Au_40 = 1.298e1; % [cm2/g]
mu_Au_50 = 7.256; % [cm2/g]
mu_Au_45 = (mu_Au_40+mu_Au_50)/2;
mu_Au_60 = 4.528; % [cm2/g]
mu_Au_150 = 1.860; % [cm2/g]
% mu_Au_100 = 5.158; % [cm2/g]
mu_Au_100 = 2.027; % mu_en
mu_Pb_100 = 5.549; % [cm2/g]
mu_Pb_150 = 2.014; % [cm2/g]
mu_Fe_60 = 1.205; % [cm2/g]
mu_Si_60 = 3.207e-01; % [cm2/g]


rho_Au = 1.932e1; % [g/cm3]
rho_Pb = 1.135e1; % [g/cm3]
rho_Fe = 7.874; % [g/cm3]
rho_Si = 2.3290; % [g/cm3]

h_G2 = -log(I_I0)/(mu_Au_50*rho_Au); % [cm]
% h = -log(I_I0)/(mu_Fe_60*rho_Fe) % [cm]
h_G2 = h_G2*1e4 % [um]

%% Grating weight

rho_Au = 1.932e1; % [g/cm3]
rho_Au = rho_Au * 1e6; % [g/m3]
rho_Si = 2.3290; % [g/cm3]
rho_Si = rho_Si * 1e6; % [g/m3]

height = 120; % [um]
height = height * 1e-6; % [m]

extra_gold = 10; % [um]
extra_gold = extra_gold * 1e-6; % [m]

dc = 0.5;

wafer_height = 200; % [um]
wafer_height = wafer_height * 1e-6; % [m]

wafer_size = 8; % [cm], 6cm with edges
wafer_size = wafer_size * 1e-2; % [m]

si_volume = wafer_size^2*(wafer_height + dc*height); % [m3]
si_weight = si_volume*rho_Si % [g]

au_volume = wafer_size^2*(dc*height + extra_gold); % [m3]
au_weight = au_volume*rho_Au % [g]

total_grating_weight = si_weight + au_weight % [g]

rho_PMMA = 1.18; % [g/cm3]
rho_PMMA = rho_PMMA * 1e6; % [g/m3]

support_height = 2; % [mm]
support_height = support_height * 1e-3; % [m]

PMMA_volume = wafer_size^2*support_height;
PMMA_weight = PMMA_volume*rho_PMMA % [g]

total_weight = total_grating_weight + PMMA_weight % [g]

%%

total_weight = 2000; % [g]

r = 0.5; % [m]
u = 2*pi*r %[m]
rpm = 50;
rps = rpm/60; % [Hz]

v = u*rps % [m/s]
% total_grating_weight = total_grating_weight * 1e-3 ; % [kg]
% F = total_grating_weight * v^2 / r % [N = kg * m/s^2]
total_weight = total_weight * 1e-3 ; % [kg]
F_total = total_weight * v^2 / r % [N]

%%

h = 120; % [um]
h = h/1e4; % [cm]


I_I0 = exp(-mu_Au_45*rho_Au*h)
% I_I0 = exp(-mu_Si_30*rho_Si*h)

%%

% Aspect ratio
aspr = 40; % aspr = 2*h/p
% possible pitch
p2 = 2*h_G2/aspr

%% Phase shift

delta_SI_60 = 0.13402e-06; % []
lambda_60 = 0.020664; % [nm]
lambda_60 = lambda_60*1e-7; % [cm]

delta_SI_50 = 0.19305E-06; % []
lambda_50 = 0.0247968; % [nm]
lambda_50 = lambda_50*1e-7; % [cm]

delta_SI_40 = 0.30180E-06; % []
lambda_40 = 0.030996; % [nm]
lambda_40 = lambda_40*1e-7; % [cm]

delta_SI_30 = 0.53707E-06; % []
lambda_30 = 0.0413281; % [nm]
lambda_30 = lambda_30*1e-7; % [cm]

delta_SI_28 = 0.61673E-06; % []
lambda_28 = 0.0442801; % [nm]
lambda_28 = lambda_28*1e-7; % [cm]

delta_SI_25 = 0.77406E-06; % []
lambda_25 = 0.0495937; % [nm]
lambda_25 = lambda_25*1e-7; % [cm]

type = pi;
% type = pi/2;

h_G1 = (type*lambda_25)/(2*pi*delta_SI_25); % [cm]
h_G1 = h_G1*1e4 % [um]

h_G1 = (type*lambda_30)/(2*pi*delta_SI_30); % [cm]
h_G1 = h_G1*1e4 % [um]

% possible pitch
p1 = 2*h_G1/aspr

%%

% b = 0.5;
% m = (100-0.5)/420;
% x = 320; % [mm]
% 
% f = b + m*x
%%

% x = [20:0.01:40]*1e4; % [um]
% 
% m_50 = 50/60e4;
% y_50 = m_50*x;
% m_40 = 40/60e4;
% y_40 = m_40*x;
% m_60 = 60/60e4;
% y_60 = m_60*x;
% 
% figure()
% plot(x*1e-4, 2*y_40, 'r', x*1e-4, 2*y_50,x*1e-4, 2*y_60, 'g')
% legend('G2: 80 um', 'G2: 100 um', 'G2: 120 um')
% xlabel('Distance source to G1 [cm]')
% ylabel('G1 thickness [um]')
% title('Thicknes of G1 depending in distance from sourve to G1')
% 
% %%
% 
% x = [0:0.01:2*pi];
% 
% y_flat = 1+sin(x);
% y_breast = 0.8+0.9*sin(x-0.4);
% y_interest = 0.7+0.8*sin(x-0.5);
% 
% figure()
% % plot(x, y_flat, 'r', x, y_breast, x, y_interest, 'g')
% % legend('Flat', 'Breast tissue', 'Lesion')
% plot(x, y_flat, 'r', x, y_breast)
% legend('Flat', 'Breast')
% xlabel('dx')
% ylabel('PSC')
% %title('Thicknes of G1 depending in distance from sourve to G1')


%%

a = 60; % [um]
a = a*1e-6; % [m]
ga = sqrt(3/8);

l = 30; % [cm]
l = l*1e-2; % [m]
lambda_p = 0.243106; % [A]
lambda_p = lambda_p*1e-1; % [nm]
lambda_p = lambda_p*1e-9; % [m]
dc = 0.5;
g0 = 4.8621; % [um]
g0 = g0*1e-6; % [m]
ex = (l*lambda_p)/(2*pi*g0*dc); % [m]
% ex = ex*1e9; % [um]

dphi = 11.42;
dmu = 0.0098;
f = 1.2;

n = (2/exp(1))/f*(ga^2/a^2)*ex^2*(dphi/dmu)^2


%
g1 = 6; % [um]
g1 = g1*1e-6; % [m]
n=[1, 3, 5];
Dn = n.*(g1/2)^2/(2*lambda_p); % [m] (pick the one that is smaller than l!)
Dn = Dn(Dn<l);
d2 = l*Dn./(l-Dn); % [m]
s_total = l+d2; % [m]
g0 = (l/d2)*((l+d2)/l)*g1/2; % [m]
g0*1e6 % [um]
%

s = (2*lambda_p*d2)/g1;
v_eff = exp(-s^2/(2*ex^2));
v_eff = 0.1
e_orig = sqrt(2)*v_eff*(l*lambda_p/g1)*(dphi/dmu);
n = (ga^2/a^2)*e_orig^2


%%

t = -200:0.01:200;       % sample freq (um)
d = -200:4:200;           % repetition frequency

w = 2;

pulse_train = pulstran(t, d, 'rectpuls', w);

figure
plot(pulse_train)


transmission_function = exp(1i*pi.*pulse_train);

figure
plot(real(transmission_function))
figure
plot(imag(transmission_function))