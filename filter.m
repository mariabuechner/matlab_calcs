%% Calc filter
clear all

%% Input
% Spectrum
energies = [20:100]; % [keV]
% photons [counts/pixel/s]
load spectrum_mean

% Angles
thetas = 0;

%% Constants
number_phase_steps = 31;
exposure_time = 1; % [s]

% attenuation coefficients
mu_sample = 0; % mu_sample(energy) [1/um] (PMMA)
mu_filter = 0; % mu_filter(energy) [1/um] (Al)
mu_g0 = 0; % mu_g0(energy) [1/um] (Au)
mu_g1 = mu_g0;
mu_g2 = mu_g0;

% Thicknesses (gratings: remaining height, since old gratings where already
% considered in counts calculation)
d0 = 80; % [um]
d1 = 1.2; % [um]
d2 = 60; % [um]

% Duty cycle
dc0 = 0.5;
dc2 = dc0;

% sample thickness
d_s = 0; % d_s(theta) [mm]
d_s = d_s * 1e3; % [um]


%% functions

% transmission
function [t] = height_to_transmission(mu, height)
t = exp(-mu*height); % mu=[1/um], height=[um]
end



% sample transmission
function [c_t] = sample_transmission(theta)
% transmission = exp(mu_sample(energy)*thickness_sample(theta))
% mu_sample(energy) (see above)
% thickness_sample(theta)
t = exp(mu_sample*thickness_sample);
c_t = 1+1./t;
end

%% sigma

% constants and fixed parameters
c_0 = 1/(2*pi*pi*number_phase_steps);

% gi visibility (visibility(energy)=[%])
% grating transmissions (ti = ti(energy))
t0 = height_to_transmission(mu_0, d0);
t1 = height_to_transmission(mu_1, d1);
t2 = height_to_transmission(mu_2, d2);
% with sinc = sin(pi*x)/pi*x
v = (4/pi) * ((dc2.*(1-t2)*sinc(dc2))./(t2+dc2.*(1-t2))) .* ((dc0.*(1-t0)*sinc(dc0))./(t0+dc0.*(1-t0)));
c_v = 1./v.^2; % c_v(energy)

% sample transmission

% c1(energy, theta)
c1 = c_0.*c_v.*c_t

%% Minimization function

% least mean square error of all signas (over all pixels)