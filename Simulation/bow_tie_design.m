%% Design the bow tie filter
% material of the bow tie, for now Al

%% angluar sampling
ang = -12:0.1:12; % [degree]

%% sample geometry
SOD = 498; % [mm]
SID = 996; % [mm]
sample_diameter = 180; % [mm]

%% spectrum, not normalized, including intensity info
E = 5:70; % unit in KeV
Spectrum = ones(size(E)); %TODO

%% Breast33_67 mass attenuation coefficient
[e mu] = breast_tissue_attenuation_coeff_Nist('Breast33_67');
e = e(9:end);
mu = mu(9:end);
mass_density = 940/1000; %g/cm^3
mu = mu*mass_density; %linear attenuation coefficient, unit in 1/cm

% fit the NIST data in order to get the right linear attenuation
% coefficient for the considered energy
% the fitting between e and log(mu) is the best
% mu = -1.562*exp(0.002561*e)+7.989*exp(-0.1006*e) 
% logMU = -1.562*exp(0.002561*E)+7.989*exp(-0.1006*E);
% MU = exp(logMU);
% figure;plot(E,MU);hold on;plot(e,mu,'r');

% or linear attenuation coefficient is obtained by interpolation
MU_intp = interp1(e,mu,E,'linear');

%% Filter attenuation coefficient: Al, the K-edge is not handled well here.
Al_density = 2.6989; % g/cm^3
[e1 mu1] = material_attenuation_coeff_Nist('Al');
mu1 = mu1*Al_density;
MU_Al = interp1(e1,mu1,E,'linear');

%% The sample thickness vs angle
sample_length_profile = zeros(size(ang));
for i = 1:numel(ang)
    sample_length_profile(i) = 2*sqrt((sample_diameter/2).^2 - (SOD*tan(abs(ang(i))/180*pi)).^2);
end

%% TODO exclude the out of sample area

%% Intensity along the profile
bow_tie_profile = zeros(size(ang));
I = zeros(size(ang));
for i = 1:numel(ang)
    for j = 1:numel(E)
        % sample attenuation and 
        I(i) = I(i) + Spectrum(j)*exp(-sample_length_profile(i)*mu(j));
    end
end