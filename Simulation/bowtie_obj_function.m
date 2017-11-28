function [ y, I, meanE, meanE_onlyfilter] = bowtie_obj_function( x, ang, E, spectrum, MU_breast, MU_Al, sample_length_profile)
%Define an objective function for optimizing the bow tie design
%   Inputs:
%       x:      the thickness of the bow tie design, change with respect to beam [cm]
%       angle:  angluar sampling
%   Outputs:
%       y:      the objective value
%       I:      the intensity profile after the bow-tie

%% Material of the bowtie filter
% material of the bow tie, for now just use 

%% Intensity along the profile
I = zeros(size(ang));
meanE = zeros(size(ang)); % the mean energy for each angle
meanE_onlyfilter = zeros(size(ang)); % the mean energy f
for i = 1:numel(ang)
    spectrum_tmp = zeros(size(E));
    spectrum_tmp1 = zeros(size(E));
    for j = 1:numel(E)
        % sample attenuation and 
        I(i) = I(i) + spectrum(j)*exp(-sample_length_profile(i)/10*MU_breast(j))*exp(-x(i)*MU_Al(j)); % convert sample profile to cm
        spectrum_tmp(j) = spectrum(j)*exp(-sample_length_profile(i)/10*MU_breast(j))*exp(-x(i)*MU_Al(j)); % after filter and sample
        spectrum_tmp1(j) = spectrum(j)*exp(-x(i)*MU_Al(j));
    end
    meanE(i) = sum(spectrum_tmp.*E)/sum(spectrum_tmp);
    meanE_onlyfilter(i) = sum(spectrum_tmp1.*E)/sum(spectrum_tmp1);
end

%% The objective function
%lambda = 0.5;
%y(1) = std(I)+lambda*sum(abs(diff(I)));
y = std(I)/mean(I);
%y = sum(abs(diff(I)));
% 
%y(2) = -mean(I);


end

