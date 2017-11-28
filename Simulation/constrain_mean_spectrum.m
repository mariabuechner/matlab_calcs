function [ y, meanE, meanE_onlyfilter, spectrum_after] = constrain_mean_spectrum( x, ang, minE, E, spectrum, MU_breast, MU_Al, sample_length_profile)
%Define an objective function for optimizing the bow tie design
%   Inputs:
%       x:      the thickness of the bow tie design, change with respect to beam [cm]
%       angle:  angluar sampling
%   Outputs:
%       y:      the objective value

%% Material of the bowtie filter
% material of the bow tie, for now just use 
meanE            = zeros(size(ang)); % the mean energy for each angle
meanE_onlyfilter = zeros(size(ang)); % the mean energy f
spectrum_after   = zeros(numel(ang),numel(E));
for i = 1:numel(ang)
    spectrum_tmp = zeros(size(E));
    spectrum_tmp1 = zeros(size(E));
    for j = 1:numel(E)
        spectrum_tmp(j) = spectrum(j)*exp(-sample_length_profile(i)/10*MU_breast(j))*exp(-x(i)*MU_Al(j)); % after filter and sample
        spectrum_tmp1(j) = spectrum(j)*exp(-x(i)*MU_Al(j));
        spectrum_after(i,j) = spectrum_tmp(j);
    end
    meanE(i) = sum(spectrum_tmp.*E)/sum(spectrum_tmp);
    meanE_onlyfilter(i) = sum(spectrum_tmp1.*E)/sum(spectrum_tmp1);
end

%% The objective function
y = minE - meanE_onlyfilter;

end

