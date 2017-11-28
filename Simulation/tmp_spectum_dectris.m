%% generate the spectrum for Dectris
load 70kvp_no_filtration.mat

%% Filter attenuation coefficient: Al, the K-edge is not handled well here.
[e1, mu1] = material_attenuation_coeff_Nist('Al'); % linear attenuation
ac_Al = interp1(e1,mu1,E,'linear');

%% Global filteration to cut low energy photons, the purpose's is to reduce the heat load on G0
T_Al = 0.4; %[cm]
for j = 1:numel(E)
    Counts(j) = Counts(j)*exp(-T_Al*ac_Al(j));
end
Counts = smooth(Counts,10);
figure;plot(E, Counts);title('Spectrum after global filter');

save('70kvp_filtered','E','Counts');