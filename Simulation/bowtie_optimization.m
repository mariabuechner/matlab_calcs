%% solve the bowtie optimization problem
% unit is in [cm]
close all;
clear all;
%% angular sampling
ang = -12:0.1:0; % [degree], initial angular sampling

%% sample geometry
SOD = 498; % [mm]
sample_diameter = 140; % [mm]

%% The sample thickness vs angle  [mm]
sample_length_profile = zeros(size(ang));
for i = 1:numel(ang)
    sample_length_profile(i) = 2*sqrt((sample_diameter/2).^2 - (SOD*tan(abs(ang(i))/180*pi)).^2);
end

%% exclude the out of sample area, here the range of the angles is changed
%sample_length_profile = real(sample_length_profile);
ids = find(real(sample_length_profile)~=0);
ang = ang(ids);
sample_length_profile = sample_length_profile(ids);

%% spectrum, not normalized, including intensity info
E = 5:0.1:70; % unit in KeV
mAs = 1;
[es, spectrum] = W_spectrum();
spectrum = interp1(es, spectrum,E,'linear'); % TODO: be careful of characteristic x-rays, cause accuracy issue
% normalize the spectrum with respect to the measured flux
spectrum = spectrum*mAs;


%% Breast33_67 mass attenuation coefficient
[e, mu] = breast_tissue_attenuation_coeff_Nist('Breast33_67'); %linear attenuation
ac_breast = interp1(e,mu,E,'linear');
tr_breast = zeros(numel(sample_length_profile),numel(E));  % transmission of breast vs. E and thickness
for i = 1:numel(sample_length_profile)
    for j = 1:numel(E)
        tr_breast(i,j) = exp(-sample_length_profile(i)/10*ac_breast(j));
    end
end
[X Y] = meshgrid(E, sample_length_profile);
figure;
surf(X, Y, tr_breast);xlabel('E (keV)');ylabel('thickness (mm)');zlabel('transmission');colormap;title('breast transmission');
%figure;plot(E,tr_breast);title('breast transmission'); %the transmission of breast

transmission_breast = exp(-sample_diameter/10.*ac_breast);
figure;plot(E,transmission_breast);title(['Transmission for ', num2str(sample_diameter),' mm breast']);

%% Filter attenuation coefficient: Al, the K-edge is not handled well here.
[e1, mu1] = material_attenuation_coeff_Nist('Al'); % linear attenuation
ac_Al = interp1(e1,mu1,E,'linear');

%% Global filteration to cut low energy photons, the purpose's is to reduce the heat load on G0
T_Al = 0.4; %[cm]
for j = 1:numel(E)
    spectrum(j) = spectrum(j)*exp(-T_Al*ac_Al(j));
end
display('Mean energy:');
mean(sum(spectrum.*E)/sum(spectrum))
figure;plot(E, spectrum);title('Spectrum after global filter');

%%generating a smooth spectrum for dectris
s = smooth(spectrum,7);
figure;plot(E,s);

%% Objective function
FitnessFunction = @(x)bowtie_obj_function(x,ang, E, spectrum, ac_breast, ac_Al, sample_length_profile);

%% Non linear constrain function
nl_cons = @(x)constrain_mean_spectrum(x,ang,40,E, spectrum, ac_breast, ac_Al, sample_length_profile);

%% filter thickness boundary
minT = 0.1; % [cm];
maxT = 0.6; %[cm];
iniT = (minT+maxT)*0.5; %[cm];

%% solver
options = optimoptions('fmincon','Algorithm','interior-point','Display','iter','TolCon',1e-14);
[x,fval,exitflag,output] = fmincon(FitnessFunction,iniT*ones(size(ang)),[],[],[],[],minT*ones(size(ang)),maxT*ones(size(ang)),nl_cons, options);

[yval, I, meanE, meanE_onlyfilter] = bowtie_obj_function(x,ang, E, spectrum, ac_breast, ac_Al, sample_length_profile);
[yy, s1,s2, spectrum_after]        = constrain_mean_spectrum(x,ang,40,E, spectrum, ac_breast, ac_Al, sample_length_profile);

[X, Y] = meshgrid(E, ang);
figure;surf(X, Y, spectrum_after);xlabel('E (keV)');ylabel('ang (degree)');zlabel('spectrum');title('after sample spectrum');

figure;
subplot(231);plot(ang, sample_length_profile);xlabel('Angle [deg]');ylabel('sample thickness [mm]');
subplot(232);plot(ang,x);xlabel('Angle [deg]');ylabel('Filter thickness [cm]');
subplot(233);plot(ang,I);xlabel('Angle [deg]');ylabel('Intensity [cps]');title('Flux after filter and smaple');
subplot(234);plot(ang, meanE);xlabel('Angle [deg]');ylabel('Mean E [keV]');title('with filter and sample');
subplot(235);plot(ang, meanE_onlyfilter);xlabel('Angle [deg]');ylabel('Mean E [keV]');title('with filter');
subplot(236);plot(ang, ones(size(ang))*sum(spectrum));xlabel('Angle [deg]');ylabel('Mean E [keV]');title('Flux without filtration');

%% genetic solver
%numberOfVariables = numel(ang);
%[x,fval,exitflag,output,population,scores] = gamultiobj(FitnessFunction,numberOfVariables,[],[],[],[],zeros(size(ang)),[]);
%[x,fval,exitflag,output,population,scores] = ga(FitnessFunction,numberOfVariables,[],[],[],[],zeros(size(ang)),[]);
% options = optimoptions(@gamultiobj,'PlotFcn',{@gaplotpareto,@gaplotscorediversity});
% [x fval] = gamultiobj(FitnessFunction,numberOfVariables,[],[],[],[],[],[],options);

