%% Calc counts from mean values of titlis and spectrum

detector_threshold = [20, 10]; % [keV]
mean_counts = [2*2.996*1e4, 2*4.553*1e4]; % [/pixel/sec] (expo times was 0.5 sec, thus double the count per second)


figure
colors = ['b', 'r'];
for ii=1:2
% Load spectrum
load spectrum
% Cut at threshold
threshold_index = find(energy==detector_threshold(ii));
energy = energy(threshold_index:end);
photons = photons(threshold_index:end);
% Normalize spectrum
photons_norm = photons/sum(photons);
% Multiply with calcs
photons = photons_norm * mean_counts(ii);

% figure(2*ii)
% plot(energy, photons_norm)
% figure(2*ii+14)
% plot(energy, photons, 'r')

plot(energy, photons, colors(ii))
hold on

save(['spectrum_', num2str(detector_threshold(ii)), 'keV'], 'energy', 'photons')
end

%% Take mean of both calcs

load spectrum_20keV
energy_20kev = energy;
photons_20keV = photons;
load spectrum_10keV
energy_10kev = energy;
photons_10keV = photons;

figure
plot(energy_10kev, photons_10keV)
hold on
plot(energy_20kev, photons_20keV)

start_energy_index = find(energy_10kev==energy_20kev(1));

mean_counts = mean([photons_10keV(start_energy_index:end), photons_20keV],2);

plot(energy_20kev, mean_counts, 'gr')

energy = energy_20kev;
photons = mean_counts;
save('spectrum_mean', 'energy', 'photons')

% csvwrite('Comet100kV_counts.csv',[energy, photons])
% csvwrite('Comet100kV_counts_10keV.csv',[energy_10kev, photons_10keV])
% csvwrite('Comet100kV_counts_20keV.csv',[energy_20kev, photons_20keV])