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