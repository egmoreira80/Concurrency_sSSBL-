function [Y,psd,f] = making_spectrum(fo, Damping, phi11, vari ,vari2, fs)
%fo frequiencia pico
%Damping entre 0.3 y 0.9
%phi11 entre -0.5 y 0.5
%lo demas son las varianzas y la frequencia de muestreo.

if nargin < 6
    fs = 300;
end
if nargin < 5
    vari2 = 0.1;
end
if nargin < 4
    vari = 0.1;
end
if nargin < 3
    phi11 = 0.5;
end
if nargin < 2
    Damping = 0.9;    
end
if nargin < 1
    fo = 10;
end

f = (0:1/fs:30);
phi2 = -(Damping^2);
phi1 = 2*Damping*cos(2*pi*fo/fs);
psd1 = 2*vari2./(1+phi11^2-2*phi11*cos(2*pi*f/fs));
psd2 = 2*vari./(1+phi1^2+phi2^2-2*phi1.*(1-phi2).*cos(2*pi*f/fs)-2*phi2*cos(4*pi*f/fs));
psd = psd1 + psd2;

model1 = arima('Constant',0,'AR',{phi1,phi2},'Variance',vari);
rng('default')
Y1 = simulate(model1, 100*fs);

model2 = arima('Constant',0,'AR',{phi11},'Variance',vari2);
Y2 = simulate(model2, 100*fs);

Y = Y1+Y2;

end

