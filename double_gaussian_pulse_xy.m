% Does a double gaussian pulse using shaped_pulse_xy from Spinach

function [rho,P]=double_gaussian_pulse_xy(spin_system, drift, controls, rho, ...
                                power1, frequency1, phase1, ...
                                power2, frequency2, phase2, ...
                                duration, npoints, method)
% Set the defaults
if ~exist('method','var'), method='expv'; end   

% Create the time grid
time_grid = duration*ones(1,npoints)/npoints;
times = linspace(0,duration,npoints);

% gaussian part
mean = duration/2;
stdev = duration/6;


% first pulse amplitudes - rotating wave
amplitude1_X = 1/erf(1.5)*exp(-((mean-times).^2)/(2*stdev^2))*power1.*cos(-2*pi * frequency1*times + phase1); 
amplitude1_Y = 1/erf(1.5)*exp(-((mean-times).^2)/(2*stdev^2))*power1.*sin(-2*pi * frequency1*times + phase1);

% second pulse amplitudes - rotating wave
amplitude2_X = 1/erf(1.5)*exp(-((mean-times).^2)/(2*stdev^2))*power2.*cos(-2*pi * frequency2*times + phase2);
amplitude2_Y = 1/erf(1.5)*exp(-((mean-times).^2)/(2*stdev^2))*power2.*sin(-2*pi * frequency2*times + phase2);

%both pulses at the same time
ampX = amplitude1_X + amplitude2_X;
ampY = amplitude1_Y + amplitude2_Y;

%execute the pulse using shaped_pulse_xy
[rho,P]=shaped_pulse_xy(spin_system,drift,controls,{ampX,ampY},time_grid,rho,method);

end