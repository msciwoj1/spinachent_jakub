% Complete set of simulations related to entanglement+ tomography with gaussian pulses. Syntax:
%
%                double_pulse_echo_diag(spin_system,parameters)
%
% Parameters:
%
%      parameters.pulse_frq  - frequencies for the two 
%                              pulses, Hz
%
%      parameters.pulse_pwr  - power levels for the two
%                              pulses, Hz
%
%      parameters.pulse_dur  - durations for the two
%                              pulses, seconds
%
%      parameters.pulse_phi  - initial phases for the two 
%                              pulses, radians
%
%      parameters.offset     - receiver offset for the time
%                              domain detection, Hz
%
%
%      parameters.npoints    - number of points in the free
%                              induction decay 
%
%      parameters.rho0       - initial state
%
%      parameters.coil       - detection state
%
%      parameters.tau        - time between the start of the first pulse
%                              and the start of the second pulse
%
%      parameters.tpoints    - number of sampling points for double pulses
%
%      parameters.echo_time  - time to sample around the ex-
%                              pected echo position
%
%      parameters.echo_npts  - number of points in the echo
%                              discretization
%
%      parameters.method     - soft puse propagation method,
%                              'expv' for Krylov propagation,
%                              'expm' for exponential propa-
%                              gation, 'evolution' for Spin-
%                              ach evolution function
%
%      parameters.assumptions - Hamiltonian generation assump-
%                               tions, use 'deer' to keep two-
%                               electron flip-flop terms and 
%                               'deer-zz' to drop them

function entanglement_diag_gauss(spin_system,parameters)

% Check consistency
%grumble(parameters);


%set the measurement kind
switch parameters.main_freq
    
    case 1
        parameters.pulse_frq(4,2) = parameters.pulse_freqs_used(2);
        parameters.pulse_frq(5,1) = parameters.pulse_freqs_used(1);
        parameters.readout_frq = parameters.pulse_freqs_used(1);
        parameters.echo_phase  = parameters.echo_phase(1);
    
    case 2
        parameters.pulse_frq(4,2) = parameters.pulse_freqs_used(1);
        parameters.pulse_frq(5,1) = parameters.pulse_freqs_used(2);
        parameters.readout_frq = parameters.pulse_freqs_used(2);
        parameters.echo_phase  = parameters.echo_phase(2);
end

switch parameters.readout(1)
    
    case 'I'
        parameters.pulse_pwr(3,1) = 0;
    
    case 'X'
        parameters.pulse_pwr(3,1) = parameters.pulse_pwr(1,1);
        parameters.pulse_phi(3,1) = pi;
    
    case 'Y'
        parameters.pulse_pwr(3,1) = parameters.pulse_pwr(1,1);
        parameters.pulse_phi(3,1) = 3*pi/2;

end

switch parameters.readout(2)
    
    case 'I'
        parameters.pulse_pwr(3,2) = 0;
    
    case 'X'
        parameters.pulse_pwr(3,2) = parameters.pulse_pwr(1,2);
        parameters.pulse_phi(3,2) = pi;
    
    case 'Y'
        parameters.pulse_pwr(3,2) = parameters.pulse_pwr(1,2);
        parameters.pulse_phi(3,2) = 3*pi/2;

end

% % ESR hole burning simulations to locate the holes
% fids=powder(spin_system,@entanglement_hole,parameters,parameters.assumptions);
% 
% % Apodisation and Fourier transform
% fid_a=apodization(fids(:,1),'crisp-1d'); spectrum_a=fftshift(fft(fid_a,parameters.zerofill));
% fid_b=apodization(fids(:,2),'crisp-1d'); spectrum_b=fftshift(fft(fid_b,parameters.zerofill));
% fid_c=apodization(fids(:,3),'crisp-1d'); spectrum_c=fftshift(fft(fid_c,parameters.zerofill));
% fid_d=apodization(fids(:,4),'crisp-1d'); spectrum_d=fftshift(fft(fid_d,parameters.zerofill));
% fid_e=apodization(fids(:,5),'crisp-1d'); spectrum_e=fftshift(fft(fid_e,parameters.zerofill));
% fid_f=apodization(fids(:,6),'crisp-1d'); spectrum_f=fftshift(fft(fid_f,parameters.zerofill));
% 
% % Plotting
% figure();
% subplot(3,2,1); plot_1d(spin_system,real(spectrum_a),parameters,'r-'); title('frequency sweep epr');
% subplot(3,2,2); plot_1d(spin_system,real(spectrum_b),parameters,'r-'); title('90 degree double pulse');
% subplot(3,2,3); plot_1d(spin_system,real(spectrum_c),parameters,'r-'); title('180 degree double pulse');
% subplot(3,2,4); plot_1d(spin_system,real(spectrum_d),parameters,'r-'); title('readout operation');
% subplot(3,2,5); plot_1d(spin_system,real(spectrum_e),parameters,'r-'); title('first readout pulse');
% subplot(3,2,6); plot_1d(spin_system,real(spectrum_f),parameters,'r-'); title('second readout pulse');


% echo simulation
echo=powder(spin_system,@entanglement_sequence_wclock_gauss,parameters,parameters.assumptions);

% Echo phasing
echo=exp(1i*parameters.echo_phase)*echo;

% Time axis for the echo plot
x_axis=1e9*linspace(-parameters.echo_time/2,...
                     parameters.echo_time/2,...
                     parameters.echo_npts+1);

% Raw echo plot
figure(); subplot(1,2,1); plot(x_axis,real(sum(echo,2))); 
title('echo, time domain'); xlabel('ns'); axis tight; grid on;
% 
% Echo Fourier transform plot
subplot(1,2,2); plot(abs(fft(sum(echo,2)))); xlabel('points');
title('echo amplitude, frequency domain'); axis tight; grid on;

% Extract echo modulation
freq = -spin_system.inter.magnet*spin('E')/(2*pi)-...
                      parameters.readout_frq-parameters.offset;
npoints = parameters.echo_npts+1;
echo_time = parameters.echo_time;

 %numfreqs = linspace(echo_time*freq+0.5,echo_time*freq+1.5,60);
 %modulation = zeros(1,parameters.readout_nsteps+1);
 %for i=numfreqs
    fourier = exp(-2i*pi/npoints*(echo_time*freq)*linspace(1,npoints,npoints));                  



    modulation =fourier*echo;%+modulation;
 %end

% modulation=fft(echo,[],1); 
% modulation=modulation(parameters.echo_window,:);
% modulation=sum(modulation,1);

% Time axis for the DEER trace
x_axis=1e6*linspace(0,parameters.readout_gap,parameters.readout_nsteps+1);
filename = ['ent_data_gauss_' char(string(parameters.main_freq)) '_' parameters.readout '_' char(string(parameters.tau*1e9)) 'ns.txt'];
%filename = ['ent_data_gauss_' char(string(parameters.main_freq)) '_' parameters.readout '_' char(string(parameters.tau*1e9)) 'ns_' char(string(spin_system.rlx.temperature)) 'K.txt'];
fileID = fopen(filename,'a+');
A = [x_axis; real(modulation); imag(modulation)];
fprintf(fileID,'%6s %12d %12s\n', 'Experiment', parameters.main_freq, parameters.readout);
fprintf(fileID,'%6.4f %12.4f %12.4f\n', A);
fprintf(fileID,'\n');
fclose(fileID);

% [real_fit,imag_fit] = fit_sines_and_cosines(modulation, x_axis, parameters.deer_freq_guess*1e-6, parameters.final_buffer*1e6);
% 
% real_cos = real_fit.a;
% real_sin = real_fit.b;
% imag_cos = imag_fit.a;
% imag_sin = imag_fit.b;
% 
% real_bounds = confint(real_fit);
% imag_bounds = confint(imag_fit);
% 
% 
% filename = ['ent_full_' char(string(parameters.tau*1e9)) 'ns.txt'];
% 
% fileID = fopen(filename,'a+');
% fprintf(fileID,'%6s %12d %12s\n', 'Experiment', parameters.main_freq, parameters.readout);
% fprintf(fileID,'%6s %12s %12s %12s\n', 'real cos', 'real sin', 'imag_cos', 'imag_sin');
% fprintf(fileID,'%6.3f %12.3f %12.3f %12.3f\n', real_cos, real_sin, imag_cos, imag_sin);
% fprintf(fileID,'%6.3f %12.3f %12.3f %12.3f\n', real_bounds(1,1), real_bounds(1,2), imag_bounds(1,1), imag_bounds(1,2));
% fprintf(fileID,'%6.3f %12.3f %12.3f %12.3f\n', real_bounds(2,1), real_bounds(2,2), imag_bounds(2,1), imag_bounds(2,2));
% fprintf(fileID,'\n');
% fclose(fileID);


% DEER trace plot
figure(); subplot(1,2,1); 
plot(x_axis,real(modulation)); title('modulation, real channel'); 
ylabel('echo amplitude'); axis tight; grid on; xlabel('time, \mus');
subplot(1,2,2); plot(x_axis,imag(modulation)); xlabel('time, \mus');
title('modulation, imag channel'); axis tight; grid on;

% % DEER trace plot
% figure(); subplot(1,2,1); 
% plot(real_fit,x_axis,real(modulation)); title('modulation, real channel'); 
% ylabel('echo amplitude'); axis tight; grid on; xlabel('time, \mus');
% subplot(1,2,2); plot(imag_fit,x_axis,imag(modulation)); xlabel('time, \mus');
% title('modulation, imag channel'); axis tight; grid on;

drawnow();



end


