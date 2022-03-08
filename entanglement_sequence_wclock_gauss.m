% Entanglement + tomography pulse sequence with phase tracking and Gaussian pulses. Syntax: 
%
%            echo=double_pulse_soft_echo(spin_system,parameters,H,R,K)
%
% where H is the Hamiltonian commutation superoperator, R is the
% relaxation superoperator and K is the chemical kinetics super-
% operator. The following parameters are required:
%
%      parameters.pulse_frq  - frequencies for the two 
%                              pulses, Hz
%
%      parameters.pulse_pwr  - power levels for the two
%                              pulses, rad/s
%
%      parameters.pulse_dur  - durations for the two
%                              pulses, seconds
%
%      parameters.pulse_phi  - initial phases for the two 
%                              pulses, radians
%
%      parameters.tau        - time between the start of the first pulse
%                              and the start of the second pulse
%
%      parameters.tpoints    - number of sampling points for double pulses
%                              
%                             
%      parameters.echo_time  - time to sample around the ex-
%                              pected echo position
%
%      parameters.echo_npts  - number of points in the echo
%                              discretization
%
%      parameters.rho0       - initial state
%
%      parameters.coil       - detection state
%
%      parameters.method     - soft pulse propagation method,
%                              'expv' for Krylov propagation,
%                              'expm' for exponential propa-
%                              gation, 'evolution' for Spin-
%                              ach evolution function
%

function ent=entanglement_sequence_wclock_gauss(spin_system,parameters,H,R,K)

% Check consistency
%grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Electron pulse operators
Ep=operator(spin_system,'L+',parameters.spins{1});
Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i;

% Frequency offsets
parameters.pulse_frq=-spin_system.inter.magnet*spin('E')/(2*pi)-...
                      parameters.pulse_frq-parameters.offset;
                  
%First pulse
[rho,P] = double_gaussian_pulse_xy(spin_system, L, {Ex,Ey}, parameters.rho0, parameters.pulse_pwr(1,1), parameters.pulse_frq(1,1), parameters.pulse_phi(1,1), ...
                                                         parameters.pulse_pwr(1,2), parameters.pulse_frq(1,2), parameters.pulse_phi(1,2), parameters.pulse_dur(1), ...
                                                         parameters.tpoints,  parameters.method);                                       
% Evolution
evol_time = parameters.tau - parameters.pulse_dur(1);
rho=evolution(spin_system,L,[],rho,evol_time,1,'final');

total_elapsed = parameters.tau;

% Second pulse 

phase1 = -2*pi*parameters.pulse_frq(1,1)*total_elapsed;
phase2 = -2*pi*parameters.pulse_frq(1,2)*total_elapsed;

[rho,P] = double_gaussian_pulse_xy(spin_system, L, {Ex,Ey}, rho, parameters.pulse_pwr(2,1), parameters.pulse_frq(2,1), parameters.pulse_phi(2,1)+phase1, ...
                                                         parameters.pulse_pwr(2,2), parameters.pulse_frq(2,2), parameters.pulse_phi(2,2)+phase2, parameters.pulse_dur(2), ...
                                                         parameters.tpoints,  parameters.method);

                                                     % Evolution
evol_time = parameters.tau -  parameters.pulse_dur(2);
rho=evolution(spin_system,L,[],rho,evol_time,1,'final');

total_elapsed = 2*parameters.tau;

% Third pulse - spin operation

phase1 = -2*pi*parameters.pulse_frq(1,1)*total_elapsed;
phase2 = -2*pi*parameters.pulse_frq(1,2)*total_elapsed;

[rho,P] = double_gaussian_pulse_xy(spin_system, L, {Ex,Ey}, rho, parameters.pulse_pwr(3,1), parameters.pulse_frq(3,1), parameters.pulse_phi(3,1)+phase1, ...
                                                         parameters.pulse_pwr(3,2), parameters.pulse_frq(3,2), parameters.pulse_phi(3,2)+phase2, parameters.pulse_dur(3), ...
                                                         parameters.tpoints,  parameters.method);

% DEER readout from here

rho = evolution(spin_system,L,[],rho,parameters.readout_buffer-parameters.pulse_dur(3),1,'final');

% Evolution
stepsize=parameters.readout_gap/parameters.readout_nsteps;
rho_stack=evolution(spin_system,L,[],rho,stepsize,parameters.readout_nsteps,'trajectory');

% First DEER readout pulse
rho_stack=shaped_pulse_af(spin_system,L,Ex,Ey,rho_stack,parameters.pulse_frq(4,2),parameters.pulse_pwr(4,2),...
                                                        parameters.pulse_dur(4),parameters.pulse_phi(4,2),...
                                                        parameters.pulse_rnk(1),parameters.method);
% Evolution
rho_stack(:,end:-1:1)=evolution(spin_system,L,[],rho_stack(:,end:-1:1),stepsize,parameters.readout_nsteps,'refocus');
rho_stack = evolution(spin_system,L,[],rho_stack,parameters.readout_buffer,1,'final');

total_elapsed = total_elapsed + 2*parameters.readout_buffer + parameters.readout_gap;
% Second DEER readout pulse

phase1 = -2*pi*parameters.pulse_frq(5,1)*total_elapsed;

rho_stack=shaped_pulse_af(spin_system,L,Ex,Ey,rho_stack,parameters.pulse_frq(5,1),parameters.pulse_pwr(5,1),...
                                                        parameters.pulse_dur(5),parameters.pulse_phi(5,1)+phase1,...
                                                        parameters.pulse_rnk(2),parameters.method);
% Evolve to the start of the echo
echo_location=parameters.readout_gap+2*parameters.readout_buffer-parameters.pulse_dur(3)+parameters.pulse_dur(4)/2+parameters.pulse_dur(5)-parameters.echo_time/2;
rho_stack=evolution(spin_system,L,[],rho_stack,echo_location,1,'final');

total_elapsed = total_elapsed + parameters.pulse_dur(5)+echo_location;
phase1 = -2*pi*parameters.pulse_frq(5,1)*total_elapsed;                                                
% Sample the echo
stepsize=parameters.echo_time/parameters.echo_npts;
ent=exp(phase1*1i)*evolution(spin_system,L,parameters.coil,rho_stack,stepsize,parameters.echo_npts,'observable');                                               
end



