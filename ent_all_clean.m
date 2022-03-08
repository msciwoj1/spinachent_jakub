% Entangling gate + tomography simulation in a two-electron system. 
% Soft pulses are simulated using the Fokker-Planck formalism.
%
% Based on the 3-pulse DEER script by Ilya Kuprov
% Calculation time: 3-5 minutes for a single simulation


function ent_all_clean()


% % Magnet field
sys.magnet=0.3451805;

% Magnet field
%sys.magnet=1.201;

% Isotopes
sys.isotopes={'E','E'};


% Zeeman interactions - g tensors input as principal values + Euler angles
inter.zeeman.eigs=cell(1,2);
inter.zeeman.euler=cell(1,2);
inter.zeeman.eigs{1}=[2.0123 1.9774 2.0322];
inter.zeeman.euler{1}=[0 0 0]*(pi/180);
inter.zeeman.eigs{2}=[2.0123 2.0322 1.9774];
inter.zeeman.euler{2}=[0 0 0]*(pi/180);

% Uncomment if dipolar coupling is to be included
% Specify the relative coordinates (Angstrom) of two electrons
inter.coordinates=cell(2,1);
inter.coordinates{1}=[0.00  0.00 0.00];
inter.coordinates{2}=[30.00 0.00 0.00];

% Uncomment if temperature needs to be specified. 
% Also requires initial state recalculation

% inter.temperature = 100;
% parameters.needs = {'aniso_eq'};

% Uncomment if isotropic coupling is to be included
% Given in frequency units
% inter.coupling.scalar = cell(2,2);
% J = 3*1e6;
% inter.coupling.scalar{1,2} = J;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'Lz','E');
parameters.coil=state(spin_system,'L+','E');
parameters.grid='rep_2ang_6400pts_sph';
parameters.method='expm';
parameters.verbose=0;

% EPR parameters
parameters.offset=-0.5e8;
parameters.sweep=2e9;
parameters.npoints=256;
parameters.zerofill=2048;
parameters.axis_units='GHz-labframe';
parameters.derivative=0;
parameters.invert_axis=1;
% Choosing 'deer' instead includes the pseudo-secular terms
% and might make the simulation much longer
parameters.assumptions='deer-zz';

% ent measurement selection 
parameters.main_freq = 1;
parameters.readout = 'YI';

% pulse parameters
parameters.pulse_freqs_used = [9.6e9 9.76e9];
parameters.pulse_rnk=[2 2];
parameters.pulse_dur=[20e-9 20e-9 20e-9 20e-9 20e-9];
parameters.pulse_phi=[pi/2 pi/2; pi/2 pi/2; pi/2 3*pi/2; pi/2 pi/2; pi/2 pi/2];
parameters.pulse_pwr=2*pi*8e6*[1.563 1.563; 3.126 3.126; 1.563 1.563; 3.126 3.126; 3.126 3.126];
parameters.pulse_frq=repmat(parameters.pulse_freqs_used,5,1);
parameters.tpoints = 100; %for double pulses sampling only

%Uncomment if using Gaussian pulses.
% duration = sigma*6, for normalization sigma = t/sqrt(2*pi)
%
%parameters.pulse_dur = parameters.pulse_dur.*(6/sqrt(2*pi));


% Echo timing parameters
parameters.tau=160e-9;
parameters.tau_nsteps = 250;
parameters.tau_stepsize = 10e-9;

parameters.time = 300e-9;

parameters.p1_p3_gap=1e-6;
parameters.p2_nsteps=100;

parameters.readout_buffer = 0e-9;
parameters.readout_gap=4e-6;
parameters.readout_nsteps=250;

parameters.echo_time=100e-9;
parameters.echo_npts=200;

parameters.phase_start = 0; 
parameters.phase_end = 2*pi;
parameters.phase_points = 200;

parameters.final_buffer = parameters.readout_buffer;% + parameters.pulse_dur(3);

% Echo integration parameters
parameters.echo_window=[38 182];
parameters.echo_phase=pi/180*[192 111];



% Simulation and plotting

% Uncomment if just one simulation needed
% entanglement_diag_new(spin_system,parameters);

% Uncomment if one simulation with Gaussian pulses needed
% entanglement_diag_gauss(spin_system,parameters);

% Uncomment to loop through multiple taus and entire tomography procedure
% for ddd = linspace(750,900,31)
%
%           parameters.tau = ddd*1e-9;
%    for i=[1 2]
%        for j='IXY'
%            for k='IXY'
%                parameters.main_freq = i;
%                parameters.readout = [k j];
%                entanglement_diag_new(spin_system,parameters);
%            end
%        end
%   end
%   end


% Uncomment to run a loop with temperature.
% Recalculation of the initial state is required
% for temp = [0.1 ]
%     inter.temperature = temp;
%     spin_system=create(sys,inter);
%     spin_system=basis(spin_system,bas);
%     
%     parameters.rho0=state(spin_system,'Lz','E');
%     parameters.coil=state(spin_system,'L+','E');
%     for ddd = [150 155 160 165 170 175 180 185 190 195 200]
%         parameters.tau = ddd*1e-9;
%         for i=[1 2]
%             for j='IXY'
%                 for k='IXY'
%                     parameters.main_freq = i;
%                     parameters.readout = [k j];
%                     entanglement_diag_gauss(spin_system,parameters);
%                 end
%             end
%         end
%     end
%  end
 end