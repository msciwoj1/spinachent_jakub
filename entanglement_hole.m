% Hole diagnostics of entanglement + tomography
%
%      parameters.pulse_frq  - two frequencies used, Hz
%
%      parameters.pulse_pwr  - power levels for the 
%                              pulses, rad/s
%
%      parameters.pulse_dur  - durations for the
%                              pulses, seconds
%
%      parameters.pulse_phi  - initial phases for the  
%                              pulses, radians
%
%      parameters.tpoints    - number of points in the pulses
%
%      parameters.offset     - receiver offset for the time
%                              domain detection, Hz
%
%      parameters.sweep      - sweep width for time domain
%                              detection, Hz
%
%      parameters.npoints    - number of points in the free
%                              induction decay 
%
%      parameters.rho0       - initial state
%
%      parameters.coil       - detection state
%


function fids=entanglement_hole(spin_system,parameters,H,R,K)

% Check consistency
%grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Pulse operators
Ep=operator(spin_system,'L+',parameters.spins{1});
Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i;

% Frequency offsets
parameters.pulse_frq=-spin_system.inter.magnet*spin('E')/(2*pi)-...
                      parameters.pulse_frq-parameters.offset;

% Soft pulses

[rho1,P] = double_pulse_xy(spin_system, L, {Ex,Ey}, parameters.rho0, parameters.pulse_pwr(1,1), parameters.pulse_frq(1,1), parameters.pulse_phi(1,1), ...
                                                         parameters.pulse_pwr(1,2), parameters.pulse_frq(1,2), parameters.pulse_phi(1,2), parameters.pulse_dur(1), ...
                                                         parameters.tpoints,  parameters.method);  
                                                     
[rho2,P] = double_pulse_xy(spin_system, L, {Ex,Ey}, parameters.rho0, parameters.pulse_pwr(2,1), parameters.pulse_frq(2,1), parameters.pulse_phi(2,1), ...
                                                         parameters.pulse_pwr(2,2), parameters.pulse_frq(2,2), parameters.pulse_phi(2,2), parameters.pulse_dur(2), ...
                                                         parameters.tpoints,  parameters.method);

[rho3,P] = double_pulse_xy(spin_system, L, {Ex,Ey},  parameters.rho0, parameters.pulse_pwr(3,1), parameters.pulse_frq(3,1), parameters.pulse_phi(3,1), ...
                                                         parameters.pulse_pwr(3,2), parameters.pulse_frq(3,2), parameters.pulse_phi(3,2), parameters.pulse_dur(3), ...
                                                         parameters.tpoints,  parameters.method);                                                     

rho4     = shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,parameters.pulse_frq(4,2),parameters.pulse_pwr(4,2),...
                                                        parameters.pulse_dur(4),parameters.pulse_phi(4,2),...
                                                        parameters.pulse_rnk(1),parameters.method);
rho5     = shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,parameters.pulse_frq(5,1),parameters.pulse_pwr(5,1),...
                                                        parameters.pulse_dur(5),parameters.pulse_phi(5,1),...
                                                        parameters.pulse_rnk(2),parameters.method);                                                  
                                                     
                                                     
% Hard pulse
parameters.rho0=step(spin_system,Ey,[parameters.rho0 rho1 rho2 rho3 rho4 rho5],pi/2);
                                  
% Acquisition
fids=acquire(spin_system,parameters,H,R,K);

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only available in Liouville space.');
end
if ~isfield(parameters,'pulse_frq')
    error('pulse frequencies must be specified in parameters.pulse_frq field.');
end
if (~isnumeric(parameters.pulse_frq))||(~isreal(parameters.pulse_frq))||...
   (numel(parameters.pulse_frq)~=3)
    error('parameters.pulse_frq must have three real elements.');
end
if ~isfield(parameters,'pulse_pwr')
    error('pulse powers must be specified in parameters.pulse_pwr field.');
end
if (~isnumeric(parameters.pulse_pwr))||(~isreal(parameters.pulse_pwr))||...
   (numel(parameters.pulse_pwr)~=3)||any(parameters.pulse_pwr<=0)
    error('parameters.pulse_pwr must have three positive real elements.');
end
if ~isfield(parameters,'pulse_dur')
    error('pulse durations must be specified in parameters.pulse_dur field.');
end
if (~isnumeric(parameters.pulse_dur))||(~isreal(parameters.pulse_dur))||...
   (numel(parameters.pulse_dur)~=3)||any(parameters.pulse_dur<=0)
    error('parameters.pulse_dur must have three positive real elements.');
end
if ~isfield(parameters,'pulse_phi')
    error('pulse phases must be specified in parameters.pulse_phi field.');
end
if (~isnumeric(parameters.pulse_phi))||(~isreal(parameters.pulse_phi))||...
   (numel(parameters.pulse_phi)~=3)
    error('parameters.pulse_phi must have three real elements.');
end
if ~isfield(parameters,'pulse_rnk')
    error('pulse grid ranks must be specified in parameters.pulse_rnk field.');
end
if (~isnumeric(parameters.pulse_rnk))||(~isreal(parameters.pulse_rnk))||...
   (numel(parameters.pulse_rnk)~=3)||any(mod(parameters.pulse_rnk,1)~=0)
    error('parameters.pulse_rnk must have three integer real elements.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'method')
    error('shaped pulse simulation method must be specified in parameters.method field.');
end
if ~isfield(parameters,'sweep')
    error('width of the detection window must be specified in parameters.sweep field.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
   (~isscalar(parameters.sweep))||(parameters.sweep<=0)
    error('parameters.sweep must be a positive real scalar.');
end
if ~isfield(parameters,'npoints')
    error('number of points in the FID must be specified in parameters.npoints field.');
end
if (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
   (~isscalar(parameters.npoints))||(parameters.npoints<1)||...
   (mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a positive real integer.');
end
end