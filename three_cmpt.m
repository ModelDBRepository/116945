%
%
% solve the 2 compartment lgmd model
%   for a constant current injection
%
% The model contains
%
% 1) I_Na, I_DR in the axon compartment with variables m_naa, h_nai, m_dra
%
% 2) I_Ca, I_AHP in the axon compartment with variables m_caa, i_ahp, ca_c
%
%
% the state variable vector is:
%
% y = [va vd h_nai m_dra ca_c]
%
% where va = axon membrane potential, vd = dendritic membrane potential,
% h_nai = sodium conductance inactivation variable, m_dra = delayed
% rectifier activation variable, ca_c = calcium concentration. 
% All other variables are set to their instantaneous steady-state value.
%
%
%  usage   three_cmpt(T,[t_s t_e Iv cmpt],param_struct)     e.g.,  lgmd_mod(500,[50 450 1],p_s);
%
%  where param_struct is a structure containing the parameters of the
%  model. Current injection parameters are determined as follows:
%
%          T = duration of simulation (ms)
%          t_s, t_e = start and end time of current stimulation,
%          respectively
%          Iv = current injection value (muA/cm2)
%          cmpt = compartment (1 = axon 2 = dendrite 3 = mid compartment)*
%
% * YES, the wiring *IS* weird --  2 <-> 3 <-> 1 is the scheme.
%
% The parameter structure (p_s) is organized as follows:
%
%          p_s.v_ss = steady-state vm for initialization purposes 
%          p_s.... = other variables of the model
%
%  figure produces Vs and Vd vs. time
%

function [t, y] = three_cmpt(ps_mod)

%initial steady state values
v_ss = ps_mod.v_ss;
cacss = ca_c_ss(v_ss, ps_mod.alphaCa, ps_mod.tauCa, ps_mod.VCa, ps_mod.gCa, ps_mod.gCaV12, ps_mod.gCaSl);
%cacss = 0; % force start [ca] = 0
y0 = [v_ss v_ss v_ss m_naa_ss(v_ss) h_nai_ss(v_ss) m_dra_ss(v_ss) cacss m_ih_ss(v_ss) m_caa_ss(v_ss, ps_mod.gCaV12, ps_mod.gCaSl), 0, 0, 0];

% an interesting plot
if ( 0 == 1 )
  v = -100:1:100;
  for i=1:length(v)
    cass_a(i) = ca_c_ss(v(i), ps_mod.alphaCa, ps_mod.tauCa, ps_mod.VCa, ps_mod.gCa, ps_mod.gCaV12, ps_mod.gCaSl);
    cass_b(i) = ca_c_ss(v(i), 2*ps_mod.alphaCa, ps_mod.tauCa, ps_mod.VCa, ps_mod.gCa, ps_mod.gCaV12, ps_mod.gCaSl);
    cass_c(i) = ca_c_ss(v(i), .5*ps_mod.alphaCa, ps_mod.tauCa, ps_mod.VCa, ps_mod.gCa, ps_mod.gCaV12, ps_mod.gCaSl);
    cass_d(i) = ca_c_ss(v(i), ps_mod.alphaCa, ps_mod.tauCa, ps_mod.VCa, ps_mod.gCa, ps_mod.gCaV12, 2*ps_mod.gCaSl);
    cass_e(i) = ca_c_ss(v(i), ps_mod.alphaCa, ps_mod.tauCa, ps_mod.VCa, ps_mod.gCa, ps_mod.gCaV12, .5*ps_mod.gCaSl);
  
    % and "m_ahp(ca)"
    mahp_a(i) = cass_a(i)/(cass_a(i) + ps_mod.kD);
    mahp_b(i) = cass_a(i)/(cass_a(i) + 2*ps_mod.kD);
    mahp_c(i) = cass_a(i)/(cass_a(i) + 4*ps_mod.kD);
    mahp_d(i) = cass_b(i)/(cass_b(i) + ps_mod.kD);
    mahp_e(i) = cass_c(i)/(cass_c(i) + ps_mod.kD);
  end
  subplot(2,1,1);
  plot (v,cass_a, 'k');
  hold on;
  plot (v,cass_b, 'rx');
  plot (v,cass_c, 'gx');
  plot (v,cass_d, 'm');
  plot (v,cass_e, 'b');

  subplot(2,1,2);
  plot (v,mahp_a, 'k');
  hold on;
  plot (v,mahp_b, 'rx');
  plot (v,mahp_c, 'gx');
  plot (v,mahp_d, 'm');
  plot (v,mahp_e, 'b');
  pause
end

% Setup the current vectors -- if necessary
if (isfield(ps_mod, 'I_of_t') == 0)
  ps_mod.I_of_t = zeros(3,length(0:ps_mod.dt:ps_mod.duration));
  ps_mod.I_of_t(3,round(ps_mod.I_inj_start/ps_mod.dt):round(ps_mod.I_inj_end/ps_mod.dt)) = ps_mod.I_inj_nA(3)*.001/ps_mod.Am; % convert to uA/cm2
  ps_mod.I_of_t(2,round(ps_mod.I_inj_start/ps_mod.dt):round(ps_mod.I_inj_end/ps_mod.dt)) = ps_mod.I_inj_nA(2)*.001/ps_mod.Ad; % convert to uA/cm2
  ps_mod.I_of_t(1,round(ps_mod.I_inj_start/ps_mod.dt):round(ps_mod.I_inj_end/ps_mod.dt)) = ps_mod.I_inj_nA(1)*.001/ps_mod.Aa;
end

if (0 == 1)
  figure;
  plot(0:ps_mod.dt:ps_mod.duration, ps_mod.I_of_t(1,:), 'b-'); 
  hold on;
  plot(0:ps_mod.dt:ps_mod.duration, ps_mod.I_of_t(2,:), 'r-'); 
  plot(0:ps_mod.dt:ps_mod.duration, ps_mod.I_of_t(3,:), 'k-'); 
  pause;
end


% Synaptic noise
g_syn_of_t_exc_spont = get_synaptic_noise(ps_mod.n_syns_per_facet_exc,...
                    ps_mod.spontaneous_firing_rate_exc,ps_mod.n_minis_avg_exc, ps_mod.tau_syn_exc, ...
                    ps_mod.g_max_spont_exc, ps_mod.duration, ps_mod.dt);
g_syn_of_t_inh_spont = get_synaptic_noise(ps_mod.n_syns_per_facet_inh,...
                    ps_mod.spontaneous_firing_rate_inh,ps_mod.n_minis_avg_inh, ps_mod.tau_syn_inh, ...
                    ps_mod.g_max_spont_inh, ps_mod.duration, ps_mod.dt);
ps_mod.g_syn_of_t_exc = ps_mod.g_syn_of_t_exc + g_syn_of_t_exc_spont;
ps_mod.g_syn_of_t_inh = ps_mod.g_syn_of_t_inh + g_syn_of_t_inh_spont;

% For debug
if ( 0 == 1 )
  figure;
	plot(0:ps_mod.dt:ps_mod.duration, ps_mod.g_syn_of_t_exc', 'r'); 
  hold on; 
end


% Note: ode45 does not converge, probably because the equations are stiff.
% ode15s and ode23 seem to give similar solutions
disp('starting simulation ...');
tic;
[t,y] = ode23(@red,0:ps_mod.dt:ps_mod.duration,y0,[],ps_mod); %options = []
lensec = toc;
disp([num2str(lensec) ' seconds elapsed; ' num2str(1000*(lensec/(ps_mod.duration/ps_mod.dt))) ' ms/step']);

% save?
if (isfield(ps_mod, 'savefile') ~= 0)
  save(ps_mod.savefile, 't', 'y', 'ps_mod');
end 

return

%
% reduced wang model
%

function dy = red(t,y,ps_mod)

Cm = ps_mod.Cm; %capacitance, muF/cm2

dy = zeros(12,1);

% --- current injection
Ia = ps_mod.I_of_t(1,1+round(t/ps_mod.dt));
Id = ps_mod.I_of_t(2,1+round(t/ps_mod.dt));
Im = ps_mod.I_of_t(3,1+round(t/ps_mod.dt));
 
% ---- axon compartment v - leak, Na HH, Kdr HH
IL_a = ps_mod.gLa*(y(1)-ps_mod.VL);
INa = ps_mod.gNa*y(4)^3*y(5)*(y(1)-ps_mod.VNa);
IKDR = ps_mod.gKDR*y(6)^4*(y(1)-ps_mod.VK);

dy(1) = (-IL_a - INa -IKDR - (ps_mod.gc_ma)*(y(1)-y(3)) + Ia)/ps_mod.Cm;

% --- dendritic compartment v - syn, inj, leak, H
IL_d = ps_mod.gLd*(y(2)-ps_mod.VL);
% NO I_H for now 
IH = 0;

dy(2) = (-IL_d -IH - (ps_mod.gc_md)*(y(2)-y(3)) + Id)/ps_mod.Cm;

% --- Middle compartment v - leak, 
IL_m = ps_mod.gLm*(y(3)-ps_mod.VL);
dy(3) = (-IL_m -(ps_mod.gc_am)*(y(3)-y(1)) - (ps_mod.gc_dm)*(y(3)-y(2)) + Im)/ps_mod.Cm;

% --- synaptic current 
Isyn_exc = ps_mod.g_syn_of_t_exc(1+round(t/ps_mod.dt))*(y(ps_mod.exc_syn_cmpt)-ps_mod.erev_syn_exc);
dy(ps_mod.exc_syn_cmpt) = dy(ps_mod.exc_syn_cmpt) - Isyn_exc/ps_mod.Cm;

Isyn_inh = ps_mod.g_syn_of_t_inh(1+round(t/ps_mod.dt))*(y(ps_mod.inh_syn_cmpt)-ps_mod.erev_syn_inh);
dy(ps_mod.inh_syn_cmpt) = dy(ps_mod.inh_syn_cmpt) - Isyn_inh/ps_mod.Cm;

% --- Calcium compartment
ICa = ps_mod.gCa*y(9)*(y(ps_mod.ca_cmpt)-ps_mod.VCa);
ca_av = ps_mod.sf_bapta*y(7);
IKAHP = ps_mod.gKAHP*(ca_av/(ca_av+ps_mod.kD))*(y(ps_mod.ca_cmpt)-ps_mod.VK);
dy(ps_mod.ca_cmpt) = dy(ps_mod.ca_cmpt) + (-ICa - IKAHP)/ps_mod.Cm;

%m_naa
phi_m = 4; %temperature factor
dy(4) = phi_m*(m_naa_ss(y(1))-y(4))/(ps_mod.tau_naa_sf/(a_naa(y(1)) + b_naa(y(1))));

%h_nai
phi_n = 4; %temperature factor
dy(5) = phi_m*(h_nai_ss(y(1))-y(5))/(ps_mod.tau_nai_sf/(a_nai(y(1)) + b_nai(y(1))));

%m_dra
phi_h = 4; %temperature factor
dy(6) = phi_h*(m_dra_ss(y(1))-y(6))/(ps_mod.tau_kdr_sf/(a_dra(y(1)) + b_dra(y(1))));

%calcium concentration (FREELY AVAILABLE!)
dy(7) = -ps_mod.alphaCa*ICa - y(7)/ps_mod.tauCa ;

% IH gate
dy(8) = phi_m*(m_ih_ss(y(2))-y(8))/tau_ih(y(2));

% Ca gate
dy(9) = phi_m*(m_caa_ss(y(ps_mod.ca_cmpt), ps_mod.gCaV12, ps_mod.gCaSl)-y(9))/tau_caa(y(ps_mod.ca_cmpt), ps_mod.gCaTau);

% bound bapta - DO NOT USE THIS
dy(10) = 0;

% Ca bound to AHP gate; DO NOT USE
dy(11) = 0;

% Ca bound to internal buffers - these should have a low concentration and get overwhelmed, but still delay onset of ICa - UNUSED
dy(12) = 0;



return

%steady-state sodium activation
function val = m_naa_ss(v)
  val = a_naa(v)./(a_naa(v) + b_naa(v));

% sodium activation forward rate constant    
function val = a_naa(v)
  if ( v == -33 )
    %a_naa is not defined there
    val = 1;
  else
    val= -0.1*(v+33)./(exp(-0.1*(v+33))-1);
  end
    
% sodium activation backward rate constant    
function val = b_naa(v)
  val = 4*exp(-(v+60)/12);

%steady-state sodium inactivation
function val = h_nai_ss(v)
    val = a_nai(v)./(a_nai(v) + b_nai(v));

%sodium inactivation variable forward rate constant
function val = a_nai(v)
    val = 0.07*exp(-0.1*(v+50));
    
%sodium inactivation variable backward rate constant
function val = b_nai(v) 
    val = 1/(exp(-0.1*(v+20)) + 1);
    
%steady-state delayed-rectifier activation    
function val = m_dra_ss(v)
    val = a_dra(v)./(a_dra(v) + b_dra(v));

%delayed rectifier activation variable forward rate constant
function val = a_dra(v)
  if ( v == -34 )
    val = 0.1;
  else
    val = -0.01*(v+34)./(exp(-0.1*(v+34))-1); 
  end;

%delayed rectifier activation variable backward rate constant
function val = b_dra(v)
  val = 0.125*exp(-(v+44)/25);

%steady-state calcium conductance activation
function val = m_caa_ss(v, v12, sl)
    val = 1./(1+exp((v-v12)/sl)); % Voltage 
    
% very fast, but not instantaneous
function val = tau_caa(v, tau)
  %val = tau; 
  val = (0.25)*tau + (0.75*tau)./(1+exp((v+10)/-5)); % Voltage 

%steady-state calcium concentration
function val = ca_c_ss(v, alphaCa, tauCa, VCa, gCa, Cav12, CaSl)
  %steady-state calcium current
  ICa_ss = gCa*m_caa_ss(v, Cav12, CaSl)*(v-VCa);
    
  %steady-state Ca concentration
  val = -1*tauCa*alphaCa*ICa_ss;
    
% steady-state Ih activation
function val = m_ih_ss(v)
  val = 1./(1+exp((v+92)/27.4));

% time constant for ih
function val = tau_ih(v)
  val = 1./(.0018*exp((v+86.6)/-47.5) + .0097*exp((v+85.2)/17.7));    
