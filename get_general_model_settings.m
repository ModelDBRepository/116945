%
% returns default model settings -- three cmpt model
%  use_set:
%    1 = Nature Neuroscience
%    0 = Biological Cybernetics; this is the DEFAULT (i.e., values below)
%
function ps_mod = get_general_model_settings(use_set);
	% --- No feed forward inhibition
  if (use_set == 0)
		% geometry
		ps_mod.rad_m = .0004; % 5 um = .0005 cm
		ps_mod.rad_a = .0015; % 10 um
		ps_mod.len_m = .02; % 200 um ; from morphology data
		ps_mod.len_a = .1; % 1 mm = .1 cm
		ps_mod.Ad = 0.0005; % the one area that is fixed - approx 100,000 um 2 or .001 cm 2
		ps_mod.len_d = .03; % 300 um; morpho data

		% calculate variables that depend on geometry
		ps_mod = compute_geometry(ps_mod);

		% Other passive properties
		ps_mod.Cm=1.5000; % uF / cm2
		ps_mod.v_ss=-65; % mV
		ps_mod.gLd=.11; % mS / cm2 [ inverse of kOhm cm2 ] -- without iH, gLd is .222 (1000/4500)
		ps_mod.gLa=.11;
		ps_mod.gLm=.11;
		ps_mod.VL=-75;
		
		% conductances - set gXX to 0 for all to disable
		ps_mod.VNa=70;
		ps_mod.gNa=90; % 45 orig 5x9  %%% 22
		ps_mod.tau_naa_sf = .25; % speed up activation time constant
		ps_mod.tau_nai_sf = .3; % speed up activation time constant

		ps_mod.VK=-80;
		ps_mod.gKDR=22; % 18 orig 2x9
		ps_mod.tau_kdr_sf = .35; % speed up activation time constant
		ps_mod.gKAHP=50; % 5 default - 50 since Am = .1 Aa
	 
		ps_mod.VIh=-35; % gIh = g_rest_total * (1/rin_change_after_zd) * (1/p_open(rest)); -35 usually
		ps_mod.gIh=.4; % mS / cm2 - assumes ~30% open at rest, contributing to half of rm; - ~0.4
		
		% calcium stuff
		ps_mod.ca_cmpt = 3; % 3 - mid; 2 - dend ; 1 - ax -- this is tricky, because in a 2cmpt model, hard to decide;
												% ideally, we would have dend--mid--ax and mid would have this current
		ps_mod.gCa=1;
		ps_mod.gCaV12=-25; %x-25
		ps_mod.gCaSl=-3; %x-5
		ps_mod.VCa=90;
		ps_mod.gCaTau = 10; % in ms
		% ps_mod.alphaCa=0.06; % uM (ms uA)^-1 cm2 -- scalign variable for rate of influx %x.012
		ps_mod.alphaCa=0.12; % uM (ms uA)^-1 cm2 -- scalign variable for rate of influx %x.012
		ps_mod.tauCa=130; % ms - rate of efflux ; 131 in paper
		ps_mod.kD=35; % uM - binding constant for k ahp -- originally 15 -- lowering this increases ca sensitivity

		% variables for simulating bapta . . . 
		ps_mod.cBAPTA = 0; % uM BAPTA concentration
		ps_mod.tau_diff_CaBAPTA = 1; % diffusion of BAPTA*Ca away from Ca influx site modeled with this time constant
		ps_mod.kon_BAPTA = .4; % uM-1 ms-1; from Meinrenken et al. 2002
		ps_mod.koff_BAPTA = 2; % ms-1; from Meinrenken et al. 2002
		ps_mod.sf_bapta = 1; % alternate bapta implementation -- fraction of [Ca] available to KAHP channel
		 
		ps_mod.bapta_base = ps_mod.tauCa;
		ps_mod.bapta_vis = 20;
		ps_mod.bapta_cur = 20;

		% excitatory synapse stuff
		ps_mod.exc_syn_cmpt = 2; % dendrite duh!
		ps_mod.erev_syn_exc = 0;
		ps_mod.tau_syn_exc = 0.3;
		ps_mod.n_syns_exc = 15000; % DEPRECATE
		ps_mod.n_syns_per_facet_exc = 1; % number of synapses -- 'contact points'  -- per ommatidium ...
		ps_mod.spontaneous_firing_rate_exc = .05; % in Hz - for individual vesicular release
		ps_mod.g_max_syn_exc = 0.94; % mS/cm2 - for ONE synapse - .47 is you assume old data (which was 47 nS/synapse or .47mS/cm2) %%% .235
		ps_mod.n_minis_avg_exc = 1; % number of minis avg per event - lambda of poisson distro
		ps_mod.g_max_spont_exc = .0047; % mini size; .0047 if 1/100 old data (that is, scaled to SA and THEN divided by 100) %%% .0282
		ps_mod.tvnf_exc = [0 0 10 50]; % delay = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4)
		ps_mod.tsdvnf_exc = [0 20 10 10]; % delay_jitter = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4)
		ps_mod.velsf_exc = [0 1 1 1]; % velocity scaling factor = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4) - turned OFF

		% inhibitory synapse stuff
		ps_mod.inh_syn_cmpt = 3; % ca cmpt
		ps_mod.erev_syn_inh = -75;
		ps_mod.tau_syn_inh = 2;
		ps_mod.n_syns_inh = 1000; 
		ps_mod.spontaneous_firing_rate_inh = 2; % in Hz - for individual vesicular release %%% .5
		ps_mod.g_max_syn_inh = 0; % mS/cm2 - ORIGINALLY .6
		ps_mod.n_minis_avg_inh = 1; % number of minis avg per event - lambda of poisson distro
		% O: ps_mod.g_max_spont_inh = .1; % mini size; .0047 if 1/100 old data (that is, scaled to SA and THEN divided by 100) %%% 0
		ps_mod.g_max_spont_inh = 0; % mini size; .0047 if 1/100 old data (that is, scaled to SA and THEN divided by 100) %%% 0
		ps_mod.n_syns_per_facet_inh = 1;
		ps_mod.tvnf_inh = [0 20 10 50]; % delay = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4)
		ps_mod.tsdvnf_inh = [0 20 10 10]; % delay_jitter = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4)
		ps_mod.velsf_inh = [0 1 1 1]; % velocity scaling factor = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4) -- NO EFFECT

		% generic inputs vectors
		ps_mod.I_inj_nA = [0 0 0]; % nA -- [ax dend]
		ps_mod.I_inj_start = 1; % ms
		ps_mod.I_inj_end = 1; % ms

	% --- NATURE NEUROSCEINCE PAPER PARAMETERS	 -- FEED FORWARD INHIBTION
	elseif (use_set == 1)
    % NOTE - these are IGNORED!
		ps_mod.rad_m = .0004; % 5 um = .0005 cm
		ps_mod.rad_a = .0015; % 10 um
		ps_mod.len_m = .02; % 200 um ; from morphology data
		ps_mod.len_a = .1; % 1 mm = .1 cm
		ps_mod.Ad = 0.0005; % the one area that is fixed - approx 100,000 um 2 or .001 cm 2
		ps_mod.len_d = .03; % 300 um; morpho data

		% derive remaining geometric variables
		ps_mod.rad_d = ps_mod.Ad/(2*pi*ps_mod.len_d);
		ps_mod.Ad = 0.0005; 
		ps_mod.Aa = 0.00095;
		ps_mod.Am = 0.00005;
		
		% passive properties - old wang model
		ps_mod.p_am=ps_mod.Am/(ps_mod.Am + ps_mod.Aa); % axonal area / total area
		ps_mod.p_ma=ps_mod.Aa/(ps_mod.Aa + ps_mod.Am); % axonal area / total area
		ps_mod.p_dm=ps_mod.Am/(ps_mod.Am + ps_mod.Ad); % axonal area / total area
		ps_mod.p_md=ps_mod.Ad/(ps_mod.Ad + ps_mod.Am); % axonal area / total area

		
		% intercompartmental coupling - new approach
		ps_mod.rax = .060; % 60 ohm cm --> .06 kOhm cm sicne we want to invert and get mS
		rd = (ps_mod.rax * ps_mod.len_d)/(2*pi*ps_mod.rad_d*ps_mod.rad_d);
		ra = (ps_mod.rax * ps_mod.len_a)/(2*pi*ps_mod.rad_a*ps_mod.rad_a);
		rm = (ps_mod.rax * ps_mod.len_m)/(2*pi*ps_mod.rad_m*ps_mod.rad_m);
		ps_mod.gc_dm = 16.12;
		ps_mod.gc_am = 12.3;
		ps_mod.gc_md = 1.62;
		ps_mod.gc_ma = .66;
		% calculate variables that depend on geometry

		% Other passive properties
		ps_mod.Cm=1.5000; % uF / cm2
		ps_mod.v_ss=-65; % mV
		ps_mod.gLd=.11; % mS / cm2 [ inverse of kOhm cm2 ] -- without iH, gLd is .222 (1000/4500)
		ps_mod.gLa=.11; %
		ps_mod.gLm=.11; %
		ps_mod.VL=-75; %
		
		% conductances - set gXX to 0 for all to disable
		ps_mod.VNa=70;
		ps_mod.gNa=22; % 45 orig 5x9
		ps_mod.tau_naa_sf = .25; % speed up activation time constant
		ps_mod.tau_nai_sf = .3; % speed up activation time constant

		ps_mod.VK=-80;
		ps_mod.gKDR=22; % 18 orig 2x9
		ps_mod.tau_kdr_sf = .35; % speed up activation time constant
		ps_mod.gKAHP=50; % 5 default - 50 since Am = .1 Aa
	 
		ps_mod.VIh=-35; % gIh = g_rest_total * (1/rin_change_after_zd) * (1/p_open(rest)); -35 usually
		ps_mod.gIh=.4; % mS / cm2 - assumes ~30% open at rest, contributing to half of rm; - ~0.4
		
		% calcium stuff
		ps_mod.ca_cmpt = 3; % 3 - mid; 2 - dend ; 1 - ax -- this is tricky, because in a 2cmpt model, hard to decide;
												% ideally, we would have dend--mid--ax and mid would have this current
		ps_mod.gCa=1;
		ps_mod.gCaV12=-25; %x-25
		ps_mod.gCaSl=-3; %x-5
		ps_mod.VCa=90;
		ps_mod.gCaTau = 10; % in ms
		ps_mod.alphaCa=0.12; % uM (ms uA)^-1 cm2 -- scalign variable for rate of influx %x.012
		ps_mod.tauCa=132; % ms - rate of efflux ; 131 in paper
		ps_mod.kD=35; % uM - binding constant for k ahp -- originally 15 -- lowering this increases ca sensitivity

		% variables for simulating bapta . . . 
		ps_mod.cBAPTA = 0; % uM BAPTA concentration
		ps_mod.tau_diff_CaBAPTA = 1; % diffusion of BAPTA*Ca away from Ca influx site modeled with this time constant
		ps_mod.kon_BAPTA = .4; % uM-1 ms-1; from Meinrenken et al. 2002
		ps_mod.koff_BAPTA = 2; % ms-1; from Meinrenken et al. 2002
		ps_mod.sf_bapta = 1; % alternate bapta implementation -- fraction of [Ca] available to KAHP channel
		 
		ps_mod.bapta_base = ps_mod.tauCa;
		ps_mod.bapta_vis = 20;
		ps_mod.bapta_cur = 20;

		% excitatory synapse stuff
		ps_mod.exc_syn_cmpt = 2; % dendrite duh!
		ps_mod.erev_syn_exc = 0;
		ps_mod.tau_syn_exc = 0.3;
		ps_mod.n_syns_per_facet_exc = 5; % number of synapses -- 'contact points'  -- per ommatidium ...
		ps_mod.spontaneous_firing_rate_exc = .05; % in Hz - for individual vesicular release
		ps_mod.g_max_syn_exc = .47; % mS/cm2 - for ONE synapse - .47 is you assume old data (which was 47 nS/synapse or .47mS/cm2)
		ps_mod.n_minis_avg_exc = 1; % number of minis avg per event - lambda of poisson distro
		ps_mod.g_max_spont_exc = .047; % mini size; .0047 if 1/100 old data (that is, scaled to SA and THEN divided by 100)
		ps_mod.tvnf_exc = [50 0 10 50]; % delay = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4)
		ps_mod.tsdvnf_exc = [50 5 10 5]; % delay_jitter = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4)
		ps_mod.velsf_exc = [0 1 1 1]; % velocity scaling factor = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4) - turned OFF

		% inhibitory synapse stuff
		ps_mod.inh_syn_cmpt = 3; % ca cmpt
		ps_mod.erev_syn_inh = -75;
		ps_mod.tau_syn_inh = 2;
		ps_mod.n_syns_per_facet_inh = 1;
		ps_mod.g_max_syn_inh = 6; % mS/cm2 - ORIGINALLY .6 10 works for loom
		ps_mod.spontaneous_firing_rate_inh = .5; % in Hz - for individual vesicular release
		ps_mod.n_minis_avg_inh = 1; % number of minis avg per event - lambda of poisson distro
		ps_mod.g_max_spont_inh = 0.05; % mini size; .0047 if 1/100 old data (that is, scaled to SA and THEN divided by 100)
		ps_mod.tvnf_inh = [120 20 100 80]; % delay = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4)
		ps_mod.tsdvnf_inh = [5 20 10 1]; % delay_jitter = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4)
		ps_mod.velsf_inh = [0 1 1 1]; % velocity scaling factor = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4) -- NO EFFECT

		% generic inputs vectors
		ps_mod.I_inj_nA = [0 0 0]; % nA -- [ax dend]
		ps_mod.I_inj_start = 1; % ms
		ps_mod.I_inj_end = 1; % ms

		ps_mod.g_max_syn_exc = .47*.5; % mS/cm2 - for ONE synapse - .47 is you assume old data (which was 47 nS/synapse or .47mS/cm2)
		ps_mod.g_max_spont_exc = .047*.6; % mini size; .0047 if 1/100 old data (that is, scaled to SA and THEN divided by 100)
		ps_mod.g_max_syn_inh = 6*.3; % mS/cm2 - ORIGINALLY .6 10 works for loom
		ps_mod.g_max_spont_inh = .05*.2; % mini size; .0047 if 1/100 old data (that is, scaled to SA and THEN divided by 100)

    ps_mod.tsdvnf_exc = [0 1 1 5]; % I did this to not rely on variable input
  end

%
% geometry-dependent variables
%
function ps_mod = compute_geometry(ps_mod)
	% derive remaining geometric variables
	ps_mod.rad_d = ps_mod.Ad/(2*pi*ps_mod.len_d);
  ps_mod.Aa=2*pi*ps_mod.rad_a*ps_mod.len_a;
  ps_mod.Am=2*pi*ps_mod.rad_m*ps_mod.len_m;
  
  % passive properties - old wang model
  ps_mod.p_am=ps_mod.Am/(ps_mod.Am + ps_mod.Aa); % axonal area / total area
  ps_mod.p_ma=ps_mod.Aa/(ps_mod.Aa + ps_mod.Am); % axonal area / total area
  ps_mod.p_dm=ps_mod.Am/(ps_mod.Am + ps_mod.Ad); % axonal area / total area
  ps_mod.p_md=ps_mod.Ad/(ps_mod.Ad + ps_mod.Am); % axonal area / total area

  
	% intercompartmental coupling - new approach
	ps_mod.rax = .060; % 60 ohm cm --> .06 kOhm cm sicne we want to invert and get mS
	rd = (ps_mod.rax * ps_mod.len_d)/(2*pi*ps_mod.rad_d*ps_mod.rad_d);
	ra = (ps_mod.rax * ps_mod.len_a)/(2*pi*ps_mod.rad_a*ps_mod.rad_a);
	rm = (ps_mod.rax * ps_mod.len_m)/(2*pi*ps_mod.rad_m*ps_mod.rad_m);
	ps_mod.gc_dm = (1/(rm + rd))*(1/ps_mod.Am); % gc_source,target
	ps_mod.gc_am = (1/(rm + ra))*(1/ps_mod.Am);
	ps_mod.gc_md = (1/(rm + rd))*(1/ps_mod.Ad);
	ps_mod.gc_ma = (1/(rm + ra))*(1/ps_mod.Aa);

