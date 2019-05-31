%
% Generates the simulations for paper as mat files to be run by parallel processor
%
function generate_simulations(n_sims)
  if (exist('n_sims', 'var') == 0) 
	  n_sims = 10;
  end

	disp('Disregard any "vector too big" warnings.');

  % call this for nature neuroscience
  parpath = 'nn_par_mat/';
  save_path = 'nn_par_out/';
  natureneuro_sims(parpath, save_path, n_sims);

%
% generates for nature neuroscience
%
function natureneuro_sims(parpath, save_path, n_sims) 
  % ----------------------------------------------------------
  % 'Realistic' visual simulations -- repeat each n_sims times
  for n=1:n_sims
    % --- ones for which you have data
    % regular box
    ps_mod = get_default_psmod(1);
		ps_mod.g_syn_of_t_inh = get_synaptic_pattern(1, 1, 40 , ps_mod.tau_syn_inh, ps_mod.g_max_syn_inh,  ps_mod.duration, ...
														ps_mod.dt, 'uniform_synmap.mat', ps_mod.n_syns_per_facet_inh, ps_mod.tvnf_inh, ...
														ps_mod.tsdvnf_inh, ps_mod.velsf_inh);
		ps_mod.g_syn_of_t_exc = get_synaptic_pattern(1, 1, 40 , ps_mod.tau_syn_exc, ps_mod.g_max_syn_exc,  ps_mod.duration, ...
														ps_mod.dt, 'uniform_synmap.mat', ps_mod.n_syns_per_facet_exc, ps_mod.tvnf_exc, ...
														ps_mod.tsdvnf_exc, ps_mod.velsf_exc);
    ps_mod.savefile = [save_path '/realistic_nobapta_transbox_' num2str(n) '.mat'];
    params(1).value = ps_mod ; 
		generate(parpath, 'three_cmpt', params);

    ps_mod.tauCa = 20; 
    ps_mod.savefile = [save_path '/realistic_bapta_transbox_' num2str(n) '.mat'];
    params(1).value = ps_mod ; 
		generate(parpath, 'three_cmpt', params);
    
		% loom boxes
    L = [10 30 50];
    for l=1:length(L)
			ps_mod = get_default_psmod(1);
			ps_mod.g_syn_of_t_inh = get_synaptic_pattern(3, 1, L(l) , ps_mod.tau_syn_inh, ps_mod.g_max_syn_inh,  ps_mod.duration, ...
															ps_mod.dt, 'uniform_synmap.mat', ps_mod.n_syns_per_facet_inh, ps_mod.tvnf_inh, ...
															ps_mod.tsdvnf_inh, ps_mod.velsf_inh);
			ps_mod.g_syn_of_t_exc = get_synaptic_pattern(3, 1, L(l), ps_mod.tau_syn_exc, ps_mod.g_max_syn_exc,  ps_mod.duration, ...
															ps_mod.dt, 'uniform_synmap.mat', ps_mod.n_syns_per_facet_exc, ps_mod.tvnf_exc, ...
															ps_mod.tsdvnf_exc, ps_mod.velsf_exc);
													
			ps_mod.savefile = [save_path '/realistic_nobapta_loom_lv_' num2str(L(l)) '_' num2str(n) '.mat'];
			params(1).value = ps_mod ; generate(parpath, 'three_cmpt', params);
			ps_mod.tauCa = 20; 
			ps_mod.savefile = [save_path '/realistic_bapta_loom_lv_' num2str(L(l)) '_' num2str(n) '.mat'];
			params(1).value = ps_mod ; generate(parpath, 'three_cmpt', params);
    end
	end


%
% returns defaults
%
function ps_mod = get_default_psmod(use_set)
  % ----------------------------------------
  % model settings
 
  ps_mod = get_general_model_settings(use_set);

  % other stuff
  ps_mod.dt = .1; % ms
  ps_mod.duration = 2500; % ms 
  
  % current injection parameters
  ps_mod.g_syn_of_t_exc = zeros(1,length(0:ps_mod.dt:ps_mod.duration));
  ps_mod.g_syn_of_t_inh = zeros(1,length(0:ps_mod.dt:ps_mod.duration));
  ps_mod.I_of_t = zeros(3,length(0:ps_mod.dt:ps_mod.duration));
	ps_mod.I_inj_nA = [0 0 0]; % nA -- [ax dend mid]
	 
