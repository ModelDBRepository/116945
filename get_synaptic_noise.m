%
% Returns a spontaneous synaptic noise trace as a conductance
%
%  n_syns_per_input - how many synapses per facet (based on size of synaptic
%                     map vector).
%  single_syn_spontaneous_freq_hz - how often do single syns fire?
%  poss_lambda_single_syn - lambda parameter for spontaneous activity poiss pdf
%  syn_tau - time of peak for alpha function
%  syn_gmax - max conductance of single synapse
%  duration - of simulation
%  dt - time step size
%  synaptic_map_path - which map to use?  not specified -> use default (see code)
%
function g_syn_of_t = get_synaptic_noise(n_syns_per_input, single_syn_spontaneous_freq_hz, ...
                      poiss_lambda_single_syn, syn_tau, syn_gmax, duration, dt, ...
											synaptic_map_path)
  % synaptic map -- to determine how many synapses get triggered
  if (exist('synaptic_map_path') == 0)
    synaptic_map_path = 'uniform_synmap.mat'; % The synaptic mapping [cmpt frac az el]
  end
  load(synaptic_map_path);
	n_syns = ceil(length(synmap)*n_syns_per_input);

  % --- we derive our synapse timing by using a uniform distribution, with a frequency of
	% ---  n_syns * freq_single_syn -- that is, we want (n_syns*(duration/1000)*single_syn_freq_in_hz)
	timing = unifrnd(0,duration,1,round(n_syns*(duration/1000)*single_syn_spontaneous_freq_hz));
	
  % --- we want to have a poisson probability distribution for the number of vesicles released
  poiss_scale = poissrnd(poiss_lambda_single_syn,1,length(timing));

	% --- compute the conductance-as-function-of-time
	g_syn_of_t = zeros(1,round(duration/dt)+1);
	for t=0:dt:duration
	  active_idx = find (timing >= t-syn_tau*20 & timing <= t);
		for a = 1:length(active_idx)
		  T = t-timing(active_idx(a));
      S = poiss_scale(active_idx(a));
		  g_syn_of_t(1+round(t/dt)) = g_syn_of_t(1+round(t/dt)) + S*syn_gmax*(T/syn_tau)*exp(1-T/syn_tau);
		end
	end



