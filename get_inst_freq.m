%
% Computes instantaneous firing frequency and returns it along with SD
%
function [inst_freq_vals] = get_inst_freq(time_vals, spike_idx)

  % The heart of the matter
	inst_freq_vals=zeros(1,length(time_vals));
	spike_times = time_vals(spike_idx);
	if (length(spike_times) > 1)
		for i=1:length(inst_freq_vals)
			pre_spike = find(spike_times < time_vals(i));
			pre_spike = max(pre_spike);
			post_spike = find(spike_times > time_vals(i));
			post_spike = min(post_spike);

			% Are we AT a spike?
			if (length(find(spike_times == time_vals(i))) == 1)
				if (length(pre_spike) > 0 & length(post_spike) > 0)
					inst_freq_vals(i) = 0.5/(time_vals(i)-spike_times(pre_spike)) ...
					+ 0.5/(spike_times(post_spike)-time_vals(i)); 
				% At first spike?
				elseif (length(pre_spike) == 0)
					inst_freq_vals(i) = 1/(spike_times(post_spike)-time_vals(i)); 
				% At last spike?
				elseif (length(post_spike) == 0)
					inst_freq_vals(i) = 1/(time_vals(i)-spike_times(pre_spike)); 
				end
			% Not at a spike
			else
				% Before first spike or after last?
				if (length(pre_spike) == 0 | length(post_spike) == 0) 
					inst_freq_vals(i) = 0; 
				else
					inst_freq_vals(i) = 1/(spike_times(post_spike)-spike_times(pre_spike)); 
				end
			end
		end
	end       
