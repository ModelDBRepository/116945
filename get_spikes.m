	%
	% Extracts spike idx, based on threshold + peak method
	% 
	function spike_idx = get_spikes(thresh, data)
		greater = find(data > thresh);
		diffs = diff(data);
    greater = greater(find(greater < length(diffs)));
	
  spike_idx = [];
  greater = greater(find(greater < length(diffs)));
	
	for g=1:length(greater)-1
	  if (diffs(greater(g)) > 0 & diffs(greater(g)+1) < 0) 
		  spike_idx(length(spike_idx)+1) = greater(g);
		end
	end
	
