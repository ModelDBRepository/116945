% This will return a conductance-as-function-of-time for treating dendrites
%  as single compartment
%
%   class: 1 = 10x10 box
%          2 = 10x80 bar
%          3 = looming
%          4 = accelarating box
%          5 = loomlike box
%   velocity - in deg/s -- for class = 4, specify start and end velocity in 2-element vector
%   velocity - in deg/s -- for class = 5, specify the angle
%   direction: 1 = AP, -1 = PA, 2=DV, -2=VD
%   dt: time step for *simulation*
%   synaptic_map_path: which map to use? uniform is current dflt
%   n_syns_per_input: how many synapses *per facet*?
%
% The following variables alter timing to account for feedforward inhibition and the observation
%  that larger inputs are more synchronized.  
%
%   dvn: parameters of delay-versus-nfacets function - delay = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4)
%   sdvn: parameters of sd of delay-versus-nsyns function - delay_jitter = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4)
%         note that the actual delay is jittered by drawing from a normal distro with the above sd;
%         you read it correctly -- dvn is based on the number of facets, while sdvn is based on the
%         number of synapses
%   vel_sf: this accounts for the observation that there is speed modulation speed_sf = v(1)./(1+exp((nsyns-v(2))/v(3)))+v(4)
%
function g_syn_of_t = get_synaptic_pattern(class, direction, velocity, syn_tau, syn_gmax, ...
                      duration, dt, synaptic_map_path, n_syns_per_input, dvn, sdvn, vel_sf)
  % --- Definitions
  % synaptic map -- to determine how many synapses get triggered
  if (exist('synaptic_map_path') == 0 | synaptic_map_path == -1)
    synaptic_map_path = 'uniform_synmap.mat'; % The synaptic mapping [cmpt frac az el]
  end
	dt_image = 5; % screen inter-step time; monitor specific 
  half_angle = 2.5; % acceptance half-angle in degrees; add

  % for translating stimuli:
  edge_sf = [1 0]; % apply 1 to the leading edge, 2 to the trailing edge [1 0], for instance, means leading only

  % --- Preliminaries
  load(synaptic_map_path);
	n_syns = length(synmap);
	
  % --- derive the timing pattern
	switch (class)
	  case 1 % 10x10 box
			vert_lims = [-5-half_angle 5+half_angle];
			hor_lims = [85-half_angle 95+half_angle];
			bar_thickness = 10+(2*half_angle); % in dimension parallel to that of travel
		  switch(direction)
			  case 1 % AP
          [synpos, timing, inst_vel] = pattern_translating_bar(synmap, 'v', velocity, [50 130], dt_image, bar_thickness, vert_lims, edge_sf);
				case -1 % PA
					[synpos, timing, inst_vel] = pattern_translating_bar(synmap, 'v', velocity, [130 50], dt_image, bar_thickness, vert_lims, edge_sf);
				case 2 % DV
					[synpos, timing, inst_vel] = pattern_translating_bar(synmap, 'h', velocity, [40 -40], dt_image, bar_thickness, hor_lims, edge_sf);
				case -2 % VD
					[synpos, timing, inst_vel] = pattern_translating_bar(synmap, 'h', velocity, [-40 40], dt_image, bar_thickness, hor_lims, edge_sf);
		  end;
	  case 1.1 % 10x10 box, but with bounds that are used for the SFA 'review'/model paper
			vert_lims = [-5-half_angle 5+half_angle];
			hor_lims = [85-half_angle 95+half_angle];
			bar_thickness = 10+(2*half_angle); % in dimension parallel to that of travel
		  switch(direction)
			  case 1 % AP
          [synpos, timing, inst_vel] = pattern_translating_bar(synmap, 'v', velocity, [60 120], dt_image, bar_thickness, vert_lims, edge_sf);
				case -1 % PA
					[synpos, timing, inst_vel] = pattern_translating_bar(synmap, 'v', velocity, [120 60], dt_image, bar_thickness, vert_lims, edge_sf);
				case 2 % DV
					[synpos, timing, inst_vel] = pattern_translating_bar(synmap, 'h', velocity, [30 -30], dt_image, bar_thickness, hor_lims, edge_sf);
				case -2 % VD
					[synpos, timing, inst_vel] = pattern_translating_bar(synmap, 'h', velocity, [-30 30], dt_image, bar_thickness, hor_lims, edge_sf);
		  end;
	  case 2 % 10x80 bar
			vert_lims = [-40-half_angle 40+half_angle];
			hor_lims = [50-half_angle 130+half_angle];
			bar_thickness = 10+(2*half_angle); % in dimension parallel to that of travel
		  switch(direction)
			  case 1 % AP
          [synpos, timing, inst_vel] = pattern_translating_bar(synmap, 'v', velocity, [50 130], dt_image, bar_thickness, vert_lims, edge_sf);
				case -1 % PA
					[synpos, timing, inst_vel] = pattern_translating_bar(synmap, 'v', velocity, [130 50], dt_image, bar_thickness, vert_lims, edge_sf);
				case 2 % DV
					[synpos, timing, inst_vel] = pattern_translating_bar(synmap, 'h', velocity, [40 -40], dt_image, bar_thickness, hor_lims, edge_sf);
				case -2 % VD
					[synpos, timing, inst_vel] = pattern_translating_bar(synmap, 'h', velocity, [-40 40], dt_image, bar_thickness, hor_lims, edge_sf);
		  end;
    case 3 % loom l/v = velocity
      % calculate so that you finish with 80 deg, then go for 100 ms more
      end_time = velocity/(tand(80/2));
      start_time = duration - (end_time + 100);
      theta_start = 2*atand(velocity/start_time);
      [synpos, timing, inst_vel] = pattern_circular_looming_stimulus(synmap, velocity, [0 90], theta_start, dt_image);
	  case 4 % 10x10 box, accelarating linearly
			vert_lims = [-5-half_angle 5+half_angle];
			hor_lims = [85-half_angle 95+half_angle];
			bar_thickness = 10+(2*half_angle); % in dimension parallel to that of travel
		  switch(direction)
			  case 1 % AP
          [synpos, timing, inst_vel] = pattern_accelarating_bar(synmap, 'v', velocity, [60 120], 5, bar_thickness, vert_lims, edge_sf);
				case -1 % PA
					[synpos, timing, inst_vel] = pattern_accelarating_bar(synmap, 'v', velocity, [120 60], 5, bar_thickness, vert_lims, edge_sf);
				case 2 % DV
					[synpos, timing, inst_vel] = pattern_accelarating_bar(synmap, 'h', velocity, [30 -30], 5, bar_thickness, hor_lims, edge_sf);
				case -2 % VD
					[synpos, timing, inst_vel] = pattern_accelarating_bar(synmap, 'h', velocity, [-30 30], 5, bar_thickness, hor_lims, edge_sf);
		  end;
    case 5 % 10x10 box, loomlike linearly
			vert_lims = [-5-half_angle 5+half_angle];
			hor_lims = [85-half_angle 95+half_angle];
			bar_thickness = 10+(2*half_angle); % in dimension parallel to that of travel
		  switch(direction)
			  case 1 % AP
          [synpos, timing, inst_vel] = pattern_loomlike_bar(synmap, 'v', velocity, [60 120], 5, bar_thickness, vert_lims, edge_sf);
				case -1 % PA
					[synpos, timing, inst_vel] = pattern_loomlike_bar(synmap, 'v', velocity, [120 60], 5, bar_thickness, vert_lims, edge_sf);
				case 2 % DV
					[synpos, timing, inst_vel] = pattern_loomlike_bar(synmap, 'h', velocity, [30 -30], 5, bar_thickness, hor_lims, edge_sf);
				case -2 % VD
					[synpos, timing, inst_vel] = pattern_loomlike_bar(synmap, 'h', velocity, [-30 30], 5, bar_thickness, hor_lims, edge_sf);
		  end;
  end;

	% --- the preceding logic returns a timing vector, which tells you the timing of ommatidial
	% --- activation.  however, response at the level of the LGMD is non-linear; specifically,
	% --- synaptic resposne increases superlinearly and occurs faster for larger luminance changes;
	% --- this is true for both inhibitory and excitatory inputs.  We will do this by making the 
	% --- assumption that intensity change is approximated by the fraction of synapses activated per timestep

  % --- for each element of timing, calculate the scaling factor based on velocity ...
	% --- this must be done BEFORE any other timing alterations as it depends on the OMMATIDIAL
	% --- transition time, which is not the same as the timing vetor after considering delay
	vsf = ones(size(timing));
	tv = inst_vel.t_vec;
	for t=0:dt:duration
	  idx = find (timing >= t & timing < t + dt_image);
		ivi = min(find(tv >= t-1));
		if (length(ivi) == 0)  
		  %disp(['warning: ivi not found; t = ' num2str(t)]) ; 
			continue; 
	  end
		vsf(idx) = vel_sf(1)./(1+exp((inst_vel.vel_vec(ivi)-vel_sf(2))/vel_sf(3)))+vel_sf(4);
%if(length(idx) > 0) ;		disp(['vsf: ' num2str(vsf(idx(1))) ' vel: ' num2str(inst_vel.vel_vec(ivi)) ' t: ' num2str(t) ' ivi: ' num2str(ivi)]); end
	end

	% --- first, introduce mean delay via dvn based on number of facets hit in a given timestep
	n_timing = timing; 
	for t=0:dt_image:duration
	  idx = find(timing >= t & timing < t + dt_image); 
	  nfacets = length(idx);
		n_timing(idx) = timing(idx) + dvn(1)./(1+exp((nfacets-dvn(2))/dvn(3)))+dvn(4);
	end
	timing = n_timing;


	% --- Now account for the fact that you do not necessarily have the same number
	% --- of synapses as facests
  n_timing = zeros(1,floor(length(synpos)*n_syns_per_input));
	cf = 1;
	for t=0:dt_image:duration
	  idx = find(timing >= t & timing < t + dt_image); 
	  nfacets = length(idx);
		nsyns = nfacets*n_syns_per_input;

    % For a non-integer number of synapses, scale based on a normal pdf
	  if (nsyns ~= round(nsyns))
		  rem = nsyns - floor(nsyns);
			if (rem > rand) 
			  nsyns = floor(nsyns) + 1;
			else
			  nsyns = floor(nsyns);
			end
    end
	  if (nsyns ~= round(nsyns)) ; disp (['wtf:' num2str(nsyns) ' rem: ' num2str(rem)]);  end

    % Now make nsyns in n_timing have the time
    n_timing(cf:min(length(n_timing),cf+nsyns)) = t;
		if (length(idx) > 0)
			n_vsf(cf:min(length(n_timing),cf+nsyns)) = mean (vsf(idx));
		else
			n_vsf(cf:min(length(n_timing),cf+nsyns)) = 1;
		end
		cf = cf + nsyns;
		if (cf > length(n_timing))
		  disp('Warning: vector too big');
			cf = length(n_timing);
		end
	end
	timing = n_timing(find(n_timing > 0));
	vsf = n_vsf;

  % --- it is still unknown if there is an actual increase in input strength from stronger stimuli;
	% --- what is certain is that jitter declines.  do jitter here. *still based on number of facets*
	n_timing = timing; 
	for t=0:dt_image:duration
	  idx = find(timing >= t & timing < t + dt_image); 
	  nsyns = length(idx);
		jitter = normrnd(0, sdvn(1)./(1+exp(nsyns-sdvn(2))/sdvn(3))+sdvn(4), length(idx), 1);
		if (length(idx) > 0) ; n_timing(idx) = timing(idx) + jitter'; end
	end
	timing = n_timing;

	% --- compute the conductance-as-function-of-time
	g_syn_of_t = zeros(1,round(duration/dt)+1);
	for t=0:dt:duration
	  active_idx = find (timing >= t-syn_tau*20 & timing <= t);
		for a = 1:length(active_idx)
		  T = t-timing(active_idx(a));
			S = vsf(active_idx(a));
		  g_syn_of_t(1+round(t/dt)) = g_syn_of_t(1+round(t/dt)) + S*syn_gmax*(T/syn_tau)*exp(1-T/syn_tau);
		end
	end
	
	% --- sanity check
	if (0 == 1)
	  plot(0:dt:duration, g_syn_of_t);
  end
	
% 
% This should return timing and synapses for loomlike bars of CONSTANT 
%  ANGULAR SIZE AND CHANGING ANGULAR VELOCITY.  Each position is activated
%  TWICE - once  for the ON and once for the OFF.
%
% RETURNS:
%  synpos - [cmpt_idx cmpt_fraction]
%  timing - in ms
%
% ASSIGNED:
%  synmap - [cmpt cmpt_frac az el]
%  direction - 'h' 'v' for horizontal, vertical (DIRECTION OF MOTION, not bar
%              orientation!) - default is vertical
%  velocity - in degrees / second -- l/v ; negative --> receding
%  limits - [start end] where to start; if horizontal motion, give azimuths; 
%           if vertical, elevations. 
%  dt - temporal resolution, in ms
%  thick - thickness of the bar, in degrees (i.e., dimension along axis PARALLEL
%          to motion vector)
%  length_limits - edges of the bar, in degrees for dimension along axis PERPENDICULAR to
%           motion vector
%  
%
function [synpos, timing, inst_vel] = pattern_loomlike_bar(synmap, direction, velocity, ...
                                                 limits, dt, thick, length_limits, ...
                                                 edge_sf)
  % equivalent disc RADIUS at which to start
  theta_start = 5;
 
  % Preliminary
  timing = [];

  % Orientation-dependent component
  if (strcmp('h', direction))
    positional_array = synmap(:,4); 
    length_array = synmap(:,3); 
  else
    positional_array = synmap(:,3); 
    length_array = synmap(:,4); 
  end
  used = zeros(1,length(positional_array));

	% Introduce variation into the eye -- vary each ommatidial view by sd of 5 degrees
  eye_deformation = normrnd(0,10,length(positional_array),1);
  positional_array = positional_array + eye_deformation;
  
  % Determine valid based on length limits
  length_valids = find(length_array > min(length_limits) & length_array < max(length_limits));

  % determine total time based on the limits and the accelaration
  % loom time
  t_o = velocity/tand(theta_start); 
  t_f = velocity/tand(diff(limits) + theta_start);
  if (velocity < 0) ; t_o = -1*t_o; t_f = -1*t_f ; end
  total_time = 5*(1+floor(abs(t_f-t_o)/5));

  % instantaneous velocity vector
  inst_vel.vel_vec = zeros(length(0:dt:total_time),1);
	ivi = 1;

  % loop thru time, and add accordingly
  position = limits(1); % LEADING edge
  for time=0:dt:total_time
    last_position = position;  
    if (velocity < 0) % 'receding'
      theta_t = atand(velocity/(t_f+time));
      theta_tp1 = atand(velocity/(t_f+(time+dt)));
    else % normal
      theta_t = atand(velocity/(t_o-time));
      theta_tp1 = atand(velocity/(t_o-(time+dt)));
    end
    dp = abs(theta_tp1-theta_t);
    
    % change in position/time is velocity!
    int_vel.vel_vec = dp;

    if (time+dt > t_o) ; dp = 0; end

    stimulated_leading = [];
    stimulated_trailing = [];
    
    % Increasing position value (e.g., rightward)
    if (limits(1) < limits(2))
      position = position + dp;
      % Leading edge
      if (position < limits(2))
        stimulated_leading = find(positional_array >= last_position & positional_array < position);
      end
      % Trailing edge
      if (position-thick > limits(1))
        stimulated_trailing = find(positional_array >= last_position-thick & positional_array < position-thick);
      end
    % For decreasing position value
    else 
      position = position - dp;
      % Leading edge
      if (position < limits(1))
        stimulated_leading = find(positional_array < last_position & positional_array >= position);
      end
      % Trailing edge
      if (position+thick > limits(2))
        stimulated_trailing = find(positional_array < last_position+thick & positional_array >= position+thick);
      end
    end
 
    % Constrain by length
    stimulated_trailing = intersect(stimulated_trailing,length_valids);
    stimulated_leading = intersect(stimulated_leading,length_valids);
    if (edge_sf(1) > 0)
      stimulated_leading = stimulated_leading(1:round(length(stimulated_leading)*edge_sf(1)));
    else
      stimulated_leading = [];
    end
    if (edge_sf(2) > 0)
      stimulated_trailing = stimulated_trailing(1:round(length(stimulated_trailing)*edge_sf(2)));
    else
      stimulated_trailing = [];
    end
      
    % And append returned arrays if needbe
    if (length(stimulated_leading) > 0)
      for s=1:length(stimulated_leading)
        synpos(length(timing) + 1, :) = [synmap(stimulated_leading(s),1) synmap(stimulated_leading(s),2)];
        timing(length(timing) + 1) = time;
      end
    end
    if (length(stimulated_trailing) > 0)
      for s=1:length(stimulated_trailing)
        synpos(length(timing) + 1, :) = [synmap(stimulated_trailing(s),1) synmap(stimulated_trailing(s),2)];
        timing(length(timing) + 1) = time;
      end
    end
  end  

  % Scale instantaneous velocity by dt -- in s
  inst_vel.t_vec = 0:dt:total_time;
  inst_vel.vel_vec = abs([inst_vel.vel_vec/(.001*dt)])';
	
% 
% This should return timing and synapses for accelarating bars of CONSTANT 
%  ANGULAR SIZE AND CHANGING ANGULAR VELOCITY.  Each position is activated
%  TWICE - once  for the ON and once for the OFF.
%
% RETURNS:
%  synpos - [cmpt_idx cmpt_fraction]
%  timing - in ms
%
% ASSIGNED:
%  synmap - [cmpt cmpt_frac az el]
%  direction - 'h' 'v' for horizontal, vertical (DIRECTION OF MOTION, not bar
%              orientation!) - default is vertical
%  velocity - in degrees / second (1): start (2): end
%  limits - [start end] where to start; if horizontal motion, give azimuths; 
%           if vertical, elevations. Bar will 'edge in' from one side and edge
%           out the other (i.e., it does not start wholly formed).
%  dt - temporal resolution, in ms
%  thick - thickness of the bar, in degrees (i.e., dimension along axis PARALLEL
%          to motion vector)
%  length_limits - edges of the bar, in degrees for dimension along axis PERPENDICULAR to
%           motion vector
%  
%
function [synpos, timing, inst_vel] = pattern_accelarating_bar(synmap, direction, velocity, ...
                                                    limits, dt, thick, length_limits, ...
                                                    edge_sf)
  % Preliminary
  timing = [];

  % Orientation-dependent component
  if (strcmp('h', direction))
    positional_array = synmap(:,4); 
    length_array = synmap(:,3); 
  else
    positional_array = synmap(:,3); 
    length_array = synmap(:,4); 
  end
  used = zeros(1,length(positional_array));

	% Introduce variation into the eye -- vary each ommatidial view by sd of 5 degrees
  eye_deformation = normrnd(0,10,length(positional_array),1);
  positional_array = positional_array + eye_deformation;
  
  % Determine valid based on length limits
  length_valids = find(length_array > min(length_limits) & length_array < max(length_limits));
  
  % determine total time based on the limits and the accelaration
  % total time = total_distance/mean_velocity [assuming LINEAR accelaration]
  mean_velocity = mean(velocity);
  total_time = abs(abs(diff(limits)))/(mean_velocity/1000);

  % instantaneous velocity vector
  inst_vel.vel_vec = zeros(length(0:dt:total_time),1);
	ivi = 1;
 
  % loop thru time, and add accordingly
  position = limits(1); % LEADING edge
  dv_dt = diff(velocity)/1000/total_time; % again, convert to deg/ms2
  v_o = velocity(1)/1000; % convert to d/ms
  for time=0:dt:total_time
    last_position = position;  
    vel_t = dv_dt*time + v_o;
    dp = vel_t*dt;    
    
    % change in position/time is velocity!
    int_vel.vel_vec = dp;

    stimulated_leading = [];
    stimulated_trailing = [];
    
    % Increasing position value (e.g., rightward)
    if (limits(1) < limits(2))
      position = position + dp;
      % Leading edge
      if (position < limits(2))
        stimulated_leading = find(positional_array >= last_position & positional_array < position);
      end
      % Trailing edge
      if (position-thick > limits(1))
        stimulated_trailing = find(positional_array >= last_position-thick & positional_array < position-thick);
      end
    % For decreasing position value
    else 
      position = position - dp;
      % Leading edge
      if (position < limits(1))
        stimulated_leading = find(positional_array < last_position & positional_array >= position);
      end
      % Trailing edge
      if (position+thick > limits(2))
        stimulated_trailing = find(positional_array < last_position+thick & positional_array >= position+thick);
      end
    end
 
    % Constrain by length
    stimulated_trailing = intersect(stimulated_trailing,length_valids);
    stimulated_leading = intersect(stimulated_leading,length_valids);
    if (edge_sf(1) > 0)
      stimulated_leading = stimulated_leading(1:round(length(stimulated_leading)*edge_sf(1)));
    else
      stimulated_leading = [];
    end
    if (edge_sf(2) > 0)
      stimulated_trailing = stimulated_trailing(1:round(length(stimulated_trailing)*edge_sf(2)));
    else
      stimulated_trailing = [];
    end
      
    % And append returned arrays if needbe
    if (length(stimulated_leading) > 0)
      for s=1:length(stimulated_leading)
        synpos(length(timing) + 1, :) = [synmap(stimulated_leading(s),1) synmap(stimulated_leading(s),2)];
        timing(length(timing) + 1) = time;
      end
    end
    if (length(stimulated_trailing) > 0)
      for s=1:length(stimulated_trailing)
        synpos(length(timing) + 1, :) = [synmap(stimulated_trailing(s),1) synmap(stimulated_trailing(s),2)];
        timing(length(timing) + 1) = time;
      end
    end
  end  

  % Scale instantaneous velocity by dt -- in s
  inst_vel.t_vec = 0:dt:total_time;
  inst_vel.vel_vec = abs([inst_vel.vel_vec/(.001*dt)])';

	
% 
% This should return timing and synapses for translating bars of CONSTANT 
%  ANGULAR SIZE AND ANGULAR VELOCITY.  Each position is activated TWICE - once
%  for the ON and once for the OFF.
%
% RETURNS:
%  synpos - [cmpt_idx cmpt_fraction]
%  timing - in ms
%
% ASSIGNED:
%  synmap - [cmpt cmpt_frac az el]
%  direction - 'h' 'v' for horizontal, vertical (DIRECTION OF MOTION, not bar
%              orientation!) - default is vertical
%  velocity - in degrees / second
%  limits - [start end] where to start; if horizontal motion, give azimuths; 
%           if vertical, elevations. Bar will 'edge in' from one side and edge
%           out the other (i.e., it does not start wholly formed).
%  dt - temporal resolution, in ms
%  thick - thickness of the bar, in degrees (i.e., dimension along axis PARALLEL
%          to motion vector)
%  length_limits - edges of the bar, in degrees for dimension along axis PERPENDICULAR to
%           motion vector
%  edge_sf - scales the input strength of each edge
%  
%
function [synpos, timing, inst_vel] = pattern_translating_bar(synmap, direction, velocity, ...
                                                    limits, dt, thick, length_limits, ...
                                                    edge_sf )
  % Preliminary
  timing = [];

  % Orientation-dependent component
  if (strcmp('h', direction))
    positional_array = synmap(:,4); 
    length_array = synmap(:,3); 
  else
    positional_array = synmap(:,3); 
    length_array = synmap(:,4); 
  end
  used = zeros(1,length(positional_array));

	% Introduce variation into the eye -- vary each ommatidial view by sd of 5 degrees
  eye_deformation = normrnd(0,10,length(positional_array),1);
  positional_array = positional_array + eye_deformation;
  
  % Determine valid based on length limits
  length_valids = find(length_array > min(length_limits) & length_array < max(length_limits));
  
  % determine total time ... limits + thickness scaled by velocity
  total_time = abs(abs(diff(limits)) + thick)/(velocity/1000);

  % loop thru time, and add accordingly
  position = limits(1); % LEADING edge
  dp = (velocity/1000)*dt;
  for time=0:dt:total_time
    last_position = position;

    stimulated_leading = [];
    stimulated_trailing = [];
    
    % Increasing position value (e.g., rightward)
    if (limits(1) < limits(2))
      position = position + dp;
      % Leading edge
      if (position < limits(2))
        stimulated_leading = find(positional_array >= last_position & positional_array < position);
      end
      % Trailing edge
      if (position-thick > limits(1))
        stimulated_trailing = find(positional_array >= last_position-thick & positional_array < position-thick);
      end
    % For decreasing position value
    else 
      position = position - dp;
      % Leading edge
      if (position < limits(1))
        stimulated_leading = find(positional_array < last_position & positional_array >= position);
      end
      % Trailing edge
      if (position+thick > limits(2))
        stimulated_trailing = find(positional_array < last_position+thick & positional_array >= position+thick);
      end
    end
 
    % Constrain by length
    stimulated_trailing = intersect(stimulated_trailing,length_valids);
    stimulated_leading = intersect(stimulated_leading,length_valids);
    if (edge_sf(1) > 0)
      stimulated_leading = stimulated_leading(1:round(length(stimulated_leading)*edge_sf(1)));
    else
      stimulated_leading = [];
    end
    if (edge_sf(2) > 0)
      stimulated_trailing = stimulated_trailing(1:round(length(stimulated_trailing)*edge_sf(2)));
    else
      stimulated_trailing = [];
    end
    
    % And append returned arrays if needbe
    if (length(stimulated_leading) > 0)
      for s=1:length(stimulated_leading)
        synpos(length(timing) + 1, :) = [synmap(stimulated_leading(s),1) synmap(stimulated_leading(s),2)];
        timing(length(timing) + 1) = time;
      end
    end
    if (length(stimulated_trailing) > 0)
      for s=1:length(stimulated_trailing)
        synpos(length(timing) + 1, :) = [synmap(stimulated_trailing(s),1) synmap(stimulated_trailing(s),2)];
        timing(length(timing) + 1) = time;
      end
    end
  end  

  % Scale instantaneous velocity by dt -- in s
  inst_vel.t_vec = 0:dt:total_time;
  inst_vel.vel_vec = abs(velocity*ones(size(inst_vel.t_vec)));

% 
% This should return timing and synapses for looming stimulus. Each position
%  is activated only once.  The stimulus is square.
%  synpos - [cmpt_idx cmpt_fraction]
%  timing - in ms
%  synmap - [cmpt cmpt_frac az el]
%  l_over_v - the l/v parameter
%  center - [el az] where to start
%  theta_start - size at start; zero would take infinite time. in degrees.
%  dt - temporal resolution, in ms
%
function [synpos, timing] = pattern_square_looming_stimulus(synmap, l_over_v, center, theta_start, dt)
  % Preliminary
  timing = [];
  theta_final = 80; % Degrees - final size, as in screen

  % Start, end time - compute it
  start_time = l_over_v/(tand(theta_start/2));
  end_time = l_over_v/(tand(theta_final/2));

  % GEnerate a center-centric coordinate set
  centered_el = synmap(:,4)-center(1);
  centered_az = synmap(:,3)-center(2);

  % determine total time ...
  total_time = start_time - end_time;
 
  % loop thru time, and add accordingly
  theta_0 = theta_start; % Theta at START of step
  for time=0:dt:total_time
    theta_f = 2*atand(l_over_v/(total_time + end_time - time));

    % Find any points between theta_0 and theta_f object size
    stimulated_1 = find(centered_el < theta_f & centered_el >= (-1*theta_f) & centered_az > theta_0 & centered_az <= theta_f);
    stimulated_2 = find(centered_el < theta_f & centered_el >= (-1*theta_f) & centered_az < -1*theta_0 & centered_az >= -1*theta_f);
    stimulated_3 = find(centered_az > -1*theta_0 & centered_az < theta_0 & centered_el >= -1*theta_f & centered_el < -1*theta_0);
    stimulated_4 = find(centered_az > -1*theta_0 & centered_az < theta_0 & centered_el <= theta_f & centered_el > theta_0);

    stimulated = unique([stimulated_1' stimulated_2' stimulated_3' stimulated_4']);
    
    % And append returned arrays if needbe
    if (length(stimulated) > 0)
      for s=1:length(stimulated)
        synpos(length(timing) + 1, :) = [synmap(stimulated(s),1) synmap(stimulated(s),2)];
        timing(length(timing) + 1) = time;
      end
    end
    theta_0 = theta_f;
  end  

% 
% This should return timing and synapses for looming stimulus. Each position
%  is activated only once.  The stimulus is square.
%  synpos - [cmpt_idx cmpt_fraction]
%  timing - in ms
%  inst_vel - the instantaneous velocity of the stimulus
%
%  synmap - [cmpt cmpt_frac az el]
%  l_over_v - the l/v parameter
%  center - [el az] where to start
%  theta_start - size at start; zero would take infinite time. in degrees.
%  dt - temporal resolution, in ms
%
function [synpos, timing, inst_vel] = pattern_circular_looming_stimulus(synmap, l_over_v, center, theta_start, dt)
  % Preliminary
  timing = [];
  theta_final = 80; % Degrees - final size, as in screen

  recede = 0;
  if (l_over_v < 0) % recede
    recede = 1;
    l_over_v = -1*l_over_v;
    theta_start = -1*theta_start;
  end

  % Start, end time - compute it
  start_time = l_over_v/(tand(theta_start/2));
  end_time = l_over_v/(tand(theta_final/2));

  % GEnerate a center-centric coordinate set
  centered_el = synmap(:,4)-center(1);
  centered_az = synmap(:,3)-center(2);

  % determine total time ...
  total_time = start_time - end_time;
 
  % For each point, determine the angular deflection from the center 
  for s=1:length(centered_el)
    alpha(s) = acosd(cosd(centered_el(s))*cosd(centered_az(s)));
  end

  % instantaneous velocity vector
  inst_vel.vel_vec = zeros(length(0:dt:total_time),1);
	ivi = 1;
  
  % loop thru time, and add accordingly
  if (recede == 0)
		theta_0 = theta_start; % Theta at START of step
		for time=0:dt:total_time
			theta_f = 2*atand(l_over_v/(total_time + end_time - time));

			% Find any points between theta_0 and theta_f object size
			stimulated = find(alpha < theta_f & alpha > theta_0);
			
			% And append returned arrays if needbe
			if (length(stimulated) > 0)
				for s=1:length(stimulated)
					synpos(length(timing) + 1, :) = [synmap(stimulated(s),1) synmap(stimulated(s),2)];
					timing(length(timing) + 1) = time;
				end
			end
			inst_vel.vel_vec(ivi) = theta_f - theta_0; ivi = ivi + 1;
			theta_0 = theta_f;
		end  
  else
		theta_0 = theta_final;
		for time=total_time:-1*dt:0
			theta_f = 2*atand(l_over_v/(total_time + end_time - time));
			% Find any points between theta_0 and theta_f object size
			stimulated = find(alpha > theta_f & alpha < theta_0);
			% And append returned arrays if needbe
			if (length(stimulated) > 0)
				for s=1:length(stimulated)
					synpos(length(timing) + 1, :) = [synmap(stimulated(s),1) synmap(stimulated(s),2)];
					timing(length(timing) + 1) = total_time-time;
				end
			end
			inst_vel.vel_vec(ivi) = theta_f - theta_0; ivi = ivi + 1;
			theta_0 = theta_f;
		end 
  end

  % Scale instantaneous velocity by dt -- in s
  inst_vel.t_vec = 0:dt:total_time;
  inst_vel.vel_vec = abs([inst_vel.vel_vec/(.001*dt)])';

