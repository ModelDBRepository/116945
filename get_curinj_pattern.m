%
% S Peron 5/12/08
%
% This will return a current injection pattern - setting t_inj if you do assign it,
%  otherwise using I_inj_end and I_inj_start.
%
% pattern_type: 1 - fixed
%               2 - linear change (negative for DECREASE)
%               3 - loomlike (negative for RECEDE)
%
% params - if pattern_Type ==3, params(1) = start angle; params(2) = end angle ;  params(3) = l/v
% 
function ps_mod = get_curinj_pattern(ps_mod, pattern_type, I_base, I_max, params, t_inj)
	ps_mod.I_inj_nA = [0 0 0]; % Don't let other stuff mess with current vectors
  if (exist('t_inj') ~= 0)
	  ps_mod.I_inj_start= t_inj(1);
	  ps_mod.I_inj_end = t_inj(2);
  end

  ps_mod.I_of_t = zeros(3,length(0:ps_mod.dt:ps_mod.duration));
  ps_mod.I_of_t(2,:) = I_base*ones(1,length(0:ps_mod.dt:ps_mod.duration));
  switch (pattern_type)
    case 1 % fixed
      ps_mod.I_of_t(2,round(ps_mod.I_inj_start/ps_mod.dt):round(ps_mod.I_inj_end/ps_mod.dt)) = I_max;

    case -2 % linear
			inj_vec = linspace(I_max,I_base,length(round(ps_mod.I_inj_start/ps_mod.dt):round(ps_mod.I_inj_end/ps_mod.dt))); % linear
			ps_mod.I_of_t(2,round(ps_mod.I_inj_start/ps_mod.dt):round(ps_mod.I_inj_end/ps_mod.dt)) = inj_vec;

    case 2 % linear
			inj_vec = linspace(I_base,I_max,length(round(ps_mod.I_inj_start/ps_mod.dt):round(ps_mod.I_inj_end/ps_mod.dt))); % linear
			ps_mod.I_of_t(2,round(ps_mod.I_inj_start/ps_mod.dt):round(ps_mod.I_inj_end/ps_mod.dt)) = inj_vec;

    case -3 % loomlike -- RECEDE
      theta_0 = params(1);
      theta_f = params(2);
      l_over_v = params(3);
      t_f = -1*l_over_v/tand(theta_f);
      t_0 = -1*l_over_v/tand(theta_0);

      t_vec = t_f:-1*ps_mod.dt:t_0;
      inj_vec = -1*atand(l_over_v./t_vec);
      inj_vec = I_max * inj_vec/max(inj_vec); % scale to max
      inj_vec(find(inj_vec < I_base)) = I_base ; % scale to base
			ps_mod.I_of_t(2,round(ps_mod.I_inj_start/ps_mod.dt):round(ps_mod.I_inj_start/ps_mod.dt)+length(inj_vec)-1) = inj_vec';

    case 3 % loomlike
      theta_0 = params(1);
      theta_f = params(2);
      l_over_v = params(3);
      t_f = -1*l_over_v/tand(theta_f);
      t_0 = -1*l_over_v/tand(theta_0);

      t_vec = t_0:ps_mod.dt:t_f;
      inj_vec = -1*atand(l_over_v./t_vec);
      inj_vec = I_max * inj_vec/max(inj_vec); % scale to max
      inj_vec(find(inj_vec < I_base)) = I_base ; % scale to base
			ps_mod.I_of_t(2,round(ps_mod.I_inj_start/ps_mod.dt):round(ps_mod.I_inj_start/ps_mod.dt)+length(inj_vec)-1) = inj_vec';
  end

  % Morphological tweak
  ps_mod.I_of_t(2,:) = ps_mod.I_of_t(2,:)*.001/ps_mod.Ad; % convert to uA/cm2
