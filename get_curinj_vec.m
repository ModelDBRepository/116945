%
% returns the curinj vector for a particular l/v
%
function [I_of_t max_time]= get_curinj_vec(l_over_v, I_max, I_base)
	theta_0 = 2;
	theta_f = 62;
	t_f = -1*l_over_v/tand(theta_f);
	t_0 = -1*l_over_v/tand(theta_0);

	t_vec = t_0:0.1:t_f;
	inj_vec = -1*atand(l_over_v./t_vec);
	inj_vec = I_max * inj_vec/max(inj_vec); % scale to max
	I_of_t = I_base*ones(1,20001);
	I_of_t(5000:5000+length(inj_vec)-1) = inj_vec';

	[irr max_idx] = max(I_of_t);
	max_time = max_idx/10;

