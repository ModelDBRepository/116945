%
%  This generates the simulation results for nature neuroscience figs
%
function plot_nn_sims()
	% directory 
	rootdir = 'nn_par_out/';
  
	% definitions
  gauss = 0.1*normpdf(-100:.1:100,0,20); % our gaussian for convolutions
  ctrl_color = [0 0 0];
  ctrl_patch = [.8 .8 .8];
  bapta_color = [1 0 0];
  bapta_patch = [.7 .5 .5];

  disp (' *************** STARTING ANALYSIS **************** ');
	disp([' *** DIRECTORY: ' rootdir ' ***']);

  % prelimiaries -- process data if needbe
	fname_roots = {'realistic_nobapta_transbox', ...
	          'realistic_bapta_transbox', ...
						'realistic_nobapta_loom_lv_10', ...
						'realistic_bapta_loom_lv_10', ...
						'realistic_nobapta_loom_lv_30', ...
						'realistic_bapta_loom_lv_30', ...
						'realistic_nobapta_loom_lv_50', ...
						'realistic_bapta_loom_lv_50'};

	% class 0 - transl ; 1 - loom
	class = [0 0 1 1 1 1 1 1]; 

  % -------------------------------
	% data processing
  for f=1:length(fname_roots)
	  if (exist([rootdir fname_roots{f} '_avg.mat']) == 0)
			num_trials = length(dir([rootdir fname_roots{f} '_*.mat']));
			gauss_means = zeros(num_trials,25001);
			tpeak = zeros(num_trials,1);
			fss = zeros(num_trials,1);
			fmax = zeros(num_trials,1);

			for n=1:num_trials % assume 10 for each
				load([rootdir fname_roots{f} '_' num2str(n) '.mat']);
				spike_idx = get_spikes(-40,y(:,1));
				inst_freq = 1000*get_inst_freq(t, spike_idx);
				gauss_mean_ifc = conv(inst_freq, gauss);
				gauss_means(n,:) = gauss_mean_ifc((length(gauss)-1)/2:length(t) + (length(gauss)-1)/2 -1);
				if (class(f) == 0)
					[fss(n) fssse fmax(n) fmaxse] = trans_resp_props(gauss_means(n,:), gauss_means(n,:));
				else
					[fss(n) fmax(n) tpeak(n)] = loom_resp_props(gauss_means(n,:), gauss_means(n,:));
				end
			end

			gauss_mean = mean(gauss_means,1);
			gauss_sd = std(gauss_means,1);
			gauss_se = gauss_sd/sqrt(num_trials);
			t_peak.mu = mean(tpeak);
			t_peak.sd = std(tpeak);
			t_peak.se = std(tpeak)/sqrt(num_trials);
			f_max.mu = mean(fmax);
			f_max.sd = std(fmax);
			f_max.se = std(fmax)/sqrt(num_trials);
			f_ss.mu = mean(fss);
			f_ss.sd = std(fss);
			f_ss.se = std(fss)/sqrt(num_trials);

			% write to file ...
			t_peak_raw = tpeak;
			save([rootdir fname_roots{f} '_avg.mat'], 'gauss_mean', 'gauss_sd', 'gauss_se', 'num_trials', 't', 'f_ss', 'f_max', 't_peak', 't_peak_raw', 'fmax', 'fss');
			disp(['saved: par_out/' fname_roots{f} '_avg.mat']);
		else
		  disp('No processing - file found');
		end
	end

  % -------------------------------
	% Plot 1: Response to translation
	figure; 
  tvec = 0:1:2500;
	stimsize = [linspace(0, 10, 250) linspace(10,10,2001) linspace(10,0,250)];
		
	% No bapta then bapta
	subplot('position', [.1 .15 .75 .2]);
	hold on;
	load([rootdir fname_roots{1} '_avg.mat']);
  [fss_ctrl fssse_ctrl fmax_ctrl fmaxse_ctrl] = trans_resp_props(gauss_mean, gauss_se);
	fss_ctrl_raw = fss;
	fmax_ctrl_raw = fmax;
  plot_err_poly (t+250, gauss_mean, gauss_se, ctrl_color, ctrl_patch);
	disp(['ctrl fss(se): ' num2str(f_ss.mu) '(' num2str(f_ss.se) ') ctrl fmax(se): ' num2str(f_max.mu) '(' num2str(f_max.se) ')' ]);
	load([rootdir fname_roots{2} '_avg.mat']);
  fss_drug_raw = fss;
  fmax_drug_raw = fmax;
  [fss_bapta fssse_bapta fmax_bapta fmaxse_bapta] = trans_resp_props(gauss_mean, gauss_se);
  plot_err_poly (t+250, gauss_mean, gauss_se, bapta_color, bapta_patch);
	disp(['bapta fss(se): ' num2str(f_ss.mu) '(' num2str(f_ss.se) ') bapta fmax(se): ' num2str(f_max.mu) '(' num2str(f_max.se) ')' ]);
	axis([0 2500 0 100]);
	set (gca,'TickDir', 'out');

	% fmax, fss
  disp(['ctrl fss of AVG: ' num2str(fss_ctrl) ' se: ' num2str(fssse_ctrl) ' fmax: ' num2str(fmax_ctrl) ' se: ' num2str(fmaxse_ctrl)]);
  disp(['bapta fss of AVG: ' num2str(fss_bapta) ' se: ' num2str(fssse_bapta) ' fmax: ' num2str(fmax_bapta) ' se: ' num2str(fmaxse_bapta)]);
  draw_bar_group (40, 0, 10, 40, [ctrl_color; bapta_color], [fmax_ctrl fmax_bapta], 1, [fmaxse_ctrl fmaxse_bapta])
  draw_bar_group (150, 0, 10, 40, [ctrl_color; bapta_color], [fss_ctrl fss_bapta], 1, [fssse_ctrl fssse_bapta])
	disp(['RS effect of simulated bapta on fmax: ' num2str(ranksum(fmax_ctrl_raw, fmax_drug_raw))]);
	disp(['RS effect of simulated bapta on fss: ' num2str(ranksum(fss_ctrl_raw, fss_drug_raw))]);
	disp(['% chg fmax: ' num2str(100*(mean(fmax_drug_raw)-mean(fmax_ctrl_raw))/mean(fmax_ctrl_raw))]);
	disp(['% chg fss: ' num2str(100*(mean(fss_drug_raw)-mean(fss_ctrl_raw))/mean(fss_ctrl_raw))]);

  % Stimulus
	subplot('position', [.1 .1 .75 .05]);
  plot(tvec, stimsize, 'Color', [.5 .5 .5], 'LineWidth', 2);
	axis([0 2500 -1 11]);
	set (gca,'TickDir', 'out');


  % -------------------------------
	% Plot 2: Response to looming
	figure;

	% l/v = 10 No bapta then bapta
	disp('l/v 10 loom');
  l_over_v = 10; 
	tvec = (500:-0.1:0) ; 
	stimsize = [2*atand(l_over_v./tvec) ones(1,1000)*(180)]; 
	idx_80 = min(find(stimsize > 80));
	t_80 = tvec(idx_80);
	stimsize(find(stimsize > 80)) = 80;
	tvec = -500:0.1:100;
		
	subplot('position', [.1 .75 .5 .2]);
	hold on;
	load([rootdir fname_roots{3} '_avg.mat']);
  fss_ctrl_raw = fss;
  fmax_ctrl_raw = fmax;
  tpeak_ctrl_raw(1,:) = t_peak_raw-2400+t_80;
	disp(['ctrl fss(se): ' num2str(f_ss.mu) '(' num2str(f_ss.se) ') ctrl fmax(se): ' num2str(f_max.mu) '(' num2str(f_max.se) ')' ]);
  plot_err_poly (t-2000+t_80, gauss_mean, gauss_se, ctrl_color, ctrl_patch);
  [fss_ctrl fmax_ctrl tpeak_ctrl fssse_ctrl fmaxse_ctrl] = loom_resp_props(gauss_mean, gauss_se);
  disp(['ctrl fss of AVG: ' num2str(fss_ctrl) ' se: ' num2str(fssse_ctrl) ' fmax: ' num2str(fmax_ctrl) ' se: ' num2str(fmaxse_ctrl)]);

	load([rootdir fname_roots{4} '_avg.mat']);
  fss_drug_raw = fss;
  fmax_drug_raw = fmax;
	tpeak_bapta_raw(1,:) = t_peak_raw-2400+t_80;
	disp(['bapta fss(se): ' num2str(f_ss.mu) '(' num2str(f_ss.se) ') bapta fmax(se): ' num2str(f_max.mu) '(' num2str(f_max.se) ')' ]);
  plot_err_poly (t-2000+t_80, gauss_mean, gauss_se, bapta_color, bapta_patch);
  [fss_bapta fmax_bapta tpeak_bapta fssse_bapta fmaxse_bapta] = loom_resp_props(gauss_mean, gauss_se);
  disp(['bapta fss of AVG: ' num2str(fss_bapta) ' se: ' num2str(fssse_bapta) ' fmax: ' num2str(fmax_bapta) ' se: ' num2str(fmaxse_bapta)]);
	axis([-100 500 0 400]);
	set (gca,'TickDir', 'out');

  draw_bar_group (10, 0, 5, 20, [ctrl_color; bapta_color], [fmax_ctrl fmax_bapta], 1, [fmaxse_ctrl fmaxse_bapta])
  draw_bar_group (60, 0, 5, 20, [ctrl_color; bapta_color], [fss_ctrl fss_bapta], 1, [fssse_ctrl fssse_bapta])

  % Stimulus
	subplot('position', [.1 .7 .5 .05]);
  plot(tvec, stimsize, 'Color', [.5 .5 .5], 'LineWidth', 2);
	axis([-500 100 -1 91]);
	set (gca,'TickDir', 'out');
	disp(['RS effect of simulated bapta on fmax: ' num2str(ranksum(fmax_ctrl_raw, fmax_drug_raw))]);
	disp(['RS effect of simulated bapta on fss: ' num2str(ranksum(fss_ctrl_raw, fss_drug_raw))]);
	disp(['% chg fmax: ' num2str(100*(mean(fmax_drug_raw)-mean(fmax_ctrl_raw))/mean(fmax_ctrl_raw))]);
	disp(['% chg fss: ' num2str(100*(mean(fss_drug_raw)-mean(fss_ctrl_raw))/mean(fss_ctrl_raw))]);


	% l/v = 30 No bapta then bapta
	disp('l/v 30 loom');
  l_over_v = 30; 
	tvec = (500:-0.1:0) ; 
	stimsize = [2*atand(l_over_v./tvec) ones(1,1000)*(180)]; 
	idx_80 = min(find(stimsize > 80));
	t_80 = tvec(idx_80);
	stimsize(find(stimsize > 80)) = 80;
	tvec = -500:0.1:100;
		
	subplot('position', [.1 .45 .5 .2]);
	hold on;
	load([rootdir fname_roots{5} '_avg.mat']);
  fss_ctrl_raw = fss;
  fmax_ctrl_raw = fmax;
	tpeak_ctrl_raw(2,:) = t_peak_raw-2400+t_80;
	disp(['ctrl fss(se): ' num2str(f_ss.mu) '(' num2str(f_ss.se) ') ctrl fmax(se): ' num2str(f_max.mu) '(' num2str(f_max.se) ')' ]);
  plot_err_poly (t-2000+t_80, gauss_mean, gauss_se, ctrl_color, ctrl_patch);
  [fss_ctrl fmax_ctrl tpeak_ctrl fssse_ctrl fmaxse_ctrl] = loom_resp_props(gauss_mean, gauss_se);
  disp(['ctrl fss of AVG: ' num2str(fss_ctrl) ' se: ' num2str(fssse_ctrl) ' fmax: ' num2str(fmax_ctrl) ' se: ' num2str(fmaxse_ctrl)]);

	load([rootdir fname_roots{6} '_avg.mat']);
  fss_drug_raw = fss;
  fmax_drug_raw = fmax;
	tpeak_bapta_raw(2,:) = t_peak_raw-2400+t_80;
	disp(['bapta fss(se): ' num2str(f_ss.mu) '(' num2str(f_ss.se) ') bapta fmax(se): ' num2str(f_max.mu) '(' num2str(f_max.se) ')' ]);
  plot_err_poly (t-2000+t_80, gauss_mean, gauss_se, bapta_color, bapta_patch);
  [fss_bapta fmax_bapta tpeak_bapta fssse_bapta fmaxse_bapta] = loom_resp_props(gauss_mean, gauss_se);
  disp(['bapta fss of AVG: ' num2str(fss_bapta) ' se: ' num2str(fssse_bapta) ' fmax: ' num2str(fmax_bapta) ' se: ' num2str(fmaxse_bapta)]);
	axis([-100 500 0 400]);
	set (gca,'TickDir', 'out');

  draw_bar_group (10, 0, 5, 20, [ctrl_color; bapta_color], [fmax_ctrl fmax_bapta], 1, [fmaxse_ctrl fmaxse_bapta])
  draw_bar_group (60, 0, 5, 20, [ctrl_color; bapta_color], [fss_ctrl fss_bapta], 1, [fssse_ctrl fssse_bapta])

  % Stimulus
	subplot('position', [.1 .4 .5 .05]);
  plot(tvec, stimsize, 'Color', [.5 .5 .5], 'LineWidth', 2);
	axis([-500 100 -1 91]);
	set (gca,'TickDir', 'out');
	disp(['RS effect of simulated bapta on fmax: ' num2str(ranksum(fmax_ctrl_raw, fmax_drug_raw))]);
	disp(['RS effect of simulated bapta on fss: ' num2str(ranksum(fss_ctrl_raw, fss_drug_raw))]);
	disp(['% chg fmax: ' num2str(100*(mean(fmax_drug_raw)-mean(fmax_ctrl_raw))/mean(fmax_ctrl_raw))]);
	disp(['% chg fss: ' num2str(100*(mean(fss_drug_raw)-mean(fss_ctrl_raw))/mean(fss_ctrl_raw))]);

 

	% l/v = 50 No bapta then bapta
	disp('l/v 50 loom');
  l_over_v = 50; 
	tvec = (500:-0.1:0) ; 
	stimsize = [2*atand(l_over_v./tvec) ones(1,1000)*(180)]; 
	idx_80 = min(find(stimsize > 80));
	t_80 = tvec(idx_80);
	stimsize(find(stimsize > 80)) = 80;
	tvec = -500:0.1:100;
		
	subplot('position', [.1 .15 .5 .2]);
	hold on;
	load([rootdir fname_roots{7} '_avg.mat']);
  fss_ctrl_raw = fss;
  fmax_ctrl_raw = fmax;
	tpeak_ctrl_raw(3,:) = t_peak_raw-2400+t_80;
	disp(['ctrl fss(se): ' num2str(f_ss.mu) '(' num2str(f_ss.se) ') ctrl fmax(se): ' num2str(f_max.mu) '(' num2str(f_max.se) ')' ]);
  plot_err_poly (t-2000+t_80, gauss_mean, gauss_se, ctrl_color, ctrl_patch);
  [fss_ctrl fmax_ctrl tpeak_ctrl fssse_ctrl fmaxse_ctrl] = loom_resp_props(gauss_mean, gauss_se);
  disp(['ctrl fss of AVG: ' num2str(fss_ctrl) ' se: ' num2str(fssse_ctrl) ' fmax: ' num2str(fmax_ctrl) ' se: ' num2str(fmaxse_ctrl)]);

  load([rootdir fname_roots{8} '_avg.mat']);
  fss_drug_raw = fss;
  fmax_drug_raw = fmax;
	tpeak_bapta_raw(3,:) = t_peak_raw-2400+t_80;
	disp(['bapta fss(se): ' num2str(f_ss.mu) '(' num2str(f_ss.se) ') bapta fmax(se): ' num2str(f_max.mu) '(' num2str(f_max.se) ')' ]);
  plot_err_poly (t-2000+t_80, gauss_mean, gauss_se, bapta_color, bapta_patch);
  [fss_bapta fmax_bapta tpeak_bapta fssse_bapta fmaxse_bapta] = loom_resp_props(gauss_mean, gauss_se);
  disp(['bapta fss of AVG: ' num2str(fss_bapta) ' se: ' num2str(fssse_bapta) ' fmax: ' num2str(fmax_bapta) ' se: ' num2str(fmaxse_bapta)]);
	axis([-100 500 0 400]);
	set (gca,'TickDir', 'out');

  draw_bar_group (10, 0, 5, 20, [ctrl_color; bapta_color], [fmax_ctrl fmax_bapta], 1, [fmaxse_ctrl fmaxse_bapta])
  draw_bar_group (60, 0, 5, 20, [ctrl_color; bapta_color], [fss_ctrl fss_bapta], 1, [fssse_ctrl fssse_bapta])

  % Stimulus
	subplot('position', [.1 .1 .5 .05]);
  plot(tvec, stimsize, 'Color', [.5 .5 .5], 'LineWidth', 2);
	axis([-500 100 -1 91]);
	set (gca,'TickDir', 'out');
	disp(['RS effect of simulated bapta on fmax: ' num2str(ranksum(fmax_ctrl_raw, fmax_drug_raw))]);
	disp(['RS effect of simulated bapta on fss: ' num2str(ranksum(fss_ctrl_raw, fss_drug_raw))]);
	disp(['% chg fmax: ' num2str(100*(mean(fmax_drug_raw)-mean(fmax_ctrl_raw))/mean(fmax_ctrl_raw))]);
	disp(['% chg fss: ' num2str(100*(mean(fss_drug_raw)-mean(fss_ctrl_raw))/mean(fss_ctrl_raw))]);


  % ---- l/v v. ttc
  tpeak_ctrl_raw = abs(tpeak_ctrl_raw);
  tpeak_bapta_raw = abs(tpeak_bapta_raw);

	l_over_vs = [10 30 50];
  mean_tpeak_ctrl = mean(tpeak_ctrl_raw,2);
  se_tpeak_ctrl = std(tpeak_ctrl_raw')/sqrt(length(tpeak_ctrl_raw));
  mean_tpeak_bapta = mean(tpeak_bapta_raw,2);
  se_tpeak_bapta = std(tpeak_bapta_raw')/sqrt(length(tpeak_bapta_raw));
  

  subplot('position', [.7 .5 .25 .25]);
	hold on;
	plot(l_over_vs, mean_tpeak_ctrl, '^', 'MarkerEdgeColor', ctrl_color, 'MarkerFaceColor', ctrl_color);  
	plot(l_over_vs+1, mean_tpeak_bapta, '^', 'MarkerEdgeColor', bapta_color, 'MarkerFaceColor', bapta_color);  
	for l=1:length(l_over_vs)
	  plot([l_over_vs(l) l_over_vs(l)], [se_tpeak_ctrl(l) -1*se_tpeak_ctrl(l)]+mean_tpeak_ctrl(l), 'Color', ctrl_color);
	  plot([l_over_vs(l) l_over_vs(l)]+1, mean_tpeak_bapta(l)+[se_tpeak_bapta(l) -1*se_tpeak_bapta(l)], 'Color', bapta_color);
	  rs = ranksum(tpeak_ctrl_raw(l,:), tpeak_bapta_raw(l,:));
		ntp = length(tpeak_ctrl_raw(l,:));
		disp([num2str(l_over_vs(l)) ' mean tpeak pre (se): ' num2str(mean_tpeak_ctrl(l)) '(' num2str(se_tpeak_ctrl(l)) ') post (se): ' ...
		      num2str(mean_tpeak_bapta(l)) '(' num2str(se_tpeak_bapta(l)) ') p=' num2str(rs) ' n=' num2str(ntp)]);
	end
	axis ([5 55 0 200]);
	set (gca,'TickDir', 'out');
  
  % Correlation coefficients please
	corr_ctrl = corr(l_over_vs', mean_tpeak_ctrl);
	disp(['rho_ctrl: ' num2str(corr_ctrl)]);
	corr_bapta = corr(l_over_vs', mean_tpeak_bapta);
	disp(['rho_bapta: ' num2str(corr_bapta)]);

  % Fit line and plot to the poitns
  fit_pre = polyfit(l_over_vs, mean_tpeak_ctrl', 1);
  fit_post = polyfit(l_over_vs, mean_tpeak_bapta', 1);
  
  plot([0 60], fit_pre(1)*[0 60] + fit_pre(2), ':', 'Color', ctrl_color);
  plot([0 60], fit_post(1)*[0 60] + fit_post(2), ':', 'Color', bapta_color);

  % labels
  ylabel('t_c_o_l_l_i_s_i_o_n - t_p_e_a_k (ms)');
  xlabel('l/v (ms)');
  
  % ---  theta_thresh (threshold angle) v. delta (intercept)
  subplot('position', [.7 .05 .25 .25]);
  axis([20 60 0 30]);
  hold on;

  for d=1:length(tpeak_ctrl_raw)
    ind_fit_pre(d,:) = polyfit(l_over_vs, tpeak_ctrl_raw(:,d)', 1);
    ind_fit_post(d,:) = polyfit(l_over_vs, tpeak_bapta_raw(:,d)', 1);
  end
  
  % Compute individual ANIMAL (eventually) deltas and theta-thresholds
  ind_delta_pre = abs(ind_fit_pre(:,2));
  ind_delta_post = abs(ind_fit_post(:,2));

  disp(['delta pre: ' num2str(mean(ind_delta_pre)) ' se: ' num2str(std(ind_delta_pre)/sqrt(length(ind_delta_pre)))]);
  disp(['delta post: ' num2str(mean(ind_delta_post)) ' se: ' num2str(std(ind_delta_post)/sqrt(length(ind_delta_post)))]);
  
  ind_theta_thresh_pre = 2*atand(1./ind_fit_pre(:,1));
  ind_theta_thresh_post = 2*atand(1./ind_fit_post(:,1));

  disp(['tht thresh pre: ' num2str(mean(ind_theta_thresh_pre)) ' se: ' num2str(std(ind_theta_thresh_pre)/sqrt(length(ind_theta_thresh_pre)))]);
  disp(['tht thresh post: ' num2str(mean(ind_theta_thresh_post)) ' se: ' num2str(std(ind_theta_thresh_post)/sqrt(length(ind_theta_thresh_post)))]);

  % Do a statistical test on the significance of the difference in deltas
  p_val = kruskalwallis([ind_delta_pre;ind_delta_post]', [], 'off');
  disp(' ');
  disp(['p-value for hypothesis that there is a drug effect on deltas (kruskallwallis): ' num2str(p_val)]);
  p_val = ranksum(ind_delta_pre, ind_delta_post);
  disp(['p-value for hypothesis that there is a drug effect on deltas (ranksum): ' num2str(p_val)]);

  % Do a statistical test on the significance of the difference in theta
  % thresh
  p_val = kruskalwallis([ind_theta_thresh_pre;ind_theta_thresh_post]', [], 'off');
  disp(' ');
  disp(['p-value for hypothesis that there is a drug effect on theta thresh (kruskallwallis): ' num2str(p_val)]);
  p_val = ranksum(ind_theta_thresh_pre, ind_theta_thresh_post);
  disp(['p-value for hypothesis that there is a drug effect on theta thresh (ranksum): ' num2str(p_val)]);
  
  sf = sqrt(length(ind_theta_thresh_pre));

  plot(mean(ind_delta_pre), mean(ind_theta_thresh_pre), 'o', 'Color', ctrl_color, 'MarkerFaceColor', ctrl_color);
  plot([mean(ind_delta_pre) mean(ind_delta_pre)], ...
    [mean(ind_theta_thresh_pre)-std(ind_theta_thresh_pre)/sf mean(ind_theta_thresh_pre)+std(ind_theta_thresh_pre)/sf], ...
    '-', 'Color', ctrl_color);
  plot([mean(ind_delta_pre)-std(ind_delta_pre)/sf mean(ind_delta_pre)+std(ind_delta_pre)/sf], ...
    [mean(ind_theta_thresh_pre) mean(ind_theta_thresh_pre)], ...
    '-', 'Color', ctrl_color);
  plot(mean(ind_delta_post), mean(ind_theta_thresh_post), 'o', 'Color', bapta_color, 'MarkerFaceColor', bapta_color);
  plot([mean(ind_delta_post) mean(ind_delta_post)], ...
    [mean(ind_theta_thresh_post)-std(ind_theta_thresh_post)/sf mean(ind_theta_thresh_post)+std(ind_theta_thresh_post)/sf], ...
    '-', 'Color', bapta_color);
  plot([mean(ind_delta_post)-std(ind_delta_post)/sf mean(ind_delta_post)+std(ind_delta_post)/sf], ...
    [mean(ind_theta_thresh_post) mean(ind_theta_thresh_post)], ...
    '-', 'Color', bapta_color);

  % labels
  ylabel('theta_t_h_r_e_s_h (deg)');
  xlabel('delta (ms)');
  

% 
% Computes f_ss and f_max (for translation) as per figure for data in manuscript
%
function [fss fssse fmax fmaxse] = trans_resp_props(gauss_mean, gauss_se)
  min_freq = 2; % Considered minimal frequency

  [fmax i_max] = max(gauss_mean(1:2500));

  fmaxse = gauss_se(i_max);

  % Find the 'dip'
  d_freq = diff(gauss_mean);
  candidates = find(d_freq > 0);
	i_start = candidates(min(find(candidates > i_max + 2)));
	gauss_mean(i_start);
  i_end = max(find(gauss_mean > min_freq));
	fss = mean(gauss_mean(i_start:i_end));

	fssse = mean(gauss_se(i_start:i_end));

  % Plot the points where the fmean is cutoff?
	if ( 0 == 1)
		plot([i_start/10 i_start/10]+250, [0 100], 'k:');
		plot([i_end/10 i_end/10]+250, [0 100], 'k:');
  end


% 
% Computes t_peak, f_ss and f_max (for loom) as per figure for data in manuscript
% In this case, f_ss is over the last 500 ms of approach (i.e., 500 to 0 in ttc)
%
function [fss fmax tpeak fssse fmaxse] = loom_resp_props(gauss_mean, gauss_se)
  min_freq = 2; % Considered minimal frequency

  % fmax and tpeak
  [fmax i_max] = max(gauss_mean);
	tpeak = i_max/10;

  % fss
	i_start = 19000;
	i_end = 24000;
	fss = mean(gauss_mean(i_start:i_end));
	fssse = mean(gauss_se(i_start:i_end));
	fmaxse = gauss_se(i_max);




%
% For plotting a nice overlay with SDs etc. - pre AND post data should be given
% if possible
%
function plot_err_poly (time, v_mean, v_sd, color, patch_color)
  % First plot the SDs
  [x_err_poly, y_err_poly] = get_sem_poly(time, v_mean, v_sd);
  patch(x_err_poly,y_err_poly, patch_color, 'EdgeColor', 'none');
 
  hold on;
  % And finally the main line
  plot(time, v_mean, 'Color', color, 'Linewidth', 2);
  
%
% Returns your data as a polygon for bounding y at each x value with +/- y_off
%  for real nice SEM/SD plots
%
function [ret_x, ret_y] = get_sem_poly(x, y, y_off);
  ret_x = zeros(2*length(x),1);
  ret_y = zeros(2*length(x),1);
  
  l = length(x);
  
  for i=1:length(ret_x)
    if (i < l)
      ret_x(i) = x(i);
      ret_y(i) = y(i) + y_off(i);
    elseif (i == l || i == l+1)
      ret_x(i) = x(l)+2;
      ret_y(i) = y(l) + y_off(l);
      if ( i == l + 1) ; ret_y(i) = y(l) - y_off(l); end
    else
      ret_x(i) = x(2*l - i + 1);
      ret_y(i) = y(2*l - i + 1) - y_off(2*l - i + 1);
    end
  end

%
% Draws some wicked bars
%   x,y_offsets are the bottom-left coordinates of the group
%   bar_spacing and width are the x(or y)-unit based spaces between bars and widths
%     thereof
%   bar_colors is a n_bars x 3 matrix corresponding to bar_data, which stores the y-values
%      or x-values, depending on orientation
%   orientation: 1 - vertical bars; 2 - horizontal bars
%   line_data: sd values
%
function draw_bar_group (x_offset, y_offset, bar_spacing, bar_size, ...
                         bar_colors, bar_data, orientation, line_data)
  x_base = x_offset;
  y_base = y_offset;


  % vertical
  if (orientation == 1)
    for b=1:length(bar_data)
      % The bar
      x = [x_base x_base x_base+bar_size x_base+bar_size x_base];
      y = [y_base y_base+bar_data(b) y_base+bar_data(b) y_base y_base];

      patch(x,y,bar_colors(b,:), 'EdgeColor', 'none');
      % error/whatever lines
      if (line_data(b) > 0)
        plot([x_base+(bar_size/2) x_base+(bar_size/2)], ...
             [y_base+bar_data(b) y_base+bar_data(b)+line_data(b)], ...
             'Color', bar_colors(b,:));
      end
    
      % increment x_base
      x_base = x_base + bar_size + bar_spacing;
    end
  % horizontal orientation
  else
    for b=1:length(bar_data)
      x = [x_base x_base x_base+bar_data(b) x_base+bar_data(b) x_base];
      y = [y_base y_base+bar_size y_base+bar_size y_base y_base];

      patch(x,y,bar_colors(b,:), 'EdgeColor', 'none');

      % error/whatever lines
      if (line_data(b) > 0)
        plot([x_base+bar_data(b) x_base+bar_data(b)+line_data(b)], ...
             [y_base+(bar_size/2) y_base+(bar_size/2)], ...
             'Color', bar_colors(b,:));
      end

      % increment x_base
      y_base = y_base + bar_size + bar_spacing;
    end
  end

  
