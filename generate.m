%
% This should be called within a for-loop to add its ONE LINE execution to the
%  parallel system queue, so that listeners can pick it up.
%
function generate(parpath, function_name, function_params)
  lock_path = [parpath '/lock.m'];

  % lock the file -- OR WAIT
  while(exist(lock_path, 'file') == 2)
    disp('Lock exists -- waiting 1 s');
    pause(1);
  end
  fid = fopen(lock_path, 'w');
  fclose(fid);

  % Generate the entry .mat file
  flist = dir([parpath '/par_*.mat']);
  next_idx = length(flist) + 1;
  newfname = [parpath '/par_' sprintf('%.6d', next_idx) '.mat'];
  save(newfname, 'function_name', 'function_params');
	disp(['Generating: ' newfname]);

  % Unlock
  delete(lock_path);
