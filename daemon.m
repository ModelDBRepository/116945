%
% This should be called by a client process -- looks to the specified queue
%
function daemon(parpath)
  % --- some additions to the path -- include directory where your scripts are 
	% if the daemon resides elsewhere
  % path(path,'somedirectory');
   
  % --- paths
  lock_path = [parpath '/lock.m'];
  stop_path = [parpath '/stop_par.m'];

  disp(['Using path: ' parpath]);

  % If stop_par.m exists, this will terminate ...
  while(exist(stop_path, 'file') ~= 2)
		% lock the file -- OR WAIT
		while(exist(lock_path, 'file') == 2)
			disp('Lock exists -- waiting 1 s');
			pause(1);
		end
		fid = fopen(lock_path, 'w');
		fclose(fid);

		% Read the LAST .mat file - do this so that length of flist is the amount left
		%  in Q letting generate work
		flist = dir([parpath '/par_*.mat']);
    if (length(flist) > 0)
  		my_file = flist(length(flist)).name;
		
		  load([parpath '/' my_file]);
		
      movefile([parpath '/' my_file],[parpath '/' my_file '.tmp']);
    end

		% Unlock
		delete(lock_path);

    % Execute ...
    if (length(flist) > 0)
      exec_str = [function_name '('];
      if (length(function_params) > 0)
        exec_str = [exec_str 'function_params(1).value'];
      end
      for p=2:length(function_params)
        exec_str = [exec_str ',function_params(' num2str(p) ').value'];
      end
      exec_str = [exec_str ')'];
      
	  	disp(['Executing ' exec_str]);
      eval(exec_str);
      delete([parpath '/' my_file '.tmp']);
    else
      disp ('No file found; waiting 10 s');
      pause(10);
    end
  end
  disp('Exit signal');
	exit;
