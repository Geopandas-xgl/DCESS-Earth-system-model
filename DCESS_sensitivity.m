%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%    U isotope sensitivity        %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function MassBalance_dZn_sens
clear all

k_A = [0.5,0.9,1,1.1,1.5];

parsSens.sensruns = length(k_A);

pool = gcp("nocreate");

if isempty(pool)
    parpool
end

out = cell(1,parsSens.sensruns);
successID = false(1,parsSens.sensruns);
errorInfo = cell(1,parsSens.sensruns);
workerID = zeros(1,parsSens.sensruns);

tStart = tic;

parfor N = 1: parsSens.sensruns
    task = getCurrentTask();
    if isempty(task)
        workerID(N) = 0;
    else
        workerID(N) = task.ID;
    end

    % fprintf("Run # %d of %d \n", N, pars.sensruns) ;
    try
        % Call your initialisation + simulation function
        out{N} = Thilda_M(k_A(N));
        successID(N) = true;   % mark as success
    catch ME
        % If error occurs, record the exception, leave result empty / NaN / whatever
        errorInfo{N} = ME;
        successID(N) = false;  % mark as failed
        fprintf('Error in run %d: %s\n', N, ME.message);
        % Optionally continue to next iteration — but in parfor, that's implicit after catch
    end
end
parsSens.elapsed = toc(tStart);

parsSens.k_A = k_A;
disp(unique(workerID))
out = out(successID);
run = [out{:}];

parsSens.dur = seconds(parsSens.elapsed);
parsSens.dur.Format = "hh:mm:ss";
fprintf("Done: time = %s\n", string(parsSens.dur));
%%
if ~exist("Data","dir")
    mkdir("Data")
end

save("Data\senanal_results.mat","run","parsSens",'-mat') ;
%% 
DCESS_plot
