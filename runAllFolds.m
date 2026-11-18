function summary = runAllFolds(varargin)
% runAllFolds  The model-comparison sweep: every model on every
% leave-one-experiment-out fold.
%
%   summary = runAllFolds()                       % the full sweep, as run on 2026-11-07
%   summary = runAllFolds('Models', {'A','D'})    % a subset of the models
%   summary = runAllFolds('OutDir', 'rerun')      % write results/ and logs/ under rerun/
%                                                 % instead of the project folder
%
% 6 models x 3 folds = 18 fits, run one after another by fitModel, all with
% the same seed. The folds hold out each of the three experiments in turn:
%
%   fold  train   test
%    1    5, 6     4
%    2    4, 6     5
%    3    4, 5     6
%
% Each fit writes its own results/<runTag>.mat, logs/<runTag>.txt and
% run_log.csv row (see fitModel). This function adds a summary table,
% results/sweep_summary_<yyyyMMdd_HHmm>.csv, with lsqnonlin's exit flag and
% how much the resnorm was still moving at the end of each fit. A failing
% fit is reported and skipped; the sweep carries on.
%
% The 2026-11-07 sweep took about 36 minutes on a local parallel pool.
% Afterwards run makeFigures to score the held-out experiments and draw
% the figures.

    setup_paths();

    ip = inputParser;
    ip.addParameter('Models', {'0', '0B', 'A', 'B', 'C', 'D'}, @(x) iscellstr(x) || isstring(x));
    ip.addParameter('Seed',   20260903, @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('OutDir', fileparts(mfilename('fullpath')), @(x) ischar(x) || isstring(x));
    ip.parse(varargin{:});
    models = cellstr(ip.Results.Models);
    seed   = ip.Results.Seed;
    outDir = char(ip.Results.OutDir);

    folds = [ struct('train', [5 6], 'test', 4)
              struct('train', [4 6], 'test', 5)
              struct('train', [4 5], 'test', 6) ];

    comment = sprintf('runAllFolds sweep %s: models %s, LOO folds over exp 4/5/6, seed %d', ...
        string(datetime('now', 'Format', 'yyyy-MM-dd')), strjoin(models, ','), seed);

    % Draw nothing on screen during the sweep.
    oldVis = get(groot, 'defaultFigureVisible');
    set(groot, 'defaultFigureVisible', 'off');
    restoreVis = onCleanup(@() set(groot, 'defaultFigureVisible', oldVis));

    nRuns = numel(models) * numel(folds);
    fprintf('\n===== runAllFolds: %d runs (%d models x %d folds), seed %d =====\n', ...
        nRuns, numel(models), numel(folds), seed);
    sweepStart = tic;

    rows = cell(nRuns, 1);
    n    = 0;
    for mi = 1:numel(models)
        for fi = 1:numel(folds)
            n = n + 1;
            fprintf('\n----- [%d/%d] model %s | train %s | test %d -----\n', ...
                n, nRuns, models{mi}, mat2str(folds(fi).train), folds(fi).test);

            runStart = tic;
            runTag = ""; AIC = NaN; BIC = NaN; errMsg = "";
            conv = struct('exitflag', NaN, 'iterations', NaN, 'relChange10', NaN, 'resnorm', NaN);
            try
                r = fitModel(models{mi}, folds(fi).train, folds(fi).test, seed, ...
                             'Comment', comment, 'OutDir', outDir);
                runTag = string(r.runTag);  AIC = r.AIC;  BIC = r.BIC;
                conv = struct('exitflag', r.sol.exitflag, 'iterations', r.sol.output.iterations, ...
                              'relChange10', r.sol.relResnormChange10, 'resnorm', r.sol.resnorm);
                status = "ok";
            catch ME
                diary off;
                status = "FAILED";
                errMsg = string(regexprep(ME.message, '\s+', ' '));
                fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
            end
            minutes = toc(runStart) / 60;
            fprintf('----- [%d/%d] %s in %.1f min -----\n', n, nRuns, status, minutes);

            rows{n} = table(string(models{mi}), string(mat2str(folds(fi).train)), ...
                folds(fi).test, status, runTag, AIC, BIC, ...
                conv.exitflag, conv.iterations, conv.relChange10, conv.resnorm, ...
                minutes, errMsg, ...
                'VariableNames', {'Model', 'TrainExps', 'TestExp', 'Status', ...
                                  'RunTag', 'AIC', 'BIC', 'ExitFlag', 'Iterations', ...
                                  'RelResnormChange10', 'Resnorm', 'Minutes', 'Error'});
            close all;
        end
    end

    summary = vertcat(rows{:});
    nFailed = sum(summary.Status == "FAILED");
    fprintf('\n===== runAllFolds finished: %d/%d ok, %.1f min total =====\n', ...
        nRuns - nFailed, nRuns, toc(sweepStart) / 60);
    disp(summary);

    summaryPath = fullfile(outDir, 'results', sprintf('sweep_summary_%s.csv', ...
        string(datetime('now', 'Format', 'yyyyMMdd_HHmm'))));
    writetable(summary, summaryPath);
    fprintf('Summary written to %s\n', summaryPath);
end
