function makeFigures(varargin)
% makeFigures  All figures of the sweep, from the saved fits.
%
%   makeFigures()                                    % every figure
%   makeFigures('Figures', {'fits', 'heldout'})      % only some groups
%   makeFigures('ResultsDir', 'rerun/results', 'OutDir', 'rerun/figures')
%
% Nothing is re-fitted: every figure is drawn from results/<runTag>.mat
% (one per model and fold, written by runAllFolds / fitModel).
%
% FIGURE GROUPS, and where they go under OutDir (default figures/):
%   'data'         data/data_overview.png         raw data of experiments 4, 5, 6
%   'fits'         fit_train/<runTag>_pfit_train_exp<N>.png    fit on each training experiment
%                  fit_held_out/<runTag>_pfit_heldout_exp<N>.png  prediction of the held-out one
%   'sensitivity'  sensitivity/<runTag>_psens_exp<N>_{param,ic}_{mean,peak}.png
%                  local sensitivity heatmaps per training experiment. The
%                  sensitivity solve (~2-5 min per experiment) is only run for
%                  experiments without a saved results/sensitivity/<runTag>_exp<N>.mat
%                  (or for all with 'RecomputeSensitivity', true).
%   'correlation'  corr/<runTag>_pcorrelation.png  parameter correlations from the
%                  Fisher information of the kinetic parameters (needs the
%                  sensitivity files above)
%   'heldout'      heldout/heldout_rmse.png, heldout_rmse_shared.png,
%                  heldout_fractions.png, heldout_fraction_bars.png
%                  held-out scores of all fits (computeHeldOutMetrics, which
%                  also saves results/heldout_metrics.mat)
%   'reduced'      corr_reduced/<runTag>_reduced[_rgratio[_kdfixed]]_pcorrelation.png,
%                  modelD_forest_reduced.png/.pdf
%                  correlations and confidence intervals of reduced
%                  parameterisations (reduceFIM), models 0, 0B and D

    root = fileparts(mfilename('fullpath'));
    ip = inputParser;
    ip.addParameter('Figures', {'data','fits','sensitivity','correlation','heldout','reduced'}, ...
        @(x) iscellstr(x) || isstring(x));
    ip.addParameter('ResultsDir', fullfile(root, 'results'), @(x) ischar(x) || isstring(x));
    ip.addParameter('OutDir',     fullfile(root, 'figures'), @(x) ischar(x) || isstring(x));
    ip.addParameter('RecomputeSensitivity', false, @(x) islogical(x) && isscalar(x));
    ip.parse(varargin{:});
    opt        = ip.Results;
    groups     = cellstr(opt.Figures);
    resultsDir = char(opt.ResultsDir);
    outDir     = char(opt.OutDir);
    setup_paths();

    oldVis = get(groot, 'defaultFigureVisible');
    set(groot, 'defaultFigureVisible', 'off');
    restoreVis = onCleanup(@() set(groot, 'defaultFigureVisible', oldVis));

    % --- the runs: every model<Name>_train<..>_test<..>_<stamp>.mat ----------------
    files = dir(fullfile(resultsDir, 'model*_train*_test*_*.mat'));
    runTags = sort(erase({files.name}, '.mat'));
    fprintf('\n===== makeFigures: %d runs in %s =====\n', numel(runTags), resultsDir);

    %% Data overview
    if ismember('data', groups)
        prettyDataOverview('OutDir', fullfile(outDir, 'data'));
        close all;
    end

    %% Fits on the training and the held-out experiments
    if ismember('fits', groups)
        for i = 1:numel(runTags)
            run = loadRun(runTags{i}, resultsDir);
            ms  = getModelSpec(run.modelName);
            for k = 1:numel(run.trainExps)
                prettyFit(run.expStructs{k}, run.pFinal, ms.full, ...
                    'Anchor', run.Dcell{k}, 'RunTag', runTags{i}, ...
                    'OutDir', fullfile(outDir, 'fit_train'), ...
                    'FileSuffix', sprintf('train_exp%d', run.trainExps(k)), ...
                    'TMax', plotTimeLimit(run.trainExps(k)), ...
                    'ShowFractions', ms.usesFractions);
            end
            eTest = quiet(@() loadExperimentData(run.testExp));
            Dt    = quiet(@() prepareFitData(eTest, run.testExp, run.cfg));
            prettyFit(eTest, run.pFinal, ms.full, ...
                'Anchor', Dt, 'RunTag', runTags{i}, ...
                'OutDir', fullfile(outDir, 'fit_held_out'), ...
                'FileSuffix', sprintf('heldout_exp%d', run.testExp), ...
                'TMax', plotTimeLimit(run.testExp), ...
                'ShowFractions', ms.usesFractions);
            close all;
        end
    end

    %% Sensitivity heatmaps (and the Fisher information files they come with)
    if ismember('sensitivity', groups)
        for i = 1:numel(runTags)
            run = loadRun(runTags{i}, resultsDir);
            for k = 1:numel(run.trainExps)
                f = sensFile(resultsDir, runTags{i}, run.trainExps(k));
                if opt.RecomputeSensitivity || exist(f, 'file') ~= 2
                    SA = computeSensitivities(run, k, resultsDir);
                else
                    SA = load(f);
                end
                prettySensitivity(SA, runTags{i}, fullfile(outDir, 'sensitivity'));
                close all;
            end
        end
    end

    %% Parameter correlations from the kinetic Fisher information
    if ismember('correlation', groups)
        for i = 1:numel(runTags)
            run = loadRun(runTags{i}, resultsDir);
            F   = runFIM(resultsDir, run);
            prettyCorrelation(F, run.P29.name, runTags{i}, fullfile(outDir, 'corr'));
            close all;
        end
    end

    %% Held-out scores
    if ismember('heldout', groups)
        M = computeHeldOutMetrics(runTags, resultsDir);
        hDir = fullfile(outDir, 'heldout');
        prettyHeldOutRMSE(M, hDir);
        prettyHeldOutRMSEShared(M, hDir);
        prettyHeldOutFractions(M, hDir);
        prettyHeldOutFractionBars(M, hDir);
        close all;
    end

    %% Reduced parameterisations (models 0, 0B and D)
    if ismember('reduced', groups)
        rDir  = fullfile(outDir, 'corr_reduced');
        ratio = {'mu_rgmx', {'mu_rgman', 'mu_rgxyl'}};
        % {model, file suffix, reduceFIM options}
        variants = {
            '0',  'reduced',                 {}
            '0B', 'reduced',                 {}
            'D',  'reduced',                 {}
            'D',  'reduced_rgratio',         {'Tie', ratio}
            'D',  'reduced_rgratio_kdfixed', {'Tie', ratio, 'Fix', {'kd'}}
        };
        Dfolds = {};
        for i = 1:numel(runTags)
            run = loadRun(runTags{i}, resultsDir);
            F   = runFIM(resultsDir, run);
            for v = find(strcmp(variants(:,1), run.modelName))'
                [Fr, names] = reduceFIM(F, run.P29.name, variants{v,3}{:});
                prettyCorrelation(Fr, names, [runTags{i} '_' variants{v,2}], rDir, ...
                    'FontSize', 9, 'Annotate', true);
                close all;
            end
            if strcmp(run.modelName, 'D')
                Dfolds(end+1, :) = {run.testExp, F, run.P29.name}; %#ok<AGROW>
            end
        end
        % forest plot, folds in the order of the held-out experiment (4, 5, 6)
        [~, o] = sort(cell2mat(Dfolds(:,1)));
        prettyForestReduced(Dfolds(o,2), Dfolds{1,3}, rDir);
        close all;
    end
end


function f = sensFile(resultsDir, runTag, expNum)
    f = fullfile(resultsDir, 'sensitivity', sprintf('%s_exp%d.mat', runTag, expNum));
end


function F = runFIM(resultsDir, run)
% Kinetic Fisher information of a run: the sum over its training experiments.
    F = zeros(height(run.P29));
    for k = 1:numel(run.trainExps)
        f = sensFile(resultsDir, run.runTag, run.trainExps(k));
        assert(exist(f, 'file') == 2, 'makeFigures:noSensitivity', ...
            'No %s -- run makeFigures(''Figures'', {''sensitivity''}) first.', f);
        S = load(f, 'fim');
        assert(isequal(cellstr(S.fim.kinetic.paramNames(:)), cellstr(run.P29.name(:))), ...
            'makeFigures:fimMismatch', 'Fisher information parameters in %s do not match the run.', f);
        F = F + S.fim.kinetic.F;
    end
end


function out = quiet(fn) %#ok<INUSD> -- fn is called inside evalc
% Run fn() without its console output (data loading is chatty by design).
    out = [];
    evalc('out = fn();');
end
