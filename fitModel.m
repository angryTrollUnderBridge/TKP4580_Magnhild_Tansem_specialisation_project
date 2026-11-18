function out = fitModel(modelName, trainExps, testExp, seed, varargin)
% fitModel  Fit one model to the training experiments of one fold.
%
%   out = fitModel('D', [5 6], 4, 20260903)
%   out = fitModel(..., 'Comment', 'text for run_log.csv', 'OutDir', 'myRuns')
%
% modelName  '0', '0B', 'A', 'B', 'C' or 'D' (see getParamConfig / src/models)
% trainExps  the two training experiments, e.g. [5 6]
% testExp    the held-out experiment, e.g. 4 (only recorded here; scored by
%            makeFigures / computeHeldOutMetrics)
% seed       random seed of the genetic algorithm
%
% Writes, under OutDir (default: the project folder):
%   results/<runTag>.mat   everything needed to reproduce figures and scores
%   results/run_log.csv    one appended row (AIC, BIC, convergence, ...)
%   logs/<runTag>.txt      the console output of the run
% with runTag = model<Name>_train<Exps>_test<Exp>_<yyyyMMdd_HHmm>.
%
% THE FIT, IN TWO STAGES
%   Stage 1  kinetic parameters (everything except the CO2 yields), fitted to
%            OD, FCM strain fractions and HPLC sugars of both training
%            experiments at once, in log10 parameter space:
%              1a geneticSearch      global search (GA)
%              1b paramLeastSquares  local refinement (lsqnonlin), started
%                                    from the GA result
%            The model integrated here (modelSpec.kinetic) has no CO2 state
%            and no air-flow input: CO2 does not feed back into biomass or
%            sugars, and leaving it out makes each solve far cheaper.
%   Stage 2  CO2 yields, with the kinetics fixed: the CO2 balance is then
%            linear in 1/Yco2, so it is one bounded linear least-squares
%            problem (estimateCO2Yields).

    ip = inputParser;
    ip.addParameter('Comment', '', @(x) ischar(x) || isstring(x));
    ip.addParameter('OutDir',  fileparts(mfilename('fullpath')), @(x) ischar(x) || isstring(x));
    ip.parse(varargin{:});
    comment = char(ip.Results.Comment);
    outDir  = char(ip.Results.OutDir);

    %% ===== Fit-data settings (as used for every run of the sweep) ========
    % FCM samples closer together than this (h) are technical replicates of
    % one burst.
    cfg.burstGap = 0.2;
    % The initial condition is taken at the first FCM burst that has a
    % measured OD and an HPLC sample within this many hours (prepareFitData).
    cfg.anchorMaxHplcGap = 0.5;
    % Experiments whose FCM fractions are converted to biomass with a smoothed
    % biomass curve instead of the OD on the same row (exp 4: its OD probe
    % stops before FCM sampling ends).
    cfg.smoothBiomassExps = 4;
    %% =====================================================================

    resultsDir = fullfile(outDir, 'results');
    logsDir    = fullfile(outDir, 'logs');
    if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
    if ~exist(logsDir, 'dir'),    mkdir(logsDir);    end

    runTag = sprintf('model%s_train%s_test%d_%s', modelName, ...
        strjoin(string(trainExps), ''), testExp, string(datetime('now', 'Format', 'yyyyMMdd_HHmm')));

    diary(fullfile(logsDir, [runTag '.txt']));
    diaryOff = onCleanup(@() diary('off'));
    fprintf('=== RUN: %s ===\n', runTag);
    fprintf('Comment: %s\n', comment);
    fprintf('Seed: %d\n', seed);
    fprintf('Started: %s\n\n', string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));

    %% --- Training data ----------------------------------------------------
    nExp = numel(trainExps);
    expStructs     = cell(1, nExp);
    Dcell          = cell(1, nExp);
    for k = 1:nExp
        expStructs{k} = loadExperimentData(trainExps(k));
        Dcell{k}      = prepareFitData(expStructs{k}, trainExps(k), cfg);
    end

    %% --- Model and parameters ------------------------------------------------
    P         = getParamConfig(modelName);
    modelSpec = getModelSpec(modelName);

    % Single-population models (0/0B) have no strain resolution, so their FCM
    % fraction data is not fitted. Everything else is identical to A-D.
    if ~modelSpec.usesFractions
        for k = 1:nExp
            Dcell{k}.expstates(:, 2:4) = NaN;
        end
        fprintf('Model %s: FCM fraction columns excluded from the fit.\n', modelName);
    end

    expstates_cell = cellfun(@(D) D.expstates, Dcell, 'UniformOutput', false);
    exptime_cell   = cellfun(@(D) D.tSim,      Dcell, 'UniformOutput', false);
    qair_cell      = cellfun(@(D) D.qairFcn,   Dcell, 'UniformOutput', false);
    initCond_cell  = cellfun(@(D) D.x0,        Dcell, 'UniformOutput', false);

    % Stage 1 fits the kinetic rows only (P29 -- named after model D's 29).
    co2Idx = find(endsWith(P.name, 'co2'));
    kinIdx = setdiff((1:height(P))', co2Idx);
    P29    = P(kinIdx, :);
    pLowLog  = log10(P29.low);
    pHighLog = log10(P29.high);

    % The model functions take the full parameter vector; the CO2 yields get
    % a placeholder that no stage-1 equation reads.
    yco2Nominal = sqrt(P.low(co2Idx) .* P.high(co2Idx));
    stage1Model = @(t, y, pk, qair) modelSpec.kinetic(t, y, expandParams(pk, kinIdx, co2Idx, yco2Nominal));

    %% --- Stage 1: kinetic parameters ---------------------------------------
    sol_ga = geneticSearch(pLowLog, pHighLog, qair_cell, initCond_cell, expstates_cell, ...
                           exptime_cell, stage1Model, seed);
    sol    = paramLeastSquares(sol_ga.p_log, pLowLog, pHighLog, qair_cell, initCond_cell, ...
                               expstates_cell, exptime_cell, stage1Model);

    %% --- Stage 2: CO2 yields -------------------------------------------------
    fprintf('\n=== Stage 2: CO2 yields ===\n');
    co2sol = estimateCO2Yields(expandParams(sol.p, kinIdx, co2Idx, yco2Nominal), ...
                               Dcell, P, modelSpec.co2Basis);
    fprintf('CO2 stage: lsqlin exitflag = %d, resnorm = %.4g\n', co2sol.exitflag, co2sol.resnorm);
    for k = 1:nExp
        fprintf('  exp %d: CO2 RMS = %.4g (N = %d)\n', Dcell{k}.expNum, co2sol.rms(k), co2sol.N(k));
    end

    %% --- Result ---------------------------------------------------------------
    pFinal = expandParams(sol.p, kinIdx, co2Idx, co2sol.Yco2);
    ciFull = nan(height(P), 2);
    ciFull(kinIdx, :) = sol.ci_linear;   % no CI for the linear CO2 stage
    disp(table(P.name, pFinal, ciFull(:,1), ciFull(:,2), P.unit, ...
        'VariableNames', {'Parameter','Estimate','CI_lo','CI_hi','Unit'}));

    [AIC, BIC] = logRunResult(sol, P29, modelName, runTag, seed, comment, resultsDir);

    save(fullfile(resultsDir, [runTag '.mat']), 'sol', 'co2sol', 'pFinal', 'runTag', ...
        'trainExps', 'testExp', 'modelName', 'AIC', 'BIC', 'initCond_cell', 'Dcell', ...
        'sol_ga', 'P', 'P29', 'expStructs', 'seed', 'cfg');

    out = struct('runTag', runTag, 'AIC', AIC, 'BIC', BIC, 'sol', sol);
end
