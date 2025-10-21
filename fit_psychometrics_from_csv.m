% fit_psychometrics_from_csv.m
% Minimal pipeline to fit cumulative-Gaussian PFs via anamax
% Conditions: Active, Passive, Control, Sensory Attenuation (SA)
% Output: table with PSE and JND (Weber fraction and McFadden's R² available but commented)

clear; clc;

%% ----------------------- USER SETUP --------------------------------------
dataPath = '';        % e.g., 'D:\myproject\data\' (include trailing separator) or '' to use current folder

% Array of subjects (file name stems)
subjects = {'S01','S02'};  % set your subjects here

% Number of blocks per condition (adjust if needed)
maxActiveBlocks  = 5;    % expects <SUBJ>_1.csv ... <SUBJ>_5.csv
maxPassiveBlocks = 5;    % expects <SUBJ>_p_1.csv ... <SUBJ>_p_5.csv

% Column indices per condition (match your recorded CSV layout)
% Active files: response in col2, intensity in col3, latency col7, duration col8 (as before)
idx.Active.intensity = 3; idx.Active.response = 2; idx.Active.latency = 7; idx.Active.duration = 8;

% Passive files: response in col2, intensity in col3 (no latency/duration used)
idx.Passive.intensity = 3; idx.Passive.response = 2; idx.Passive.latency = NaN; idx.Passive.duration = NaN;

% Control files: <SUBJ>_C.csv with intensity in col1, response in col2
idx.Control.intensity = 1; idx.Control.response = 2; idx.Control.latency = NaN; idx.Control.duration = NaN;

% SA files: <SUBJ>_SA.csv with intensity in col1, response in col2
idx.SA.intensity = 1; idx.SA.response = 2; idx.SA.latency = NaN; idx.SA.duration = NaN;

% Trial quality filtering (set false to disable)
use_nan_duration_reject = true;   % drops rows with NaN duration when duration exists
use_latency_duration_sd = true;   % drops rows beyond ±2 SD for latency/duration when both exist

% anamax options (as in your previous script)
anamaxOptions = {'Fit','Boot10','Over','Num','Silent'};

% Optional CSV export (leave '' to skip)
output_csv = '';  % e.g., 'PSE_JND_results.csv'
%% ------------------------------------------------------------------------

results = table();

for s = 1:numel(subjects)
    subj = subjects{s};

    % -------- LOAD & CONCAT: ACTIVE --------
    [I_active, R_active, L_active, D_active] = load_blocks( ...
        subj, @(subj,blk) sprintf('%s_%d.csv', subj, blk), ...
        maxActiveBlocks, dataPath, idx.Active);

    % -------- LOAD & CONCAT: PASSIVE -------
    [I_pass, R_pass, L_pass, D_pass] = load_blocks( ...
        subj, @(subj,blk) sprintf('%s_p_%d.csv', subj, blk), ...
        maxPassiveBlocks, dataPath, idx.Passive);

    % -------- LOAD: CONTROL ----------------
    [I_ctrl, R_ctrl, L_ctrl, D_ctrl] = load_single( ...
        fullfile(dataPath, sprintf('%s_C.csv', subj)), idx.Control);

    % -------- LOAD: SA ---------------------
    [I_sa, R_sa, L_sa, D_sa] = load_single( ...
        fullfile(dataPath, sprintf('%s_SA.csv', subj)), idx.SA);

    % -------- FILTERING (Active uses latency/duration if present) --------
    [I_active, R_active] = basic_filter(I_active, R_active, L_active, D_active, ...
        use_nan_duration_reject, use_latency_duration_sd);
    % Passive/control/SA typically lack latency/duration; filter is no-op
    [I_pass,  R_pass ]  = basic_filter(I_pass,  R_pass,  L_pass,  D_pass,  false, false);
    [I_ctrl,  R_ctrl ]  = basic_filter(I_ctrl,  R_ctrl,  L_ctrl,  D_ctrl,  false, false);
    [I_sa,    R_sa   ]  = basic_filter(I_sa,    R_sa,    L_sa,    D_sa,    false, false);

    % -------- FITS --------------------------------------------------------
    rows = {};
    rows = [rows; fit_one(subj, 'Active',  I_active, R_active, anamaxOptions)];
    rows = [rows; fit_one(subj, 'Passive', I_pass,   R_pass,   anamaxOptions)];
    rows = [rows; fit_one(subj, 'Control', I_ctrl,   R_ctrl,   anamaxOptions)];
    rows = [rows; fit_one(subj, 'SA',      I_sa,     R_sa,     anamaxOptions)];

    % Append to results table
    T = cell2table(rows, 'VariableNames', {'Subject','Condition','PSE','JND' ...
        % ,'WeberFraction','R2'  % <- uncomment both here and in fit_one() to include these columns
    });
    results = [results; T]; %#ok<AGROW>
end

disp(results);
if ~isempty(output_csv)
    writetable(results, output_csv);
end

%% ============================ FUNCTIONS =================================
function rows = fit_one(subj, condName, I, R, anamaxOptions)
    if isempty(I) || numel(I) < 8 || numel(unique(I)) < 3
        rows = {subj, condName, NaN, NaN ...
            % ,NaN,NaN   % WeberFraction, McFadden's R² (uncomment if enabling)
        };
        return;
    end

    Rmat = [I(:), R(:)];
    DEP  = '2';  % responses column in Rmat
    IDEP = '1';  % intensities column in Rmat

    try
        % PJJ assumed as in your previous script:
        %   PJJ(1,1) = PSE
        %   PJJ(3,1) = JND
        %   PJJ(5,1) = Weber fraction (if anamax provides it)
        [PJJ, ~, ~, ~, ~, ~] = anamax(Rmat, DEP, IDEP, anamaxOptions{:});
        PSE = PJJ(1,1);
        JND = PJJ(3,1);

        % ---- Optional metrics (kept commented by default) ----------------
        % To include Weber fraction and McFadden's R² in outputs:
        % 1) Uncomment these lines, and
        % 2) Uncomment the corresponding columns in the table headers above.
        %
        % WFR = PJJ(5,1);
        %
        % p_model = 0.5 * (1 + erf( (I(:) - PSE) ./ (JND * sqrt(2)) ));
        % p_model = clamp01(p_model);
        % p_null  = clamp01( mean(R(:)) ) * ones(size(R(:)));
        % logL_model = binom_loglike(R(:), p_model);
        % logL_null  = binom_loglike(R(:), p_null);
        % R2_mcfadden = 1 - (logL_model / logL_null);

        rows = {subj, condName, PSE, JND ...
            % ,WFR, R2_mcfadden
        };

    catch
        rows = {subj, condName, NaN, NaN ...
            % ,NaN,NaN
        };
    end
end

function [I, R, L, D] = load_blocks(subj, namefun, maxBlocks, dataPath, idxs)
    I = []; R = []; L = []; D = [];
    for b = 1:maxBlocks
        file = fullfile(dataPath, namefun(subj,b));
        if ~exist(file,'file'), continue; end
        T = readtable(file, 'VariableNamingRule','preserve');
        Ii = T{:, idxs.intensity}; if iscell(Ii), Ii = cell2mat(Ii); end
        Ri = T{:, idxs.response }; if iscell(Ri), Ri = cell2mat(Ri); end
        I = [I; Ii]; %#ok<AGROW>
        R = [R; Ri]; %#ok<AGROW>
        if ~isnan(idxs.latency)  && idxs.latency  <= width(T)
            Li = T{:, idxs.latency};  if iscell(Li), Li = cell2mat(Li); end
        else
            Li = nan(size(Ii));
        end
        if ~isnan(idxs.duration) && idxs.duration <= width(T)
            Di = T{:, idxs.duration}; if iscell(Di), Di = cell2mat(Di); end
        else
            Di = nan(size(Ii));
        end
        L = [L; Li]; %#ok<AGROW>
        D = [D; Di]; %#ok<AGROW>
    end
end

function [I, R, L, D] = load_single(file, idxs)
    I = []; R = []; L = []; D = [];
    if ~exist(file,'file'), return; end
    T = readtable(file, 'VariableNamingRule','preserve');
    I = T{:, idxs.intensity}; if iscell(I), I = cell2mat(I); end
    R = T{:, idxs.response }; if iscell(R), R = cell2mat(R); end
    L = nan(size(I)); D = nan(size(I));
end

function [Iout, Rout] = basic_filter(I, R, L, D, dropNaNDur, dropSD)
    keep = true(size(I));
    if dropNaNDur && any(~isnan(D))
        keep = keep & ~isnan(D);
    end
    if dropSD && any(~isnan(L)) && any(~isnan(D))
        lm = nanmean(L); ls = nanstd(L);
        dm = nanmean(D); ds = nanstd(D);
        out = (L < lm - 2*ls) | (L > lm + 2*ls) | (D < dm - 2*ds) | (D > dm + 2*ds);
        keep = keep & ~out;
    end
    Iout = I(keep); Rout = R(keep);
end

function y = clamp01(x)
    eps_val = 1e-10;
    y = min(max(x, eps_val), 1 - eps_val);
end

function ll = binom_loglike(resp01, p)
    ll = sum(resp01 .* log(p) + (1 - resp01) .* log(1 - p));
end
