classdef ICCPermutations < handle & mlwong.ICC
    %% line1
    %  line2
    %  
    %  Created 28-Jun-2025 23:33:56 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlwong/src/+mlwong.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    
    methods
        function this = ICCPermutations(varargin)
            this = this@mlwong.ICC(varargin{:});
        end
    end

    methods (Static)
        function analysis()
            %% PET Test-Retest Reproducibility Analysis - Small Sample Robust Methods
            % N=5 subjects with permutation-based inference

            import mlwong.ICCPermutations.make_readable_methods
            import mlwong.ICCPermutations.load_pet_data
            import mlwong.ICCPermutations.find_best_method

            clear; close all;

            %% Define analysis parameters
            regions = mlwong.ICC.REGIONS;

            % Define models with clear interpretable names
            models = struct();

            % Arterial input models
            models.arterial = {
                'twotcm',    % Two-tissue compartment model
                'logan', % Logan graphical with arterial input
                'ma2'      % Multilinear analysis (Ichise)
                };

            % Reference tissue models
            models.reference = {
                'mrtm2',   % Multilinear reference tissue model
                'logan_ref', % Logan reference tissue
                'srtm2'    % Simplified reference tissue model
                };

            % Define model characteristics for interpretation
            model_properties = struct();

            % Arterial models
            model_properties.('twotcm').params = 'Vt, Vs, K1, k2, k3, k4, BPnd';
            model_properties.('twotcm').benefits = 'Gold standard; Full kinetic characterization';
            model_properties.('twotcm').risks = 'Sensitive to noise; Requires arterial sampling; Parameter identifiability issues';

            model_properties.('logan').params = 'Vt, Intercept';
            model_properties.('logan').benefits = 'Robust; Fewer parameters; Visual QC possible';
            model_properties.('logan').risks = 'Bias with noise; Requires equilibrium time (t*)';

            model_properties.('ma2').params = 'Vt, Vs';
            model_properties.('ma2').benefits = 'No t* selection; Computationally efficient';
            model_properties.('ma2').risks = 'Less intuitive; Still requires arterial input';

            % Reference tissue models
            model_properties.('mrtm2').params = 'BPnd, k2'', R1';
            model_properties.('mrtm2').benefits = 'No arterial sampling; Uses fixed k2prime';
            model_properties.('mrtm2').risks = 'Assumes reference region has no specific binding; k2prime dependency';

            model_properties.('logan_ref').params = 'BPnd, DVR, Intercept, k2''';
            model_properties.('logan_ref').benefits = 'Simple; Robust for reversible tracers';
            model_properties.('logan_ref').risks = 'Noise-induced bias; Requires equilibrium';

            model_properties.('srtm2').params = 'BPnd, R1, k2''';
            model_properties.('srtm2').benefits = 'Models delivery differences; No arterial sampling';
            model_properties.('srtm2').risks = 'Assumes single tissue compartment in reference; Computational complexity';

            % Reference regions
            ref_regions = {'pons', 'cbm'};

            %% Initialize storage
            n_subjects = 5;
            n_permutations = 10000; % For permutation testing, prefer 10000

            % Prepare all method combinations
            all_methods = {};
            method_types = {}; % Track if arterial or reference
            reference_used = {}; % Track which reference region
            measure_used = {}; % Most representative measure for method

            for i = 1:length(models.arterial)
                all_methods{end+1} = models.arterial{i};
                method_types{end+1} = 'arterial';
                reference_used{end+1} = 'none';
                measure_used{end+1} = 'Vt';
            end

            for ref = 1:length(ref_regions)
                for i = 1:length(models.reference)
                    all_methods{end+1} = sprintf('%s_%s', models.reference{i}, ref_regions{ref});
                    method_types{end+1} = 'reference';
                    reference_used{end+1} = ref_regions{ref};
                    measure_used{end+1} = 'BPnd';
                end
            end

            readable_methods = make_readable_methods(all_methods);

            n_methods = length(all_methods);
            n_regions = length(regions);

            %% Main analysis with small sample considerations

            % Initialize results
            icc_values = NaN(n_regions, n_methods);
            icc_lower_ci = NaN(n_regions, n_methods);
            icc_upper_ci = NaN(n_regions, n_methods);
            icc_perm_pvalue = NaN(n_regions, n_methods);

            % Additional small-sample metrics
            spearman_rho = NaN(n_regions, n_methods);
            kendall_tau = NaN(n_regions, n_methods);
            individual_differences = cell(n_regions, n_methods);

            %% Analysis loop with permutation testing
            fprintf('Starting analysis with permutation testing...\n');

            for r = 1:n_regions
                fprintf('Processing region: %s\n', regions{r});

                for m = 1:n_methods
                    try
                        % Assuming data stored in structured format
                        test_data = load_pet_data(regions{r}, all_methods{m}, measure_used{m}, 'test');
                        retest_data = load_pet_data(regions{r}, all_methods{m}, measure_used{m}, 'retest');

                        % Store individual differences for transparency
                        individual_differences{r,m} = table(...
                            (1:n_subjects)', test_data, retest_data, ...
                            'VariableNames', {'Subject', 'Test', 'Retest'});

                        %% Calculate ICC with original data
                        M = [test_data, retest_data];
                        [icc_obs, LB, UB, ~, ~, ~, p_param] = ICC(M, 'A-1', 0.05);

                        icc_values(r, m) = icc_obs;

                        %% Permutation-based confidence intervals and p-value
                        % More robust for small samples
                        icc_perm_dist = zeros(n_permutations, 1);

                        for perm = 1:n_permutations
                            % Randomly swap test/retest within subjects
                            swap_idx = rand(n_subjects, 1) > 0.5;
                            M_perm = M;
                            M_perm(swap_idx, :) = M_perm(swap_idx, [2 1]);

                            % Calculate ICC for permuted data
                            icc_perm_dist(perm) = ICC(M_perm, 'A-1', 0.05);
                        end

                        % Permutation-based 95% CI
                        icc_lower_ci(r, m) = prctile(icc_perm_dist, 2.5);
                        icc_upper_ci(r, m) = prctile(icc_perm_dist, 97.5);

                        % Permutation p-value (test against null of ICC=0)
                        % Create null distribution by permuting across all subjects
                        icc_null_dist = zeros(n_permutations, 1);
                        all_values = M(:);

                        for perm = 1:n_permutations
                            % Random pairing under null hypothesis
                            perm_values = all_values(randperm(length(all_values)));
                            M_null = reshape(perm_values, n_subjects, 2);
                            icc_null_dist(perm) = ICC(M_null, 'A-1', 0.05);
                        end

                        icc_perm_pvalue(r, m) = mean(abs(icc_null_dist) >= abs(icc_obs));

                        %% Non-parametric correlation measures
                        % More robust for small samples
                        spearman_rho(r, m) = corr(test_data, retest_data, 'Type', 'Spearman');
                        kendall_tau(r, m) = corr(test_data, retest_data, 'Type', 'Kendall');

                    catch ME
                        fprintf('Error in %s - %s: %s\n', regions{r}, all_methods{m}, ME.message);
                    end
                end
            end

            %% Create comprehensive visualization with model interpretation
            figure('Position', [50 50 1400 900]);

            % 1. ICC heatmap with significance markers
            subplot(2,3,1);
            imagesc(icc_values);
            colormap(magma(256));
            colorbar;
            caxis([0 1]);
            title('ICC Values (* = p<0.05 by permutation)');

            % Add significance stars
            hold on;
            [sig_r, sig_m] = find(icc_perm_pvalue < 0.05);
            plot(sig_m, sig_r, 'g*', 'MarkerSize', 12);

            set(gca, 'XTick', 1:n_methods, 'XTickLabel', strrep(readable_methods, '_', ' '), ...
                'XTickLabelRotation', 45);
            set(gca, 'YTick', 1:n_regions, 'YTickLabel', regions);

            % 2. Model comparison with individual points
            subplot(2,3,2);
            model_categories = [];
            for m = 1:n_methods
                if strcmp(method_types{m}, 'arterial')
                    model_categories(m) = 1;
                elseif contains(all_methods{m}, 'pons')
                    model_categories(m) = 2;
                else
                    model_categories(m) = 3;
                end
            end

            % Show individual region points
            hold on;
            for cat = 1:3
                cat_methods = find(model_categories == cat);
                cat_data = icc_values(:, cat_methods); %#ok<*FNDSB>

                x_positions = cat + (rand(size(cat_data)) - 0.5) * 0.3;
                plot(x_positions(:), cat_data(:), 'o', 'MarkerSize', 6);
            end

            % Add median lines
            for cat = 1:3
                cat_methods = find(model_categories == cat);
                cat_data = icc_values(:, cat_methods);
                median_val = nanmedian(cat_data(:));
                plot([cat-0.4 cat+0.4], [median_val median_val], 'r-', 'LineWidth', 2);
            end

            set(gca, 'XTick', 1:3, 'XTickLabel', {'Arterial', 'Ref: Pons', 'Ref: Cerebellum'});
            ylabel('ICC');
            title('Model Type Comparison (red = median)');
            ylim([-0.1 1.1]);
            grid on;

            % 3. Individual subject variability
            subplot(2,3,3);
            % Select best performing method for each category
            best_arterial = find_best_method(icc_values, method_types, 'arterial');
            best_ref = find_best_method(icc_values, method_types, 'reference');

            % Plot example individual differences
            region_idx = find(strcmp(regions, 'Frontal')); % Example region

            if ~isempty(region_idx)
                data_art = individual_differences{region_idx, best_arterial};
                data_ref = individual_differences{region_idx, best_ref};

                plot(1:n_subjects, data_art.Test, 'bo-', 'LineWidth', 2);
                hold on;
                plot(1:n_subjects, data_art.Retest, 'b^--', 'LineWidth', 2);
                plot(1:n_subjects, data_ref.Test, 'ro-', 'LineWidth', 2);
                plot(1:n_subjects, data_ref.Retest, 'r^--', 'LineWidth', 2);

                xlabel('Subject');
                ylabel('Parameter Value');
                title('Best Individual Test-Retest: Frontal');
                legend({[readable_methods{best_arterial} ' Vt Test'], [readable_methods{best_arterial} ' Vt Retest'], ...
                    [readable_methods{best_ref} ' BPnd Test'], [readable_methods{best_ref} ' BPnd Retest']}, ...
                    'Location', 'best', 'Interpreter', 'none');
            end

            % 4. Correlation comparison (Spearman vs ICC)
            subplot(2,3,4);
            scatter(icc_values(:), spearman_rho(:), 30, 'filled');
            hold on;
            plot([0 1], [0 1], 'k--');
            xlabel('ICC');
            ylabel('Spearman \rho');
            title('ICC vs Spearman Correlation');
            text(0.1, 0.9, sprintf('r = %.3f', corr(icc_values(:), spearman_rho(:), 'rows', 'complete')));
            grid on;

            % 5. Model characteristics table
            subplot(2,3,5);
            axis off;
            model_names = unique([models.arterial, models.reference], 'stable');
            % model_names = make_readable_methods(model_names);
            text_y = 1;
            text(0.5, text_y, 'Model Characteristics', 'FontWeight', 'bold', 'FontSize', 12, ...
                'HorizontalAlignment', 'center');

            for i = 1:length(model_names)
                model = model_names{i};
                if isfield(model_properties, model)
                    text_y = text_y - 0.0618;
                    text(0, text_y, [model ':'], 'FontWeight', 'bold');
                    text_y = text_y - 0.033;
                    text(0.05, text_y, ['Parameters available: ' model_properties.(model).params], ...
                        'FontSize', 8, 'Color', [0.1 0.1 0.1]);
                    text_y = text_y - 0.033;
                    text(0.05, text_y, ['Benefits: ' model_properties.(model).benefits], ...
                        'FontSize', 8, 'Color', [0 0.5 0]);
                    text_y = text_y - 0.033;
                    text(0.05, text_y, ['Risks: ' model_properties.(model).risks], ...
                        'FontSize', 8, 'Color', [0.8 0 0]);
                end
            end

            % 6. Bootstrap confidence interval width comparison
            subplot(2,3,6);
            ci_width = icc_upper_ci - icc_lower_ci;
            mean_ci_width = nanmean(ci_width, 1);
            std_ci_width = nanstd(ci_width, 0, 1);

            bar(1:n_methods, mean_ci_width);
            hold on;
            errorbar(1:n_methods, mean_ci_width, std_ci_width, 'k', 'LineStyle', 'none');
            set(gca, 'XTick', 1:n_methods, 'XTickLabel', strrep(readable_methods, '_', ' '), ...
                'XTickLabelRotation', 45);
            ylabel('95% CI Width');
            title('Uncertainty by Method (Permutation-based)');

            %% Generate detailed report accounting for small sample
            fprintf('\n=== DETAILED REPORT FOR SMALL SAMPLE ANALYSIS ===\n');
            fprintf('N = %d subjects\n', n_subjects);
            fprintf('Permutation iterations = %d\n', n_permutations);

            fprintf('\n--- Model Performance Summary ---\n');
            for model = unique([models.arterial, models.reference], 'stable')
                model_indices = find(contains(all_methods, model));
                if ~isempty(model_indices)
                    model_icc = icc_values(:, model_indices);
                    fprintf('\n%s:\n', model{1});
                    fprintf('  Median ICC = %.3f [Range: %.3f - %.3f]\n', ...
                        nanmedian(model_icc(:)), nanmin(model_icc(:)), nanmax(model_icc(:)));
                    fprintf('  Regions with ICC > 0.75: %d/%d\n', ...
                        sum(model_icc(:) > 0.75), sum(~isnan(model_icc(:))));

                    % Report specific concerns for small samples
                    if strcmp(model{1}, 'twotcm')
                        fprintf('  ⚠ Small sample warning: With %d subjects, %d parameters may be unstable\n', ...
                            n_subjects, length(model_properties.('twotcm').params));
                    end
                end
            end

            %% Save results with full transparency
            results_struct = struct();
            results_struct.icc_values = icc_values;
            results_struct.icc_permutation_ci = [icc_lower_ci, icc_upper_ci];
            results_struct.permutation_pvalues = icc_perm_pvalue;
            results_struct.spearman_rho = spearman_rho;
            results_struct.kendall_tau = kendall_tau;
            results_struct.individual_data = individual_differences;
            results_struct.n_subjects = n_subjects;
            results_struct.n_permutations = n_permutations;
            results_struct.regions = regions;
            results_struct.methods = all_methods;
            results_struct.model_properties = model_properties;

            save('PET_reproducibility_small_sample_results.mat', 'results_struct');

            %% Export key results table for manuscript
            % Create table with method recommendations
            recommendation_table = table();
            row = 1;

            for m = 1:n_methods
                % Calculate summary statistics across regions
                method_icc = icc_values(:, m);
                valid_icc = method_icc(~isnan(method_icc));

                if ~isempty(valid_icc)
                    recommendation_table.Method{row} = readable_methods{m};
                    recommendation_table.Type{row} = method_types{m};
                    recommendation_table.Reference{row} = reference_used{m};
                    recommendation_table.Median_ICC(row) = median(valid_icc);
                    recommendation_table.IQR_ICC(row) = iqr(valid_icc);
                    recommendation_table.Prob_Good_ICC(row) = sum(valid_icc > 0.75) / length(valid_icc);
                    recommendation_table.Mean_CI_Width(row) = mean(ci_width(~isnan(method_icc), m));

                    % Small sample specific recommendation
                    if recommendation_table.Median_ICC(row) > 0.8 && recommendation_table.Mean_CI_Width(row) < 0.3
                        recommendation_table.Recommendation{row} = 'Recommended';
                    elseif recommendation_table.Median_ICC(row) > 0.7
                        recommendation_table.Recommendation{row} = 'Acceptable';
                    else
                        recommendation_table.Recommendation{row} = 'Use with caution';
                    end

                    row = row + 1;
                end
            end

            writetable(recommendation_table, 'Method_Recommendations_N5.xlsx');
            
        end

        %% Helper functions

        function readable_methods = make_readable_methods(all_methods)
            % Replace multiple substrings
            replacements = { ...
                'twotcm', '2TCM'; ...
                'logan', 'Logan'; ...
                'ma2', 'Ichise MA2'; ...
                'mrtm2', 'Ichise MRTM2'; ...
                'srtm2', 'SRTM2'};
            readable_methods = all_methods;
            for i = 1:size(replacements, 1)
                readable_methods = cellfun(@(x) strrep(x, replacements{i,1}, replacements{i,2}), ...
                    readable_methods, ...
                    'UniformOutput', false);
            end
        end

        function data = load_pet_data(region, model, measure, scan_type)
            % Example: Load from structured .mat files
            % filename = sprintf('PET_data_%s_%s_%s.mat', region, method, scan_type);
            % data_struct = load(filename);
            % data = data_struct.values;

            if strcmpi(scan_type, 'test'); scan_type = 'PET1'; end
            if strcmpi(scan_type, 'retest'); scan_type = 'PET2'; end

            persistent data_struct;
            if isempty(data_struct)
                EVAT = mlwong.Eisai_VAT;
                data_struct = EVAT.struct_for_anat_regions(EVAT.data_struct);
            end

            data = [];
            for eva = ["EVA109", "EVA112", "EVA114", "EVA115", "EVA118"]
                tbl = data_struct.(eva).(scan_type).(model);
                data = [data, tbl{region, measure}]; %#ok<*AGROW>
            end
            data = ascol(data);
        end

        function best_idx = find_best_method(icc_values, method_types, type)
            type_indices = find(strcmp(method_types, type));
            mean_icc = nanmean(icc_values(:, type_indices), 1);
            [~, best_within_type] = max(mean_icc);
            best_idx = type_indices(best_within_type);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
