classdef ICC < handle
    %% line1
    %  line2
    %  
    %  Created 28-Jun-2025 23:33:56 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlwong/src/+mlwong.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    
    properties (Constant)
        REGIONS = { ...
            'Frontal', 'Temporal', 'Parietal', 'Occipital', ...
            'CaudateNucl', 'NuclAccumb', 'Putamen', ...
            'Hippocampus', 'Thalamus', ...
            'White'};  
        % 'Cingulate', 'Insula', 'Amygdala', 'Pallidum', ...
    end

    methods
        function this = ICC(varargin)
        end
    end

    methods (Static)
        function analysis()
            %% PET Test-Retest Reproducibility Analysis
            % Using MATLAB Statistics Toolbox ICC function

            import mlwong.ICC.load_pet_data

            clear; close all;

            %% Define analysis parameters
            regions = mlwong.ICC.REGIONS;

            % Define all methods with clear naming
            arterial_methods = {'twotcm_Vt', 'logan_Vt', 'ma2_Vt'};
            reference_methods_pons = {'mrtm2_pons_BPnd', 'logan_ref_pons_BPnd', 'srtm2_pons_BPnd'};
            reference_methods_cereb = {'mrtm2_cbm_BPnd', 'logan_ref_cbm_BPnd', 'srtm2_cbm_BPnd'};

            all_methods = [arterial_methods, reference_methods_pons, reference_methods_cereb];

            % Initialize results storage
            n_regions = length(regions);
            n_methods = length(all_methods);
            n_subjects = 5;

            % Main results matrices
            icc_values = NaN(n_regions, n_methods);
            icc_lower_ci = NaN(n_regions, n_methods);
            icc_upper_ci = NaN(n_regions, n_methods);
            icc_pvalues = NaN(n_regions, n_methods);

            % Supplementary metrics
            abs_diff_pct = NaN(n_regions, n_methods);
            wsCV = NaN(n_regions, n_methods);
            mean_test = NaN(n_regions, n_methods);
            mean_retest = NaN(n_regions, n_methods);

            %% Main analysis loop
            for r = 1:n_regions
                for m = 1:n_methods

                    % Load your data here (replace with actual data loading)
                    % Example structure - adjust to match your data format
                    try
                        % Assuming data stored in structured format
                        test_data = load_pet_data(regions{r}, all_methods{m}, 'test');
                        retest_data = load_pet_data(regions{r}, all_methods{m}, 'retest');

                        % Ensure column vectors
                        test_data = test_data(:);
                        retest_data = retest_data(:);

                        % Create data matrix for ICC
                        M = [test_data, retest_data];

                        % Check for valid data
                        if any(isnan(M(:))) || any(isinf(M(:)))
                            fprintf('Warning: Invalid data for %s - %s\n', regions{r}, all_methods{m});
                            continue;
                        end

                        %% Calculate ICC using Statistics Toolbox
                        % ICC(2,1): Two-way mixed, single measurement, absolute agreement
                        [r_icc, LB, UB, F, df1, df2, p] = ICC(M, 'A-1', 0.05);

                        % Store ICC results
                        icc_values(r, m) = r_icc;
                        icc_lower_ci(r, m) = LB;
                        icc_upper_ci(r, m) = UB;
                        icc_pvalues(r, m) = p;

                        %% Calculate supplementary metrics
                        % Mean values
                        mean_test(r, m) = mean(test_data);
                        mean_retest(r, m) = mean(retest_data);

                        % Absolute percentage difference
                        diff = test_data - retest_data;
                        avg = (test_data + retest_data) / 2;
                        abs_diff_pct(r, m) = mean(abs(diff) ./ avg) * 100;

                        % Within-subject coefficient of variation
                        wsCV(r, m) = std(diff) / (sqrt(2) * mean(avg)) * 100;

                    catch ME
                        fprintf('Error processing %s - %s: %s\n', regions{r}, all_methods{m}, ME.message);
                    end
                end
            end

            %% Create comprehensive visualizations
            % 1. Heatmap of ICC values
            figure('Position', [100 100 1200 800]);
            subplot(2,2,1);
            imagesc(icc_values);
            colorbar;
            caxis([0 1]);
            colormap(parula);
            title('ICC Values by Region and Method');
            set(gca, 'XTick', 1:n_methods, 'XTickLabel', strrep(all_methods, '_', ' '), ...
                'XTickLabelRotation', 45);
            set(gca, 'YTick', 1:n_regions, 'YTickLabel', regions);

            % Add text annotations for values
            for i = 1:n_regions
                for j = 1:n_methods
                    if ~isnan(icc_values(i,j))
                        text(j, i, sprintf('%.2f', icc_values(i,j)), ...
                            'HorizontalAlignment', 'center', 'FontSize', 8);
                    end
                end
            end

            % 2. Method comparison plot
            % subplot(2,2,2);
            % method_categories = [ones(1,length(arterial_methods)), ...
            %     2*ones(1,length(reference_methods_pons)), ...
            %     3*ones(1,length(reference_methods_cereb))];
            % 
            % boxplot(icc_values', method_categories);
            % set(gca, 'XTickLabel', {'Arterial', 'Ref: Pons', 'Ref: Cereb'});
            % ylabel('ICC');
            % title('ICC Distribution by Method Type');
            % ylim([0 1]);

            % 3. Regional average ICC
            subplot(2,2,3);
            regional_mean_icc = nanmean(icc_values, 2);
            regional_sem_icc = nanstd(icc_values, 0, 2) / sqrt(sum(~isnan(icc_values), 2));

            bar(regional_mean_icc);
            hold on;
            errorbar(1:n_regions, regional_mean_icc, regional_sem_icc, 'k', 'LineStyle', 'none');
            set(gca, 'XTick', 1:n_regions, 'XTickLabel', regions, 'XTickLabelRotation', 45);
            ylabel('Mean ICC ± SEM');
            title('Average ICC by Region');
            ylim([0 1]);
            line([0 n_regions+1], [0.75 0.75], 'Color', 'r', 'LineStyle', '--');
            text(n_regions/2, 0.77, 'Good reliability threshold', 'HorizontalAlignment', 'center');

            % 4. Scatter plot of test vs retest for best method
            subplot(2,2,4);
            [best_icc, best_idx] = max(icc_values(:));
            [best_region, best_method] = ind2sub(size(icc_values), best_idx);

            % Plot for best performing region/method combination
            test_best = load_pet_data(regions{best_region}, all_methods{best_method}, 'test');
            retest_best = load_pet_data(regions{best_region}, all_methods{best_method}, 'retest');

            scatter(test_best, retest_best, 100, 'filled');
            hold on;
            % Identity line
            plot_range = [min([test_best; retest_best]), max([test_best; retest_best])];
            plot(plot_range, plot_range, 'k--');
            xlabel('Test');
            ylabel('Retest');
            title(sprintf('Best ICC: %s - %s (ICC=%.3f)', regions{best_region}, ...
                strrep(all_methods{best_method}, '_', ' '), best_icc));
            axis equal;

            %% Generate summary report
            % Create results table
            fprintf('\n=== SUMMARY REPORT ===\n');
            fprintf('\nTop 5 Most Reproducible Combinations:\n');
            [sorted_icc, sorted_idx] = sort(icc_values(:), 'descend', 'MissingPlacement', 'last');
            for i = 1:min(5, sum(~isnan(sorted_icc)))
                [r_idx, m_idx] = ind2sub(size(icc_values), sorted_idx(i));
                fprintf('%d. %s - %s: ICC = %.3f [%.3f, %.3f]\n', ...
                    i, regions{r_idx}, all_methods{m_idx}, ...
                    icc_values(r_idx, m_idx), icc_lower_ci(r_idx, m_idx), icc_upper_ci(r_idx, m_idx));
            end

            % Method comparison summary
            fprintf('\n\nMethod Category Summary:\n');
            arterial_icc = icc_values(:, 1:length(arterial_methods));
            pons_icc = icc_values(:, length(arterial_methods)+1:length(arterial_methods)+length(reference_methods_pons));
            cereb_icc = icc_values(:, end-length(reference_methods_cereb)+1:end);

            fprintf('Arterial methods: Mean ICC = %.3f ± %.3f\n', ...
                nanmean(arterial_icc(:)), nanstd(arterial_icc(:)));
            fprintf('Reference (Pons): Mean ICC = %.3f ± %.3f\n', ...
                nanmean(pons_icc(:)), nanstd(pons_icc(:)));
            fprintf('Reference (Cerebellum): Mean ICC = %.3f ± %.3f\n', ...
                nanmean(cereb_icc(:)), nanstd(cereb_icc(:)));

            %% Export results
            % Save to Excel for easy sharing
            results_table = array2table(icc_values, 'VariableNames', all_methods, 'RowNames', regions);
            writetable(results_table, 'ICC_results.xlsx', 'WriteRowNames', true);

            % Save MATLAB workspace
            save('PET_reproducibility_results.mat');

            %% Create publication-ready figure
            figure('Position', [100 100 800 600]);
            % Select key regions and methods for cleaner visualization
            key_regions = mlwong.ICC.REGIONS;
            key_methods = {'twotcm_Vt', 'logan_Vt', 'twotcm_Vt', 'mrtm2_cbm_BPnd', 'logan_ref_cbm_BPnd', 'srtm2_cbv_BPnd'};

            [~, r_idx] = ismember(key_regions, regions);
            [~, m_idx] = ismember(key_methods, all_methods);

            selected_icc = icc_values(r_idx, m_idx);
            selected_ci_lower = icc_lower_ci(r_idx, m_idx);
            selected_ci_upper = icc_upper_ci(r_idx, m_idx);

            % Create grouped bar plot with error bars
            b = bar(selected_icc');
            hold on;

            % Add error bars
            ngroups = size(selected_icc, 2);
            nbars = size(selected_icc, 1);
            groupwidth = min(0.8, nbars/(nbars + 1.5));

            for i = 1:nbars
                x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
                errorbar(x, selected_icc(i,:), ...
                    selected_icc(i,:) - selected_ci_lower(i,:), ...
                    selected_ci_upper(i,:) - selected_icc(i,:), ...
                    'k', 'LineStyle', 'none');
            end

            % Formatting
            set(gca, 'XTickLabel', strrep(key_methods, '_', ' '));
            ylabel('ICC (95% CI)');
            legend(key_regions, 'Location', 'best');
            ylim([0 1]);
            line([0 ngroups+1], [0.75 0.75], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
            title('Test-Retest Reproducibility of PET Quantification Methods');

            % Save figure
            print('ICC_publication_figure', '-dpng', '-r300');
        end

        function data = load_pet_data(region, method, scan_type)
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

            ss = strsplit(method, '_');
            measure = ss{end};
            model = regexprep(method, "_" + measure + "$", "");
            data = [];
            for eva = ["EVA109", "EVA112", "EVA114", "EVA115", "EVA118"]
                tbl = data_struct.(eva).(scan_type).(model);
                data = [data, tbl{region, measure}]; %#ok<*AGROW>
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
