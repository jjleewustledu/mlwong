classdef Test_Eisai_VAT < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 12-Jun-2025 22:25:49 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlwong/test/+mlwong_unittest.
    %  Developed on Matlab 25.1.0.2943329 (R2025a) for MACA64.  Copyright 2025 John J. Lee.
    
    properties
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlwong.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end

        function test_final_calculation(this)
            c = this.testObj.data_cell;

            fluid = 'WB';  % 'PL'
            
            times = {};
            kBq = {};
            for idx = 1:length(c)
                fc = c{idx}.final_calculation;
                fc = fc(contains(fc.Barcode, fluid), :);
                dt = mean([fc.drawStart, fc.drawFinish], 2);
                min = minutes(dt - dt(1));
                times{idx} = min; %#ok<*AGROW>
                kBq{idx} = fc.kBq;
            end

            % plot subjects with same colors

            figure;
            hold on;

            colors = parula(6);
            h = zeros(6, 1);  % Store handles for legend

            h(1) = plot(times{1}, kBq{1}, ':', 'Color', colors(1, :), 'LineWidth', 3);
            for subIdx = 2:6                
                scan1 = 2*subIdx - 2;
                scan2 = 2*subIdx - 1;

                % Store handle from first plot for legend
                h(subIdx) = plot(times{scan1}, kBq{scan1}, '-', 'Color', colors(subIdx, :), 'LineWidth', 1.5);
                plot(times{scan2}, kBq{scan2}, '--', 'Color', colors(subIdx, :), 'LineWidth', 1.5);
            end

            hold off;
            xlabel('Time (min)');
            ylabel('Activity (kBq)');
            title("Gamma counting for " + fluid);
            grid on;
            legend(h, ...
                ["EVA105", arrayfun(@(x) sprintf("EVA%i", x), this.testObj.EVA_ids(2:end))], ...
                'Location', 'best');
            fontsize(scale=1.618)
        end

        function test_tacs(this)
            c = this.testObj.data_cell;

            region = "CaudateNucl";  % "Putamen";
            regions = region + ["_l", "_r"];
            
            standardTimes = this.testObj.suvs_times_mid/60;  % min
            kBq = {};
            for idx = 1:length(c)
                S = c{idx}.suvs;

                % Extract all specified variables
                regionVars = S(:, regions);  % Extract columns by name
                regiaonAvg = mean(table2array(regionVars), 2);  % Convert to array and average

                times{idx} = standardTimes; %#ok<*AGROW>
                kBq{idx} = regiaonAvg;
            end

            % plot subjects with same colors

            figure;
            hold on;

            colors = parula(6);
            h = zeros(6, 1);  % Store handles for legend

            h(1) = plot(times{1}, kBq{1}, ':', 'Color', colors(1, :), 'LineWidth', 3);
            for subIdx = 2:6                
                scan1 = 2*subIdx - 2;
                scan2 = 2*subIdx - 1;

                % Store handle from first plot for legend
                h(subIdx) = plot(times{scan1}, kBq{scan1}, '-', 'Color', colors(subIdx, :), 'LineWidth', 1.5);
                plot(times{scan2}, kBq{scan2}, '--', 'Color', colors(subIdx, :), 'LineWidth', 1.5);
            end

            hold off;
            xlabel('Time (min)');
            ylabel('Activity (kBq/mL)');
            title("PET TAC for " + region);
            grid on;
            legend(h, ...
                ["EVA105", arrayfun(@(x) sprintf("EVA%i", x), this.testObj.EVA_ids(2:end))], ...
                'Location', 'best');
            fontsize(scale=1.618)
        end
    end
    
    methods (TestClassSetup)
        function setupEisai_VAT(this)
            import mlwong.*
            this.testObj_ = Eisai_VAT();
        end
    end
    
    methods (TestMethodSetup)
        function setupEisai_VATTest(this)
            this.testObj = this.testObj_;
            this.addTeardown(@this.cleanTestMethod)
        end
    end
    
    properties (Access = private)
        testObj_
    end
    
    methods (Access = private)
        function cleanTestMethod(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
