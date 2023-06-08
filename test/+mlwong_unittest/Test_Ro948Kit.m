classdef Test_Ro948Kit < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 06-Jun-2023 18:35:41 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlwong/test/+mlwong_unittest.
    %  Developed on Matlab 9.14.0.2254940 (R2023a) Update 2 for MACI64.  Copyright 2023 John J. Lee.
    
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
        function test_ctor(this)
            cd(fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01384/derivatives/sub-R21-004/ses-20230330/twilite'))
            obj = mlwong.Ro948Kit();
            disp(obj)
        end
        function test_call(this)
            cd(fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01384/derivatives/sub-R21-004/ses-20230330/twilite'))
            obj = mlwong.Ro948Kit();
            obj.call();
        end
        function test_deconv(this)
            cd(fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01384/derivatives/sub-R21-004/ses-20230330/twilite'));
            inveff = 68/0.27; % counts/s to Bq/mL
            idx0 = 1262; % TOI - 30 s
            idx_toss = 422; % trunc at nominal 5 min
            %idx_5min = 1592; % nominal 5 min blood draw, 11:42:01

            crv = mlswisstrace.CrvData.createFromFilename('R21_004_03302023_D1.crv');            
            disp(crv)

            % deconvBayes
            M_ = crv.timetable().Coincidence(idx0:end-idx_toss)*inveff;
            cath = mlswisstrace.Catheter_DT20190930( ...
                'Measurement', M_, ...
                'hct', 42.2, ...
                'tracer', '18F'); % t0 reflects rigid extension + Luer valve + cath in Twilite cradle
            M = zeros(size(crv.timetable().Coincidence));
            M(idx0:end-idx_toss) = cath.deconvBayes( ...
                't0_forced', 47, 'decay_correct', true, 'xlim', [0 347]);

            % write csv
            toi = datetime(2023,3,30,11,37,1, TimeZone="local");
            select = toi <= crv.time & crv.time <= toi + seconds(301);
            Time = seconds(crv.time(select) - toi);
            wbKBq_mL = ascol(M(select));
            T = table(Time, wbKBq_mL);            
            writetable(T, "R21_004_03302023_D1_deconv.csv");

            % plot
            figure
            plot(T, "Time", "wbKBq_mL", LineStyle="none", Marker="+", MarkerSize=6)
            xlabel("time (s)", FontSize=14)
            ylabel("wb activity (kBq/mL)", FontSize=14)
            title("From R21_004_03302023_D1.crv", FontSize=18, Interpreter="none")    
        end
        function test_Hill(this)
            cd(fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01384/derivatives/sub-R21-004/ses-20230330/chemistry'));
            T = readtable("R21_004_03302023_frac_intact.csv");
            T.Properties.VariableNames = ["Time", "FractionIntact"];
            T.Time = T.Time*60;
            hill = mlwong.Hill(T, dt=1, visualize_anneal=false);
            T1 = hill.solve();
            plot(hill)
            disp(hill)
            disp(hill.results.ks)
            disp(hill.results.sse)
            writetable(T1, "R21_004_03302023_Hill.csv")
        end
        function test_pow(this)
            cd(fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01384/derivatives/sub-R21-004/ses-20230330/chemistry'));
            T = readtable("R21_004_03302023_minimal.csv");
            toi = datetime(2023,3,30,11,37,1);
            draw = (T.drawFinish - T.drawStart)/2 + T.drawStart;
            draw.Year = toi.Year;
            draw.Month = toi.Month;
            draw.Day = toi.Day;
            Time = seconds(draw - toi);
            T = addvars(T, Time, Before=1, NewVariableNames="Time");
            PlasmaOverWb = T.plasmaKBq_mL ./ T.wbKBq_mL;
            T = addvars(T, PlasmaOverWb, NewVariableNames="PlasmaOverWb");
            U = T(T.Time < 30*60, :);
            
            % regresion for early times
            mdl = fitlm(U, "PlasmaOverWb~Time");
            disp(mdl)

            Time1 = ascol(0:30*60);
            Pow1 = feval(mdl, Time1);
            hold('on')
            plot(T, "Time", "PlasmaOverWb", LineStyle="none", Marker="o", MarkerSize=9)
            plot(Time1, Pow1, LineWidth=2)
            xlabel("time (s)", FontSize=14)
            ylabel("plasma activity / wb activity", FontSize=14)
            title("From R21_004_03302023_minimal.csv", FontSize=18, Interpreter="none")            
            
            V = table(Time1, Pow1, VariableNames=["Time", "PlasmaOverWb"]);
            writetable(V, "R21_004_03302023_pow.csv")
        end
    end
    
    methods (TestClassSetup)
        function setupRo948Kit(this)
            import mlwong.*
            %this.testObj_ = Ro948Kit();
        end
    end
    
    methods (TestMethodSetup)
        function setupRo948KitTest(this)
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
