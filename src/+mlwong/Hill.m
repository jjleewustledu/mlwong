classdef Hill < handle
    %% Fits the fraction of intact tracer in plasma.
    %  See also http://www.turkupetcentre.net/petanalysis/input_parent_fitting_hill.html
    %  
    %  Created 07-Jun-2023 16:42:57 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlwong/src/+mlwong.
    %  Developed on Matlab 9.14.0.2254940 (R2023a) Update 2 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        dt
        interp_method
        T
        Tout
        visualize_anneal
    end

    properties (Dependent)
        product
        results
 		timeInterpolants
    end

    methods %% GET, SET
        function g = get.product(this)
            g = this.product_;
        end
        function g = get.results(this)
            g = this.product_;
        end
        function g = get.timeInterpolants(this)
            if isempty(this.timeInterpolants_)
                secs = this.T.Time(end);
                assert(isnumeric(secs))
                g = 0:this.dt:secs-1;
                return
            end
            g = this.timeInterpolants_;
        end
        function set.timeInterpolants(this, s)
            this.timeInterpolants_ = s;
        end
    end

    methods
        function this = Hill(T, opts)
            arguments
                T table % ["Time", "FractionIntact"]
                opts.dt double = 1 % sampling interval (sec)
                opts.d_min double = 0.9
                opts.interp_method {mustBeTextScalar} = "makima"
                opts.visualize_anneal logical = false;
            end
            T = T(~isnan(T.FractionIntact), :);
            if max(T.FractionIntact) > 1
                % ensure fraction
                T.FractionIntact = T.FractionIntact/100;
            end
            this.T = T;
            this.dt = opts.dt;
            this.visualize_anneal = opts.visualize_anneal;

            % a, b, c, d, e; http://www.turkupetcentre.net/petanalysis/input_parent_fitting_hill.html
            e_max = T.Time(end)/4;
            seps = 0.00001;
            this.ks_lower = [seps, 1, seps, opts.d_min, 0];
            this.ks0 = [0.03, 1.5, 2*seps, 0.975, 60]; 
            this.ks_upper = [0.06447, 4, inf, 1, e_max];
            this.ks_names = {'a' 'b' 'c' 'd' 'e'};
            this.interp_method = opts.interp_method;
            this.solved_ = false;
        end
        function h = plot(this, varargin)
            T_pcnt_ = this.T;
            T_pcnt_.FractionIntact = 100*T_pcnt_.FractionIntact;
            Tout_pcnt_ = this.Tout;
            Tout_pcnt_.FractionIntact = 100*Tout_pcnt_.FractionIntact;

            h = figure;
            hold("on")
            plot(T_pcnt_, "Time", "FractionIntact", LineStyle="none", Marker="o", MarkerSize=9)
            plot(Tout_pcnt_, "Time", "FractionIntact", LineStyle="--", Color="#D95319", LineWidth=2)
            hold("off")
            xlabel("time (s)", FontSize=14, FontWeight="bold")
            ylabel("fraction of parent tracer (%)", FontSize=14, FontWeight="bold")
            if this.solved_
                title("Extended Hill Function", FontSize=18, FontWeight="bold")
                annotation('textbox', [.45 .35 .8 .5], 'String', sprintfModel(this), 'FitBoxToText', 'on', 'FontSize', 14, 'LineStyle', 'none')
                legend(["HPLC measurements", "Hill function fit"])
            else
                title("Interpolation " + this.interp_method, FontSize=18, FontWeight="bold")
                legend(["HPLC measurements", "interp1"])
            end
        end
        function Tout1 = interp1(this)
            T_ = this.T;
            Time = ascol(this.timeInterpolants);
            FractionIntact = interp1(T_.Time, T_.FractionIntact, Time, this.interp_method);
            this.Tout = table(Time, FractionIntact);
            Tout1 = this.Tout;

            FractionIntact = makima(T_.Time, T_.FractionIntact, Time);
            T1 = table(ascol(Time), ascol(FractionIntact), variableNames=["Time", "FractionIntact"]);
            h = figure;
            plot(T1, "Time", "FractionIntact");
            xlabel("Time");
            ylabel("Fraction Intact");
            title("Fraction Intact by Makima (replacing Hill)");
        end
        function Tout1 = solve(this)
            %% options
            options_fmincon = optimoptions('fmincon', ...
                'FunctionTolerance', 1e-12, ...
                'OptimalityTolerance', 1e-12, ...
                'TolCon', 1e-14, ...
                'TolX', 1e-14);
            if this.visualize_anneal
                options = optimoptions('simulannealbnd', ...
                    'AnnealingFcn', 'annealingboltz', ...
                    'FunctionTolerance', 1e-10, ...
                    'HybridFcn', {@fmincon, options_fmincon}, ...
                    'InitialTemperature', 20, ...
                    'MaxFunEvals', 50000, ...
                    'ReannealInterval', 200, ...
                    'TemperatureFcn', 'temperatureexp', ...
                    'Display', 'diagnose', ...
                    'PlotFcns', {@saplotbestx,@saplotbestf,@saplotx,@saplotf,@saplotstopping,@saplottemperature});
            else
                options = optimoptions('simulannealbnd', ...
                    'AnnealingFcn', 'annealingboltz', ...
                    'FunctionTolerance', 1e-10, ...
                    'HybridFcn', {@fmincon, options_fmincon}, ...
                    'InitialTemperature', 20, ...
                    'MaxFunEvals', 50000, ...
                    'ReannealInterval', 200, ...
                    'TemperatureFcn', 'temperatureexp');
            end

            %%

            times_sampled = double(this.T.Time);
            measurement = double(this.T.FractionIntact);
 			[ks_,loss,exitflag,output] = simulannealbnd( ...
                @(ks__) this.loss_function(ks__, times_sampled, measurement), ...
                this.ks0, this.ks_lower, this.ks_upper, options);
            this.product_ = struct('ks0', this.ks0, 'ks', ks_, 'loss', loss, 'exitflag', exitflag, 'output', output);            

            Time = ascol(this.timeInterpolants);
            FractionIntact = ascol(mlwong.Hill.solution(ks_, Time));
            this.Tout = table(Time, FractionIntact);
            Tout1 = this.Tout;
            this.solved_ = true;
        end
        function s = sprintfModel(this)
            s = sprintf('extended Hill function:\n');
            ks = this.results.ks;
            kl = this.ks_lower;
            ku = this.ks_upper;
            for ky = 1:length(ks)
                s = [s sprintf('\t%s = %g in [%g, %g]\n', this.ks_names{ky}, ks(ky), kl(ky), ku(ky))]; %#ok<AGROW>
            end

            times_sampled = double(this.T.Time);
            measurement = double(this.T.FractionIntact);
            s = [s sprintf('\tloss = %g\n', this.loss_function(ks, times_sampled, measurement))];
            s = [s sprintf('\tloss = %g\n', this.results.loss)];
        end
    end

    methods (Static)
        function loss = loss_function(ks, times_sampled, measurement)
            estimation = mlwong.Hill.solution(ks, times_sampled); % \in [0 1] 
            measurement1 = measurement/max(measurement); % \in [0 1] 
            positive = measurement1 > 0;
            eoverm = estimation(positive)./measurement1(positive);            
            loss = mean(abs(1 - eoverm));
        end
        function fp = solution(ks, t)
            a = ks(1);
            b = ks(2);
            c = ks(3);
            d = ks(4);
            e = ks(5);
            
            fp = nan(size(t));
            fp(t <= e) = d;
            rational = ((d - a)*(t - e).^b)./(c + (t - e).^b);
            fp(t > e) = d - rational(t > e);
        end
    end

    %% PRIVATE

    properties (Access = private)
        ks0
        ks_lower
        ks_names
        ks_upper
        product_
        solved_
        timeInterpolants_
    end

    methods (Access = private)
    end

    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
