classdef TZ3108 < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 12-Jan-2024 23:29:16 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlwong/src/+mlwong.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2024 John J. Lee.
   
    properties
        N_ensemble

        products
        ptac
        rois
        signals
        tacs
    end

    methods
        function this = build_multinest(this)
            pwd0 = pushd(this.kmdata_.filepath);

            this.ptac = this.kmdata_.imaging_metab_corr_ptac();
            this.rois = lower(this.kmdata_.regions);
            % this.rois = { ...
            %     'frontal', 'hippocampus', 'occipital', ...
            %     'striatum', ...
            %     'white', 'cerebellum', 'midbrain'};
            % this.rois = { ...
            %     'wholebrain', ...
            %     'prefrontal', 'parietal', 'occipital', ...
            %     'thalamus', ...
            %     'corticalwhitematter', 'cerebellum', 'midbrain'};
            ROIs = this.rois';
            Nr = length(ROIs);

            k1 = NaN(Nr, 1);
            k2 = NaN(Nr, 1);
            k3 = NaN(Nr, 1);
            k4 = NaN(Nr, 1);
            VP = NaN(Nr, 1);
            VN_plus_VS = NaN(Nr, 1);
            V_star = NaN(Nr, 1);
            BP_ND = NaN(Nr, 1);
            BP_P = NaN(Nr, 1);
            loss = NaN(Nr, 1);
            for r = 1:Nr
                this.tacs{r} = this.kmdata_.imaging_tac(ROIs{r});
                ichise = mlnest.Ichise.create( ...
                    input_func=this.ptac, ...
                    scanner=this.tacs{r}, ...
                    model_kind="vasc", ...
                    tracer="tz3108");
                nest = mlnest.MultiNest(context=ichise);
                tic
                nest = nest.solve( ...
                    signal_model=@ichise.signalmodel, ...
                    verbose=false, ...
                    Nlive=100, ...
                    Nmcmc=0); % 25 -> 600
                toc
                
                fqfn = fullfile(this.kmdata_.filepath, sprintf("%s_%s.mat", stackstr(), ROIs{r}));
                nest.saveas(fqfn);

                %nest.plot_posteriors(singles=true);
                %saveFigures(sprintf("%s_r=%i_Nlive=200_Nmcmc=400", stackstr(), r), closeFigure=false);
                %close('all');

                this.signals{r} = ichise.simulate(nest.product);
                p = nest.product;
                p.mean_ks(1:4) = 1./p.mean_ks(1:4); 
                p.std_ks(1:4) = 1./p.std_ks(1:4); 
                p.V_star = p.mean_ks(5) + p.mean_ks(6);
                p.BP_ND = p.mean_ks(3)/p.mean_ks(4);
                p.BP_P = p.mean_ks(1)*p.mean_ks(3)/(p.mean_ks(2)*p.mean_ks(4));
                fprintf("======== %s ========\n", ROIs{r})
                disp(p)                
                this.products{r} = p;

                k1(r) = p.mean_ks(1)*60;
                k2(r) = p.mean_ks(2)*60;
                k3(r) = p.mean_ks(3)*60;
                k4(r) = p.mean_ks(4)*60;
                VP(r) = p.mean_ks(5);
                VN_plus_VS(r) = p.mean_ks(6);
                V_star(r) = p.V_star;
                BP_ND(r) = p.BP_ND;  
                BP_P(r) = p.BP_P;
                loss(r) = p.loss;
            end

            T = table(ROIs, k1, k2, k3, k4, VP, VN_plus_VS, V_star, BP_ND, BP_P, loss);
            disp(T)
            writetable(T, sprintf("%s_%s.csv", this.kmdata_.fqfp, stackstr()));
        end
        function h = plot_multinest(this)
            pwd0 = pushd(this.kmdata_.filepath);

            h = figure;
            % ptac
            ptac = this.ptac/1e3;
            ptac.json_metadata.timesMid = ptac.json_metadata.timesMid/60; % -> min
            plot_over_figure(ptac, '--', LineWidth=3); 
            hold on;

            % wholebrain
            % sig = this.signals{1}; % fitted
            % sig = sig/1e3; % -> kBq
            % sig.json_metadata.timesMid = sig.json_metadata.timesMid/60; % -> min
            % plot_over_figure(sig, ':', LineWidth=5);
            % tac = this.tacs{1}; % measured
            % tac = tac/1e3; % -> kBq
            % tac.json_metadata.timesMid = tac.json_metadata.timesMid/60; % -> min
            % plot_over_figure(tac, 's', MarkerSize=18, MarkerEdgeColor=[.2,.2,.2]);

            Nr = length(this.tacs);
            for r = 1:Nr
                sig = this.signals{r}; % fitted
                sig = sig/1e3; % -> kBq
                sig.json_metadata.timesMid = sig.json_metadata.timesMid/60; % -> min
                plot_over_figure(sig, '-', LineWidth=3);
            end
            for r = 1:Nr
                tac = this.tacs{r}; % measured
                tac = tac/1e3; % -> kBq
                tac.json_metadata.timesMid = tac.json_metadata.timesMid/60; % -> min
                plot_over_figure(tac, 'o', MarkerSize=12, MarkerEdgeColor=[.5,.5,.5]);
            end            
            hold off;
            ylabel("activity (kBq/cm^3)")
            xlabel("time (min)")
            %legend([{'ptac', 'wholebrain', ''}, this.rois(2:end), repmat({''}, [1, Nr-1])])
            legend([{'ptac'}, this.rois, repmat({''}, [1, Nr])])
            fontsize(scale=2);
            set(h, position=[100,100,1000,700]) % position and size of window on display
            saveFigures(stackstr(), closeFigure=false);
            %close('all')

            popd(pwd0);
        end

        function this = build_simul_anneal(this)
            T_ptac = this.kmdata_.table_metab_corr_ptac();
            T_tac = this.kmdata_.table_tac('WholeBrain');
            mdl = mlkinetics.Ichise2002Model.create();
            mdl.set_times_sampled(T_tac.Time); % sec
            mdl.set_artery_interpolated(T_ptac.pTAC);
            mdl.Data = struct( ...
                'measurement_sampled', asrow(T_tac.TAC), ...
                'times', asrow(T_tac.Time), ...
                'taus', this.kmdata_.taus);
            mdl.build_model( ...
                map=mdl.preferredMap_unknown, ...
                measurement=asrow(T_tac.TAC));

            sol = cell(this.N_ensemble, 1);
            loss = NaN(this.N_ensemble, 1);
            for e = 1:this.N_ensemble
                sol{e} = mdl.solver.solve(@mlkinetics.Ichise2002Model.loss_function);
                loss(e) = sol{e}.loss();
            end

            this.solver_ = mdl.solver.find_best_soln(sol, loss);
        end
        function h = plot_simul_anneal(this, opts)
            %% PLOT 0:length() -> this.dispersedAif();
            %       this.TimesSampled -> this.Measurement;
            %       this.TimesSampled -> this.model.sampled();
            %       
            %  also plot xs{i} -> ys{i};           

            arguments
                this mlwong.TZ3108
                opts.activityUnits {mustBeTextScalar} = "Bq/mL"
                opts.timeUnits {mustBeTextScalar} = "s"
                opts.tag {mustBeTextScalar} = stackstr()
                opts.xlim {mustBeNumeric} = [-10 10600] % sec
                opts.ylim {mustBeNumeric} = []
                opts.xs cell = {} % additional xs to plot
                opts.ys cell = {} % additional ys to plot
                opts.legends cell = {} % of additional xs, ys
                opts.colorArt {mustBeTextScalar} = "#A2142F" % maroon
                opts.colorMeas {mustBeTextScalar} = "#0072BD" % navy
                opts.colorModel {mustBeTextScalar} = "0072BD" % navy
                opts.colors cell = {} % consider #EDB120 ~ mustard
                opts.zoomArt double = 1
                opts.zoomMeas double = 1 
                opts.zoomModel double = 1 
                opts.zooms cell = {}
            end
            assert(length(opts.xs) == length(opts.ys))
            if ~isempty(opts.xs) && isempty(opts.colors)
                opts.colors = repmat("#EDB120", size(opts.xs));
            end

            solver = this.solver_;
            solver.zoom = struct( ...
                'zoomArt', opts.zoomArt, ...
                'zoomMeas', opts.zoomMeas, ...
                'zoomModel', opts.zoomModel, ...
                'zooms', opts.zooms);
            
            % var notations
            %ad = mlaif.AifData.instance();
            %tBuffer = ad.tBuffer;
            TS = solver.TimesSampled;
            TSInt = 0:TS(end);
            ArtInt = solver.ArteryInterpolated;  
            Meas = solver.Measurement;

            [qs,A] = solver.model.sampled(solver.ks, solver.Data, ArtInt, TS);
            Fitted = solver.rescaleModelEstimate(A*qs);
            
            % build legends
            legendCell = {};
            legendCell = [legendCell, sprintf('Arterial x%i', opts.zoomArt)];
            legendCell = [legendCell, sprintf('Measurement x%i', opts.zoomMeas)];
            legendCell = [legendCell, sprintf('Model x%i', opts.zoomModel)];
            legendCell = [legendCell, opts.legends]; % of additional xs, ys

            % plotting implementation
            [TSInt,ArtInt] = solver.match_lengths(TSInt,ArtInt);
            [TS, Meas] = solver.match_lengths(TS, Meas);
            [TS1, Fitted] = solver.match_lengths(TS, Fitted);
            h = figure;
            hold("on");
            plot(TSInt, opts.zoomArt*ArtInt, '-', LineWidth=2, Color=opts.colorArt)
            plot(TS, opts.zoomMeas*Meas, 'o', MarkerSize=12, Color=opts.colorMeas)
            plot(TS1, opts.zoomModel*Fitted, '--', LineWidth=2, Color=opts.colorMeas)
            for ci = 1:length(opts.xs)
                [xs,ys] = solver.match_lengths(opts.xs{ci}, opts.ys{ci});
                plot(xs, opts.zooms{ci}*ys, ':', LineWidth=2, Color=opts.colors{ci})
            end
            legend(legendCell);
            if ~isempty(opts.xlim); xlim(opts.xlim); end
            if ~isempty(opts.ylim); ylim(opts.ylim); end
            xlabel(sprintf('times (%s)', opts.timeUnits))
            ylabel(sprintf('activity / (%s)', opts.activityUnits))
            annotation('textbox', [.5 .1 .3 .8], 'String', sprintfModel(solver), 'FitBoxToText', 'on', 'LineStyle', 'none') % 'FontSize', 10, 
            opts.tag = strrep(opts.tag, "_", " ");
            title([stackstr(use_spaces=true)+";"; string(opts.tag); ""], Interpreter="none")
            hold("off");
            set(h, position=[100,100,1000,700]) % position and size of window on display
        end
    end

    methods (Static)
        function this = create(kmd_fqfn, opts)
            arguments
                kmd_fqfn {mustBeFile}
                opts.N_ensemble double {mustBeInteger} = 10
            end

            this = mlwong.TZ3108();
            this.N_ensemble = opts.N_ensemble;
            this.kmdata_ = mlpmod.KMData.create(kmd_fqfn);
        end
    end

    %% PRIVATE

    properties (Access = private)
        kmdata_
        solver_
    end

    methods (Access = private)
        function this = TZ3108()
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
