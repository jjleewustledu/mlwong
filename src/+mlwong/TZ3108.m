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

    properties (Constant)
        home = "/Users/jjlee/Documents/PapersMine/PapersInProgress/Will_Tu/minus-TZ3108";
        bids_home = "/Users/jjlee/Singularity/TZ3108/derivatives"
    end

    properties (Dependent)
        proc_verified_kmdata_xls
        proc_aif_kmdata_xls
    end

    methods  %% GET
        function g = get.proc_verified_kmdata_xls(this)
            g = mglob(fullfile(this.home, "sub-*/ses-*/chemistry/sub-*_ses-*_proc-verified-pkin.kmData.xls"));
        end
        function g = get.proc_aif_kmdata_xls(this)
            g = mglob(fullfile(this.home, "sub-*/ses-*/chemistry/sub-*_ses-*_proc-aif-*-pkin.kmData.xls"));
        end
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
            this.pmod_ = mlpmod.Pmod(filename=kmd_fqfn);
        end

        function create_all_niftis()
            this = mlwong.TZ3108();
            pwd0 = pushd(this.home);
            
            globbed = [mglob(this.proc_verified_kmdata_xls), mglob(this.proc_aif_kmdata_xls)];
            for idx = 1:length(globbed)
                p = mlpmod.Pmod(filename=globbed(idx));
                [aif_,tacs_] = p.build_nifti();
                aif_.save();
                tacs_.save();
            end

            popd(pwd0);
        end

        function [t,h] = create_population_aif()
            %% Normalizes measured AIFs with time average of PET, \int^T dt whole_brain(t)/T,
            %  then interpolates to 1 sec, then obtains mean across measurements.
            %  Saves table to mlwong.TZ3108.home as stackstr() + ".mat".
            %  
            globbed = mglob( ...
                fullfile(mlwong.TZ3108.home, "sub-*/ses-*/chemistry/sub-*_ses-*_proc-verified-pkin.kmData.xls"));
            Ng = length(globbed);
            assert(Ng > 0)
            ps = cell(1, Ng);
            aifs = cell(1, Ng);
            suv = [];
            for idx = 1:Ng
                ps{idx} = mlpmod.Pmod(filename=globbed(idx)); 
                aifs{idx} = ps{idx}.aif_as_suv(interp_method="pchip", T=150); 
                suv = [suv, aifs{idx}.suv]; %#ok<AGROW>
            end
            suv = mean(suv, 2);
            time_min = aifs{1}.time_min;
            t = table(time_min, suv);
            save(fullfile(mlwong.TZ3108.home, stackstr(use_underscores=true) + ".mat"), "t");

            h = figure; 
            hold("on");
            for idx = 1:Ng
                plot(ps{idx}.aif_as_suv, "time_min", "suv", LineStyle="-", Marker="o", MarkerSize=8); 
            end
            plot(t, "time_min", "suv", LineStyle="-", LineWidth=2, Color="k");
            legend([ ...
                "n.h. primate 1, (-), Yun122", ...
                "n.h. primate 2, (-), N_2O ", ...
                "n.h. primate 2, (-)", ...
                "n.h. primate 3, (+)", ...
                "n.h. primate 4, (-)", ...
                "cohort"], Interpreter="none");
            xlim([-10, 160]);
            xlabel("time (min)");
            ylabel("AIF SUV (g/cm^3");
            title("Cohort-based AIF from metabolite-corrected plasma AIFs");
            hold("off");
            fontsize(scale=1.25);
        end

        function [htiled,h] = plot_age_dependence(fns, opts)
            arguments
                fns {mustBeText}
                opts.lwidth double = 2
                opts.fscale double = 1.5
            end

            p_bases = cell(1, length(fns));
            t_bases = cell(1, length(fns));
            suv_bases = cell(1, length(fns));
            for fidx = 1:length(fns)
                p_bases{fidx} = mlpmod.Pmod(filename=fns{fidx});
                t_bases{fidx} = p_bases{fidx}.tacs;
                suv_base_ = table2array(t_bases{fidx}(:, 4:end));
                suv_bases{fidx} = suv_base_ * p_bases{fidx}.weight / p_bases{fidx}.dose / 37;
            end

            % remove blocking
            is_blocking = contains(fns, "ses-20140724") | contains(fns, "ses-20160309");
            p_bases = p_bases(~is_blocking);
            t_bases = t_bases(~is_blocking);
            suv_bases = suv_bases(~is_blocking);            

            % sort by age
            ages = cellfun(@(x) x.age_acquisition, p_bases);  % numeric
            [sorted_ages, sorting_idx] = sort(ages);
            sorted_legend = sorted_ages + " y";
            p_bases = p_bases(sorting_idx);
            t_bases = t_bases(sorting_idx);
            suv_bases = suv_bases(sorting_idx);

            % cmap = cbrewer2('qual', 'Set2', 3, 'pchip');
            cmap = viridis(length(ages));
            % region_names = string(t_base.Properties.VariableNames(4:end));
            % region_names = strrep(region_names, "tac_", "");
            % region_names = strrep(region_names, "_kBq_cc_", "");
            region_names = [ ...                
                "whole brain", ...
                "prefrontal", ...
                "base frontal", ...
                "anterior cingulate", ...
                "caudate", ...
                "putamen", ...
                "insula", ...
                "lateral temporal", ...
                "medial temporal", ...
                "amygdala", ...
                "hippocampus", ...
                "thalamus", ...
                "hypothalamus", ...
                "parietal", ...
                "posterior cingulate", ...
                "occipital", ...
                "pons", ...
                "midbrain", ...
                "corpus callosum", ...
                "cortical white matter", ...
                "cerebellum"];

            suv_timesMid = t_bases{1}{:, 2} / 60;  % min
            times = 0:1:suv_timesMid(end);  % + suv_taus(end)/2);
            the_max = max( ...
                cellfun(@(x) max(x, [], "all"), suv_bases), [], "all");


            %% whole brain

            h = figure;   
            r = 1;

            legend_ = sorted_legend;
            for fidx = 1:length(fns)
                try
                    suv_base = suv_bases{fidx};
                    if any(isnan(suv_base(:, r)), "all")
                        legend_(fidx) = [];
                        continue
                    end
                    suv_timesMid = t_bases{fidx}{:, 2} / 60;  % min
                    interp_base = interp1(suv_timesMid, suv_base(:, r), times, "spline");
                    interp_base(interp_base < 0 ) = 0;

                    hold("on")
                    plot(times, interp_base, LineStyle="-", LineWidth=opts.lwidth, Color=cmap(fidx, :));
                    hold("off")
                catch ME
                    handwarning(ME)
                end
            end
            legend(legend_);
            xlabel("Time (min)");
            ylabel("SUV (g/cm^3)");
            ylim([0, the_max]);
            title(region_names(r));

            fontsize(scale=opts.fscale);
            saveFigure2(h, "whole_brain_age")


            %% cerebellum

            h = figure;   
            r = 21;

            legend_ = sorted_legend;
            for fidx = 1:length(fns)
                try
                    suv_base = suv_bases{fidx};
                    if any(isnan(suv_base(:, r)), "all")
                        legend_(fidx) = [];
                        continue
                    end
                    suv_timesMid = t_bases{fidx}{:, 2} / 60;  % min
                    interp_base = interp1(suv_timesMid, suv_base(:, r), times, "spline");
                    interp_base(interp_base < 0 ) = 0;

                    hold("on")
                    plot(times, interp_base, LineStyle="-", LineWidth=opts.lwidth, Color=cmap(fidx, :));
                    hold("off")
                catch ME
                    handwarning(ME)
                end
            end
            legend(legend_);
            xlabel("Time (min)");
            ylabel("SUV (g/cm^3)");
            ylim([0, the_max]);
            title(region_names(r));

            fontsize(scale=opts.fscale);
            saveFigure2(h, "cerebellum_age")


            %% regions

            htiled = figure;            
            tiledlayout(4, 5)
            for r = 2:21
                nexttile

                legend_ = sorted_legend;
                for fidx = 1:length(fns)
                    try
                        suv_base = suv_bases{fidx};
                        if any(isnan(suv_base(:, r)), "all")
                            legend_(fidx) = [];
                            continue
                        end
                        suv_timesMid = t_bases{fidx}{:, 2} / 60;  % min
                        interp_base = interp1(suv_timesMid, suv_base(:, r), times, "spline");
                        interp_base(interp_base < 0 ) = 0;

                        hold("on")
                        plot(times, interp_base, LineStyle="-", LineWidth=opts.lwidth, Color=cmap(fidx, :));
                        hold("off")
                    catch ME
                        handwarning(ME)
                    end
                end
                % legend(sorted_legend);
                if r == 2 || r == 7 || r == 12 || r == 17 
                    ylabel("SUV (g/cm^3)");
                else
                    set(gca, yticklabels=[])
                end
                if r == 17 || r == 18 || r == 19 || r == 20 || r == 21
                    xlabel("Time (min)");
                else
                    set(gca, xticklabels=[])
                end
                ylim([0, the_max]);
                title(region_names(r));
            end
            fontsize(scale=opts.fscale);
            saveFigure2(htiled, "regions_age")
        end

        function [htiled,h] = plot_blocking(fn_base, fn_block, fn_block2, opts)
            arguments
                fn_base {mustBeFile}
                fn_block {mustBeFile}
                fn_block2 {mustBeFile}
                opts.region_key {mustBeTextScalar} = "cerebellum"
                opts.legend {mustBeText} = ["", "baseline", "", "blocked by YUN122", "", "blocked by SA4503"]
                opts.lwidth double = 4
                opts.mark {mustBeTextScalar} = "+"
                opts.msize double = 12
                opts.fscale double = 1.5
            end

            p_base = mlpmod.Pmod(filename=fn_base);
            t_base = p_base.tacs;
            suv_base = table2array(t_base(:, 4:end));
            suv_base = suv_base * p_base.weight / p_base.dose / 37;

            p_block = mlpmod.Pmod(filename=fn_block);
            t_block = p_block.tacs;
            suv_block = table2array(t_block(:, 4:end));
            suv_block = suv_block * p_block.dose / p_block.weight / 37;

            p_block2 = mlpmod.Pmod(filename=fn_block2);
            t_block2 = p_block2.tacs;
            suv_block2 = table2array(t_block2(:, 4:end));
            suv_block2 = suv_block2 * p_block2.weight / p_block2.dose / 37;

            % cmap = cbrewer2('qual', 'Set2', 3, 'pchip');
            cmap = cividis(3);
            % region_names = string(t_base.Properties.VariableNames(4:end));
            % region_names = strrep(region_names, "tac_", "");
            % region_names = strrep(region_names, "_kBq_cc_", "");
            % region_index = contains(region_names, opts.region_key, IgnoreCase=true);
            region_names = [ ...                
                "whole brain", ...
                "prefrontal", ...
                "base frontal", ...
                "anterior cingulate", ...
                "caudate", ...
                "putamen", ...
                "insula", ...
                "lateral temporal", ...
                "medial temporal", ...
                "amygdala", ...
                "hippocampus", ...
                "thalamus", ...
                "hypothalamus", ...
                "parietal", ...
                "posterior cingulate", ...
                "occipital", ...
                "pons", ...
                "midbrain", ...
                "corpus callosum", ...
                "cortical white matter", ...
                "cerebellum"];

            % suv_diff = suv_block - suv_base;
            suv_taus = t_base{:, 3};
            suv_timesMid = t_base{:, 2} / 60;  % min
            times = 0:1:suv_timesMid(end);  % + suv_taus(end)/2);
            the_max = max([suv_base, suv_block, suv_block2], [], "all");


            %% whole brain

            h = figure;   
            r = 1;

            interp_base = interp1(suv_timesMid, suv_base(:, r), times, "spline");
            interp_base(interp_base < 0 ) = 0;
            interp_block = interp1(suv_timesMid, suv_block(:, r), times, "spline");
            interp_block(interp_block < 0 ) = 0;
            interp_block2 = interp1(suv_timesMid, suv_block2(:, r), times, "spline");
            interp_block2(interp_block2 < 0 ) = 0;

            hold("on")
            plot(suv_timesMid, suv_base(:, r), LineStyle="none", Marker=opts.mark, MarkerSize=opts.msize, Color=cmap(1, :));
            plot(times, interp_base, LineStyle="-", LineWidth=opts.lwidth, Color=cmap(1, :));
            plot(suv_timesMid, suv_block(:, r), LineStyle="none", Marker=opts.mark, MarkerSize=opts.msize, Color=cmap(2, :));
            plot(times, interp_block, LineStyle="-", LineWidth=opts.lwidth, Color=cmap(2, :));
            plot(suv_timesMid, suv_block2(:, r), LineStyle="none", Marker=opts.mark, MarkerSize=opts.msize, Color=cmap(3, :));
            plot(times, interp_block2, LineStyle="-", LineWidth=opts.lwidth, Color=cmap(3, :));
            hold("off")
            legend(opts.legend);
            xlabel("Time (min)");
            ylabel("SUV (g/cm^3)");
            ylim([0, the_max]);
            title(region_names(r));

            fontsize(scale=opts.fscale);
            saveFigure2(h, "whole_brain_including_blocking")


            %% regions

            htiled = figure;            
            tiledlayout(4, 5)
            for r = 2:21
                nexttile

                interp_base = interp1(suv_timesMid, suv_base(:, r), times, "spline");
                interp_base(interp_base < 0 ) = 0;
                interp_block = interp1(suv_timesMid, suv_block(:, r), times, "spline");
                interp_block(interp_block < 0 ) = 0;
                interp_block2 = interp1(suv_timesMid, suv_block2(:, r), times, "spline");
                interp_block2(interp_block2 < 0 ) = 0;

                hold("on")
                plot(suv_timesMid, suv_base(:, r), LineStyle="none", Marker=opts.mark, MarkerSize=opts.msize, Color=cmap(1, :));
                plot(times, interp_base, LineStyle="-", LineWidth=opts.lwidth, Color=cmap(1, :));
                plot(suv_timesMid, suv_block(:, r), LineStyle="none", Marker=opts.mark, MarkerSize=opts.msize, Color=cmap(2, :));
                plot(times, interp_block, LineStyle="-", LineWidth=opts.lwidth, Color=cmap(2, :));
                plot(suv_timesMid, suv_block2(:, r), LineStyle="none", Marker=opts.mark, MarkerSize=opts.msize, Color=cmap(3, :));
                plot(times, interp_block2, LineStyle="-", LineWidth=opts.lwidth, Color=cmap(3, :));
                hold("off")
                % legend(opts.legend);
                if r == 2 || r == 7 || r == 12 || r == 17 
                    ylabel("SUV (g/cm^3)");
                else
                    set(gca, yticklabels=[])
                end
                if r == 17 || r == 18 || r == 19 || r == 20 || r == 21
                    xlabel("Time (min)");
                else
                    set(gca, xticklabels=[])
                end
                ylim([0, the_max]);
                title(region_names(r));
            end
            fontsize(scale=opts.fscale);
            saveFigure2(htiled, "regions_including_blocking")
        end

        function s = rename_animals(s)
            given_names = ["bud", "cheech", "lou", "ollie"];
            anon_names = "nhprimate" + (1:4);
            for idx = 1:length(given_names)
                s = strrep(s, given_names(idx), anon_names(idx));
            end
        end
    end

    %% PRIVATE

    properties (Access = private)
        kmdata_
        pmod_
        solver_
    end

    methods (Access = private)
        function this = TZ3108()
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
