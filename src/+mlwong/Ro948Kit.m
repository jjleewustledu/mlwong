classdef Ro948Kit < handle & mlsystem.IHandle
    %% Has resemblance to a builder design pattern.
    %  All time activity curves are decay-corrected to the time of drawing blood.
    %  
    %  Created 06-Jun-2023 18:35:41 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlwong/src/+mlwong.
    %  Developed on Matlab 9.14.0.2254940 (R2023a) Update 2 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        Delta
        do_Hill
        half_life
        hct
        Hill_d_min
        N_average_minimal
        ptac_units
        rescaling_twilite
        sesdir
        t0_forced
        timeCliff
        toi
        tracer

        color_ocean = ""
        color_scarlet = ""
        color_orange = ""
        color_grape = ""
        color_green = ""
        color_sky = ""
        color_maroon = ""
    end

    properties (Constant)
        do_decay_correct = false  % Rick does external decay correction
    end

    properties (Dependent)
        chemdir
        crv
        crv_phant
        fileprefix
        inveff % counts/s to Bq/mL
        ses
        sesid
        sub
        subid
        twildir

        T_frac_intact
        T_hill
        T_metab_corr_ptac
        T_minimal
        T_parent_frac
        T_pow
        T_total_ptac
        T_total_ptac_no_twil
        T_total_ptac_twil
        T_twil
    end

    methods %% GET,SET
        function g = get.chemdir(this)
            g = fullfile(this.sesdir, "chemistry");
        end
        function g = get.crv(this)
            g = this.crv_;
        end
        function g = get.crv_phant(this)
            g = this.crv_phant_;
        end
        function g = get.fileprefix(this)
            if ~isemptytext(this.fileprefix_)
                g = this.fileprefix_;
                return
            end
            this.fileprefix_ = sprintf("%s_%s", this.subid, this.sesid);
            g = this.fileprefix_;
        end
        function g = get.inveff(this)
            g = 68/0.27; 
        end
        function g = get.ses(this)
            ss = pathsplit(this.sesdir, filesep);
            g = ss{end};
        end
        function g = get.sesid(this)
            re = regexp(this.ses, "ses-(?<sesid>\d+)", "names");
            assert(~isempty(re))
            g = re.sesid;
        end
        function g = get.sub(this)
            ss = pathsplit(this.sesdir, filesep);
            g = ss{end-1};
        end
        function g = get.subid(this)
            re = regexp(this.ses, "sub-(?<subid>\S+)", "names");
            assert(~isempty(re))
            g = re.subid;
        end
        function g = get.twildir(this)
            g = fullfile(this.sesdir, "twilite");
            ensuredir(g);
        end

        function g = get.T_frac_intact(~)
            g = readtable(this.fileprefix+"_frac_intact.csv");
        end
        function g = get.T_hill(this)
            if ~isempty(this.T_hill_)
                g = this.T_hill_;
                return
            end

            T_hill_csv = fullfile(this.chemdir, this.fileprefix+"_Hill.csv");
            if isfile(T_hill_csv)
                this.T_hill_ = readtable(T_hill_csv);
                g = this.T_hill_;
                return
            end

            [~,this.T_hill_] = this.build_Hill();
            g = this.T_hill_;
        end
        function g = get.T_metab_corr_ptac(this)
            %% plasma time activity curves, w/ metab corrections

            if ~isempty(this.T_metab_corr_ptac_)
                g = this.T_metab_corr_ptac_;
                return
            end

            Time = this.T_total_ptac.Sec;
            PF = makima(this.T_parent_frac.Sec, this.T_parent_frac.PF, Time);
            mc_pTAC = this.T_total_ptac.pTAC .* PF/100;

            DateTime = this.toi + seconds(Time);
            Min = Time/60;
            Sec = Time; 
            g = table(DateTime, Min, Sec, mc_pTAC);
            this.T_metab_corr_ptac_ = g;
        end
        function g = get.T_minimal(this)
            if ~isempty(this.T_minimal_)
                g = this.T_minimal_;
                return
            end
            T_csv = fullfile(this.chemdir, this.fileprefix+"_minimal.csv");
            assert(isfile(T_csv))
            T = readtable(T_csv);
            assert(length(T.Properties.VariableNames) == 4, stackstr())
            T.Properties.VariableNames = ["drawStart", "drawFinish", "wbKBq_mL", "plasmaKBq_mL"];
            if isnumeric(T.drawStart) || isnumeric(T.drawFinish)
                if T.drawFinish(end) > 600
                    T.drawStart = this.toi + seconds(T.drawStart);
                    T.drawFinish = this.toi + seconds(T.drawFinish);
                else
                    T.drawStart = this.toi + minutes(T.drawStart);
                    T.drawFinish = this.toi + minutes(T.drawFinish);
                end
            end
            T.drawStart.TimeZone = "local";
            T.drawFinish.TimeZone = "local";
            draw = duration(T.drawFinish - T.drawStart)/2 + T.drawStart;
            %doi = datetime(this.toi.Year, this.toi.Month, this.toi.Day, 0,0,0);
            Time = seconds(draw - this.toi);
            this.T_minimal_ = addvars(T, Time, Before=1, NewVariableNames="Time");

            % decay-correct table minimal from syringe measurements
            if this.do_decay_correct
                this.T_minimal_.wbKBq_mL = this.T_minimal_.wbKBq_mL .* 2.^(Time/this.half_life);
                this.T_minimal_.plasmaKBq_mL = this.T_minimal_.plasmaKBq_mL .* 2.^(Time/this.half_life);
            end

            g = this.T_minimal_;
        end
        function g = get.T_parent_frac(this)
            %% parent fraction from extended Hill function

            if ~isempty(this.T_parent_frac_)
                g = this.T_parent_frac_;
                return
            end

            DateTime = this.toi + seconds(this.T_hill.Time);
            Min = this.T_hill.Time/60;
            Sec = this.T_hill.Time; 
            PF = this.T_hill.FractionIntact*100; % percent
            g = table(DateTime, Min, Sec, PF);
            this.T_parent_frac_ = g;
        end
        function g = get.T_pow(this)
            if ~isempty(this.T_pow_)
                g = this.T_pow_;
                return
            end
            T_pow_csv = fullfile(this.chemdir, this.fileprefix+"_pow.csv");
            if isfile(T_pow_csv)
                this.T_pow_ = readtable(T_pow_csv);
                g = this.T_pow_;
                return
            end
            [~,this.T_pow_] = this.build_pow();
            g = this.T_pow_;
        end
        function g = get.T_total_ptac(this)
            %% plasma time activity curves, w/o metab corrections

            if ~isempty(this.T_total_ptac_)
                g = this.T_total_ptac_;
                return
            end

            if ~isempty(this.crv)
                g = this.T_total_ptac_twil;
            else
                g = this.T_total_ptac_no_twil;
            end
        end
        function g = get.T_total_ptac_twil(this)
            %% plasma time activity curves, w/o metab corrections

            Time = [this.T_twil.Time(1:end-1); this.T_minimal.Time];
            N_twil = length(this.T_twil.wbKBq_mL);
            twil_plasmaKBq_mL = this.T_twil.wbKBq_mL(1:end-1).*this.T_pow.PlasmaOverWb(1:N_twil-1);
            pTAC = [twil_plasmaKBq_mL; this.T_minimal.plasmaKBq_mL];

            DateTime = this.toi + seconds(Time);
            Min = Time/60;
            Sec = Time; 
            pTAC = this.convert_ptac_units(pTAC);
            g = table(DateTime, Min, Sec, pTAC);
            this.T_total_ptac_ = g;
        end
        function g = get.T_total_ptac_no_twil(this)
            %% plasma time activity curves, w/o metab corrections

            Time = this.T_minimal.Time;
            DateTime = this.toi + seconds(Time);
            Min = Time/60;
            Sec = Time; 
            pTAC = this.T_minimal.plasmaKBq_mL;
            pTAC = this.convert_ptac_units(pTAC); 
            g = table(DateTime, Min, Sec, pTAC);
            this.T_total_ptac_ = g;
        end
        function g = get.T_twil(this)
            if isempty(this.crv)
                g = [];
                return
            end

            [~,fp] = myfileparts(this.crv.filename);
            deconv_csv = fullfile(this.twildir, fp+"_deconv.csv");
            T = readtable(deconv_csv); 

            % rescale twilite
            first_T_minimal = this.T_minimal.wbKBq_mL(1:this.N_average_minimal);
            first_T_minimal = mean(first_T_minimal);
            last_T_wbKBq_mL = T.wbKBq_mL(T.wbKBq_mL > 0);
            last_T_wbKBq_mL = last_T_wbKBq_mL(end);
            rescaling = first_T_minimal/last_T_wbKBq_mL;
            assert(isfinite(rescaling), stackstr())
            T.wbKBq_mL = T.wbKBq_mL * rescaling;
            fprintf("%s.rescaling->%g\n", stackstr(), rescaling)
            this.T_twil_ = T;
            g = this.T_twil_;
            this.rescaling_twilite = rescaling;
        end
    end

    methods

        function [h,T] = build_deconv(this, hct, t0_forced, opts)
            arguments
                this mlwong.Ro948Kit
                hct {mustBeScalarOrEmpty} = this.hct
                t0_forced {mustBeScalarOrEmpty} = this.t0_forced
                opts.timeCliff double = this.timeCliff
            end

            pwd0 = pushd(this.twildir);
            
            idx0 = this.index_toi_minus_30sec(); % TOI - 30 s
            idxF = this.index_toi_plus_time_cliff_min(); % nominal TOI + 5 min blood draw 
            toi_ = this.toi;
            select = toi_ <= this.crv.time & this.crv.time <= toi_ + seconds(this.timeCliff); % 5 min + 1, col

            % deconvBayes
            M_ = this.crv.timetable().Coincidence(idx0:idxF)*this.inveff; % col
            cath = mlswisstrace.Catheter_DT20190930( ...
                Measurement=M_, hct=hct, tracer=this.tracer); 
            M = zeros(size(this.crv.timetable().Coincidence)); % col
            % t0 reflects rigid extension + Luer valve + cath in Twilite cradle
            xlim = [0 this.timeCliff+t0_forced];
            M(idx0:idxF) = ascol( ...
                cath.deconvBayes( ...
                t0_forced=t0_forced, xlim=xlim));
            fqfp = fullfile(this.twildir, this.fileprefix+"_deconvBayes");
            %if ~isfile(fqfp+".fig")
                saveFigure2(gcf, fqfp);
            %end

            % disperse Measurement
            M_select = this.disperse(M(select));

            % write csv
            Time = seconds(this.crv.time(select) - toi_);
            wbKBq_mL = ascol(M_select)/1e3; % cps -> kcps
            T = table(Time, wbKBq_mL);            
            [~,fp] = myfileparts(this.crv.filename);
            writetable(T, fullfile(this.twildir, fp+"_deconv.csv"));

            % plot wb specific activity alone
            h = figure;
            plot(T, "Time", "wbKBq_mL", LineStyle="none", Marker="+", MarkerSize=6)
            xlabel("time (s)", FontSize=14, FontWeight='bold')
            ylabel("wb activity (kcps)", FontSize=14, FontWeight='bold')
            title(sprintf("From %s", this.crv.filename), FontSize=16, Interpreter="none")   
            fqfp = fullfile(this.twildir, fp+"_deconv");
            %if ~isfile(fqfp+".fig")
                saveFigure2(h, fqfp)
            %end

            popd(pwd0);
        end


        function [h,T] = build_deconv1(this, hct, t0_forced, opts)
            arguments
                this mlwong.Ro948Kit
                hct {mustBeScalarOrEmpty} = this.hct
                t0_forced {mustBeScalarOrEmpty} = this.t0_forced
                opts.timeCliff double = this.timeCliff
            end

            pwd0 = pushd(this.twildir);

            popd(pwd0);
        end

        function [h,V] = build_pow(this, opts)
            %% build plasma over whole-blood fractions

            arguments
                this mlwong.Ro948Kit
                opts.visible_interval double = 30 % visible from start of syringe samples (min)
            end

            pwd0 = pushd(this.chemdir);

            T = readtable(this.fileprefix+"_minimal.csv");
            T.Properties.VariableNames = ["drawStart", "drawFinish", "wbKBq_mL", "plasmaKBq_mL"];
            if isnumeric(T.drawStart) || isnumeric(T.drawFinish)
                drawStart_ = this.toi + minutes(T.drawStart);
                drawFinish_ = this.toi + minutes(T.drawFinish);
                T.drawStart = drawStart_;
                T.drawFinish = drawFinish_;
            end

            T.drawFinish.TimeZone = "local";
            T.drawStart.TimeZone = "local";
            draw = duration(T.drawFinish - T.drawStart)/2 + T.drawStart;
            %date_toi = datetime(this.toi.Year, this.toi.Month, this.toi.Day, 0,0,0);
            Time = seconds(draw - this.toi);
            T = addvars(T, Time, Before=1, NewVariableNames="Time");
            PlasmaOverWb = T.plasmaKBq_mL ./ T.wbKBq_mL;
            T = addvars(T, PlasmaOverWb, NewVariableNames="PlasmaOverWb");
            U = T(T.Time < opts.visible_interval*60, :);
            
            % regresion for early times
            mdl = fitlm(U, "PlasmaOverWb~Time");
            disp(mdl)

            Time1 = ascol(0:opts.visible_interval*60);
            Pow1 = feval(mdl, Time1);
            h = figure;
            hold('on')
            plot(T, "Time", "PlasmaOverWb", LineStyle="none", Marker="o", MarkerSize=9)
            plot(Time1, Pow1, LineWidth=2)
            xlabel("time (s)", FontSize=14, FontWeight="bold")
            ylabel("plasma activity / wb activity", FontSize=14, FontWeight="bold")
            title(sprintf("From "+this.fileprefix+"_minimal.csv"), FontSize=16, Interpreter="none")            
            saveFigure2(h, fullfile(this.chemdir, this.fileprefix + "_pow"))
            
            V = table(Time1, Pow1, VariableNames=["Time", "PlasmaOverWb"]);
            writetable(V, this.fileprefix+"_pow.csv")

            popd(pwd0);
        end        
        function [h,T1] = build_Hill(this)
            pwd0 = pushd(this.chemdir);

            T = readtable(this.fileprefix+"_frac_intact.csv");
            T.Properties.VariableNames = ["Time", "FractionIntact"];
            if max(T.Time) <= 190 
                % ensure sec
                T.Time = T.Time*60; 
            end
            if max(T.FractionIntact) > 1
                T.FractionIntact = 0.01*T.FractionIntact;
            end
            if this.do_Hill
                hill = mlwong.Hill(T, dt=1, d_min=this.Hill_d_min, visualize_anneal=false);
                T1 = hill.solve();
                h = plot(hill);
                disp(hill)
                disp(hill.results.ks)
                disp(hill.results.loss)
            else
                timeInterp = 0:1:T.Time(end)-1;
                FractionIntact = makima(T.Time, T.FractionIntact, timeInterp);
                T1 = table(ascol(timeInterp), ascol(FractionIntact), variableNames=["Time", "FractionIntact"]);
                h = plot(T1, "Time", "FractionIntact");
                xlabel("Time");
                ylabel("Fraction Intact");
                title("Fraction Intact by Makima (replacing Hill)");
            end
            saveFigure2(h, this.fileprefix+"_Hill");
            writetable(T1, this.fileprefix+"_Hill.csv")

            popd(pwd0);
        end
        
        function h = build_ptacs(this)
            if ~isempty(this.crv)
                h = this.build_ptacs_twil();
            else
                h = this.build_ptacs_no_twil();
            end
        end
        function h = build_ptacs_no_twil(this)
            %% total ptac & m.c. ptac

            ptac_2 = this.T_total_ptac;
            mc_ptac_2 = this.T_metab_corr_ptac;

            h = figure;
            hold("on")
            plot(ptac_2, "Min", "pTAC", LineStyle="none", Marker="x", MarkerEdgeColor="#7E2F8E", MarkerSize=10)
            plot(mc_ptac_2, "Min", "mc_pTAC", LineStyle="none", Marker="o", MarkerEdgeColor="#0072BD", MarkerSize=10)
            hold("off")
            ylim([0 this.convert_ptac_units(18.5)])
            xlabel("Time (min)", FontSize=14, FontWeight="bold")
            ylabel("Plasma time activity curves "+this.ptac_units, FontSize=14, FontWeight="bold")
            title("Plasma time activity curves", FontSize=18, Interpreter="none")    
            legend(["Syringe total", "Syringe metab. corr."]) 
            saveFigure2(h, fullfile(this.chemdir, this.fileprefix+"_ptacs_no_twil"));
        end
        function h = build_ptacs_twil(this)
            %% total ptac & m.c. ptac

            Nt = length(this.T_twil.Time) - 1;
            ptac_1 = this.T_total_ptac(1:Nt,:);
            ptac_2 = this.T_total_ptac(Nt+1:end,:);
            mc_ptac_1 = this.T_metab_corr_ptac(1:Nt,:);
            mc_ptac_2 = this.T_metab_corr_ptac(Nt+1:end,:);

            h = figure;
            hold("on")
            plot(ptac_1, "Min", "pTAC", LineStyle="--", Color="#7E2F8E", LineWidth=3, Marker="none")
            plot(ptac_2, "Min", "pTAC", LineStyle="none", Marker="x", MarkerEdgeColor="#7E2F8E", MarkerSize=10)
            plot(mc_ptac_1, "Min", "mc_pTAC", LineStyle="-", Color="#0072BD", LineWidth=2, Marker="none")
            plot(mc_ptac_2, "Min", "mc_pTAC", LineStyle="none", Marker="o", MarkerEdgeColor="#0072BD", MarkerSize=10)
            hold("off")
            ylim([0 this.convert_ptac_units(18.5)])
            xlabel("Time (min)", FontSize=14, FontWeight="bold")
            ylabel("Plasma time activity curves "+this.ptac_units, FontSize=14, FontWeight="bold")
            title("Plasma time activity curves", FontSize=18, Interpreter="none")    
            legend(["Twilite total", "Syringe total", "Twilite metab. corr.", "Syringe metab. corr."]) 
            saveFigure2(h, fullfile(this.chemdir, this.fileprefix+"_ptacs_twil"));
        end
        function call(this)
            if ~isempty(this.crv)
                this.build_deconv();
                this.build_pow();
            end
            this.build_Hill();
            this.build_ptacs();
            this.plot_crvs();

            %this.T_total_ptac % disp head & tail

            writetable(this.T_total_ptac, ...
                fullfile(this.chemdir, this.fileprefix + "_total_ptac.csv"));
            writetable(this.T_metab_corr_ptac, ...
                fullfile(this.chemdir, this.fileprefix + "_metab_corr_ptac.csv"));
            writetable(this.T_parent_frac, ...
                fullfile(this.chemdir, this.fileprefix + "_parent_frac.csv"));
        end
        function M1 = disperse(this, M)
            %% disperse Measurement, which has been trimmed to interval of measurement

            if ~isrow(M)
                M = asrow(M);
                times = 0:length(M)-1;
                AUC = trapz(exp(-times*this.Delta));
                M__ = conv(M, exp(-times*this.Delta))/AUC;
                M1 = M__(1:length(M));
                M1 = ascol(M1);
                return
            end
    
            % M is row
            times = 0:length(M)-1;
            AUC = trapz(exp(-times*this.Delta));
            M__ = conv(M, exp(-times*this.Delta))/AUC;
            M1 = M__(1:length(M));
        end
        function [h,h1,h2,h3] = plot(this)

            h = this.plot_crvs();
            [~,fp] = myfileparts(this.crv_.filename);
            fqfp = fullfile(this.twildir, fp);
            %if ~isfile(fqfp+".fig")
                saveFigure2(h, fqfp)
            %end

            T = readtable(fullfile( ...
                this.chemdir, this.fileprefix + "_parent_frac.csv"));
            h1 = figure;
            plot(T, "Min", "PF");
            xlabel("Time (Min)")
            ylabel("Parent Fraction (%)")
            title(this.fileprefix + "_parent_frac.csv", Interpreter="none")
            saveFigure2(h1, fullfile( ...
                this.chemdir, this.fileprefix + "_parent_frac"))

            T = readtable(fullfile( ...
                this.chemdir, this.fileprefix + "_total_ptac.csv"));
            h2 = figure;
            plot(T, "Min", "pTAC");
            xlabel("Time (Min)")
            ylabel("Activity "+this.ptac_units)
            title(this.fileprefix + "_total_ptac.csv", Interpreter="none")
            saveFigure2(h2, fullfile( ...
                this.chemdir, this.fileprefix + "_total_ptac"))

            T = readtable(fullfile( ...
                this.chemdir, this.fileprefix + "_metab_corr_ptac.csv"));
            h3 = figure;
            plot(T, "Min", "mc_pTAC");
            xlabel("Time (Min)")
            ylabel("Activity "+this.ptac_units)
            title(this.fileprefix + "_metab_corr_ptac.csv", Interpreter="none")
            saveFigure2(h3, fullfile( ...
                this.chemdir, this.fileprefix + "_metab_corr_ptac"))
        end
        function h = plot_crvs(this)
            h = figure;
            this.crv_.plotAll()
            x = [0.25, 0.25];
            y = [0.5, 0.65];
            x1 = [0.7, 0.8];
            y1 = [0.5, 0.65];
            annotation('textarrow', x, y, String='time of injection')
            annotation('textarrow', x1, y1, String='disconnetion of Twilite from radial artery')
            saveFigure2(h, fullfile(this.chemdir, this.fileprefix + "_crvs"))

            if ~isempty(this.crv_phant_)
                figure
                this.crv_phant_.plotAll()
            end
        end
        
        function T = frac_intact(this)
            %% parent fraction from HPLC

            DateTime = this.toi + seconds(this.T_frac_intact.Time);
            Min = this.T_frac_intact.Time;
            Sec = this.T_frac_intact.Time*60; 
            PF = this.T_frac_intact.FractionIntact*100; % percent
            T = table(DateTime, Min, Sec, PF);
            this.T_frac_intact_ = T;
        end
                
        function this = Ro948Kit(opts)
            arguments
                opts.toi datetime = NaT
                opts.sesdir {mustBeFolder} = pwd
                opts.crv = []
                opts.crv_phant = []
                opts.fileprefix {mustBeTextScalar} = "" % preferred
                opts.hct {mustBeScalarOrEmpty} = 43.75;
                opts.t0_forced {mustBeScalarOrEmpty} = 40;
                opts.Hill_d_min double = 0.95
                opts.do_Hill = true
                opts.ptac_units = "uCi/mL" % "uCi/mL", "kBq/mL"
                opts.half_life = 109.77
                opts.tracer = "RO948"
                opts.timeCliff double = 300
                opts.N_average_minimal double = 1
                opts.Delta double = 1/60
            end
            this.toi = opts.toi;
            this.sesdir = opts.sesdir;
            this.crv_ = opts.crv;
            this.crv_phant_ = opts.crv_phant;
            this.fileprefix_ = opts.fileprefix;
            this.hct = opts.hct;
            this.t0_forced = opts.t0_forced;
            this.Hill_d_min = opts.Hill_d_min;
            this.do_Hill = opts.do_Hill;
            this.ptac_units = opts.ptac_units;
            this.half_life = opts.half_life;
            this.tracer = opts.tracer;
            this.timeCliff = opts.timeCliff;
            this.N_average_minimal = opts.N_average_minimal;
            this.Delta = opts.Delta;
        end
        function initialize(this)
            this.T_frac_intact = readtable(fullfile(this.chemdir, this.fileprefix + "_frac_intact.csv"));
            this.T_frac_intact.Properties.VariableNames = ["Time", "FractionIntact"];
            % scrub nans
            this.T_frac_intact = this.T_frac_intact(~isnan(this.T_frac_intact.FractionIntact), :);            
        end

        %% UTILITES

        function ptac = convert_ptac_units(this, ptac)
            switch convertStringsToChars(this.ptac_units)
                case 'kBq/mL'
                case {'uCi/mL','\mu Ci/mL','muCi/mL'}
                    ptac = ptac/37;
                otherwise
                    error("mlwong:ValueError", "%s: this.ptac_units->%s", stackstr(), this.ptac_units)
            end
        end
        function idx = index_toi_minus_30sec(this)
            dur = duration(this.toi - this.datetime_crv_init(this.crv)) - seconds(30);
            idx = round(seconds(dur) + 1);
        end
        function idx = index_toi_plus_time_cliff_min(this)
            dur = duration(this.toi - this.datetime_crv_init(this.crv)) + seconds(this.timeCliff);
            idx = round(seconds(dur) + 1);
        end
    end

    methods (Static)
        function dt = datetime_crv_init(crv)
            dt = crv.timetable().Time(1);
            dt.TimeZone = "local";
        end
        function hct = estimated_hct(opts)
            arguments
                opts.sex char = ''
            end
            switch lower(opts.sex)                
                case 'f'
                    hct = mean([36, 48]);
                case 'm'
                    hct = mean([41, 50]);
                otherwise 
                    hct = 43.75;
            end
        end  
    end

    %% PROTECTED

    properties (Access = protected)
        crv_
        crv_phant_
        fileprefix_

        T_frac_intact_
        T_hill_
        T_metab_corr_ptac_
        T_minimal_
        T_parent_frac_  
        T_pow_
        T_total_ptac_
        T_twil_      
    end

    methods (Access = protected)
        function that = copyElement(this)
            that = copyElement@matlab.mixin.Copyable(this);
            if ~isempty(this.crv_)
                that.crv_ = copy(this.crv_); end
            if ~isempty(this.crv_phant_)
                that.crv_phant_ = copy(this.crv_phant_); end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
