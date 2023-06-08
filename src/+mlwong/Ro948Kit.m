classdef Ro948Kit < handle
    %% line1
    %  line2
    %  
    %  Created 06-Jun-2023 18:35:41 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlwong/src/+mlwong.
    %  Developed on Matlab 9.14.0.2254940 (R2023a) Update 2 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        fileprefix
        toi
        T_frac_intact
        T_hill
        T_minimal
        T_pow
        T_twil
    end

    methods
        function call(this)
            this.total_ptac();
            this.parent_frac();
            this.metab_corr_ptac();
            this.frac_intact();
            this.plot()
        end
        function [h1,h2,h3] = plot(this)

            %% total ptac & m.c. ptac
            ptac_twil = this.T_total_ptac_(1:300,:);
            ptac_syringe = this.T_total_ptac_(301:end,:);
            mc_ptac_twil = this.T_metab_corr_ptac_(1:300,:);
            mc_ptac_syringe = this.T_metab_corr_ptac_(301:end,:);

            h1 = figure;
            hold("on")
            plot(ptac_twil, "Min", "pTAC", LineStyle="--", LineWidth=2, Marker="none")
            plot(ptac_syringe, "Min", "pTAC", LineStyle="none", Marker="x", MarkerSize=9)
            plot(mc_ptac_twil, "Min", "mc_pTAC", LineStyle="-", LineWidth=2, Marker="none")
            plot(mc_ptac_syringe, "Min", "mc_pTAC", LineStyle="none", Marker="o", MarkerSize=9)
            hold("off")
            ylim([0 0.5]) % \mu Ci/mL
            xlabel("Time (min)", FontSize=14)
            ylabel("Plasma time activity curves ({\mu}Ci/mL)", FontSize=14)
            title("Plasma time activity curves", FontSize=18, Interpreter="none")    
            legend(["Twilite total", "Syringe total", "Twilite metab. corr.", "Syringe metab. corr."])
            savefig(h1, this.fileprefix+"_ptac.fig")

            %% total ptac & m.c. ptac, 0-299 sec
            h2 = figure;
            hold("on")
            plot(ptac_twil, "Sec", "pTAC", LineStyle="--", LineWidth=2, Marker="none")
            plot(mc_ptac_twil, "Sec", "mc_pTAC", LineStyle="-", LineWidth=2, Marker="none")
            hold("off")
            xlim([0 300])
            xlabel("Time (sec)", FontSize=14)
            ylabel("Plasma time activity curves ({\mu}Ci/mL)", FontSize=14)
            title("Plasma time activity curves (early times)", FontSize=18, Interpreter="none")    
            legend(["Twilite total", "Twilite metab. corr."])
            savefig(h1, this.fileprefix+"_ptac_earlytimes.fig")
            
            %% parent frac
            h3 = figure;
            hold("on")
            plot(this.T_frac_intact_, "Min", "PF", LineStyle="none", Marker="o", MarkerSize=9)
            plot(this.T_parent_frac_, "Min", "PF", LineStyle="--", LineWidth=2, Marker="none")
            hold("off")
            xlabel("Time (min)", FontSize=14)
            ylabel("Parent fraction (%)", FontSize=14)
            title("Parent fraction", FontSize=18, Interpreter="none")    
            legend(["HPLC measurement", "Hill function fit"])
            savefig(h3, this.fileprefix+"_parent_frac.fig")
        end
        function T = total_ptac(this)
            %% plasma time activity curves, w/o metab corrections

            Time = [this.T_twil.Time(1:end-1); this.T_minimal.Time];
            N_twil = length(this.T_twil.wbKBq_mL);
            twil_plasmaKBq_mL = this.T_twil.wbKBq_mL(1:end-1).*this.T_pow.PlasmaOverWb(1:N_twil-1);
            pTAC = [twil_plasmaKBq_mL; this.T_minimal.plasmaKBq_mL];

            DateTime = this.toi + seconds(Time);
            Min = Time/60;
            Sec = Time; 
            pTAC = pTAC/37; % \mu Ci/mL
            T = table(DateTime, Min, Sec, pTAC);
            this.T_total_ptac_ = T;
            writetable(T, sprintf("%s_total_ptac.csv", this.fileprefix))
        end
        function T = metab_corr_ptac(this)
            %% plasma time activity curves, w/ metab corrections

            Time = this.T_total_ptac_.Sec;
            PF = interp1(this.T_parent_frac_.Sec, this.T_parent_frac_.PF, Time);
            mc_pTAC = this.T_total_ptac_.pTAC .* PF/100;

            DateTime = this.toi + seconds(Time);
            Min = Time/60;
            Sec = Time; 
            T = table(DateTime, Min, Sec, mc_pTAC);
            this.T_metab_corr_ptac_ = T;
            writetable(T, sprintf("%s_metab_corr_ptac.csv", this.fileprefix))
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
        function T = parent_frac(this)
            %% parent fraction from extended Hill function

            DateTime = this.toi + seconds(this.T_hill.Time);
            Min = this.T_hill.Time/60;
            Sec = this.T_hill.Time; 
            PF = this.T_hill.FractionIntact*100; % percent
            T = table(DateTime, Min, Sec, PF);
            this.T_parent_frac_ = T;
            writetable(T, sprintf("%s_parent_frac.csv", this.fileprefix))
        end
        function this = Ro948Kit(varargin)

            twildir = fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01384/derivatives/sub-R21-004/ses-20230330/twilite');
            chemdir = fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01384/derivatives/sub-R21-004/ses-20230330/chemistry');

            this.toi = datetime(2023,3,30,11,37,1);
            this.fileprefix = "R21_004_03302023";
            this.T_frac_intact = readtable(fullfile(chemdir, "R21_004_03302023_frac_intact.csv"));
            this.T_frac_intact.Properties.VariableNames = ["Time", "FractionIntact"];
            this.T_hill = readtable(fullfile(chemdir, "R21_004_03302023_Hill.csv"));
            this.T_minimal = readtable(fullfile(chemdir, "R21_004_03302023_minimal.csv"));
            this.T_pow = readtable(fullfile(chemdir, "R21_004_03302023_pow.csv"));
            this.T_twil = readtable(fullfile(twildir, "R21_004_03302023_D1_deconv.csv"));

            % scrub nans
            this.T_frac_intact = this.T_frac_intact(~isnan(this.T_frac_intact.FractionIntact), :);
            % augment minimal
            draw = (this.T_minimal.drawFinish - this.T_minimal.drawStart)/2 + this.T_minimal.drawStart;
            draw.Year = this.toi.Year;
            draw.Month = this.toi.Month;
            draw.Day = this.toi.Day;
            Time = seconds(draw - this.toi);
            this.T_minimal = addvars(this.T_minimal, Time, Before=1, NewVariableNames="Time");

            % rescale twilite
            rescaling = this.T_minimal.wbKBq_mL(1)/this.T_twil.wbKBq_mL(end);
            this.T_twil.wbKBq_mL = this.T_twil.wbKBq_mL * rescaling;
            disp(rescaling)
        end
    end

    %% PRIVATE

    properties (Access = private)
        T_total_ptac_
        T_metab_corr_ptac_
        T_frac_intact_
        T_parent_frac_
        
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
