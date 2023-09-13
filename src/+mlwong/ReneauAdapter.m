classdef ReneauAdapter < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 14-Jul-2023 18:24:02 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlwong/src/+mlwong.
    %  Developed on Matlab 9.14.0.2286388 (R2023a) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        sesdir
        final_calc
        fractions
        fileprefix

        T_minimal
        T_frac_intact
    end

    methods
        function this = ReneauAdapter(opts)
            arguments
                opts.sesdir {mustBeFolder} = pwd
                opts.final_calc {mustBeFile} % Rick's csv file
                opts.fractions {mustBeFile} % Rick's csv file
                opts.fileprefix {mustBeTextScalar} % BIDs adherent fileprefix
            end

            this.sesdir = opts.sesdir;
            this.final_calc = opts.final_calc;
            this.fractions = opts.fractions;
            this.fileprefix = opts.fileprefix;
        end

        function call(this)
            %% reads Rick's tables, adjusts formats, writes tables understood by mlwong.Ro948Kit.

            this.readtables();
            this.writetables();
        end

        function readtables(this)

            T_m = readtable(fullfile(this.sesdir, this.final_calc));
            N = sum(contains(T_m.Barcode, 'PL'));
            drawStart = T_m.drawStart(1:N);
            drawFinish = T_m.drawFinish(1:N);
            assert(all(drawStart == T_m.drawStart(N+1:end)))
            assert(all(drawFinish == T_m.drawFinish(N+1:end)))
            wbKBq_mL = T_m.kBq(N+1:end);
            plasmaKBq_mL = T_m.kBq(1:N);
            this.T_minimal = table(drawStart, drawFinish, wbKBq_mL, plasmaKBq_mL);

            T_fi = readtable(fullfile(this.sesdir, this.fractions));
            select_min = contains(T_fi.Name, 'min');
            select_dose = contains(T_fi.Name, 'dose');
            min_ = strip(strrep(T_fi.Name(select_min), 'min', ''));
            Time_min_ = [0; ascol(cellfun(@str2num, min_, UniformOutput=false))];
            Time_min_ = cell2mat(Time_min_);
            Fraction_Intact = [T_fi.Percent_Intact(select_dose); ascol(T_fi.Percent_Intact(select_min))];            
            this.T_frac_intact = table(Time_min_, Fraction_Intact);
            this.T_frac_intact = sortrows(this.T_frac_intact, 1);
        end
        function writetables(this)
            writetable(this.T_minimal, fullfile(this.sesdir, 'chemistry', this.fileprefix+"_minimal.csv"));
            writetable(this.T_frac_intact, fullfile(this.sesdir, 'chemistry', this.fileprefix+"_frac_intact.csv"));
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
