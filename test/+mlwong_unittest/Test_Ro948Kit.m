classdef Test_Ro948Kit < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 06-Jun-2023 18:35:41 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlwong/test/+mlwong_unittest.
    %  Developed on Matlab 9.14.0.2254940 (R2023a) Update 2 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        testObj
        reneau
    end
    
    methods (Test)
        function test_dynesty(this)
            pwd0 = pushd(this.sesdir_);

            opts.final_calc = "chemistry/2023-06-29_Final_Calculation.csv";
            opts.fractions = "chemistry/2023-06-29_RO948_Metabolite_Fractions.csv";
            opts.fileprefix = stackstr();
            opts.toi = "29-Jun-2023 13:10:23";
            opts.crv_filename = "twilite/R21_016_PET2_06292023_D1.crv";
            opts.hct = "41.5";
            opts.t0_forced = "47";
            opts.time_cliff = "300";
            

            opts.toi = datetime(opts.toi, TimeZone = "local");
            assert(contains(opts.crv_filename, ".crv"))
            crv_ = mlswisstrace.CrvData(opts.crv_filename);

            reneau = mlwong.ReneauAdapter( ...
                sesdir=this.sesdir_, ...
                final_calc=opts.final_calc, ...
                fractions=opts.fractions, ...
                fileprefix=opts.fileprefix);
            reneau.readtables();
            reneau.writetables();

            this = mlwong.Ro948Kit( ...
                toi=opts.toi, ...
                crv=crv_, ...
                fileprefix=opts.fileprefix, ...
                hct = str2double(opts.hct), ...
                t0_forced = str2double(opts.t0_forced), ...
                timeCliff = str2double(opts.time_cliff));
            %disp(this)
            call_before_dynesty(this)

            popd(pwd0)
        end

        function test_standalone(this)
            pwd0 = pushd(this.sesdir_);
            cnami_standalone( ...
                final_calc="chemistry/2023-06-29_Final_Calculation.csv", ...
                fractions="chemistry/2023-06-29_RO948_Metabolite_Fractions.csv", ...
                fileprefix=stackstr(), ...
                toi="29-Jun-2023 13:10:23", ...
                crv_filename="twilite/R21_016_PET2_06292023_D1.crv", ...
                hct="41.5", ...
                time_cliff="300")
            popd(pwd0)
        end

        function test_plot_crvs(this)
            this.testObj.plot_crvs();
            [pth,fp] = fileparts(this.testObj.crv.filename);
            saveFigure2(h, fullfile(pth, fp+".fig"))
        end

        function test_ReneauAdapter_readtables(this)
            this.reneau.readtables();
            disp(this.reneau)
            disp(this.reneau.T_minimal)
            disp(this.reneau.T_frac_intact)
            this.verifyEqual(this.reneau.T_minimal{1,"wbKBq_mL"}, 6.37186782444291, RelTol=1e-9);
            this.verifyEqual(this.reneau.T_frac_intact{1,"Fraction_Intact"}, 0.171872573116216, RelTol=1e-9);
        end
        function test_ReneauAdapter_writetables(this)
            this.reneau.readtables();
            this.reneau.writetables();
            this.verifyTrue(isfile(fullfile(this.sesdir_, "chemistry", "sub-R21-016_ses-20230629_minimal.csv")))
            this.verifyTrue(isfile(fullfile(this.sesdir_, "chemistry", "sub-R21-016_ses-20230629_frac_intact.csv")))            
        end

        %% self-testing

        function test_afun(this)
            import mlwong.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_ctor(this)
            obj = this.testObj;
            this.verifyEqual(obj.toi, this.toi_)
            this.verifyEqual(obj.crv, this.crv_)
            this.verifyEqual(obj.crv_phant, this.crv_phant_)
            this.verifyEqual(obj.sesdir, this.sesdir_)
            %disp(obj)
        end
        
        function test_build_deconv(this)
            obj = this.testObj;
            h = obj.build_deconv(41.5, 47);
            saveFigure2(h, fullfile(obj.twildir, obj.fileprefix+"_deconvBayes.fig"))
        end
        function test_build_pow(this)
            %% plasma over whole-blood

            obj = this.testObj;
            h = obj.build_pow();
            saveFigure2(h, fullfile(obj.chemdir, obj.fileprefix+"_pow.fig"))
        end
        function test_Hill(this)
            obj = this.testObj;
            h = obj.build_Hill();
            saveFigure2(h, fullfile(obj.chemdir, obj.fileprefix+"_Hill.fig"))
        end
        function test_build_ptacs(this)
            obj = this.testObj;
            h = obj.build_ptacs();
            saveFigure2(h, fullfile(obj.chemdir, obj.fileprefix+"_ptacs.fig"))

            fprintf("%s:  rescaling_twilite->%g\n", stackstr(), obj.rescaling_twilite)
        end
        function test_call(this)
            obj = this.testObj;
            obj.call();

            T = readtable(fullfile( ...
                obj.chemdir, "sub-R21_016_20230629_parent_frac.csv"));
            T %#ok<NOPRT> % disp head & tail
            figure
            plot(T, "Min", "PF")
            xlabel("Time (Min)")
            ylabel("Parent Fraction (%)")
            title("sub-R21_016_20230629_parent_frac.csv")

            T = readtable(fullfile( ...
                obj.chemdir, "sub-R21_016_20230629_total_ptac.csv"));
            T %#ok<NOPRT> % disp head & tail
            figure
            plot(T, "Min", "pTAC")
            xlabel("Time (Min)")
            ylabel("Activity \mu Ci/mL")
            title("sub-R21_016_20230629_total_ptac.csv")

            T = readtable(fullfile( ...
                obj.chemdir, "sub-R21_016_20230629_metab_corr_ptac.csv"));
            T %#ok<NOPRT> % disp head & tail
            figure
            plot(T, "Min", "mc_pTAC")
            xlabel("Time (Min)")
            ylabel("Activity \mu Ci/mL")
            title("sub-R21_016_20230629_metab_corr_ptac.csv")

            fprintf("%s:  rescaling_twilite->%g\n", stackstr(), obj.rescaling_twilite)
        end
        function test_plot(this)
            obj = this.testObj;
            obj.call();
            obj.plot();
        end

        %% for NHP data from Hao Jiang & co.
        function test_call_no_twil(this)
            import mlwong.*

            sesdir = "~/Tmp/(-)-TZ3108/sub-ollie/ses-20230516";
            pwd0 = pushd(sesdir);

            toi = datetime(2023,5,16,11,21,0, TimeZone="local");
            kit = Ro948Kit( ...
                toi=toi, ...
                fileprefix="sub-ollie_ses-20230516", ...
                hct=nan, ...
                t0_forced=nan, ...
                Hill_d_min=1, ...
                do_Hill=false, ...
                ptac_units="kBq/mL");
            call(kit)
            saveFigures

            popd(pwd0)
        end
    end
    
    methods (TestClassSetup)
        function setupRo948Kit(this)
            import mlwong.*

            this.sesdir_ = "/Volumes/PrecunealSSD/Singularity/CCIR_01384/derivatives/sub-R21-016/ses-20230629";
            this.pwd0_ = pushd(this.sesdir_);

            this.toi_ = datetime(2023,6,29,13,10,23, TimeZone="local");
            this.crv_ = mlswisstrace.CrvData( ...
                fullfile(this.sesdir_, "twilite/R21_016_PET2_06292023_D1.crv"));
            this.testObj_ = Ro948Kit( ...
                toi=this.toi_, ...
                crv=this.crv_, ...
                fileprefix="sub-R21-016_ses-20230629", ...
                hct = 41.5, ...
                t0_forced = 47);

            this.reneau = ReneauAdapter( ...
                sesdir=this.sesdir_, ...
                final_calc="chemistry/2023-06-29_Final_Calculation.csv", ...
                fractions="chemistry/2023-06-29_RO948_Metabolite_Fractions.csv", ...
                fileprefix="sub-R21-016_ses-20230629");   

            this.addTeardown(@this.cleanTestClass)
        end
    end
    
    methods (TestMethodSetup)
        function setupRo948KitTest(this)
            this.testObj = copy(this.testObj_);
        end
    end
    
    properties (Access = private)
        toi_
        crv_
        crv_phant_

        testObj_
        sesdir_
        pwd0_
    end
    
    methods (Access = private)
        function cleanTestClass(this)
            popd(this.pwd0_);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
