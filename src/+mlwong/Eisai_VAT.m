classdef Eisai_VAT < handle
    %% line1
    %  line2
    %  
    %  Created 04-Jun-2025 17:25:43 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlwong/src/+mlwong.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    
    properties (Constant)
        DERIVATIVES_HOME = "/Users/jjlee/Library/CloudStorage/Box-Box/Eisai_JJL/derivatives"
        EISAI_HOME = "/Users/jjlee/Library/CloudStorage/Box-Box/Eisai_JJL"
        EVA_ids = [105, 109, 112, 114, 115, 118];
        EVA_dates = { ...
             datetime(2025,1,31), ...
            [datetime(2024,11,1),  datetime(2024,11,22)], ...
            [datetime(2025,1,7),   datetime(2025,1,28)], ...
            [datetime(2025,1,8),   datetime(2025,2,3)], ...
            [datetime(2024,12,12), datetime(2025,2,13)], ...
            [datetime(2025,2,4),   datetime(2025,2,25)]};
        EVA_labels = "EVA" + [ ...
            "105_PET1", ...
            "109_PET1", ...
            "109_PET2", ...
            "112_PET1", ...
            "112_PET2", ...
            "114_PET1", ...
            "114_PET2", ...
            "115_PET1", ...
            "115_PET2", ...
            "118_PET1", ...
            "118_PET2"];
        EXCLUDED_REGIONS = [ ...
            "ThirdVentricl"; ...
            "FrontalHorn"; ...
            "TemporaHorn"];
    end

    properties
        regions_to_rename
        region_replacement
    end

    properties (Dependent)
        audit_home
        data_cell
        data_dict  % loads and saves to "data_dict.mat"
        data_struct
        suvs_miao_xlsx
        suvs_times_mid
    end

    methods  %% GET
        function g = get.audit_home(this)
            g = fullfile(this.EISAI_HOME, "audit_of_VAT_results");
        end

        function g = get.data_cell(this)
            if ~isempty(this.data_cell_)
                g = this.data_cell_;
                return
            end

            g = {};
            g{1} = this.data_dict("EVA105").PET1;
            idx = 2;
            for sub = "EVA" + this.EVA_ids(2:end)
                g{idx} = this.data_dict(sub).PET1; %#ok<AGROW>
                g{idx+1} = this.data_dict(sub).PET2; %#ok<AGROW>
                idx = idx + 2;
            end
        end

        function set.data_cell(this, s)
            assert(iscell(s))
            this.data_cell_ = s;
        end

        function g = get.data_dict(this)
            if ~isempty(this.data_dict_)
                g = this.data_dict_;
                return
            end

            try
                ld = load(fullfile(this.DERIVATIVES_HOME, "data_dict.mat"));
                this.data_dict_ = ld.data_dict;
                g = this.data_dict_;
            catch ME
                handwarning(ME)
                g = [];
            end
        end

        function set.data_dict(this, data_dict)
            this.data_dict_ = data_dict;
            save(fullfile(this.DERIVATIVES_HOME, "data_dict.mat"), "data_dict", "-v7.3");
        end

        function g = get.data_struct(this)
            if ~isempty(this.data_struct_)
                g = this.data_struct_;
                return
            end

            this.data_struct_ = struct();
            keys = this.data_dict.keys;
            for i = 1:length(keys)
                this.data_struct_.(keys{i}) = this.data_dict(keys{i});
            end
            g = this.data_struct_;
        end

        function set.data_struct(this, s)
            assert(isstruct(s))
            this.data_struct_ = s;
        end

        function g = get.suvs_miao_xlsx(this)
            g = fullfile(this.audit_home, "SUVR", "Eisai_SUV_HM_20250324.xlsx");
        end

        function g = get.suvs_times_mid(this)
            ds = this.data_dict("EVA105");
            frame_start = ds.PET1.suvs.start_seconds_;
            frame_end = ds.PET1.suvs.end_kBq_cc_;
            g = mean([frame_start, frame_end], 2);
        end
    end

    methods
        function this = Eisai_VAT(varargin)
            this.regions_to_rename = ["TL", "FL", "G", "OL", "PL"];
            replacement = ["Temp", "Front", "Gyrus", "Occ", "Par"];
            this.region_replacement = containers.Map(this.regions_to_rename, replacement);
        end

        function this = add_atom(this, sub, ses, field, atom, opts)
            arguments
                this mlwong.Eisai_VAT
                sub {mustBeTextScalar}
                ses {mustBeTextScalar}
                field {mustBeTextScalar}
                atom {mustBeNonempty}
                opts.do_save logical = false
            end
            
            dd = this.data_dict;
            ds = dd(sub);
            ds.(ses).(field) = atom;
            dd(sub) = ds;
            if opts.do_save
                this.data_dict = dd;
            else
                this.data_dict_ = dd;
            end
        end

        function this = load_miao_xlsx(this)

            filename = this.suvs_miao_xlsx;
            %[~, sheetNames] = xlsfinfo(filename)
            %return

            this.data_dict_ = containers.Map;
            this.data_dict_("EVA105") = struct( ...
                "PET1", struct( ...
                    "date", this.EVA_dates{1}, ...
                    "suvs", readtable(filename, 'Sheet', 'EVA_105_PET1')));

            for eidx = 2:length(this.EVA_ids)
                e = this.EVA_ids(eidx);  % 109
                eva = "EVA" + e;  % EVA109                
                sheet1 = "EVA_" + e + "_PET1";
                sheet2 = "EVA_" + e + "_PET2";

                this.data_dict_(eva) = struct( ...
                    "PET1", struct( ...
                        "date", this.EVA_dates{eidx}(1), ...
                        "suvs", readtable(filename, 'Sheet', char(sheet1)) ...
                    ), ...
                    "PET2", struct( ...
                        "date", this.EVA_dates{eidx}(2), ...
                        "suvs", readtable(filename, 'Sheet', char(sheet2)) ...
                    ) ...
                );

            end
        end
       
        function [M1,M2] = M_pair(this, model, measure)
            arguments
                this mlwong.Eisai_VAT
                model {mustBeTextScalar} = "ma2"
                measure {mustBeTextScalar} = "Vt"
            end

            largeStruct = this.data_struct;
            M1 = []; 
            M2 = [];
            Region = largeStruct.EVA109.PET1.(model).Region;
            for eva = ["EVA109", "EVA112", "EVA114", "EVA115", "EVA118"]

                M1_clean = this.clean_matrix_2(largeStruct.(eva).PET1.(model), measure);
                M2_clean = this.clean_matrix_2(largeStruct.(eva).PET2.(model), measure);

                try
                    M1 = [M1, M1_clean]; %#ok<*AGROW>
                    M2 = [M2, M2_clean];
                catch ME
                    handerror(ME)
                end
            end

            if strcmp(measure, "BPnd")
                M1(abs(M1) > 10) = nan;
                M2(abs(M2) > 10) = nan;
            end
        end

        function [M120,M90,M60,M30] = M_quartet(this, model, measure)
            arguments
                this mlwong.Eisai_VAT
                model {mustBeTextScalar} = "twotcm"
                measure {mustBeTextScalar} = "Vt"
            end

            largeStruct = this.data_struct;
            M120 = []; M90 = []; M60 = []; M30 = [];
            for eva = ["EVA109", "EVA112", "EVA114", "EVA115", "EVA118"]
                for pet = ["PET1", "PET2"]
                    meas = largeStruct.(eva).(pet).(model).(measure); 
                    dur = largeStruct.(eva).(pet).(model).Duration; 
                    Region = largeStruct.(eva).(pet).(model).Region;
                    Region = Region(dur==120);
                    M120 = [M120, meas(dur==120)]; 
                    M90 = [M90, meas(dur==90)]; 
                    M60 = [M60, meas(dur==60)]; 
                    M30 = [M30, meas(dur==30)];
                end
            end

            M120 = this.clean_matrix(M120, Region);
            M90 = this.clean_matrix(M90, Region);
            M60 = this.clean_matrix(M60, Region);
            M30 = this.clean_matrix(M30, Region);
        end

        %% collection with import_excel_to_struct()

        function data_struct = import_excel_to_struct(this, excel_files, output_struct)
            % IMPORT_EXCEL_TO_STRUCT Imports Excel files into organized MATLAB structure
            %
            % Inputs:
            %   excel_files - Cell array of Excel file paths or struct with file info
            %   output_struct - Pre-existing structure to populate (optional)
            %
            % Output:
            %   data_struct - Organized structure with imported data
            %
            % Example usage:
            %   files = {'logan.xlsx', 'ma2.xlsx', 'mrtm2.xlsx', 'logan_ref.xlsx', 'srtm2.xlsx', 'twotcm.xlsx'};
            %   data = import_excel_to_struct(files);
            %
            % More example usage and testing:
            %
            % % Method 1: Specify individual files
            % files = {'logan.xlsx', 'ma2.xlsx', 'mrtm2.xlsx', 'logan_ref.xlsx', 'srtm2.xlsx', 'twotcm.xlsx'};
            % data = import_excel_to_struct(files);
            %
            % % Method 2: Import from directory
            % data = import_excel_from_directory('/path/to/your/excel/files');
            %
            % % Display structure information
            % display_structure_info(data);
            %
            % % Access specific data
            % % Example: Access logan data for EVA001, PET1
            % if isfield(data, 'EVA001') && isfield(data.EVA001, 'PET1') && isfield(data.EVA001.PET1, 'logan')
            %     logan_data = data.EVA001.PET1.logan;
            %     disp(logan_data);
            % end

            % Initialize output structure if not provided
            if nargin < 3
                output_struct = struct();
            end
            data_struct = output_struct;

            % Valid file identifiers
            valid_file_ids = {'logan', 'ma2', 'mrtm2', 'mrtm2_cbm', 'logan_ref', 'logan_ref_cbm', 'srtm2', 'srtm2_cbm', 'twotcm'};

            % Process each Excel file
            for file_idx = 1:length(excel_files)
                current_file = excel_files{file_idx};

                % Extract file identifier from filename
                [~, filename, ~] = fileparts(current_file);
                file_id = this.extract_file_identifier(filename, valid_file_ids);

                if isempty(file_id)
                    warning('File %s does not match expected naming pattern. Skipping.', current_file);
                    continue;
                end

                fprintf('Processing file: %s (ID: %s)\n', current_file, file_id);

                % Get all sheet names from the Excel file
                sheet_names = sheetnames(current_file);

                % Process each sheet
                for sheet_idx = 1:length(sheet_names)
                    sheet_name = sheet_names{sheet_idx};

                    % Check if sheet name matches expected pattern
                    if ~this.matches_sheet_pattern(sheet_name)
                        warning('Sheet %s does not match expected pattern eva\\d{3}_pet\\d{1}. Skipping.', sheet_name);
                        continue;
                    end

                    % Extract EVA and PET identifiers from sheet name
                    [eva_id, pet_id] = this.extract_eva_pet_ids(sheet_name);

                    fprintf('  Processing sheet: %s (EVA: %s, PET: %s)\n', sheet_name, eva_id, pet_id);

                    try
                        % Read the table from the current sheet
                        table_data = readtable(current_file, 'Sheet', sheet_name, ...
                            'VariableNamingRule', 'preserve');

                        % Create the nested structure path
                        eva_field = upper(eva_id);  % Convert to uppercase for field name
                        pet_field = upper(pet_id);  % Convert to uppercase for field name

                        % Initialize structure levels if they don't exist
                        if ~isfield(data_struct, eva_field)
                            data_struct.(eva_field) = struct();
                        end

                        if ~isfield(data_struct.(eva_field), pet_field)
                            data_struct.(eva_field).(pet_field) = struct();
                        end

                        % Assign the table to the appropriate field
                        data_struct.(eva_field).(pet_field).(file_id) = table_data;

                        fprintf('    Successfully imported %d rows, %d columns\n', ...
                            height(table_data), width(table_data));

                    catch ME
                        warning('Error reading sheet %s from file %s: %s', ...
                            sheet_name, current_file, ME.message);
                        continue;
                    end
                end
            end

            fprintf('Import complete!\n');
        end

        % Alternative function for batch processing from a directory
        function data_struct = import_excel_from_directory(this, directory_path, output_struct)
            % IMPORT_EXCEL_FROM_DIRECTORY Import all Excel files from a directory
            %
            % Inputs:
            %   directory_path - Path to directory containing Excel files
            %   output_struct - Pre-existing structure to populate (optional)

            if nargin < 3
                output_struct = struct();
            end

            % Get all Excel files in the directory
            excel_files = dir(fullfile(directory_path, '*.xlsx'));

            if isempty(excel_files)
                error('No Excel files found in directory: %s', directory_path);
            end

            % Create full file paths
            file_paths = cell(length(excel_files), 1);
            for i = 1:length(excel_files)
                file_paths{i} = fullfile(directory_path, excel_files(i).name);
            end

            % Import the files
            data_struct = this.import_excel_to_struct(file_paths, output_struct);
        end

        % Utility function to display the structure
        function display_structure_info(this, data_struct)
            % DISPLAY_STRUCTURE_INFO Display information about the imported structure

            fprintf('\n=== Structure Information ===\n');
            eva_fields = fieldnames(data_struct);

            for i = 1:length(eva_fields)
                eva_field = eva_fields{i};
                fprintf('EVA Field: %s\n', eva_field);

                pet_fields = fieldnames(data_struct.(eva_field));
                for j = 1:length(pet_fields)
                    pet_field = pet_fields{j};
                    fprintf('  PET Field: %s\n', pet_field);

                    file_fields = fieldnames(data_struct.(eva_field).(pet_field));
                    for k = 1:length(file_fields)
                        file_field = file_fields{k};
                        table_data = data_struct.(eva_field).(pet_field).(file_field);
                        fprintf('    File: %s (%d rows, %d cols)\n', ...
                            file_field, height(table_data), width(table_data));
                    end
                end
            end
        end

        function write_summary_xlsx(this)
            mdl0 = ["twotcm", "ma2", "logan", "mrtm2", "logan_ref", "srtm2", "mrtm2_cbm", "logan_ref_cbm", "srtm2_cbm"];
            mdl1 = ["2TCM", "MA2", "Logan", "MRTM2 (pons)", "Logan Ref. Tiss. (pons)", "SRTM2 (pons)", "MRTM2 (cbm)", "Logan Ref. Tiss. (cbm)", "SRTM2 (cbm)"];
            mdl_map = containers.Map(mdl0, mdl1);
            measure = ["Vt", "Vt", "Vt", "BPnd", "BPnd", "BPnd", "BPnd", "BPnd", "BPnd"];
            meas_map = containers.Map(mdl0, measure);
            twotcm_duration = [repmat(30, [46, 1]); repmat(60, [46, 1]); repmat(90, [46, 1]); repmat(120, [46, 1])];

            S = this.data_struct;

            for imdl0 = 1:length(mdl0)
                mdl = mdl0(imdl0);
                meas = measure(imdl0);

                Region = S.EVA109.PET1.(mdl).Region;
                T = table(Region);
                if strcmp(mdl, "twotcm")
                    T = addvars(T, twotcm_duration, NewVariableNames="Scan_Duration");
                end
                T = this.clean_table(T);

                for eva = ["EVA109", "EVA112", "EVA114", "EVA115", "EVA118"]
                    for pet = ["PET1", "PET2"]
                        new_col = S.(eva).(pet).(mdl).(meas);
                        scan_region = S.(eva).(pet).(mdl).Region;
                        new_col = this.clean_matrix(new_col, scan_region);
                        T = addvars(T, new_col, NewVariableNames=sprintf("%s_%s_%s", eva, pet, meas));
                    end
                end

                writetable(T, "Eisai_VAT_complete.xlsx", Sheet=mdl1(imdl0));
            end
        end
    end

    methods (Static)
        function newTable = averageLRColumns(T)
            varNames = T.Properties.VariableNames;

            % Find columns ending with '_l' and '_r'
            lCols = endsWith(varNames, '_l');
            rCols = endsWith(varNames, '_r');

            % Get base names by removing only trailing '_l' and '_r'
            lBaseNames = regexprep(varNames(lCols), '_l$', '');
            rBaseNames = regexprep(varNames(rCols), '_r$', '');

            % Find matching pairs
            [matchingBases, ~, ~] = intersect(lBaseNames, rBaseNames);

            % Start with original table
            newTable = T;

            % Process each matching pair
            for i = 1:length(matchingBases)
                baseName = matchingBases{i};
                lColName = [baseName '_l'];
                rColName = [baseName '_r'];

                % Calculate average and replace
                newTable.(baseName) = (newTable.(lColName) + newTable.(rColName)) / 2;
                newTable = removevars(newTable, {lColName, rColName});
            end
        end

        function largeStruct = map2struct(myMap)
            largeStruct = struct();
            keys = myMap.keys;

            for i = 1:length(keys)
                key = keys{i};
                % Ensure key is a valid field name
                if isvarname(key)
                    largeStruct.(key) = myMap(key);
                else
                    % Handle invalid field names
                    validKey = matlab.lang.makeValidName(key);
                    largeStruct.(validKey) = myMap(key);
                end
            end
        end

        function new_S = struct_for_anat_regions(old_S)

            import mlwong.Eisai_VAT.table_for_anat_regions

            % Initialize result structure
            new_S = old_S;

            % Get all field names at level 1
            level1_fields = fieldnames(old_S);

            % Iterate through level 1 fields
            for i = 1:length(level1_fields)
                field1 = level1_fields{i};

                % Check if level 1 field exists and is a struct
                if isstruct(old_S.(field1))
                    % Get field names at level 2
                    level2_fields = fieldnames(old_S.(field1));

                    % Iterate through level 2 fields
                    for j = 1:length(level2_fields)
                        field2 = level2_fields{j};

                        % Check if level 2 field exists and is a struct
                        if isstruct(old_S.(field1).(field2))
                            % Get field names at level 3
                            level3_fields = fieldnames(old_S.(field1).(field2));

                            % Look for interesting fields at level 3
                            interesting_found = false;
                            for k = 1:length(level3_fields)
                                field3 = level3_fields{k};

                                % Check if field is one of the interesting ones
                                if strcmp(field3, 'suvs')
                                    % Initialize nested structure if needed
                                    if ~isfield(new_S, field1)
                                        new_S.(field1) = struct();
                                    end
                                    if ~isfield(new_S.(field1), field2)
                                        new_S.(field1).(field2) = struct();
                                    end

                                    suvs = old_S.(field1).(field2).(field3);
                                    suvs = mlwong.Eisai_VAT.averageLRColumns(suvs);
                                    suvs.Var1 = []; suvs.start_seconds_ = []; suvs.end_kBq_cc_ = [];
                                    suvs = rows2vars(suvs);
                                    suvs = renamevars(suvs, 'OriginalVariableNames', 'Region');
                                    suvs = table_for_anat_regions(suvs);
                                    new_S.(field1).(field2).(field3) = suvs;
                                end
                                if strcmp(field3, 'logan_ref') || ...
                                        strcmp(field3, 'mrtm2') || ...
                                        strcmp(field3, 'srtm2')
                                    % Initialize nested structure if needed
                                    if ~isfield(new_S, field1)
                                        new_S.(field1) = struct();
                                    end
                                    if ~isfield(new_S.(field1), field2)
                                        new_S.(field1).(field2) = struct();
                                    end

                                    % Extract the data maintaining hierarchy
                                    interesting = old_S.(field1).(field2).(field3);
                                    interesting = table_for_anat_regions(interesting);
                                    new_S.(field1).(field2).([field3, '_pons']) = interesting;
                                    interesting_found = true;
                                end
                                if strcmp(field3, 'logan') || ...
                                        strcmp(field3, 'ma2') || ...
                                        strcmp(field3, 'twotcm') || ...
                                        strcmp(field3, 'logan_ref_cbm') || ...
                                        strcmp(field3, 'mrtm2_cbm') || ...
                                        strcmp(field3, 'srtm2_cbm')
                                    % Initialize nested structure if needed
                                    if ~isfield(new_S, field1)
                                        new_S.(field1) = struct();
                                    end
                                    if ~isfield(new_S.(field1), field2)
                                        new_S.(field1).(field2) = struct();
                                    end

                                    % Extract the data maintaining hierarchy
                                    interesting = old_S.(field1).(field2).(field3);
                                    interesting = table_for_anat_regions(interesting);
                                    new_S.(field1).(field2).(field3) = interesting;
                                    interesting_found = true;
                                end
                            end

                            if interesting_found
                                fprintf('Found interesting fields at S.%s.%s\n', field1, field2);
                            end
                        end
                    end
                end
            end
        end

        function myMap = struct2map(largeStruct)
            fieldNames = fieldnames(largeStruct);
            fieldValues = struct2cell(largeStruct);
            myMap = containers.Map(fieldNames, fieldValues);
        end
    
        function T_suv = suv(T_kBq, subid, sesid)
            %% create SUVs from Noah's TACs
            %
            % e.g.:
            % cd("NoahGoldman")
            % subs = "EVA" + [109, 112, 114, 115, 118];
            % sess = ["baseline_tacs", "retest_tacs"];
            % for s = subs
            %     for e = sess
            %         T = readtable(fullfile("TACs", "sub-" + s, sprintf("sub-%s_ses-%s.csv", s, e)));
            %         T1 = mlwong.Eisai_VAT.suv(T, s, e);
            %         writetable(T1, fullfile("SUVs", "sub-" + s, sprintf("sub-%s_ses-%s.csv", s, strrep(e, "tacs", "suvs"))));
            %     end
            % end

            arguments
                T_kBq table
                subid {mustBeTextScalar}  % "EVA109"
                sesid {mustBeNonempty}  % 1 | "baseline" | "baseline_tacs"
            end
            if istext(sesid)
                if contains(sesid, "baseline")
                    sesid = 1;
                elseif contains(sesid, "retest")
                    sesid = 2;
                else
                    error("mlwong:ValueError", stackstr());
                end
            end

            DHOME = mlwong.Eisai_VAT.DERIVATIVES_HOME;

            % dosing data
            T_dosing = readtable(fullfile(DHOME, 'dosing.csv'));
            T_dosing.Subject = strrep(convertCharsToStrings(T_dosing.Subject), "_", "");
            T_dosing = addvars(T_dosing, ...
                T_dosing.Weight_lbs_ / 2.2, After="Weight_lbs_", NewVariableNames="Weight_kg_");
            T_dosing = T_dosing(strcmp(subid, T_dosing.Subject), :);

            % timing data
            T_timings = readtable(fullfile(DHOME, "frame_timings.xlsx"));

            % assembly
            T_kBq(:, 1) = [];
            T_suv = array2table(suv_of_col(table2array(T_kBq)), VariableNames=T_kBq.Properties.VariableNames);
            T_suv.Properties.VariableNames{1} = 'MidFrameTime_sec_';
            T_suv.MidFrameTime_sec_ = T_timings.mid_frameTime_s_;

            function val1 = suv_of_col(val0)
                wt = 1e3 * T_dosing{sesid, "Weight_kg_"};  % g
                dose = T_dosing{sesid, "InjectedDose_bq_"} / 1e3;  % kBq
                val1 = val0 * wt / dose;
            end
        end

        function modified_table = table_for_anat_regions(input_table)
            % Modifies regions in place

            % Remove uninformative columns
            if ismember('No', input_table.Properties.VariableNames)
                input_table.No = [];
            end
            if ismember('Model', input_table.Properties.VariableNames)
                input_table.Model = [];
            end

            % Remove less interesting rows
            R = input_table.Region;
            input_table(strcmp(R, "Corp_Callosum"), :) = [];
            R = input_table.Region;
            input_table(strcmp(R, "Brainstem"), :) = [];
            R = input_table.Region;
            input_table(strcmp(R, "ThirdVentricl"), :) = [];
            R = input_table.Region;
            input_table(strcmp(R, "Cerebellum"), :) = [];
            R = input_table.Region;
            input_table(strcmp(R, "Cerebellum_wm"), :) = [];
            R = input_table.Region;
            input_table(strcmp(R, "Pons"), :) = [];
            R = input_table.Region;
            input_table(strcmp(R, "FrontalHorn"), :) = [];
            R = input_table.Region;
            input_table(strcmp(R, "G_fus"), :) = [];
            R = input_table.Region;
            input_table(strcmp(R, "Presubgen_antCing"), :) = [];
            R = input_table.Region;
            input_table(strcmp(R, "S_nigra"), :) = [];
            R = input_table.Region;
            input_table(strcmp(R, "Subcall_area"), :) = [];
            R = input_table.Region;
            input_table(strcmp(R, "Subgen_antCing"), :) = [];
            R = input_table.Region;
            input_table(strcmp(R, "TemporaHorn"), :) = [];

            % Create region groups
            regions = input_table.Region;
            for i = 1:height(input_table)
                if startsWith(regions{i}, 'G_cing')
                    input_table.Region{i} = 'Cingulate';
                elseif startsWith(regions{i}, 'FL_')
                    input_table.Region{i} = 'Frontal';
                elseif startsWith(regions{i}, 'OL_')
                    input_table.Region{i} = 'Occipital';
                elseif startsWith(regions{i}, 'PL_')
                    input_table.Region{i} = 'Parietal';
                elseif startsWith(regions{i}, 'Ant_TL_')
                    input_table.Region{i} = 'Temporal';
                elseif startsWith(regions{i}, 'G_paraH_amb')
                    input_table.Region{i} = 'Temporal';
                elseif startsWith(regions{i}, 'G_sup_temp_ant')
                    input_table.Region{i} = 'Temporal';
                elseif startsWith(regions{i}, 'G_sup_temp_post')
                    input_table.Region{i} = 'Temporal';
                elseif startsWith(regions{i}, 'G_tem_midin')
                    input_table.Region{i} = 'Temporal';
                elseif startsWith(regions{i}, 'Post_TL')
                    input_table.Region{i} = 'Temporal';
                elseif startsWith(regions{i}, 'Corp_Callosum')
                    input_table.Region{i} = 'White';     
                elseif startsWith(regions{i}, 'White_matter')
                    input_table.Region{i} = 'White';                
                end
            end

            % Group by Region and calculate means
            selected_vars = ~strcmp(input_table.Properties.VariableNames, "Region");
            modified_table = grpstats(input_table, 'Region', 'mean', 'DataVars',selected_vars);

            % Clean up variable names (remove 'mean_' prefix)
            var_names = modified_table.Properties.VariableNames;
            for i = 1:length(var_names)
                if startsWith(var_names{i}, 'mean_')
                    modified_table.Properties.VariableNames{i} = var_names{i}(6:end);
                end
            end

            % Reorder columns to put Region first
            col_order = [{'Region'}, setdiff(modified_table.Properties.VariableNames, {'Region'})];
            modified_table = modified_table(:, col_order);
        end

        %% support for import_excel_to_struct()

        function file_id = extract_file_identifier(filename, valid_ids)
            % Extract file identifier from filename
            file_id = '';

            % Check each valid identifier
            for i = 1:length(valid_ids)
                if strcmpi(filename, valid_ids{i})
                    file_id = valid_ids{i};
                    break;
                end
            end
        end

        function is_match = matches_sheet_pattern(sheet_name)
            % Check if sheet name matches the pattern eva\d{3}_pet\d{1}
            pattern = '^eva\d{3}_pet\d{1}$';
            is_match = ~isempty(regexp(sheet_name, pattern, 'once'));
        end

        function [eva_id, pet_id] = extract_eva_pet_ids(sheet_name)
            % Extract EVA and PET identifiers from sheet name
            % Pattern to capture the numbers
            pattern = '^eva(\d{3})_pet(\d{1})$';
            tokens = regexp(sheet_name, pattern, 'tokens');

            if ~isempty(tokens)
                eva_id = ['EVA' tokens{1}{1}];  % e.g., 'EVA001'
                pet_id = ['PET' tokens{1}{2}];  % e.g., 'PET1'
            else
                eva_id = '';
                pet_id = '';
            end
        end
    end

    %% PRIVATE

    properties (Access = private)
        data_cell_
        data_dict_
        data_struct_
    end

    methods (Access = private)
        function M = clean_matrix(this, M, Regions)
            Regions = ascol(Regions);
            try
                assert(size(M, 1) == size(Regions, 1))
                M = M(~contains(Regions, this.EXCLUDED_REGIONS), :);
            catch ME
                handerror(ME)
            end
        end
        function M = clean_matrix_2(this, tbl, measure)
            selection = ~contains(tbl.Region, this.EXCLUDED_REGIONS);
            try
                M = tbl.(measure)(selection, :);
            catch ME
                handerror(ME)
            end
        end
        function T = clean_table(this, T)
            T = T(~contains(T.Region, this.EXCLUDED_REGIONS), :);
            for k = this.regions_to_rename
                T.Region = strrep(T.Region, k, this.region_replacement(k));
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
