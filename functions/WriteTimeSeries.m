function WriteTimeSeries(scName, opt, MPopt, moonName, orbNum, outFstr, FULLORBITS)
% Evaluate a given magnetic field model and optional magnetopause model at spacecraft locations and
% save to disk.
%
% Currently implemented for Galileo, Juno, and Cassini.
%
% Datasets:
%
%   * Cassini MAG: https://doi.org/10.17189/5rhj-sm88, volume CO-E/SW/J/S-MAG-4-SUMM-1MINAVG-V2.1
%   * Galileo MAG: https://doi.org/10.17189/1519667, volume GO-J-MAG-3-RDR-HIGHRES-V1.0 (moons)
%       https://doi.org/10.17189/1519668, volume GO-J-MAG-3-RDR-MAGSPHERIC-SURVEY-V1.0 (full 
%       orbits)
%   * Juno MAG: https://doi.org/10.17189/1519711, volume JNO-J-3-FGM-CAL-V1.0
%
% Parameters
% ----------
% scName : string, default="Galileo"
%   Spacecraft name for which magnetic field data will be compared against implemented models. A
%   directory must exist with this name in the ``MAG`` directory within the top-level P\lanetMag
%   directory. This directory will be searched for body-specific data files, and each of these
%   files will be loaded.
% opt : int, default=0
%   Magnetic field model option ID number to evaluate, as defined in GetModelOpts.
% MPopt : int, default=-1
%   Magnetopause model option ID number to evaluate, as defined in GetModelOpts. Passing -1 skips
%   magnetopause modeling.
% moonName : char, 1xC, default='Io'
%   Name of moon for which ``scName`` has flyby data. Flyby data is not used if ``FULLORBITS`` is
%   true, but a name must be passed for kernel loading.
% orbNum : int, default=-1
%   Orbit number (or year for Cassini) to use for data comparison. Accepts the following special
%   options:
%
%   - ``-1``: Use all orbits/years for which data is available.
%   - ``-2``: Use all orbits for which flyby data is available for ``moonName``.
%
% outFstr : char, 1xE, default='modelDescrip'
%   String to append to PDS file names to make output file name. Passing 'modelDescrip' adds the
%   magModelDescrip file name code string.
% FULLORBITS : bool, default=0
%   Whether to use full-orbit MAG data files or flyby-only files. Only applies to Galileo data.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('scName', 'var'); scName = "Galileo"; end
    if ~exist('opt', 'var'); opt = 0; end
    if ~exist('MPopt', 'var'); MPopt = -1; end
    if ~exist('moonName', 'var'); moonName = 'Io'; end
    if ~exist('orbNum', 'var'); orbNum = -1; end
    if ~exist('outFstr', 'var'); outFstr = 'modelDescrip'; end
    if ~exist('FULLORBITS', 'var'); FULLORBITS = 0; end

    cspice_kclear;
    SetPlotDefaults();
    
    sc = char(scName);
    parentName = LoadSpice(moonName, sc);
    [~, ~, ~, a_AU, ~, ~, ~, ~, ~, ~] = GetBodyParams(moonName);

    if strcmp(scName, "Cassini")
        formatSpec = '%19s%11f%11f%11f%11f%8f%7f%7f%5f%3d%[^\n\r]';
        delim = ',';
        xtn = '.TAB';
        nHeadLines = 0;
        inDir = fullfile('MAG', sc);

        if orbNum == -1 || orbNum == -2
            yearRange = 4:17;
        else
            yearRange = orbNum:orbNum;
        end
        nFiles = length(yearRange);
        fileList = strings(1,nFiles);
        for i=1:nFiles
            fileList(i) = fullfile(inDir, [num2str(2000+yearRange(i)) '_FGM_KRTP_1M']);
        end
            
    elseif strcmp(scName, "Juno")
        formatSpec = '%d%d%d%d%d%d%f%f%f%f%f%f%f%f%[^\n\r]';
        delim = '';
        xtn = '.sts';
        inDir = fullfile('MAG', sc, 'Jupiter');

        if strcmp(moonName, 'Io')
            dayRange2023 = 0:0;
            dayRange2024 = 0:0;
            nFiles2023 = length(dayRange2023);
            nFiles2024 = length(dayRange2024);
            nFiles = nFiles2023 + nFiles2024;
            fileList = strings(1,nFiles);
            for i=1:nFiles2023
                fileList(i) = fullfile(inDir, ...
                    ['fgm_jno_l3_2023' sprintf('%03d', dayRange2023(i)) 'pc_r1s_v01']);
            end
            for i=(nFiles2023+1):nFiles2024
                fileList(i) = fullfile(inDir, ...
                    ['fgm_jno_l3_2024' sprintf('%03d', dayRange2024(i)) 'pc_r1s_v01']);
            end
        elseif strcmp(moonName, 'Europa')
            dayRange = 252:290;
            nFiles = length(dayRange);
            fileList = strings(1,nFiles);
            for i=1:nFiles
                fileList(i) = fullfile(inDir, ...
                    ['fgm_jno_l3_2022' sprintf('%03d', dayRange(i)) 'pc_r1s_v01']);
            end
            fileList(dayRange==272) = fullfile(inDir, 'fgm_jno_l3_2022272pc_r1s_v02');
        elseif strcmp(moonName, 'Ganymede')
            dayRange = 133:179;
            nFiles = length(dayRange);
            fileList = strings(1,nFiles);
            for i=1:nFiles
                fileList(i) = fullfile(inDir, ...
                    ['fgm_jno_l3_2021' sprintf('%03d', dayRange(i)) 'pc_r1s_v01']);
            end
            fileList(dayRange==159) = fullfile(inDir, 'fgm_jno_l3_2021159pc_r1s_v02');
        else
            
        end
        
    elseif strcmp(scName, "Galileo")
        sc = 'Galileo Orbiter';
        xtn = '.TAB';
        nHeadLines = 0;
        delim = '';

        switch moonName
            case 'Io'
                fbList = [0, 24, 27, 31, 32];
                fbCode = 'IO';
            case 'Europa'
                fbList = [4, 11, 12, 14, 15, 19, 26];
                fbCode = 'EUR';
            case 'Ganymede'
                fbList = [1, 2, 7, 8, 28, 29];
                fbCode = 'GAN';
            case 'Callisto'
                fbList = [3, 9, 10, 30];
                fbCode = 'CALL';
        end
        
        if orbNum == -1
            orbList = linspace(0, 35, 36);
            orbList(orbList == 5) = []; % Solar conjunction on orbit 5 means no data at all

        elseif orbNum == -2
            orbList = fbList;

        else
            orbList = orbNum:orbNum;
        end

        nFiles = length(orbList);
        fileList = strings(1,nFiles);
        if FULLORBITS
            formatSpec = '%23s%24s%10f%10f%10f%10f%7f%7f%7f%f%[^\n\r]';
            inDir = fullfile('MAG', char(scName), 'Jupiter');
            for i=1:nFiles
                fileList(i) = fullfile(inDir, ['ORB' sprintf('%02d', orbList(i)) '_SYS3']);
            end
        else
            formatSpec = '%23s%10f%10f%10f%10f%7f%7f%7f%f%[^\n\r]';
            inDir = fullfile('MAG', char(scName), moonName);
            for i=1:nFiles
                fileList(i) = fullfile(inDir, ['ORB' sprintf('%02d_', orbList(i)) fbCode '_SYS3']);
            end
        end
        
    else
        error(['Spacecraft "' char(scName) '" not implemented in WriteTimeSeries.'])
    end

    %% Finished building spacecraft-specific parts, now evaluate and save

    [MagModel, CsheetModel, MPmodel, magModelDescrip, fEnd] = GetModelOpts(parentName, opt, MPopt);
    DO_MPAUSE = MPopt ~= -1;
    disp(['Evaluating ' magModelDescrip ' field model for ' num2str(nFiles) ' files.'])

    for i=1:nFiles
        thisFile = [char(fileList(i)) xtn];
        disp([thisFile ', ' num2str(i) ' of ' num2str(nFiles) '.'])
        fileID = fopen(thisFile,'r');
        if strcmp(scName, "Juno")
            [~, fileStr, ~] = fileparts(thisFile);
            nHeadLines = 0;
            line = textscan(fileID, '%s', 1, 'Delimiter', '\n');
            headLine = pad(line{1}{1}, 9);
            while ~strcmp([headLine(1:4) headLine(6:8)], fileStr(12:18)) || strcmp(headLine(9), 'T')
                line = textscan(fileID, '%s', 1, 'Delimiter', '\n');
                headLine = pad(line{1}{1}, 9);
                nHeadLines = nHeadLines + 1;
            end
                  
            frewind(fileID);
            magData = textscan(fileID, formatSpec, Inf, 'Delimiter', delim, 'TextType', 'char', ...
                'EndOfLine', '\r\n', 'HeaderLines', nHeadLines);
            yyyy = magData{1}';
            doy =  magData{2}';
            h =    magData{3}';
            m =    magData{4}';
            s =    magData{5}';
            ms =   magData{6}';
            npts = length(yyyy);
            t_UTC = strings(npts, 1);
            for j=1:npts
                t_UTC(j) = sprintf('%d-%03d//%02d:%02d:%02d.%03d', yyyy(j), doy(j), h(j), m(j), s(j), ms(j));
            end
            t_UTC = char(t_UTC);
        else
            magData = textscan(fileID, formatSpec, Inf, 'Delimiter', delim, 'TextType', 'char', ...
                'EndOfLine', '\r\n', 'HeaderLines', nHeadLines);
            t_UTC = magData{1}';
        end
        fclose(fileID);
    
        ets = cspice_str2et(t_UTC);
        t_h = ets / 3600;
        npts = length(ets);
        [r_km, theta, phi, xyz_km, S3coords] = GetPosSpice(sc, parentName, t_h);
    
        if contains(magModelDescrip, 'KS2005')
            [Bvec, Mdip_nT, Odip_km] = MagFldJupiterKS2005(r_km, theta, phi, ets);
        else
            [Bvec, Mdip_nT, Odip_km] = MagFldParent(parentName, r_km, theta, phi, MagModel, ...
                CsheetModel);
        end
        if DO_MPAUSE && ~(strcmp(MPmodel, 'None') || contains(magModelDescrip, 'KS2005'))
    
            nSW_pcc = 4 / a_AU^2 * ones(1,npts);
            vSW_kms = 400  * ones(1,npts);
            [mpBvec, OUTSIDE_MP] = MpauseFld(nSW_pcc, vSW_kms, ets, xyz_km, Mdip_nT, Odip_km, ...
                parentName, S3coords, MPmodel);
            Bvec = Bvec + mpBvec;
            Bvec(:,OUTSIDE_MP) = 0;
    
        end
    
        BxS3_nT = Bvec(1,:)';
        ByS3_nT = Bvec(2,:)';
        BzS3_nT = Bvec(3,:)';
    
        % Save data to disk
        if strcmp(outFstr, 'modelDescrip')
            outFstr = fEnd;
        end
        outFname = fullfile([char(fileList(i)) outFstr '.mat']);
        save(outFname, 't_UTC', 'BxS3_nT', 'ByS3_nT', 'BzS3_nT');
        disp(['Data saved to file: ' outFname])
    end
end
