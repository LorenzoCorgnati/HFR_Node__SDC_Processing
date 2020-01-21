%% radial_netCDF_aggregation_v212.m
% This function accesses the THREDDS catalog of the HFR networks via
% OpenDAP and creates HFR radial aggregated netCDF datasets compliant to the
% European standard data model (that integrates CMEMS-INSTAC and SDC CF extension
% requirements) for distribution on the SeaDataNet infrastructure.
% The v212 version creates aggregated dataset according the the v2.1.2
% version of the European standard data model.

% INPUT:
%         networkData: cell array containing information about the network
%                      (metadata)
%         networkFields: field names of the cell array containing
%                       information about the network.
%         stationData: cell array containing information about the station
%                      (metadata)
%         stationFields: field names of the cell array containing
%                       information about the station.
%         timeSpan: temporal interval for aggregation (expressed in months)

% OUTPUT:
%         rnA_err: error flag (0 = correct, 1 = error)
%         ncFileNoPath: filename of the generated nc file, without the full path
%         ncFilesize: size of the generated nc file.
%         tStart: starting time of the dataset (in datenum format).
%         tEnd: ending time of the dataset (in datenum format).
%         dataID: SDN local CDI id.

% Author: Lorenzo Corgnati
% Date: December 4, 2019

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [rnA_err, ncFileNoPath, ncFilesize, tStart, tEnd, dataID] = radial_netCDF_aggregation_v212(networkData,networkFields,stationData,stationFields,timeSpan)

disp(['[' datestr(now) '] - - ' 'radial_netCDF_aggregation_v212.m started.']);

rnA_err = 0;

% Set the output variables in case of function exit due to absence of data
ncFileNoPath = '';
ncFilesize = 0;
tStart = 0;
tEnd = 0;
dataID = '';

warning('off', 'all');

%% Set scale_factor and add_offset values

scaleFactor = 0.001;
addOffset = 0;

%%

%% Retrieve the origin coordinates

try
    site_lonIndex = find(not(cellfun('isempty', strfind(stationFields, 'site_lon'))));
    site_latIndex = find(not(cellfun('isempty', strfind(stationFields, 'site_lat'))));
    Origin = [stationData{site_latIndex},stationData{site_lonIndex}];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

%%

%% Retrieve EDIOS and EDMO codes, site code and platform code

try
    network_idIndex = find(not(cellfun('isempty', strfind(networkFields, 'network_id'))));
    EDIOS_Series_ID = networkData{network_idIndex};
    station_idIndex = find(not(cellfun('isempty', strfind(stationFields, 'station_id'))));
    siteCode = stationData{station_idIndex};
    EDIOS_Platform_ID = siteCode;
    EDMO_codeIndex = find(not(cellfun('isempty', strfind(stationFields, 'EDMO_code'))));
    EDMO_code = stationData{EDMO_codeIndex};
    site_code = EDIOS_Series_ID;
    platform_code = [EDIOS_Series_ID '-' EDIOS_Platform_ID];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

%%

%% Set non physical dimensions
try
    maxSite_dim = 1;
    maxInst_dim = length(EDMO_code);
    refMax_dim = 1;
    string20_dim = 20;
    string50_dim = 50;
    string80_dim = 80;
    string250_dim = 250;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

%%

%% Read the aggregated radial file from THREDDS catalog via OpenDAP

try
    % Find the SDC_OpenDAP_data_url field from station data
    SDC_OpenDAP_data_urlIndex = find(not(cellfun('isempty', strfind(stationFields, 'SDC_OpenDAP_data_url'))));
    SDC_OpenDAP_data_url = stationData{SDC_OpenDAP_data_urlIndex};
    
    % Retrieve manufacturer info
    sensorATT = ncreadatt(SDC_OpenDAP_data_url,'/','sensor');
    
    % Read time and convert it to Matlab time
    nc.time = ncread_cf_time(SDC_OpenDAP_data_url,'TIME');
    
    % Select the data range according to the aggregation time span
    t0 = datetime(datevec(now));
    tEnd = datenum(dateshift(t0,'start','month'));
    t1 = t0 - calmonths(timeSpan);
    tStart = datenum(dateshift(t1,'start','month'));
    iTime = find(nc.time>=tStart & nc.time<tEnd);
    nc.time = nc.time(iTime);
    
    % Read variables
    if (~isempty(nc.time))
        % Coordinate variables
        nc.bear = ncread(SDC_OpenDAP_data_url,'BEAR');
        nc.rnge = ncread(SDC_OpenDAP_data_url,'RNGE');
        nc.depth = ncread(SDC_OpenDAP_data_url,'DEPH');
        
        nc.latitude = ncread(SDC_OpenDAP_data_url,'LATITUDE');
        nc.longitude = ncread(SDC_OpenDAP_data_url,'LONGITUDE');
        
        % Data variables
        nc.rdva = ncread(SDC_OpenDAP_data_url,'RDVA',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.drva = ncread(SDC_OpenDAP_data_url,'DRVA',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.ewct = ncread(SDC_OpenDAP_data_url,'EWCT',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.nsct = ncread(SDC_OpenDAP_data_url,'NSCT',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.espc = ncread(SDC_OpenDAP_data_url,'ESPC',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.etmp = ncread(SDC_OpenDAP_data_url,'ETMP',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.maxv = ncread(SDC_OpenDAP_data_url,'MAXV',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.minv = ncread(SDC_OpenDAP_data_url,'MINV',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.ersc = ncread(SDC_OpenDAP_data_url,'ERSC',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.ertc = ncread(SDC_OpenDAP_data_url,'ERTC',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.xdst = ncread(SDC_OpenDAP_data_url,'XDST');
        nc.ydst = ncread(SDC_OpenDAP_data_url,'YDST');
        nc.sprc = ncread(SDC_OpenDAP_data_url,'SPRC',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.narx = ncread(SDC_OpenDAP_data_url,'NARX',min(iTime),length(iTime));
        nc.natx = ncread(SDC_OpenDAP_data_url,'NATX',min(iTime),length(iTime));
        nc.sltr = ncread(SDC_OpenDAP_data_url,'SLTR',[1,min(iTime)],[maxSite_dim,length(iTime)]);
        nc.slnr = ncread(SDC_OpenDAP_data_url,'SLNR',[1,min(iTime)],[maxSite_dim,length(iTime)]);
        nc.sltt = ncread(SDC_OpenDAP_data_url,'SLTT',[1,min(iTime)],[maxSite_dim,length(iTime)]);
        nc.slnt = ncread(SDC_OpenDAP_data_url,'SLNT',[1,min(iTime)],[maxSite_dim,length(iTime)]);
        nc.scdr = ncread(SDC_OpenDAP_data_url,'SCDR',[1,1,min(iTime)],[length(siteCode),maxSite_dim,length(iTime)]);
        nc.scdt = ncread(SDC_OpenDAP_data_url,'SCDT',[1,1,min(iTime)],[length(siteCode),maxSite_dim,length(iTime)]);
        
        % QC variables
        nc.time_seadatanet_qc = ncread(SDC_OpenDAP_data_url,'TIME_QC',min(iTime),length(iTime));
        nc.position_seadatanet_qc = ncread(SDC_OpenDAP_data_url,'POSITION_QC',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.depth_seadatanet_qc = ncread(SDC_OpenDAP_data_url,'DEPH_QC',min(iTime),length(iTime));
        
        nc.qcflag = ncread(SDC_OpenDAP_data_url,'QCflag',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.owtr_qc = ncread(SDC_OpenDAP_data_url,'OWTR_QC',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.mdfl_qc = ncread(SDC_OpenDAP_data_url,'MDFL_QC',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.vart_qc = ncread(SDC_OpenDAP_data_url,'VART_QC',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.cspd_qc = ncread(SDC_OpenDAP_data_url,'CSPD_QC',[1,1,1,min(iTime)],[length(nc.rnge),length(nc.bear),length(nc.depth),length(iTime)]);
        nc.avrb_qc = ncread(SDC_OpenDAP_data_url,'AVRB_QC',min(iTime),length(iTime));
        nc.rdct_qc = ncread(SDC_OpenDAP_data_url,'RDCT_QC',min(iTime),length(iTime));
    else
        disp(['[' datestr(now) '] - - No data available for the selected period.']);
        disp(['[' datestr(now) '] - - ' 'radial_netCDF_aggregation_v211.m successfully executed.']);
        rnA_err = 1;
        return
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
    return
end
%%

%% Map QC variables to the SDC schema

try
    nc.time_seadatanet_qc(nc.time_seadatanet_qc==1) = int8(49);
    nc.position_seadatanet_qc(nc.position_seadatanet_qc==1) = int8(49);
    nc.position_seadatanet_qc(isnan(nc.position_seadatanet_qc)) = int8(57);
    nc.depth_seadatanet_qc(nc.depth_seadatanet_qc==1) = int8(49);
    
    nc.qcflag(nc.qcflag==0) = int8(48);
    nc.qcflag(nc.qcflag==1) = int8(49);
    nc.qcflag(nc.qcflag==2) = int8(50);
    nc.qcflag(nc.qcflag==3) = int8(51);
    nc.qcflag(nc.qcflag==4) = int8(52);
    nc.qcflag(nc.qcflag==8) = int8(56);
    nc.qcflag(isnan(nc.qcflag)) = int8(57);
    
    nc.owtr_qc(nc.owtr_qc==0) = int8(48);
    nc.owtr_qc(nc.owtr_qc==1) = int8(49);
    nc.owtr_qc(nc.owtr_qc==2) = int8(50);
    nc.owtr_qc(nc.owtr_qc==3) = int8(51);
    nc.owtr_qc(nc.owtr_qc==4) = int8(52);
    nc.owtr_qc(nc.owtr_qc==8) = int8(56);
    nc.owtr_qc(isnan(nc.owtr_qc)) = int8(57);
    
    nc.mdfl_qc(nc.mdfl_qc==0) = int8(48);
    nc.mdfl_qc(nc.mdfl_qc==1) = int8(49);
    nc.mdfl_qc(nc.mdfl_qc==2) = int8(50);
    nc.mdfl_qc(nc.mdfl_qc==3) = int8(51);
    nc.mdfl_qc(nc.mdfl_qc==4) = int8(52);
    nc.mdfl_qc(nc.mdfl_qc==8) = int8(56);
    nc.mdfl_qc(isnan(nc.mdfl_qc)) = int8(57);
    
    nc.vart_qc(nc.vart_qc==0) = int8(48);
    nc.vart_qc(nc.vart_qc==1) = int8(49);
    nc.vart_qc(nc.vart_qc==2) = int8(50);
    nc.vart_qc(nc.vart_qc==3) = int8(51);
    nc.vart_qc(nc.vart_qc==4) = int8(52);
    nc.vart_qc(nc.vart_qc==8) = int8(56);
    nc.vart_qc(isnan(nc.vart_qc)) = int8(57);
    
    nc.cspd_qc(nc.cspd_qc==0) = int8(48);
    nc.cspd_qc(nc.cspd_qc==1) = int8(49);
    nc.cspd_qc(nc.cspd_qc==2) = int8(50);
    nc.cspd_qc(nc.cspd_qc==3) = int8(51);
    nc.cspd_qc(nc.cspd_qc==4) = int8(52);
    nc.cspd_qc(nc.cspd_qc==8) = int8(56);
    nc.cspd_qc(isnan(nc.cspd_qc)) = int8(57);
    
    nc.avrb_qc(nc.avrb_qc==0) = int8(48);
    nc.avrb_qc(nc.avrb_qc==1) = int8(49);
    nc.avrb_qc(nc.avrb_qc==2) = int8(50);
    nc.avrb_qc(nc.avrb_qc==3) = int8(51);
    nc.avrb_qc(nc.avrb_qc==4) = int8(52);
    nc.avrb_qc(nc.avrb_qc==8) = int8(56);
    nc.avrb_qc(isnan(nc.avrb_qc)) = int8(57);
    
    nc.rdct_qc(nc.rdct_qc==0) = int8(48);
    nc.rdct_qc(nc.rdct_qc==1) = int8(49);
    nc.rdct_qc(nc.rdct_qc==2) = int8(50);
    nc.rdct_qc(nc.rdct_qc==3) = int8(51);
    nc.rdct_qc(nc.rdct_qc==4) = int8(52);
    nc.rdct_qc(nc.rdct_qc==8) = int8(56);
    nc.rdct_qc(isnan(nc.rdct_qc)) = int8(57);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

%%

%% Build structures containing QC tests parameters

try
    maxspd_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_velocity_threshold'))));
    Radial_QC_params.VelThr = stationData{maxspd_Index};
    var_thr_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_variance_threshold'))));
    Radial_QC_params.VarThr = stationData{var_thr_Index};
    temp_der_thr_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_temporal_derivative_threshold'))));
    Radial_QC_params.TempDerThr.threshold = stationData{temp_der_thr_Index};
    med_filt_RCLim_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_median_filter_RCLim'))));
    med_filt_AngLim_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_median_filter_AngLim'))));
    med_filt_CurLim_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_median_filter_CurLim'))));
    Radial_QC_params.MedFilt = [stationData{med_filt_RCLim_Index},stationData{med_filt_AngLim_Index},stationData{med_filt_CurLim_Index}];
    avg_rad_bear_min_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_average_radial_bearing_min'))));
    avg_rad_bear_max_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_average_radial_bearing_max'))));
    Radial_QC_params.AvgRadBear = [stationData{avg_rad_bear_min_Index},stationData{avg_rad_bear_max_Index}];
    rad_cnt_Index = find(not(cellfun('isempty', strfind(stationFields, 'radial_QC_radial_count_threshold'))));
    Radial_QC_params.RadCnt = stationData{rad_cnt_Index};
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

%%

%% Set physical dimensions
try
    time_dim = length(iTime);
    bear_dim = length(nc.bear);
    rnge_dim = length(nc.rnge);
    depth_dim = 1;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

%%

%% Retrieve temporal and spatial metadata

% Set time reference
try
    timeref = datenum(1950,1,1);
    time_units = ['days since ' datestr(timeref, 'yyyy-mm-dd') 'T' datestr(timeref, 'HH:MM:SS') 'Z'];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

% Set ADCC compliant data creation, coverage, resolution and duration times
try
    % File creation datetime
    dateCreated = [datestr(now, 'yyyy-mm-dd') 'T' datestr(now, 'HH:MM:SS') 'Z'];
    % Data coverage period
    temporal_resolutionIndex = find(not(cellfun('isempty', strfind(networkFields, 'temporal_resolution'))));
    temporal_resolution = networkData{temporal_resolutionIndex};
    coverageStart = addtodate(nc.time(1), -temporal_resolution/2, 'minute');
    timeCoverageStart = [datestr(coverageStart, 'yyyy-mm-dd') 'T' datestr(coverageStart, 'HH:MM:SS') 'Z'];
    coverageEnd = addtodate(nc.time(length(nc.time)), temporal_resolution/2, 'minute');
    timeCoverageEnd = [datestr(coverageEnd, 'yyyy-mm-dd') 'T' datestr(coverageEnd, 'HH:MM:SS') 'Z'];
    % Temporal resolution
    resolutionMinutes = minutes(temporal_resolution);
    [resH,resM,resS] = hms(resolutionMinutes);
    timeCoverageResolution = 'PT';
    if(resH~=0)
        timeCoverageResolution = [timeCoverageResolution num2str(resH) 'H'];
    end
    if(resM~=0)
        timeCoverageResolution = [timeCoverageResolution num2str(resM) 'M'];
    end
    if(resS~=0)
        timeCoverageResolution = [timeCoverageResolution num2str(resS) 'S'];
    end
    % Temporal duration
    t1=datetime(datevec(nc.time(1)));
    t2=datetime(datevec(addtodate(nc.time(length(nc.time)), temporal_resolution, 'minute')));
    durationVec = datevec(between(t1,t2));
    timeCoverageDuration = 'P';
    if(durationVec(1)~=0)
        timeCoverageDuration = [timeCoverageDuration num2str(durationVec(1)) 'Y'];
    end
    if(durationVec(2)~=0)
        timeCoverageDuration = [timeCoverageDuration num2str(durationVec(2)) 'M'];
    end
    if(durationVec(3)~=0)
        timeCoverageDuration = [timeCoverageDuration num2str(durationVec(3)) 'D'];
    end
    timeCoverageDuration = [timeCoverageDuration 'T'];
    if(durationVec(4)~=0)
        timeCoverageDuration = [timeCoverageDuration num2str(durationVec(4)) 'H'];
    end
    if(durationVec(5)~=0)
        timeCoverageDuration = [timeCoverageDuration num2str(durationVec(5)) 'M'];
    end
    if(durationVec(6)~=0)
        timeCoverageDuration = [timeCoverageDuration num2str(durationVec(6)) 'S'];
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

% Latitude and longitude resolution (degrees)
try
    lonArr = unique(nc.longitude);
    lonArr = sort(lonArr);
    latArr = unique(nc.latitude);
    latArr = sort(latArr);
    latDiff = diff(latArr);
    lonDiff = diff(lonArr);
    latRes = abs(mean(latDiff));
    lonRes = abs(mean(lonDiff));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

% Evaluate measurement vertical max and resolution
try
    transmit_central_frequencyIndex = find(not(cellfun('isempty', strfind(stationFields, 'transmit_central_frequency'))));
    txFreq = stationData{transmit_central_frequencyIndex}*1e6; % transmit frequency in Hertz
    vertMax = (3e8)/(8*pi*txFreq);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

% Set nc output file name
try
    startVec = datevec(nc.time(1));
    endVec = datevec(nc.time(length(nc.time)));
    % Check if the aggregation period is exactly 1 year
    if(strcmp(timeCoverageDuration,'P1YT'))
        time_str = sprintf('%.4d',startVec(1));
    elseif(strcmp(timeCoverageDuration,'P1MT'))
        time_str = sprintf('%.4d%.2d',startVec(1),startVec(2));
    else
        time_str = sprintf('%.4d%.2d%.2d_%.4d%.2d%.2d',startVec(1),startVec(2),startVec(3),endVec(1),endVec(2),endVec(3));
    end
    % Set path name and create the output folder, if needed
    outputPathIndex = find(not(cellfun('isempty', strfind(stationFields, 'SDC_folder_path'))));
    ncFilePath = [strtrim(stationData{outputPathIndex}) filesep stationData{station_idIndex}];
    if (exist(ncFilePath, 'dir') ~= 7)
        mkdir(ncFilePath);
    end
    ncFileNoPath = ['RV_HF_' platform_code '_' time_str '.nc'];
    ncfile = [ncFilePath filesep ncFileNoPath];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

% Set citation string and distribution string
try
    citation_statementIndex = find(not(cellfun('isempty', strfind(networkFields, 'citation_statement'))));
    citation_str = ['These data were collected and made freely available by the SeaDataNet project and the programs that contribute to it. ' networkData{citation_statementIndex}];
    distribution_str = 'These data follow SeaDataNet standards; they are public and free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.';
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

% Set naming authority
try
    institution_websiteIndex = find(not(cellfun('isempty', strfind(stationFields, 'institution_website'))));
    institution_websiteStr = stationData{institution_websiteIndex};
%     if(~isempty(strfind(institution_websiteStr,'http://')))
%         tmpStr = strrep(institution_websiteStr,'http://','');
%     elseif(~isempty(strfind(institution_websiteStr,'https://')))
%         tmpStr = strrep(institution_websiteStr,'https://','');
%     else
%         tmpStr = institution_websiteStr;
%     end
%     tmpStr = strrep(tmpStr,'www.','');
%     tmpStr = strrep(tmpStr,'/','');
%     splitStr = strsplit(tmpStr,'.');
%     naming_authorityStr = [];
%     for split_idx=length(splitStr):-1:1
%         naming_authorityStr = [naming_authorityStr splitStr{split_idx}];
%         if(split_idx~=1)
%             naming_authorityStr = [naming_authorityStr '.'];
%         end
%     end
%     naming_authorityStr= naming_authorityStr(~isspace(naming_authorityStr));
    naming_authorityStr = 'eu.eurogoos';
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

% Set collection time
try
    time_coll{1} = [datestr(startVec, 'yyyy-mm-dd') 'T' datestr(startVec, 'HH:MM:SS') 'Z'];
    time_coll{2} = [datestr(endVec, 'yyyy-mm-dd') 'T' datestr(endVec, 'HH:MM:SS') 'Z'];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

try
    % Define id and metadata resources
    dataID = [EDIOS_Series_ID '-' EDIOS_Platform_ID '_' time_coll{1} '_' time_coll{2}];
    metadata_pageIndex = find(not(cellfun('isempty', strfind(networkFields, 'metadata_page'))));
    TDS_catalog = networkData{metadata_pageIndex};
    xlink = ['<sdn_reference xlink:href="' TDS_catalog '" xlink:role="" xlink:type="URL"/>'];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

% Deletes the eventually present netCDF file with the same name
try
    delete(ncfile);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

%%

%% Creates vars with their dimensions

% Set netcdf format
try
    ncfmt = 'netcdf4_classic';
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

try
    nccreate(ncfile,'TIME',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','double',...
        'Format',ncfmt);
    
    nccreate(ncfile,'BEAR',...
        'Dimensions',{'BEAR',bear_dim},...
        'Datatype','single',...
        'Format',ncfmt);
    
    nccreate(ncfile,'RNGE',...
        'Dimensions',{'RNGE',rnge_dim},...
        'Datatype','single',...
        'Format',ncfmt);
    
    nccreate(ncfile,'DEPTH',...
        'Dimensions',{'DEPTH',depth_dim},...
        'Datatype','single',...
        'Format',ncfmt);
    
    nccreate(ncfile,'LATITUDE',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim},...
        'Datatype','single',...
        'Format',ncfmt);
    
    nccreate(ncfile,'LONGITUDE',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim},...
        'Datatype','single',...
        'Format',ncfmt);
    
    nccreate(ncfile,'crs',...
        'Datatype','int16',...
        'Format',ncfmt);
    
    nccreate(ncfile,'SDN_CRUISE',...
        'Dimensions',{'STRING50',string50_dim},...
        'Datatype','char',...
        'Format',ncfmt);
    
    nccreate(ncfile,'SDN_STATION',...
        'Dimensions',{'STRING50',string50_dim},...
        'Datatype','char',...
        'Format',ncfmt);
    
    nccreate(ncfile,'SDN_LOCAL_CDI_ID',...
        'Dimensions',{'STRING80',string80_dim},...
        'Datatype','char',...
        'Format',ncfmt);
    
    nccreate(ncfile,'SDN_EDMO_CODE',...
        'Dimensions',{'MAXINST',maxInst_dim},...
        'Datatype','int16',...
        'Format',ncfmt);
    
    nccreate(ncfile,'SDN_REFERENCES',...
        'Dimensions',{'STRING250', string250_dim},...
        'Datatype','char',...
        'Format',ncfmt);
    
    nccreate(ncfile,'SDN_XLINK',...
        'Dimensions',{'STRING250',string250_dim, 'REFMAX',refMax_dim},...
        'Datatype','char',...
        'Format',ncfmt);
    
    nccreate(ncfile,'RDVA',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue', netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'DRVA',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int32',...
        'FillValue', netcdf.getConstant('NC_FILL_INT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'EWCT',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue', netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'NSCT',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'ESPC',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue', netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'ETMP',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'MAXV',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue', netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'MINV',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'ERSC',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue', netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'ERTC',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'XDST',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim},...
        'Datatype','int32',...
        'FillValue', netcdf.getConstant('NC_FILL_INT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'YDST',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim},...
        'Datatype','int32',...
        'FillValue',netcdf.getConstant('NC_FILL_INT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SPRC',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue', netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'NARX',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'NATX',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',netcdf.getConstant('NC_FILL_BYTE'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SLTR',...
        'Dimensions',{'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','int32',...
        'FillValue',netcdf.getConstant('NC_FILL_INT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SLNR',...
        'Dimensions',{'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','int32',...
        'FillValue',netcdf.getConstant('NC_FILL_INT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SLTT',...
        'Dimensions',{'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','int32',...
        'FillValue',netcdf.getConstant('NC_FILL_INT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SLNT',...
        'Dimensions',{'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','int32',...
        'FillValue',netcdf.getConstant('NC_FILL_INT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SCDR',...
        'Dimensions',{'STRING20',string20_dim,'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','char',...
        'FillValue',netcdf.getConstant('NC_FILL_CHAR'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SCDT',...
        'Dimensions',{'STRING20',string20_dim,'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','char',...
        'FillValue',netcdf.getConstant('NC_FILL_CHAR'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'TIME_SEADATANET_QC',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'POSITION_SEADATANET_QC',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'DEPTH_SEADATANET_QC',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'QCflag',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'OWTR_QC',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'MDFL_QC',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'VART_QC',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'CSPD_QC',...
        'Dimensions',{'RNGE',rnge_dim,'BEAR',bear_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'AVRB_QC',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'RDCT_QC',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    %%
    
    %% Creates attributes for the variables
    ncwriteatt(ncfile,'TIME','long_name',char('Chronological Julian Date'));
    ncwriteatt(ncfile,'TIME','standard_name',char('time'));
    ncwriteatt(ncfile,'TIME','units',char(time_units));
    ncwriteatt(ncfile,'TIME','calendar',char('julian'));
    ncwriteatt(ncfile,'TIME','axis',char('T'));
    ncwriteatt(ncfile,'TIME','sdn_parameter_name',char('Elapsed time (since 1950-01-01T00:00:00Z)'));
    ncwriteatt(ncfile,'TIME','sdn_parameter_urn',char('SDN:P01::ELTJLD01'));
    ncwriteatt(ncfile,'TIME','sdn_uom_name',char('Days'));
    ncwriteatt(ncfile,'TIME','sdn_uom_urn',char('SDN:P06::UTAA'));
    ncwriteatt(ncfile,'TIME','ancillary_variables',char('TIME_SEADATANET_QC'));
    
    ncwriteatt(ncfile,'BEAR','long_name',char('Bearing away from instrument'));
    ncwriteatt(ncfile,'BEAR','units',char('degrees_true'));
    ncwriteatt(ncfile,'BEAR','axis',char('X'));
    ncwriteatt(ncfile,'BEAR','sdn_parameter_name',char('Bearing'));
    ncwriteatt(ncfile,'BEAR','sdn_parameter_urn',char('SDN:P01::BEARRFTR'));
    ncwriteatt(ncfile,'BEAR','sdn_uom_name',char('Degrees true'));
    ncwriteatt(ncfile,'BEAR','sdn_uom_urn',char('SDN:P06::UABB'));
    ncwriteatt(ncfile,'BEAR','ancillary_variables',char('POSITION_SEADATANET_QC'));
    
    ncwriteatt(ncfile,'RNGE','long_name',char('Range away from instrument'));
    ncwriteatt(ncfile,'RNGE','units',char('km'));
    ncwriteatt(ncfile,'RNGE','axis',char('Y'));
    ncwriteatt(ncfile,'RNGE','sdn_parameter_name',char('Range (from fixed reference point) by unspecified GPS system'));
    ncwriteatt(ncfile,'RNGE','sdn_parameter_urn',char('SDN:P01::RIFNAX01'));
    ncwriteatt(ncfile,'RNGE','sdn_uom_name',char('Kilometres'));
    ncwriteatt(ncfile,'RNGE','sdn_uom_urn',char('SDN:P06::ULKM'));
    ncwriteatt(ncfile,'RNGE','ancillary_variables',char('POSITION_SEADATANET_QC'));
    
    ncwriteatt(ncfile,'DEPTH','long_name',char('Depth'));
    ncwriteatt(ncfile,'DEPTH','standard_name',char('depth'));
    ncwriteatt(ncfile,'DEPTH','units',char('m'));
    ncwriteatt(ncfile,'DEPTH','axis',char('Z'));
    ncwriteatt(ncfile,'DEPTH','positive',char('down'));
    ncwriteatt(ncfile,'DEPTH','reference',char('sea_level'));
    ncwriteatt(ncfile,'DEPTH','sdn_parameter_name',char('Depth below surface of the water body'));
    ncwriteatt(ncfile,'DEPTH','sdn_parameter_urn',char('SDN:P01::ADEPZZ01'));
    ncwriteatt(ncfile,'DEPTH','sdn_uom_name',char('Metres'));
    ncwriteatt(ncfile,'DEPTH','sdn_uom_urn',char('SDN:P06::ULAA'));
    ncwriteatt(ncfile,'DEPTH','ancillary_variables',char('DEPTH_SEADATANET_QC'));
    
    ncwriteatt(ncfile,'LATITUDE','long_name',char('Latitude'));
    ncwriteatt(ncfile,'LATITUDE','standard_name',char('latitude'));
    ncwriteatt(ncfile,'LATITUDE','units',char('degrees_north'));
    ncwriteatt(ncfile,'LATITUDE','valid_range',single( [-90, 90] ));
    ncwriteatt(ncfile,'LATITUDE','sdn_parameter_name',char('Latitude north'));
    ncwriteatt(ncfile,'LATITUDE','sdn_parameter_urn',char('SDN:P01::ALATZZ01'));
    ncwriteatt(ncfile,'LATITUDE','sdn_uom_name',char('Degrees north'));
    ncwriteatt(ncfile,'LATITUDE','sdn_uom_urn',char('SDN:P06::DEGN'));
    ncwriteatt(ncfile,'LATITUDE','grid_mapping',char('crs'));
    ncwriteatt(ncfile,'LATITUDE','ancillary_variables',char('POSITION_SEADATANET_QC'));
    
    ncwriteatt(ncfile,'LONGITUDE','long_name',char('Longitude'));
    ncwriteatt(ncfile,'LONGITUDE','standard_name',char('longitude'));
    ncwriteatt(ncfile,'LONGITUDE','units',char('degrees_east'));
    ncwriteatt(ncfile,'LONGITUDE','valid_range',single( [-180, 180] ));
    ncwriteatt(ncfile,'LONGITUDE','sdn_parameter_name',char('Longitude east'));
    ncwriteatt(ncfile,'LONGITUDE','sdn_parameter_urn',char('SDN:P01::ALONZZ01'));
    ncwriteatt(ncfile,'LONGITUDE','sdn_uom_name',char('Degrees east'));
    ncwriteatt(ncfile,'LONGITUDE','sdn_uom_urn',char('SDN:P06::DEGE'));
    ncwriteatt(ncfile,'LONGITUDE','grid_mapping',char('crs'));
    ncwriteatt(ncfile,'LONGITUDE','ancillary_variables',char('POSITION_SEADATANET_QC'));
    
    ncwriteatt(ncfile,'crs','grid_mapping_name',char('latitude_longitude'));
    ncwriteatt(ncfile,'crs','epsg_code',char('EPSG:4326'));
    ncwriteatt(ncfile,'crs','semi_major_axis',6378137.0);
    ncwriteatt(ncfile,'crs','inverse_flattening',298.257223563);
    
    ncwriteatt(ncfile,'SDN_CRUISE','long_name',char('Grid grouping label'));
    
    ncwriteatt(ncfile,'SDN_STATION','long_name',char('Grid label'));
    
    ncwriteatt(ncfile,'SDN_LOCAL_CDI_ID','long_name',char('SeaDataNet CDI identifier'));
    ncwriteatt(ncfile,'SDN_LOCAL_CDI_ID','cf_role',char('timeseries_id'));
    
    ncwriteatt(ncfile,'SDN_EDMO_CODE','long_name',char('European Directory of Marine Organisations code for the CDI partner'));
    ncwriteatt(ncfile,'SDN_EDMO_CODE','units',char('1'));
    
    ncwriteatt(ncfile,'SDN_REFERENCES','long_name',char('Usage metadata reference'));
    
    ncwriteatt(ncfile,'SDN_XLINK','long_name',char('External resource linkages'));
    
    ncwriteatt(ncfile,'RDVA','long_name',char('Radial Sea Water Velocity Away From Instrument'));
    ncwriteatt(ncfile,'RDVA','standard_name',char('radial_sea_water_velocity_away_from_instrument'));
    ncwriteatt(ncfile,'RDVA','units',char('m s-1'));
    ncwriteatt(ncfile,'RDVA','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'RDVA','add_offset',double(addOffset));
    ncwriteatt(ncfile,'RDVA','sdn_parameter_name',char('Current speed (Eulerian) in the water body by directional range-gated radar'));
    ncwriteatt(ncfile,'RDVA','sdn_parameter_urn',char('SDN:P01::LCSAWVRD'));
    ncwriteatt(ncfile,'RDVA','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'RDVA','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'RDVA','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'RDVA','valid_range',int16([(-10-addOffset)./scaleFactor, (10-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'RDVA','ancillary_variables',char('QCflag, OWTR_QC, MDFL_QC, CSPD_QC, RDCT_QC'));
    
    ncwriteatt(ncfile,'DRVA','long_name',char('Direction of Radial Vector Away From Instrument'));
    ncwriteatt(ncfile,'DRVA','standard_name',char('direction_of_radial_vector_away_from_instrument'));
    ncwriteatt(ncfile,'DRVA','units',char('degrees_true'));
    ncwriteatt(ncfile,'DRVA','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'DRVA','add_offset',double(addOffset));
    ncwriteatt(ncfile,'DRVA','sdn_parameter_name',char('Current direction (Eulerian) in the water body by directional range-gated radar'));
    ncwriteatt(ncfile,'DRVA','sdn_parameter_urn',char('SDN:P01::LCDAWVRD'));
    ncwriteatt(ncfile,'DRVA','sdn_uom_name',char('Degrees True'));
    ncwriteatt(ncfile,'DRVA','sdn_uom_urn',char('SDN:P06::UABB'));
    ncwriteatt(ncfile,'DRVA','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'DRVA','valid_range',int32([0, (360-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'DRVA','ancillary_variables',char('QCflag, OWTR_QC, MDFL_QC, AVRB_QC, RDCT_QC'));
    
    ncwriteatt(ncfile,'EWCT','long_name',char('West-east current component'));
    ncwriteatt(ncfile,'EWCT','standard_name',char('eastward_sea_water_velocity'));
    ncwriteatt(ncfile,'EWCT','units',char('m s-1'));
    ncwriteatt(ncfile,'EWCT','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'EWCT','add_offset',double(addOffset));
    ncwriteatt(ncfile,'EWCT','sdn_parameter_name',char('Eastward current velocity in the water body'));
    ncwriteatt(ncfile,'EWCT','sdn_parameter_urn',char('SDN:P01::LCEWZZ01'));
    ncwriteatt(ncfile,'EWCT','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'EWCT','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'EWCT','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'EWCT','valid_range',int16([(-10-addOffset)./scaleFactor, (10-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'EWCT','ancillary_variables',char('QCflag, OWTR_QC, MDFL_QC, VART_QC, CSPD_QC, AVRB_QC, RDCT_QC'));
    
    ncwriteatt(ncfile,'NSCT','long_name',char('South-north current component'));
    ncwriteatt(ncfile,'NSCT','standard_name',char('northward_sea_water_velocity'));
    ncwriteatt(ncfile,'NSCT','units',char('m s-1'));
    ncwriteatt(ncfile,'NSCT','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'NSCT','add_offset',double(addOffset));
    ncwriteatt(ncfile,'NSCT','sdn_parameter_name',char('Northward current velocity in the water body'));
    ncwriteatt(ncfile,'NSCT','sdn_parameter_urn',char('SDN:P01::LCNSZZ01'));
    ncwriteatt(ncfile,'NSCT','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'NSCT','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'NSCT','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'NSCT','valid_range',int16([(-10-addOffset)./scaleFactor, (10-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'NSCT','ancillary_variables',char('QCflag, OWTR_QC, MDFL_QC, VART_QC, CSPD_QC, AVRB_QC, RDCT_QC'));
    
    ncwriteatt(ncfile,'ESPC','long_name',char('Radial Standard Deviation of Current Velocity over the Scatter Patch'));
    ncwriteatt(ncfile,'ESPC','units',char('m s-1'));
    ncwriteatt(ncfile,'ESPC','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'ESPC','add_offset',double(addOffset));
    ncwriteatt(ncfile,'ESPC','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'ESPC','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'ESPC','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'ESPC','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'ESPC','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'ESPC','valid_range',int16([(-32-addOffset)./scaleFactor, (32-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'ESPC','ancillary_variables',char('QCflag, VART_QC'));
    
    ncwriteatt(ncfile,'ETMP','long_name',char('Radial Standard Deviation of Current Velocity over Coverage Period'));
    ncwriteatt(ncfile,'ETMP','units',char('m s-1'));
    ncwriteatt(ncfile,'ETMP','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'ETMP','add_offset',double(addOffset));
    ncwriteatt(ncfile,'ETMP','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'ETMP','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'ETMP','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'ETMP','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'ETMP','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'ETMP','valid_range',int16([(-32-addOffset)./scaleFactor, (32-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'ETMP','ancillary_variables',char('QCflag, VART_QC'));
    
    ncwriteatt(ncfile,'MINV','long_name',char('Radial Sea Water Velocity Away From Instrument Minimum'));
    ncwriteatt(ncfile,'MINV','units',char('m s-1'));
    ncwriteatt(ncfile,'MINV','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'MINV','add_offset',double(addOffset));
    ncwriteatt(ncfile,'MINV','sdn_parameter_name',char('Current speed (Eulerian) in the water body by directional range-gated radar'));
    ncwriteatt(ncfile,'MINV','sdn_parameter_urn',char('SDN:P01::LCSAWVRD'));
    ncwriteatt(ncfile,'MINV','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'MINV','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'MINV','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'MINV','valid_range',int16([(-10-addOffset)./scaleFactor, (10-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'MINV','ancillary_variables',char('QCflag, MDFL_QC, VART_QC, CSPD_QC'));
    
    ncwriteatt(ncfile,'MAXV','long_name',char('Radial Sea Water Velocity Away From Instrument Maximum'));
    ncwriteatt(ncfile,'MAXV','units',char('m s-1'));
    ncwriteatt(ncfile,'MAXV','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'MAXV','add_offset',double(addOffset));
    ncwriteatt(ncfile,'MAXV','sdn_parameter_name',char('Current speed (Eulerian) in the water body by directional range-gated radar'));
    ncwriteatt(ncfile,'MAXV','sdn_parameter_urn',char('SDN:P01::LCSAWVRD'));
    ncwriteatt(ncfile,'MAXV','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'MAXV','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'MAXV','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'MAXV','valid_range',int16([(-10-addOffset)./scaleFactor, (10-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'MAXV','ancillary_variables',char('QCflag, MDFL_QC, VART_QC, CSPD_QC'));
    
    ncwriteatt(ncfile,'ERSC','long_name',char('Radial Sea Water Velocity Spatial Quality Count'));
    ncwriteatt(ncfile,'ERSC','units',char('1'));
    ncwriteatt(ncfile,'ERSC','scale_factor',double(1));
    ncwriteatt(ncfile,'ERSC','add_offset',double(0));
    ncwriteatt(ncfile,'ERSC','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'ERSC','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'ERSC','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'ERSC','sdn_uom_urn',char('SDN:P06::UUUU'));
    ncwriteatt(ncfile,'ERSC','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'ERSC','valid_range',int16([0, 127]));
    ncwriteatt(ncfile,'ERSC','ancillary_variables',char('QCflag'));
    
    ncwriteatt(ncfile,'ERTC','long_name',char('Radial Sea Water Velocity Temporal Quality Count'));
    ncwriteatt(ncfile,'ERTC','units',char('1'));
    ncwriteatt(ncfile,'ERTC','scale_factor',double(1));
    ncwriteatt(ncfile,'ERTC','add_offset',double(0));
    ncwriteatt(ncfile,'ERTC','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'ERTC','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'ERTC','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'ERTC','sdn_uom_urn',char('SDN:P06::UUUU'));
    ncwriteatt(ncfile,'ERTC','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'ERTC','valid_range',int16([0, 127]));
    ncwriteatt(ncfile,'ERTC','ancillary_variables',char('QCflag'));
    
    ncwriteatt(ncfile,'XDST','long_name',char('Eastward Distance From Instrument'));
    ncwriteatt(ncfile,'XDST','units',char('km'));
    ncwriteatt(ncfile,'XDST','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'XDST','add_offset',double(addOffset));
    ncwriteatt(ncfile,'XDST','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'XDST','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'XDST','sdn_uom_name',char('Kilometres'));
    ncwriteatt(ncfile,'XDST','sdn_uom_urn',char('SDN:P06::ULKM'));
    ncwriteatt(ncfile,'XDST','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'XDST','valid_range',int32([0, (1e3-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'XDST','ancillary_variables',char('QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC'));
    
    ncwriteatt(ncfile,'YDST','long_name',char('Northward Distance From Instrument'));
    ncwriteatt(ncfile,'YDST','units',char('km'));
    ncwriteatt(ncfile,'YDST','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'YDST','add_offset',double(addOffset));
    ncwriteatt(ncfile,'YDST','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'YDST','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'YDST','sdn_uom_name',char('Kilometres'));
    ncwriteatt(ncfile,'YDST','sdn_uom_urn',char('SDN:P06::ULKM'));
    ncwriteatt(ncfile,'YDST','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'YDST','valid_range',int32([0, (1e3-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'YDST','ancillary_variables',char('QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC'));
    
    ncwriteatt(ncfile,'SPRC','long_name',char('Radial Sea Water Velocity Cross Spectra Range Cell'));
    ncwriteatt(ncfile,'SPRC','units',char('1'));
    ncwriteatt(ncfile,'SPRC','scale_factor',double(1));
    ncwriteatt(ncfile,'SPRC','add_offset',double(0));
    ncwriteatt(ncfile,'SPRC','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'SPRC','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'SPRC','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'SPRC','sdn_uom_urn',char('SDN:P06::UUUU'));
    ncwriteatt(ncfile,'SPRC','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    ncwriteatt(ncfile,'SPRC','valid_range',int16([0, 127]));
    ncwriteatt(ncfile,'SPRC','ancillary_variables',char('QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC'));
    
    ncwriteatt(ncfile,'NARX','long_name',char('Number of Receive Antennas'));
    ncwriteatt(ncfile,'NARX','units',char('1'));
    ncwriteatt(ncfile,'NARX','valid_range',int8([0 maxSite_dim]));
    ncwriteatt(ncfile,'NARX','scale_factor',int8(1));
    ncwriteatt(ncfile,'NARX','add_offset',int8(0));
    ncwriteatt(ncfile,'NARX','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'NARX','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'NARX','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'NARX','sdn_uom_urn',char('SDN:P06::UUUU'));
    
    ncwriteatt(ncfile,'NATX','long_name',char('Number of Transmit Antennas'));
    ncwriteatt(ncfile,'NATX','units',char('1'));
    ncwriteatt(ncfile,'NATX','valid_range',int8([0 maxSite_dim]));
    ncwriteatt(ncfile,'NATX','scale_factor',int8(1));
    ncwriteatt(ncfile,'NATX','add_offset',int8(0));
    ncwriteatt(ncfile,'NATX','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'NATX','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'NATX','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'NATX','sdn_uom_urn',char('SDN:P06::UUUU'));
    
    ncwriteatt(ncfile,'SLTR','long_name',char('Receive Antenna Latitudes'));
    ncwriteatt(ncfile,'SLTR','standard_name',char('latitude'));
    ncwriteatt(ncfile,'SLTR','units','degrees_north');
    ncwriteatt(ncfile,'SLTR','valid_range',int32( [(-90-addOffset)./scaleFactor (90-addOffset)./scaleFactor] ));
    ncwriteatt(ncfile,'SLTR','coordinates',char('TIME MAXSITE'));
    ncwriteatt(ncfile,'SLTR','scale_factor',single(scaleFactor));
    ncwriteatt(ncfile,'SLTR','add_offset',single(addOffset));
    ncwriteatt(ncfile,'SLTR','sdn_parameter_name',char('Latitude north'));
    ncwriteatt(ncfile,'SLTR','sdn_parameter_urn',char('SDN:P01::ALATZZ01'));
    ncwriteatt(ncfile,'SLTR','sdn_uom_name',char('Degrees north'));
    ncwriteatt(ncfile,'SLTR','sdn_uom_urn',char('SDN:P06::DEGN'));
    
    ncwriteatt(ncfile,'SLNR','long_name',char('Receive Antenna Longitudes'));
    ncwriteatt(ncfile,'SLNR','standard_name',char('longitude'));
    ncwriteatt(ncfile,'SLNR','units','degrees_east');
    ncwriteatt(ncfile,'SLNR','valid_range',int32( [(-180-addOffset)./scaleFactor (180-addOffset)./scaleFactor] ));
    ncwriteatt(ncfile,'SLNR','coordinates',char('TIME MAXSITE'));
    ncwriteatt(ncfile,'SLNR','scale_factor',single(scaleFactor));
    ncwriteatt(ncfile,'SLNR','add_offset',single(addOffset));
    ncwriteatt(ncfile,'SLNR','sdn_parameter_name',char('Longitude east'));
    ncwriteatt(ncfile,'SLNR','sdn_parameter_urn',char('SDN:P01::ALONZZ01'));
    ncwriteatt(ncfile,'SLNR','sdn_uom_name',char('Degrees east'));
    ncwriteatt(ncfile,'SLNR','sdn_uom_urn',char('SDN:P06::DEGE'));
    
    ncwriteatt(ncfile,'SLTT','long_name',char('Transmit Antenna Latitudes'));
    ncwriteatt(ncfile,'SLTT','standard_name',char('latitude'));
    ncwriteatt(ncfile,'SLTT','units','degrees_north');
    ncwriteatt(ncfile,'SLTT','valid_range',int32( [(-90-addOffset)./scaleFactor (90-addOffset)./scaleFactor] ));
    ncwriteatt(ncfile,'SLTT','coordinates',char('TIME MAXSITE'));
    ncwriteatt(ncfile,'SLTT','scale_factor',single(scaleFactor));
    ncwriteatt(ncfile,'SLTT','add_offset',single(addOffset));
    ncwriteatt(ncfile,'SLTT','sdn_parameter_name',char('Latitude north'));
    ncwriteatt(ncfile,'SLTT','sdn_parameter_urn',char('SDN:P01::ALATZZ01'));
    ncwriteatt(ncfile,'SLTT','sdn_uom_name',char('Degrees north'));
    ncwriteatt(ncfile,'SLTT','sdn_uom_urn',char('SDN:P06::DEGN'));
    
    ncwriteatt(ncfile,'SLNT','long_name',char('Transmit Antenna Longitudes'));
    ncwriteatt(ncfile,'SLNT','standard_name',char('longitude'));
    ncwriteatt(ncfile,'SLNT','units','degrees_east');
    ncwriteatt(ncfile,'SLNT','valid_range',int32( [(-180-addOffset)./scaleFactor (180-addOffset)./scaleFactor] ));
    ncwriteatt(ncfile,'SLNT','coordinates',char('TIME MAXSITE'));
    ncwriteatt(ncfile,'SLNT','scale_factor',single(scaleFactor));
    ncwriteatt(ncfile,'SLNT','add_offset',single(addOffset));
    ncwriteatt(ncfile,'SLNT','sdn_parameter_name',char('Longitude east'));
    ncwriteatt(ncfile,'SLNT','sdn_parameter_urn',char('SDN:P01::ALONZZ01'));
    ncwriteatt(ncfile,'SLNT','sdn_uom_name',char('Degrees east'));
    ncwriteatt(ncfile,'SLNT','sdn_uom_urn',char('SDN:P06::DEGE'));
    
    ncwriteatt(ncfile,'SCDR','long_name',char('Receive Antenna Codes'));
    ncwriteatt(ncfile,'SCDR','units',char('1'));
    ncwriteatt(ncfile,'SCDR','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'SCDR','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'SCDR','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'SCDR','sdn_uom_urn',char('SDN:P06::UUUU'));
    
    ncwriteatt(ncfile,'SCDT','long_name',char('Transmit Antenna Codes'));
    ncwriteatt(ncfile,'SCDT','units',char('1'));
    ncwriteatt(ncfile,'SCDT','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'SCDT','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'SCDT','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'SCDT','sdn_uom_urn',char('SDN:P06::UUUU'));
    
    ncwriteatt(ncfile,'TIME_SEADATANET_QC','long_name',char('Time SeaDataNet Quality Flag'));
    ncwriteatt(ncfile,'TIME_SEADATANET_QC','units',char('1'));
    ncwriteatt(ncfile,'TIME_SEADATANET_QC','valid_range',int8([48 65]));
    ncwriteatt(ncfile,'TIME_SEADATANET_QC','flag_values',int8([48 49 50 51 52 53 54 55 56 57 65]));
    ncwriteatt(ncfile,'TIME_SEADATANET_QC','flag_meanings',char('no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain'));
    ncwriteatt(ncfile,'TIME_SEADATANET_QC','Conventions',char('SeaDataNet measurand qualifier flags.'));
    ncwriteatt(ncfile,'TIME_SEADATANET_QC','sdn_conventions_urn',char('SDN:L20::'));
    ncwriteatt(ncfile,'TIME_SEADATANET_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'TIME_SEADATANET_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'POSITION_SEADATANET_QC','long_name',char('Position SeaDataNet Quality Flags'));
    ncwriteatt(ncfile,'POSITION_SEADATANET_QC','units',char('1'));
    ncwriteatt(ncfile,'POSITION_SEADATANET_QC','valid_range',int8([48 65]));
    ncwriteatt(ncfile,'POSITION_SEADATANET_QC','flag_values',int8([48 49 50 51 52 53 54 55 56 57 65]));
    ncwriteatt(ncfile,'POSITION_SEADATANET_QC','flag_meanings',char('no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain'));
    ncwriteatt(ncfile,'POSITION_SEADATANET_QC','Conventions',char('SeaDataNet measurand qualifier flags.'));
    ncwriteatt(ncfile,'POSITION_SEADATANET_QC','sdn_conventions_urn',char('SDN:L20::'));
    ncwriteatt(ncfile,'POSITION_SEADATANET_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'POSITION_SEADATANET_QC','add_offset',int8(0));
    ncwriteatt(ncfile,'POSITION_SEADATANET_QC','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    
    ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','long_name',char('Depth SeaDataNet Quality Flag'));
    ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','units',char('1'));
    ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','valid_range',int8([48 65]));
    ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','flag_values',int8([48 49 50 51 52 53 54 55 56 57 65]));
    ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','flag_meanings',char('no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain'));
    ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','Conventions',char('SeaDataNet measurand qualifier flags.'));
    ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','sdn_conventions_urn',char('SDN:L20::'));
    ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'QCflag','long_name',char('Overall Quality Flags'));
    ncwriteatt(ncfile,'QCflag','units',char('1'));
    ncwriteatt(ncfile,'QCflag','valid_range',int8([48 65]));
    ncwriteatt(ncfile,'QCflag','flag_values',int8([48 49 50 51 52 53 54 55 56 57 65]));
    ncwriteatt(ncfile,'QCflag','flag_meanings',char('no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain'));
    ncwriteatt(ncfile,'QCflag','Conventions',char('SeaDataNet measurand qualifier flags.'));
    ncwriteatt(ncfile,'QCflag','sdn_conventions_urn',char('SDN:L20::'));
    ncwriteatt(ncfile,'QCflag','scale_factor',int8(1));
    ncwriteatt(ncfile,'QCflag','add_offset',int8(0));
    ncwriteatt(ncfile,'QCflag','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    
    ncwriteatt(ncfile,'OWTR_QC','long_name',char('Over-water Quality Flags'));
    ncwriteatt(ncfile,'OWTR_QC','units',char('1'));
    ncwriteatt(ncfile,'OWTR_QC','valid_range',int8([48 65]));
    ncwriteatt(ncfile,'OWTR_QC','flag_values',int8([48 49 50 51 52 53 54 55 56 57 65]));
    ncwriteatt(ncfile,'OWTR_QC','flag_meanings',char('no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain'));
    ncwriteatt(ncfile,'OWTR_QC','Conventions',char('SeaDataNet measurand qualifier flags.'));
    ncwriteatt(ncfile,'OWTR_QC','sdn_conventions_urn',char('SDN:L20::'));
    ncwriteatt(ncfile,'OWTR_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'OWTR_QC','add_offset',int8(0));
    ncwriteatt(ncfile,'OWTR_QC','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    
    ncwriteatt(ncfile,'MDFL_QC','long_name',char('Median Filter Quality Flags'));
    ncwriteatt(ncfile,'MDFL_QC','units',char('1'));
    ncwriteatt(ncfile,'MDFL_QC','valid_range',int8([48 65]));
    ncwriteatt(ncfile,'MDFL_QC','flag_values',int8([48 49 50 51 52 53 54 55 56 57 65]));
    ncwriteatt(ncfile,'MDFL_QC','flag_meanings',char('no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain'));
    ncwriteatt(ncfile,'MDFL_QC','comment', ['Threshold set to ' num2str(Radial_QC_params.MedFilt(1)) ' km, ' num2str(Radial_QC_params.MedFilt(2)) ' deg, ' num2str(Radial_QC_params.MedFilt(3)) ' m/s, ']);
    ncwriteatt(ncfile,'MDFL_QC','Conventions',char('SeaDataNet measurand qualifier flags.'));
    ncwriteatt(ncfile,'MDFL_QC','sdn_conventions_urn',char('SDN:L20::'));
    ncwriteatt(ncfile,'MDFL_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'MDFL_QC','add_offset',int8(0));
    ncwriteatt(ncfile,'MDFL_QC','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    
    ncwriteatt(ncfile,'VART_QC','long_name',char('Variance Threshold Quality Flags'));
    ncwriteatt(ncfile,'VART_QC','units',char('1'));
    ncwriteatt(ncfile,'VART_QC','valid_range',int8([48 65]));
    ncwriteatt(ncfile,'VART_QC','flag_values',int8([48 49 50 51 52 53 54 55 56 57 65]));
    ncwriteatt(ncfile,'VART_QC','flag_meanings',char('no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain'));
    if(contains(sensorATT,'codar','IgnoreCase',true))
        ncwriteatt(ncfile,'VART_QC','comment',char(['Test not applicable to Direction Finding systems. The Temporal Derivative test is applied.' ...
            'Threshold set to ' num2str(Radial_QC_params.TempDerThr.threshold) ' m/s. ']));
    elseif(contains(sensorATT,'wera','IgnoreCase',true))
        ncwriteatt(ncfile,'VART_QC','comment',char(['Threshold set to ' num2str(Radial_QC_params.VarThr.threshold) ' m2/s2. ']));
    end
    ncwriteatt(ncfile,'VART_QC','Conventions',char('SeaDataNet measurand qualifier flags.'));
    ncwriteatt(ncfile,'VART_QC','sdn_conventions_urn',char('SDN:L20::'));
    ncwriteatt(ncfile,'VART_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'VART_QC','add_offset',int8(0));
    ncwriteatt(ncfile,'VART_QC','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    
    ncwriteatt(ncfile,'CSPD_QC','long_name',char('Velocity Threshold Quality Flags'));
    ncwriteatt(ncfile,'CSPD_QC','units',char('1'));
    ncwriteatt(ncfile,'CSPD_QC','valid_range',int8([48 65]));
    ncwriteatt(ncfile,'CSPD_QC','flag_values',int8([48 49 50 51 52 53 54 55 56 57 65]));
    ncwriteatt(ncfile,'CSPD_QC','flag_meanings',char('no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain'));
    ncwriteatt(ncfile,'CSPD_QC','comment',char(['Threshold set to ' num2str(Radial_QC_params.VelThr) ' m/s.']));
    ncwriteatt(ncfile,'CSPD_QC','Conventions',char('SeaDataNet measurand qualifier flags.'));
    ncwriteatt(ncfile,'CSPD_QC','sdn_conventions_urn',char('SDN:L20::'));
    ncwriteatt(ncfile,'CSPD_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'CSPD_QC','add_offset',int8(0));
    ncwriteatt(ncfile,'CSPD_QC','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    
    ncwriteatt(ncfile,'AVRB_QC','long_name',char('Average Radial Bearing Quality Flag'));
    ncwriteatt(ncfile,'AVRB_QC','units',char('1'));
    ncwriteatt(ncfile,'AVRB_QC','valid_range',int8([48 65]));
    ncwriteatt(ncfile,'AVRB_QC','flag_values',int8([48 49 50 51 52 53 54 55 56 57 65]));
    ncwriteatt(ncfile,'AVRB_QC','flag_meanings',char('no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain'));
    ncwriteatt(ncfile,'AVRB_QC','comment',char(['Thresholds set to [' num2str(Radial_QC_params.AvgRadBear(1)) '-' num2str(Radial_QC_params.AvgRadBear(2)) '] deg.']));
    ncwriteatt(ncfile,'AVRB_QC','Conventions',char('SeaDataNet measurand qualifier flags.'));
    ncwriteatt(ncfile,'AVRB_QC','sdn_conventions_urn',char('SDN:L20::'));
    ncwriteatt(ncfile,'AVRB_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'AVRB_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'RDCT_QC','long_name',char('Radial Count Quality Flag'));
    ncwriteatt(ncfile,'RDCT_QC','units',char('1'));
    ncwriteatt(ncfile,'RDCT_QC','valid_range',int8([48 65]));
    ncwriteatt(ncfile,'RDCT_QC','flag_values',int8([48 49 50 51 52 53 54 55 56 57 65]));
    ncwriteatt(ncfile,'RDCT_QC','flag_meanings',char('no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain'));
    ncwriteatt(ncfile,'RDCT_QC','comment',char(['Thresholds set to ' num2str(Radial_QC_params.RadCnt) ' vectors.']));
    ncwriteatt(ncfile,'RDCT_QC','Conventions',char('SeaDataNet measurand qualifier flags.'));
    ncwriteatt(ncfile,'RDCT_QC','sdn_conventions_urn',char('SDN:L20::'));
    ncwriteatt(ncfile,'RDCT_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'RDCT_QC','add_offset',int8(0));
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

%%

%% Writes values in variables

try
    ncwrite(ncfile, 'TIME', nc.time - datenum(1950,1,1));
    ncwrite(ncfile, 'BEAR', nc.bear);
    ncwrite(ncfile, 'RNGE', nc.rnge);
    ncwrite(ncfile, 'DEPTH', nc.depth);
    ncwrite(ncfile, 'LATITUDE', nc.latitude);
    ncwrite(ncfile, 'LONGITUDE', nc.longitude);
    ncwrite(ncfile, 'crs', 0);
    ncwrite(ncfile, 'SDN_CRUISE', site_code');
    ncwrite(ncfile, 'SDN_STATION', platform_code');
    ncwrite(ncfile, 'SDN_EDMO_CODE', EDMO_code');
    ncwrite(ncfile, 'SDN_LOCAL_CDI_ID', dataID');
    ncwrite(ncfile, 'SDN_REFERENCES', TDS_catalog');
    ncwrite(ncfile, 'SDN_XLINK', xlink');
    ncwrite(ncfile, 'RDVA', nc.rdva);
    ncwrite(ncfile, 'DRVA', nc.drva);
    ncwrite(ncfile, 'EWCT', nc.ewct);
    ncwrite(ncfile, 'NSCT', nc.nsct);
    ncwrite(ncfile, 'ESPC', nc.espc);
    ncwrite(ncfile, 'ETMP', nc.etmp);
    ncwrite(ncfile, 'MAXV', nc.maxv);
    ncwrite(ncfile, 'MINV', nc.minv);
    ncwrite(ncfile, 'ERSC', nc.ersc);
    ncwrite(ncfile, 'ERTC', nc.ertc);
    ncwrite(ncfile, 'XDST', nc.xdst);
    ncwrite(ncfile, 'YDST', nc.ydst);
    ncwrite(ncfile, 'SPRC', nc.sprc);
    ncwrite(ncfile, 'NARX', nc.narx);
    ncwrite(ncfile, 'NATX', nc.natx);
    ncwrite(ncfile, 'SLTR', nc.sltr);
    ncwrite(ncfile, 'SLNR', nc.slnr);
    ncwrite(ncfile, 'SLTT', nc.sltt);
    ncwrite(ncfile, 'SLNT', nc.slnt);
    ncwrite(ncfile, 'SCDR', nc.scdr);
    ncwrite(ncfile, 'SCDT', nc.scdt);
    ncwrite(ncfile, 'TIME_SEADATANET_QC', nc.time_seadatanet_qc);
    ncwrite(ncfile, 'POSITION_SEADATANET_QC', nc.position_seadatanet_qc);
    ncwrite(ncfile, 'DEPTH_SEADATANET_QC', nc.depth_seadatanet_qc);
    ncwrite(ncfile, 'QCflag', nc.qcflag);
    ncwrite(ncfile, 'OWTR_QC', nc.owtr_qc);
    ncwrite(ncfile, 'MDFL_QC', nc.mdfl_qc);
    ncwrite(ncfile, 'VART_QC', nc.vart_qc);
    ncwrite(ncfile, 'CSPD_QC', nc.cspd_qc);
    ncwrite(ncfile, 'AVRB_QC', nc.avrb_qc);
    ncwrite(ncfile, 'RDCT_QC', nc.rdct_qc);
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

%%

%% Define global attributes

try
    
    % MANDATORY ATTRIBUTES
    % Discovery and Identification
    ncwriteatt(ncfile, '/', 'site_code', char(site_code));
    ncwriteatt(ncfile, '/', 'platform_code', char(platform_code));
    ncwriteatt(ncfile, '/', 'data_mode', char('R'));
    DoAIndex = find(not(cellfun('isempty', strfind(stationFields, 'DoA_estimation_method'))));
    ncwriteatt(ncfile, '/', 'DoA_estimation_method', char(stationData{DoAIndex}));
    calibration_typeIndex = find(not(cellfun('isempty', strfind(stationFields, 'calibration_type'))));
    ncwriteatt(ncfile, '/', 'calibration_type', char(stationData{calibration_typeIndex}));
    last_calibration_dateIndex = find(not(cellfun('isempty', strfind(stationFields, 'last_calibration_date'))));
    ncwriteatt(ncfile, '/', 'last_calibration_date', char(stationData{last_calibration_dateIndex}));
    calibration_linkIndex = find(not(cellfun('isempty', strfind(stationFields, 'calibration_link'))));
    ncwriteatt(ncfile, '/', 'calibration_link', char(stationData{calibration_linkIndex}));
    titleIndex = find(not(cellfun('isempty', strfind(networkFields, 'title'))));
    ncwriteatt(ncfile, '/', 'title', char(networkData{titleIndex}));
    summaryIndex = find(not(cellfun('isempty', strfind(stationFields, 'summary'))));
    ncwriteatt(ncfile, '/', 'summary', char(stationData{summaryIndex}));
    ncwriteatt(ncfile, '/', 'source', char('coastal structure'));
    ncwriteatt(ncfile, '/', 'source_platform_category_code', char('17'));
    institution_nameIndex = find(not(cellfun('isempty', strfind(stationFields, 'institution_name'))));
    ncwriteatt(ncfile, '/', 'institution', char(stationData{institution_nameIndex}));
    ncwriteatt(ncfile, '/', 'institution_edmo_code', char(num2str(EDMO_code)));
    ncwriteatt(ncfile, '/', 'data_assembly_center', char('European HFR Node'));
    ncwriteatt(ncfile, '/', 'id', char(dataID));
    % Geo-spatial-temporal
    ncwriteatt(ncfile, '/', 'data_type', char('HF radar radial data'));
%     ncwriteatt(ncfile, '/', 'feature_type', char('surface'));
    geospatial_lat_minIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lat_min'))));
    ncwriteatt(ncfile, '/', 'geospatial_lat_min', char(num2str(networkData{geospatial_lat_minIndex})));
    geospatial_lat_maxIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lat_max'))));
    ncwriteatt(ncfile, '/', 'geospatial_lat_max', char(num2str(networkData{geospatial_lat_maxIndex})));
    geospatial_lon_minIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lon_min'))));
    ncwriteatt(ncfile, '/', 'geospatial_lon_min', char(num2str(networkData{geospatial_lon_minIndex})));
    geospatial_lon_maxIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lon_max'))));
    ncwriteatt(ncfile, '/', 'geospatial_lon_max', char(num2str(networkData{geospatial_lon_maxIndex})));
    ncwriteatt(ncfile, '/', 'geospatial_vertical_max', char(num2str(vertMax)));
    ncwriteatt(ncfile, '/', 'geospatial_vertical_min', char('0'));
    ncwriteatt(ncfile, '/', 'time_coverage_start', char(timeCoverageStart));
    ncwriteatt(ncfile, '/', 'time_coverage_end', char(timeCoverageEnd));
    % Conventions used
    ncwriteatt(ncfile, '/', 'format_version', char('v2.1.2'));
    ncwriteatt(ncfile, '/', 'Conventions', char('CF-1.6, OceanSITES Manual 1.2, SeaDataNet_1.0, INSPIRE'));
    % Publication information
    ncwriteatt(ncfile, '/', 'update_interval', char('void'));
    ncwriteatt(ncfile, '/', 'citation', char(citation_str));
    ncwriteatt(ncfile, '/', 'distribution_statement', char(distribution_str));
    ncwriteatt(ncfile, '/', 'publisher_name', char('European HFR Node'));
    ncwriteatt(ncfile, '/', 'publisher_url', char('http://eurogoos.eu/'));
    ncwriteatt(ncfile, '/', 'publisher_email', char('euhfrnode@azti.es'));
    licenseIndex = find(not(cellfun('isempty', strfind(networkFields, 'license'))));
    ncwriteatt(ncfile, '/', 'license', char(networkData{licenseIndex}));
    acknowledgmentIndex = find(not(cellfun('isempty', strfind(networkFields, 'acknowledgment'))));
    ncwriteatt(ncfile, '/', 'acknowledgment', char(networkData{acknowledgmentIndex}));
    % Provenance
    ncwriteatt(ncfile, '/', 'date_created', char(dateCreated));
    ncwriteatt(ncfile, '/', 'history', char(['Data collected from ' time_coll{1} ' to ' time_coll{2} '. ' dateCreated ' netCDF file created by the European HFR Node']));
    ncwriteatt(ncfile, '/', 'date_modified', char(dateCreated));
    ncwriteatt(ncfile, '/', 'date_update', char(dateCreated));
    ncwriteatt(ncfile, '/', 'processing_level', char('2B'));
    contributor_nameIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_name'))));
    ncwriteatt(ncfile, '/', 'contributor_name', char(networkData{contributor_nameIndex}));
    contributor_roleIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_role'))));
    ncwriteatt(ncfile, '/', 'contributor_role', char(networkData{contributor_roleIndex}));
    contributor_emailIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_email'))));
    ncwriteatt(ncfile, '/', 'contributor_email', char(networkData{contributor_emailIndex}));
    
    % RECOMMENDED ATTRIBUTES
    % Discovery and Identification
    projectIndex = find(not(cellfun('isempty', strfind(networkFields, 'project'))));
    ncwriteatt(ncfile, '/', 'project', char(networkData{projectIndex}));
    ncwriteatt(ncfile, '/', 'naming_authority',char( naming_authorityStr));
    ncwriteatt(ncfile, '/', 'keywords', char('OCEAN CURRENTS, SURFACE WATER, RADAR, SCR-HF'));
    ncwriteatt(ncfile, '/', 'keywords_vocabulary', char('GCMD Science Keywords'));
    commentIndex = find(not(cellfun('isempty', strfind(networkFields, 'comment'))));
    ncwriteatt(ncfile, '/', 'comment', char(networkData{commentIndex}));
    ncwriteatt(ncfile, '/', 'data_language', char('eng'));
    ncwriteatt(ncfile, '/', 'data_character_set', char('utf8'));
    ncwriteatt(ncfile, '/', 'metadata_language', char('eng'));
    ncwriteatt(ncfile, '/', 'metadata_character_set', char('utf8'));
    ncwriteatt(ncfile, '/', 'topic_category', char('oceans'));
    network_nameIndex = find(not(cellfun('isempty', strfind(networkFields, 'network_name'))));
    ncwriteatt(ncfile, '/', 'network', char(networkData{network_nameIndex}));
    % Geo-spatial-temporal
    areaIndex = find(not(cellfun('isempty', strfind(networkFields, 'area'))));
    ncwriteatt(ncfile, '/', 'area', char(networkData{areaIndex}));
    ncwriteatt(ncfile, '/', 'geospatial_lat_units', char('degrees_north'));
    ncwriteatt(ncfile, '/', 'geospatial_lon_units', char('degrees_east'));
    ncwriteatt(ncfile, '/', 'geospatial_lat_resolution', char(num2str(latRes)));
    ncwriteatt(ncfile, '/', 'geospatial_lon_resolution', char(num2str(lonRes)));
    ncwriteatt(ncfile, '/', 'geospatial_vertical_resolution', char(num2str(vertMax)));
    ncwriteatt(ncfile, '/', 'geospatial_vertical_units', char('m'));
    ncwriteatt(ncfile, '/', 'geospatial_vertical_positive', char('down'));
    ncwriteatt(ncfile, '/', 'time_coverage_duration', char(timeCoverageDuration));
    ncwriteatt(ncfile, '/', 'time_coverage_resolution', char(timeCoverageResolution));
    ncwriteatt(ncfile, '/', 'reference_system', char('EPSG:4806'));
    ncwriteatt(ncfile, '/', 'cdm_data_type', char('Grid'));
    % Conventions used
    ncwriteatt(ncfile, '/', 'netcdf_version', char(netcdf.inqLibVers));
    ncwriteatt(ncfile, '/', 'netcdf_format', char('netcdf4_classic'));
    
    % OTHER ATTRIBUTES
    ncwriteatt(ncfile, '/', 'metadata_contact',char( 'lorenzo.corgnati@sp.ismar.cnr.it'));
    ncwriteatt(ncfile, '/', 'metadata_date_stamp', char(dateCreated));
    ncwriteatt(ncfile, '/', 'standard_name_vocabulary', char('NetCDF Climate and Forecast (CF) Metadata Convention Standard Name Table Version 1.6'));
    ncwriteatt(ncfile, '/', 'sensor', char(sensorATT));
    ncwriteatt(ncfile, '/', 'institution_reference', char(institution_websiteStr));
    ncwriteatt(ncfile, '/', 'references', char('High Frequency Radar European common data and metadata model Reference Card: all you need to know about High Frequency Radar (HFR) data harmonization at a glance. http://www.marineinsitu.eu/wp-content/uploads/2018/02/HFR_Data_Model_Reference_Card_v1.pdf'));
    ncwriteatt(ncfile, '/', 'software_name', char('HFR_Combiner'));
    ncwriteatt(ncfile, '/', 'software_version',char( 'v3.3'));
    ncwriteatt(ncfile, '/', 'date_issued', char(dateCreated));
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

%%

%% Retrieve information about the nc file
try
    ncfileInfo = dir(ncfile);
    ncFilesize = ncfileInfo.bytes/1024;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    rnA_err = 1;
end

%%

if(rnA_err==0)
    disp(['[' datestr(now) '] - - ' 'radial_netCDF_aggregation_v212.m successfully executed.']);
end

return

