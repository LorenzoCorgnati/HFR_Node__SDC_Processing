%% total_netCDF_aggregation_v212.m
% This function accesses the THREDDS catalog of the HFR networks via
% OpenDAP and creates HFR total aggregated netCDF datasets compliant to the
% European standard data model (that integrates CMEMS-INSTAC and SDC CF extension
% requirements) for distribution on the SeaDataNet infrastructure.
% The v211 version creates aggregated dataset according the the v2.1.2
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
%         tnA_err: error flag (0 = correct, 1 = error)
%         ncFileNoPath: filename of the generated nc file, without the full path
%         ncFilesize: size of the generated nc file.
%         tStart: starting time of the dataset (in datenum format).
%         tEnd: ending time of the dataset (in datenum format).
%         dataID: SDN local CDI id.

% Author: Lorenzo Corgnati
% Date: December 4, 2019

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [tnA_err, ncFileNoPath, ncFilesize, tStart, tEnd, dataID] = total_netCDF_aggregation_v212(networkData,networkFields,stationData,stationFields,timeSpan)

disp(['[' datestr(now) '] - - ' 'total_netCDF_aggregation_v212.m started.']);

tnA_err = 0;

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

%% Retrieve site codes and coordinates

try
    % Find the index of the site_lon field
    site_lonIndexC = strfind(stationFields, 'site_lon');
    site_lonIndex = find(not(cellfun('isempty', site_lonIndexC)));
    
    % Find the index of the site_lat field
    site_latIndexC = strfind(stationFields, 'site_lat');
    site_latIndex = find(not(cellfun('isempty', site_latIndexC)));
    
    % Find the index of the station_id field
    station_idIndexC = strfind(stationFields, 'station_id');
    station_idIndex = find(not(cellfun('isempty', station_idIndexC)));
    
    % Build the arrays for site coordinates and codes
    sitesLon = cell2mat(stationData(:,site_lonIndex));
    sitesLat = cell2mat(stationData(:,site_latIndex));
    sitesCodes = cell2mat(stationData(:,station_idIndex));
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

%%

%% Retrieve EDMO codes and institution names

try
    % Find the EDMO_code field from network data
    NT_EDMO_codeIndex = find(not(cellfun('isempty', strfind(networkFields, 'EDMO_code'))));
    NT_EDMO_code = networkData{NT_EDMO_codeIndex};
    NT_EDMO_code = NT_EDMO_code(NT_EDMO_code~=0);
    
    % Find the EDMO_code field from station data
    ST_EDMO_codeIndex = find(not(cellfun('isempty', strfind(stationFields, 'EDMO_code'))));
    ST_EDMO_code = cell2mat(stationData(:,ST_EDMO_codeIndex));
    ST_EDMO_code = ST_EDMO_code(ST_EDMO_code~=0);
    
    % Build the cumulative EDMO code list
    [EDMO_code,ia,ic] = unique([NT_EDMO_code; ST_EDMO_code]);
    EDMO_code = EDMO_code';
    EDMO_codeStr = sprintf('%.0d ' , EDMO_code);
    %     EDMO_codeStr = EDMO_codeStr(1:end-2);% strip final comma
    
    % Find the institution_name field from network data
    NT_institution_nameIndex = find(not(cellfun('isempty', strfind(networkFields, 'institution_name'))));
    NT_institution_name = networkData{NT_institution_nameIndex};
    
    % Find the institution_name field from station data
    ST_institution_nameIndex = find(not(cellfun('isempty', strfind(stationFields, 'institution_name'))));
    ST_institution_name = stationData(:,ST_institution_nameIndex);
    ST_institution_name(cellfun('isempty',ST_institution_name)) = [];
    
    % Build the cumulative institution name list
    institutionList = [NT_institution_name; ST_institution_name];
    institution_names = institutionList(ia);
    institution_nameStr = strjoin(institution_names,'; ');
    
    % Find the institution website field from network data
    NT_institution_websiteIndex = find(not(cellfun('isempty', strfind(networkFields, 'institution_website'))));
    NT_institution_website = networkData{NT_institution_websiteIndex};
    
    % Find the institution website field from station data
    ST_institution_websiteIndex = find(not(cellfun('isempty', strfind(stationFields, 'institution_website'))));
    ST_institution_website = stationData(:,ST_institution_websiteIndex);
    ST_institution_website(cellfun('isempty',ST_institution_website)) = [];
    
    % Build the cumulative institution website list
    websiteList = [NT_institution_website; ST_institution_website];
    institution_websites = websiteList(ia);
    institution_websiteStr = strjoin(institution_websites,'; ');
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

%%

%% Retrieve DoA, calibration types, calibration links and calibrations dates
try
    % Find the last_calibration_date field from station data
    ST_last_calibration_dateIndex = find(not(cellfun('isempty', strfind(stationFields, 'last_calibration_date'))));
    ST_last_calibration_date = datenum(stationData(:,ST_last_calibration_dateIndex));
    ST_last_calibration_date = ST_last_calibration_date(ST_last_calibration_date~=0);
    lastPatternStr = [sitesCodes(1,:) ': '];
    for lcd_idx=2:length(ST_last_calibration_date)
        lastPatternStr = [lastPatternStr datestr(ST_last_calibration_date(lcd_idx-1), 'yyyy-mm-dd') 'T' datestr(ST_last_calibration_date(lcd_idx-1), 'HH:MM:SS') 'Z; ' sitesCodes(lcd_idx,:) ': '];
    end
    lastPatternStr = [lastPatternStr datestr(ST_last_calibration_date(lcd_idx), 'yyyy-mm-dd') 'T' datestr(ST_last_calibration_date(lcd_idx), 'HH:MM:SS') 'Z'];
    
    % Find the DoA from station data
    ST_DoAIndex = find(not(cellfun('isempty', strfind(stationFields, 'DoA_estimation_method'))));
    ST_DoA = stationData(:,ST_DoAIndex);
    ST_DoA(cellfun('isempty',ST_DoA)) = [];
    %     ST_DoA = uniqueStrCell(ST_DoA);
    %     DoAStr = strjoin(ST_DoA,', ');
    DoAStr = [sitesCodes(1,:) ': '];
    for doa_idx=2:length(ST_DoA)
        DoAStr = [DoAStr ST_DoA{doa_idx} '; ' sitesCodes(doa_idx,:) ': '];
    end
    DoAStr = [DoAStr ST_DoA{doa_idx}];
    
    % Find the calibration_type from station data
    ST_calibration_typeIndex = find(not(cellfun('isempty', strfind(stationFields, 'calibration_type'))));
    ST_calibration_type = stationData(:,ST_calibration_typeIndex);
    ST_calibration_type(cellfun('isempty',ST_calibration_type)) = [];
    %     ST_calibration_type = uniqueStrCell(ST_calibration_type);
    %     calibration_typeStr = strjoin(ST_calibration_type,', ');
    calibration_typeStr = [sitesCodes(1,:) ': '];
    for ct_idx=2:length(ST_calibration_type)
        calibration_typeStr = [calibration_typeStr ST_calibration_type{ct_idx-1} '; ' sitesCodes(ct_idx,:) ': '];
    end
    calibration_typeStr = [calibration_typeStr ST_calibration_type{ct_idx}];
    
    % Find the calibration_link from station data
    ST_calibration_linkIndex = find(not(cellfun('isempty', strfind(stationFields, 'calibration_link'))));
    ST_calibration_link = stationData(:,ST_calibration_linkIndex);
    ST_calibration_link(cellfun('isempty',ST_calibration_link)) = [];
    %     ST_calibration_link = uniqueStrCell(ST_calibration_link);
    %     calibration_linkStr = strjoin(ST_calibration_link,', ');
    calibration_linkStr = [sitesCodes(1,:) ': '];
    for cl_idx=2:length(ST_calibration_link)
        calibration_linkStr = [calibration_linkStr ST_calibration_link{cl_idx-1} '; ' sitesCodes(cl_idx,:) ': '];
    end
    calibration_linkStr = [calibration_linkStr ST_calibration_link{ct_idx}];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

%%

%% Set non physical dimensions
try
    maxSite_dim = 50;
    maxInst_dim = length(EDMO_code);
    refMax_dim = 1;
    string15_dim = 15;
    string50_dim = 50;
    string80_dim = 80;
    string250_dim = 250;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

%%

%% Read the aggregated total file from THREDDS catalog via OpenDAP

try
    % Find the SDC_OpenDAP_data_url field from network data
    SDC_OpenDAP_data_urlIndex = find(not(cellfun('isempty', strfind(networkFields, 'SDC_OpenDAP_data_url'))));
    SDC_OpenDAP_data_url = networkData{SDC_OpenDAP_data_urlIndex};
    
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
        nc.latitude = ncread(SDC_OpenDAP_data_url,'LATITUDE');
        nc.longitude = ncread(SDC_OpenDAP_data_url,'LONGITUDE');
        nc.depth = ncread(SDC_OpenDAP_data_url,'DEPH');
        
        % Data variables
        nc.ewct = ncread(SDC_OpenDAP_data_url,'EWCT',[1,1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(nc.depth),length(iTime)]);
        nc.nsct = ncread(SDC_OpenDAP_data_url,'NSCT',[1,1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(nc.depth),length(iTime)]);
        nc.ewcs = ncread(SDC_OpenDAP_data_url,'EWCS',[1,1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(nc.depth),length(iTime)]);
        nc.nscs = ncread(SDC_OpenDAP_data_url,'NSCS',[1,1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(nc.depth),length(iTime)]);
        nc.ccov = ncread(SDC_OpenDAP_data_url,'CCOV',[1,1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(nc.depth),length(iTime)]);
        nc.gdop = ncread(SDC_OpenDAP_data_url,'GDOP',[1,1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(nc.depth),length(iTime)]);
        nc.narx = ncread(SDC_OpenDAP_data_url,'NARX',min(iTime),length(iTime));
        nc.natx = ncread(SDC_OpenDAP_data_url,'NATX',min(iTime),length(iTime));
        nc.sltr = ncread(SDC_OpenDAP_data_url,'SLTR',[1,min(iTime)],[maxSite_dim,length(iTime)]);
        nc.slnr = ncread(SDC_OpenDAP_data_url,'SLNR',[1,min(iTime)],[maxSite_dim,length(iTime)]);
        nc.sltt = ncread(SDC_OpenDAP_data_url,'SLTT',[1,min(iTime)],[maxSite_dim,length(iTime)]);
        nc.slnt = ncread(SDC_OpenDAP_data_url,'SLNT',[1,min(iTime)],[maxSite_dim,length(iTime)]);
        nc.scdr = ncread(SDC_OpenDAP_data_url,'SCDR',[1,1,min(iTime)],[string15_dim,maxSite_dim,length(iTime)]);
        nc.scdt = ncread(SDC_OpenDAP_data_url,'SCDT',[1,1,min(iTime)],[string15_dim,maxSite_dim,length(iTime)]);
        
        % QC variables
        nc.time_seadatanet_qc = ncread(SDC_OpenDAP_data_url,'TIME_QC',min(iTime),length(iTime));
        nc.position_seadatanet_qc = ncread(SDC_OpenDAP_data_url,'POSITION_QC',[1,1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(nc.depth),length(iTime)]);
        nc.depth_seadatanet_qc = ncread(SDC_OpenDAP_data_url,'DEPH_QC',min(iTime),length(iTime));
        nc.qcflag = ncread(SDC_OpenDAP_data_url,'QCflag',[1,1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(nc.depth),length(iTime)]);
        nc.vart_qc = ncread(SDC_OpenDAP_data_url,'VART_QC',[1,1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(nc.depth),length(iTime)]);
        nc.gdop_qc = ncread(SDC_OpenDAP_data_url,'GDOP_QC',[1,1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(nc.depth),length(iTime)]);
        nc.ddns_qc = ncread(SDC_OpenDAP_data_url,'DDNS_QC',[1,1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(nc.depth),length(iTime)]);
        nc.cspd_qc = ncread(SDC_OpenDAP_data_url,'CSPD_QC',[1,1,1,min(iTime)],[length(nc.longitude),length(nc.latitude),length(nc.depth),length(iTime)]);
    else
        disp(['[' datestr(now) '] - - No data available for the selected period.']);
        disp(['[' datestr(now) '] - - ' 'total_netCDF_aggregation_v212.m successfully executed.']);
        tnA_err = 1;
        return
    end
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
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
    
    nc.vart_qc(nc.vart_qc==0) = int8(48);
    nc.vart_qc(nc.vart_qc==1) = int8(49);
    nc.vart_qc(nc.vart_qc==2) = int8(50);
    nc.vart_qc(nc.vart_qc==3) = int8(51);
    nc.vart_qc(nc.vart_qc==4) = int8(52);
    nc.vart_qc(nc.vart_qc==8) = int8(56);
    nc.vart_qc(isnan(nc.vart_qc)) = int8(57);
    
    nc.gdop_qc(nc.gdop_qc==0) = int8(48);
    nc.gdop_qc(nc.gdop_qc==1) = int8(49);
    nc.gdop_qc(nc.gdop_qc==2) = int8(50);
    nc.gdop_qc(nc.gdop_qc==3) = int8(51);
    nc.gdop_qc(nc.gdop_qc==4) = int8(52);
    nc.gdop_qc(nc.gdop_qc==8) = int8(56);
    nc.gdop_qc(isnan(nc.gdop_qc)) = int8(57);
    
    nc.ddns_qc(nc.ddns_qc==0) = int8(48);
    nc.ddns_qc(nc.ddns_qc==1) = int8(49);
    nc.ddns_qc(nc.ddns_qc==2) = int8(50);
    nc.ddns_qc(nc.ddns_qc==3) = int8(51);
    nc.ddns_qc(nc.ddns_qc==4) = int8(52);
    nc.ddns_qc(nc.ddns_qc==8) = int8(56);
    nc.ddns_qc(isnan(nc.ddns_qc)) = int8(57);
    
    nc.cspd_qc(nc.cspd_qc==0) = int8(48);
    nc.cspd_qc(nc.cspd_qc==1) = int8(49);
    nc.cspd_qc(nc.cspd_qc==2) = int8(50);
    nc.cspd_qc(nc.cspd_qc==3) = int8(51);
    nc.cspd_qc(nc.cspd_qc==4) = int8(52);
    nc.cspd_qc(nc.cspd_qc==8) = int8(56);
    nc.cspd_qc(isnan(nc.cspd_qc)) = int8(57);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

%%

%% Build structures containing QC tests parameters

try
    GDOPthreshIndex = find(not(cellfun('isempty', strfind(networkFields, 'total_QC_GDOP_threshold'))));
    Total_QC_params.GDOPThr = networkData{GDOPthreshIndex};
    var_thr_Index = find(not(cellfun('isempty', strfind(networkFields, 'total_QC_variance_threshold'))));
    Total_QC_params.VarThr = networkData{var_thr_Index};
    temp_der_thr_Index = find(not(cellfun('isempty', strfind(networkFields, 'total_QC_temporal_derivative_threshold'))));
    Total_QC_params.TempDerThr.threshold = networkData{temp_der_thr_Index};
    maxspd_Index = find(not(cellfun('isempty', strfind(networkFields, 'total_QC_velocity_threshold'))));
    Total_QC_params.VelThr = networkData{maxspd_Index};
    dataDens_Index = find(not(cellfun('isempty', strfind(networkFields, 'total_QC_data_density_threshold'))));
    Total_QC_params.DataDensityThr = networkData{dataDens_Index};
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

%%

%% Set physical dimensions
try
    time_dim = length(iTime);
    lat_dim = length(nc.latitude);
    lon_dim = length(nc.longitude);
    depth_dim = 1;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

%%


%% Retrieve temporal and spatial metadata

try
    % Set netcdf format
    ncfmt = 'netcdf4_classic';
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

try
    timeref = datenum(1950,1,1);
    time_units = ['days since ' datestr(timeref, 'yyyy-mm-dd') 'T' datestr(timeref, 'HH:MM:SS') 'Z'];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
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
    tnA_err = 1;
end

% Set geographical resolutions
try
    % Latitude and longitude resolution (degrees)
    latDiff = diff(nc.latitude);
    lonDiff = diff(nc.longitude);
    latRes = abs(mean(latDiff));
    lonRes = abs(mean(lonDiff));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

try
    transmit_central_frequencyIndex = find(not(cellfun('isempty', strfind(stationFields, 'transmit_central_frequency'))));
    txFreq = cell2mat(stationData(:,transmit_central_frequencyIndex)).*1e6; % transmit frequency in Hertz
    vertMax = (3e8)/(8*pi*min(txFreq));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

% Set nc output file name
try
    startVec = datevec(nc.time(1));
    endVec = datevec(nc.time(length(nc.time)));
    % Check if the aggregation period is exactly 1 year
    if(timeSpan==12)
        time_str = sprintf('%.4d',startVec(1));
    elseif(timeSpan==1)
        time_str = sprintf('%.4d%.2d',startVec(1),startVec(2));
    else
        time_str = sprintf('%.4d%.2d%.2d_%.4d%.2d%.2d',startVec(1),startVec(2),startVec(3),endVec(1),endVec(2),endVec(3));
    end
    outputPathIndex = find(not(cellfun('isempty', strfind(networkFields, 'SDC_folder_path'))));
    network_idIndex = find(not(cellfun('isempty', strfind(networkFields, 'network_id'))));
    ncFilePath = strtrim(networkData{outputPathIndex});
    ncFileNoPath = ['TV_HF_' networkData{network_idIndex} '_' time_str '.nc'];
    ncfile = [ncFilePath filesep ncFileNoPath];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

% Set citation string and distribution string
try
    citation_statementIndex = find(not(cellfun('isempty', strfind(networkFields, 'citation_statement'))));
    citation_str = ['These data were collected and made freely available by the SeaDataNet project and the programs that contribute to it. ' networkData{citation_statementIndex}];
    distribution_str = 'These data follow SeaDataNet standards; they are public and free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.';
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

% Set naming authority
try
    %     institution_websiteIndex = find(not(cellfun('isempty', strfind(networkFields, 'institution_website'))));
    %     institution_websiteStr = networkData{institution_websiteIndex};
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
    tnA_err = 1;
end

% Set collection time
try
    time_coll{1} = [datestr(startVec, 'yyyy-mm-dd') 'T' datestr(startVec, 'HH:MM:SS') 'Z'];
    time_coll{2} = [datestr(endVec, 'yyyy-mm-dd') 'T' datestr(endVec, 'HH:MM:SS') 'Z'];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

try
    % Define EDIOS codes, site code, platform code, id and metadata resources
    EDIOS_Series_ID = networkData{network_idIndex};
    site_code = EDIOS_Series_ID;
    platform_code = [EDIOS_Series_ID '-Total'];
    dataID = [EDIOS_Series_ID '-Total_' time_coll{1} '_' time_coll{2}];
    metadata_pageIndex = find(not(cellfun('isempty', strfind(networkFields, 'metadata_page'))));
    TDS_catalog = networkData{metadata_pageIndex};
    xlink = ['<sdn_reference xlink:href="' TDS_catalog '" xlink:role="" xlink:type="URL"/>'];
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

% Deletes the eventually present netCDF file with the same name
try
    delete(ncfile);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

%%

%% Creates vars with their dimensions
try
    nccreate(ncfile,'TIME',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','double',...
        'Format',ncfmt);
    
    nccreate(ncfile,'LATITUDE',...
        'Dimensions',{'LATITUDE',lat_dim},...
        'Datatype','single',...
        'Format',ncfmt);
    
    nccreate(ncfile,'LONGITUDE',...
        'Dimensions',{'LONGITUDE',lon_dim},...
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
    
    nccreate(ncfile,'DEPTH',...
        'Dimensions',{'DEPTH',depth_dim},...
        'Datatype','single',...
        'Format',ncfmt);
    
    nccreate(ncfile,'EWCT',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue', netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'NSCT',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'EWCS',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'NSCS',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'CCOV',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int32',...
        'FillValue',netcdf.getConstant('NC_FILL_INT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'GDOP',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int16',...
        'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'TIME_SEADATANET_QC',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'POSITION_SEADATANET_QC',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'DEPTH_SEADATANET_QC',...
        'Dimensions',{'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'QCflag',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'VART_QC',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'GDOP_QC',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'DDNS_QC',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
        'Format',ncfmt);
    
    nccreate(ncfile,'CSPD_QC',...
        'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPTH', depth_dim, 'TIME',time_dim},...
        'Datatype','int8',...
        'FillValue',int8(57),...
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
        'Dimensions',{'STRING15',string15_dim,'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','char',...
        'FillValue',netcdf.getConstant('NC_FILL_CHAR'),...
        'Format',ncfmt);
    
    nccreate(ncfile,'SCDT',...
        'Dimensions',{'STRING15',string15_dim,'MAXSITE',maxSite_dim,'TIME',time_dim},...
        'Datatype','char',...
        'FillValue',netcdf.getConstant('NC_FILL_CHAR'),...
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
    
    ncwriteatt(ncfile,'LATITUDE','long_name',char('Latitude'));
    ncwriteatt(ncfile,'LATITUDE','standard_name',char('latitude'));
    ncwriteatt(ncfile,'LATITUDE','units',char('degrees_north'));
    ncwriteatt(ncfile,'LATITUDE','axis',char('Y'));
    ncwriteatt(ncfile,'LATITUDE','sdn_parameter_name',char('Latitude north'));
    ncwriteatt(ncfile,'LATITUDE','sdn_parameter_urn',char('SDN:P01::ALATZZ01'));
    ncwriteatt(ncfile,'LATITUDE','sdn_uom_name',char('Degrees north'));
    ncwriteatt(ncfile,'LATITUDE','sdn_uom_urn',char('SDN:P06::DEGN'));
    ncwriteatt(ncfile,'LATITUDE','grid_mapping',char('crs'));
    ncwriteatt(ncfile,'LATITUDE','ancillary_variables',char('POSITION_SEADATANET_QC'));
    
    ncwriteatt(ncfile,'LONGITUDE','long_name',char('Longitude'));
    ncwriteatt(ncfile,'LONGITUDE','standard_name',char('longitude'));
    ncwriteatt(ncfile,'LONGITUDE','units',char('degrees_east'));
    ncwriteatt(ncfile,'LONGITUDE','axis',char('X'));
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
    
    ncwriteatt(ncfile,'EWCT','long_name',char('West-east current component'));
    ncwriteatt(ncfile,'EWCT','standard_name',char('eastward_sea_water_velocity'));
    ncwriteatt(ncfile,'EWCT','units',char('m s-1'));
    ncwriteatt(ncfile,'EWCT','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'EWCT','add_offset',double(addOffset));
    ncwriteatt(ncfile,'EWCT','ioos_category',char('Currents'));
    ncwriteatt(ncfile,'EWCT','coordsys',char('geographic'));
    ncwriteatt(ncfile,'EWCT','sdn_parameter_name',char('Eastward current velocity in the water body'));
    ncwriteatt(ncfile,'EWCT','sdn_parameter_urn',char('SDN:P01::LCEWZZ01'));
    ncwriteatt(ncfile,'EWCT','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'EWCT','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'EWCT','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    %        ncwriteatt(ncfile,'EWCT','cell_methods',char('time: mean over hours time'));
    ncwriteatt(ncfile,'EWCT','valid_range',int16([(-10-addOffset)./scaleFactor, (10-addOffset)./scaleFactor]));
    %         ncwriteatt(ncfile,'EWCT','valid_min',double(-10.0));
    %         ncwriteatt(ncfile,'EWCT','valid_max',double(10.0));
    ncwriteatt(ncfile,'EWCT','ancillary_variables',char('QCflag, VART_QC, CSPD_QC, DDNS_QC, GDOP_QC'));
    
    ncwriteatt(ncfile,'NSCT','long_name',char('South-north current component'));
    ncwriteatt(ncfile,'NSCT','standard_name',char('northward_sea_water_velocity'));
    ncwriteatt(ncfile,'NSCT','units',char('m s-1'));
    ncwriteatt(ncfile,'NSCT','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'NSCT','add_offset',double(addOffset));
    ncwriteatt(ncfile,'NSCT','ioos_category',char('Currents'));
    ncwriteatt(ncfile,'NSCT','coordsys',char('geographic'));
    ncwriteatt(ncfile,'NSCT','sdn_parameter_name',char('Northward current velocity in the water body'));
    ncwriteatt(ncfile,'NSCT','sdn_parameter_urn',char('SDN:P01::LCNSZZ01'));
    ncwriteatt(ncfile,'NSCT','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'NSCT','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'NSCT','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    %        ncwriteatt(ncfile,'NSCT','cell_methods',char('time: mean over hours time'));
    ncwriteatt(ncfile,'NSCT','valid_range',int16([(-10-addOffset)./scaleFactor, (10-addOffset)./scaleFactor]));
    %         ncwriteatt(ncfile,'NSCT','valid_min',double(-10.0));
    %         ncwriteatt(ncfile,'NSCT','valid_max',double(10.0));
    ncwriteatt(ncfile,'NSCT','ancillary_variables',char('QCflag, VART_QC, CSPD_QC, DDNS_QC, GDOP_QC'));
    
    ncwriteatt(ncfile,'EWCS','long_name',char('Standard Deviation of Surface Eastward Sea Water Velocity'));
    %        ncwriteatt(ncfile,'EWCS','standard_name',char('surface_eastward_sea_water_velocity_standard_error'));
    ncwriteatt(ncfile,'EWCS','units',char('m s-1'));
    ncwriteatt(ncfile,'EWCS','valid_range',int16([(-10-addOffset)./scaleFactor, (10-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'EWCS','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    %         ncwriteatt(ncfile,'EWCS','valid_min',double(-10.0));
    %         ncwriteatt(ncfile,'EWCS','valid_max',double(10.0));
    ncwriteatt(ncfile,'EWCS','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'EWCS','add_offset',double(addOffset));
    ncwriteatt(ncfile,'EWCS','sdn_parameter_name',char('Eastward current velocity standard deviation in the water body'));
    ncwriteatt(ncfile,'EWCS','sdn_parameter_urn',char('SDN:P01::SDEWZZZZ'));
    ncwriteatt(ncfile,'EWCS','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'EWCS','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'EWCS','ancillary_variables',char('QCflag, VART_QC'));
    
    ncwriteatt(ncfile,'NSCS','long_name',char('Standard Deviation of Surface Northward Sea Water Velocity'));
    %        ncwriteatt(ncfile,'NSCS','standard_name',char('surface_northward_sea_water_velocity_standard_error'));
    ncwriteatt(ncfile,'NSCS','units',char('m s-1'));
    ncwriteatt(ncfile,'NSCS','valid_range',int16([(-10-addOffset)./scaleFactor, (10-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'NSCS','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    %         ncwriteatt(ncfile,'NSCS','valid_min',double(-10.0));
    %         ncwriteatt(ncfile,'NSCS','valid_max',double(10.0));
    ncwriteatt(ncfile,'NSCS','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'NSCS','add_offset',double(addOffset));
    ncwriteatt(ncfile,'NSCS','sdn_parameter_name',char('Northward current velocity standard deviation in the water body'));
    ncwriteatt(ncfile,'NSCS','sdn_parameter_urn',char('SDN:P01::SDNSZZZZ'));
    ncwriteatt(ncfile,'NSCS','sdn_uom_name',char('Metres per second'));
    ncwriteatt(ncfile,'NSCS','sdn_uom_urn',char('SDN:P06::UVAA'));
    ncwriteatt(ncfile,'NSCS','ancillary_variables',char('QCflag, VART_QC'));
    
    ncwriteatt(ncfile,'CCOV','long_name',char('Covariance of Surface Sea Water Velocity'));
    %         ncwriteatt(ncfile,'CCOV','standard_name',char('surface_sea_water_velocity_covariance'));
    ncwriteatt(ncfile,'CCOV','units',char('m2 s-2'));
    ncwriteatt(ncfile,'CCOV','valid_range',int32([(-10-addOffset)./(scaleFactor^2), (10-addOffset)./(scaleFactor^2)]));
    ncwriteatt(ncfile,'CCOV','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    %         ncwriteatt(ncfile,'CCOV','valid_min',double(-10.0));
    %         ncwriteatt(ncfile,'CCOV','valid_max',double(10.0));
    ncwriteatt(ncfile,'CCOV','scale_factor',double(scaleFactor^2));
    ncwriteatt(ncfile,'CCOV','add_offset',double(addOffset));
    ncwriteatt(ncfile,'CCOV','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'CCOV','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'CCOV','sdn_uom_name',char('Square metres per second squared'));
    ncwriteatt(ncfile,'CCOV','sdn_uom_urn',char('SDN:P06::SQM2'));
    ncwriteatt(ncfile,'CCOV','ancillary_variables',char('QCflag'));
    
    ncwriteatt(ncfile,'GDOP','long_name',char('Geometrical Dilution Of Precision'));
    %         ncwriteatt(ncfile,'GDOP','standard_name',char('gdop'));
    ncwriteatt(ncfile,'GDOP','units',char('1'));
    ncwriteatt(ncfile,'GDOP','valid_range',int16([(-20-addOffset)./scaleFactor, (20-addOffset)./scaleFactor]));
    ncwriteatt(ncfile,'GDOP','coordinates',char('TIME DEPTH LATITUDE LONGITUDE'));
    %         ncwriteatt(ncfile,'GDOP','valid_min',double(-20.0));
    %         ncwriteatt(ncfile,'GDOP','valid_max',double(20.0));
    ncwriteatt(ncfile,'GDOP','scale_factor',double(scaleFactor));
    ncwriteatt(ncfile,'GDOP','add_offset',double(addOffset));
    ncwriteatt(ncfile,'GDOP','comment',char(['The Geometric Dilution of Precision (GDOP) is the coefficient of the uncertainty, which relates the uncertainties in radial and velocity vectors.' ...
        ' The GDOP is a unit-less coefficient, which characterizes the effect that radar station geometry has on the measurement and position determination errors.' ...
        ' A low GDOP corresponds to an optimal geometric configuration of radar stations, and results in accurate surface current data. Essentially, GDOP is a quantitative way to relate the radial and velocity vector uncertainties.'...
        ' Setting a threshold on GDOP for total combination avoids the combination of radials with an intersection angle below a certain value.' ...
        ' GDOP is a useful metric for filtering errant velocities due to poor geometry.']));
    ncwriteatt(ncfile,'GDOP','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'GDOP','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'GDOP','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'GDOP','sdn_uom_urn',char('SDN:P06::UUUU'));
    ncwriteatt(ncfile,'GDOP','ancillary_variables',char('QCflag, GDOP_QC'));
    
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
    
    ncwriteatt(ncfile,'VART_QC','long_name',char('Variance Threshold Quality Flags'));
    ncwriteatt(ncfile,'VART_QC','units',char('1'));
    ncwriteatt(ncfile,'VART_QC','valid_range',int8([48 65]));
    ncwriteatt(ncfile,'VART_QC','flag_values',int8([48 49 50 51 52 53 54 55 56 57 65]));
    ncwriteatt(ncfile,'VART_QC','flag_meanings',char('no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain'));
    if(contains(sensorATT,'codar','IgnoreCase',true))
        ncwriteatt(ncfile,'VART_QC','comment',char(['Test not applicable to Direction Finding systems. The Temporal Derivative test is applied.' ...
            'Threshold set to ' num2str(Total_QC_params.TempDerThr.threshold) ' m/s. ']));
    elseif(contains(sensorATT,'wera','IgnoreCase',true))
        ncwriteatt(ncfile,'VART_QC','comment',char(['Threshold set to ' num2str(Total_QC_params.VarThr.threshold) ' m2/s2. ']));
    end
    ncwriteatt(ncfile,'VART_QC','Conventions',char('SeaDataNet measurand qualifier flags.'));
    ncwriteatt(ncfile,'VART_QC','sdn_conventions_urn',char('SDN:L20::'));
    ncwriteatt(ncfile,'VART_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'VART_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'GDOP_QC','long_name',char('GDOP Threshold Quality Flags'));
    ncwriteatt(ncfile,'GDOP_QC','units',char('1'));
    ncwriteatt(ncfile,'GDOP_QC','valid_range',int8([48 65]));
    ncwriteatt(ncfile,'GDOP_QC','flag_values',int8([48 49 50 51 52 53 54 55 56 57 65]));
    ncwriteatt(ncfile,'GDOP_QC','flag_meanings',char('no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain'));
    ncwriteatt(ncfile,'GDOP_QC','comment',char(['Threshold set to ' num2str(Total_QC_params.GDOPThr) '.']));
    ncwriteatt(ncfile,'GDOP_QC','Conventions',char('SeaDataNet measurand qualifier flags.'));
    ncwriteatt(ncfile,'GDOP_QC','sdn_conventions_urn',char('SDN:L20::'));
    ncwriteatt(ncfile,'GDOP_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'GDOP_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'DDNS_QC','long_name',char('Data Density Threshold Quality Flags'));
    ncwriteatt(ncfile,'DDNS_QC','units',char('1'));
    ncwriteatt(ncfile,'DDNS_QC','valid_range',int8([48 65]));
    ncwriteatt(ncfile,'DDNS_QC','flag_values',int8([48 49 50 51 52 53 54 55 56 57 65]));
    ncwriteatt(ncfile,'DDNS_QC','flag_meanings',char('no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain'));
    ncwriteatt(ncfile,'DDNS_QC','comment',char(['Threshold set to ' num2str(Total_QC_params.DataDensityThr) ' radials.']));
    ncwriteatt(ncfile,'DDNS_QC','Conventions',char('SeaDataNet measurand qualifier flags.'));
    ncwriteatt(ncfile,'DDNS_QC','sdn_conventions_urn',char('SDN:L20::'));
    ncwriteatt(ncfile,'DDNS_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'DDNS_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'CSPD_QC','long_name',char('Velocity Threshold Quality Flags'));
    ncwriteatt(ncfile,'CSPD_QC','units',char('1'));
    ncwriteatt(ncfile,'CSPD_QC','valid_range',int8([48 65]));
    ncwriteatt(ncfile,'CSPD_QC','flag_values',int8([48 49 50 51 52 53 54 55 56 57 65]));
    ncwriteatt(ncfile,'CSPD_QC','flag_meanings',char('no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain'));
    ncwriteatt(ncfile,'CSPD_QC','comment',char(['Threshold set to ' num2str(Total_QC_params.VelThr) ' m/s.']));
    ncwriteatt(ncfile,'CSPD_QC','Conventions',char('SeaDataNet measurand qualifier flags.'));
    ncwriteatt(ncfile,'CSPD_QC','sdn_conventions_urn',char('SDN:L20::'));
    ncwriteatt(ncfile,'CSPD_QC','scale_factor',int8(1));
    ncwriteatt(ncfile,'CSPD_QC','add_offset',int8(0));
    
    ncwriteatt(ncfile,'NARX','long_name',char('Number of Receive Antennas'));
    ncwriteatt(ncfile,'NARX','units',char('1'));
    ncwriteatt(ncfile,'NARX','valid_range',int8([0 maxSite_dim]));
    %         ncwriteatt(ncfile,'NARX','coordinates',char('TIME'));
    ncwriteatt(ncfile,'NARX','scale_factor',int8(1));
    ncwriteatt(ncfile,'NARX','add_offset',int8(0));
    ncwriteatt(ncfile,'NARX','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'NARX','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'NARX','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'NARX','sdn_uom_urn',char('SDN:P06::UUUU'));
    
    ncwriteatt(ncfile,'NATX','long_name',char('Number of Transmit Antennas'));
    ncwriteatt(ncfile,'NATX','units',char('1'));
    ncwriteatt(ncfile,'NATX','valid_range',int8([0 maxSite_dim]));
    %         ncwriteatt(ncfile,'NATX','coordinates',char('TIME'));
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
    %     ncwriteatt(ncfile,'SCDR','valid_range',char(''));
    ncwriteatt(ncfile,'SCDR','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'SCDR','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'SCDR','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'SCDR','sdn_uom_urn',char('SDN:P06::UUUU'));
    
    ncwriteatt(ncfile,'SCDT','long_name',char('Transmit Antenna Codes'));
    ncwriteatt(ncfile,'SCDT','units',char('1'));
    %     ncwriteatt(ncfile,'SCDT','valid_range',char(''));
    ncwriteatt(ncfile,'SCDT','sdn_parameter_name',char(''));
    ncwriteatt(ncfile,'SCDT','sdn_parameter_urn',char(''));
    ncwriteatt(ncfile,'SCDT','sdn_uom_name',char('Dimensionless'));
    ncwriteatt(ncfile,'SCDT','sdn_uom_urn',char('SDN:P06::UUUU'));
    
    %% Writes values in variables
    %         ncwrite(ncfile,'TIME',int32((TUVgrid.TimeStamp-timeref)*86400));
    ncwrite(ncfile,'TIME',nc.time-timeref);
    ncwrite(ncfile,'LATITUDE',nc.latitude);
    ncwrite(ncfile,'LONGITUDE',nc.longitude);
    ncwrite(ncfile,'crs',0);
    ncwrite(ncfile,'SDN_CRUISE',site_code');
    ncwrite(ncfile,'SDN_STATION',platform_code');
    ncwrite(ncfile,'SDN_LOCAL_CDI_ID',dataID');
    ncwrite(ncfile,'SDN_EDMO_CODE',EDMO_code');
    ncwrite(ncfile,'SDN_REFERENCES',TDS_catalog');
    ncwrite(ncfile,'SDN_XLINK',xlink');
    ncwrite(ncfile,'DEPTH',nc.depth);
    ncwrite(ncfile,'EWCT',nc.ewct);
    ncwrite(ncfile,'NSCT',nc.nsct);
    ncwrite(ncfile,'EWCS',nc.ewcs);
    ncwrite(ncfile,'NSCS',nc.nscs);
    ncwrite(ncfile,'CCOV',nc.ccov);
    ncwrite(ncfile,'GDOP',nc.gdop);
    ncwrite(ncfile,'NARX',nc.narx);
    ncwrite(ncfile,'NATX',nc.natx);
    ncwrite(ncfile,'SLTR',nc.sltr);
    ncwrite(ncfile,'SLNR',nc.slnr);
    ncwrite(ncfile,'SLTT',nc.sltt);
    ncwrite(ncfile,'SLNT',nc.slnt);
    ncwrite(ncfile,'SCDR',nc.scdr);
    ncwrite(ncfile,'SCDT',nc.scdt);
    ncwrite(ncfile,'TIME_SEADATANET_QC',nc.time_seadatanet_qc);
    ncwrite(ncfile,'POSITION_SEADATANET_QC',nc.position_seadatanet_qc);
    ncwrite(ncfile,'DEPTH_SEADATANET_QC',nc.depth_seadatanet_qc);
    ncwrite(ncfile,'QCflag',nc.qcflag);
    ncwrite(ncfile,'VART_QC',nc.vart_qc);
    ncwrite(ncfile,'GDOP_QC',nc.gdop_qc);
    ncwrite(ncfile,'DDNS_QC',nc.ddns_qc);
    ncwrite(ncfile,'CSPD_QC',nc.cspd_qc);
    
    %% Define global attributes
    
    % MANDATORY ATTRIBUTES
    % Discovery and Identification
    ncwriteatt(ncfile,'/','site_code',char(site_code));
    ncwriteatt(ncfile,'/','platform_code',char(platform_code));
    ncwriteatt(ncfile,'/','data_mode',char('R'));
    ncwriteatt(ncfile,'/','DoA_estimation_method',char(DoAStr));
    ncwriteatt(ncfile,'/','calibration_type',char(calibration_typeStr));
    ncwriteatt(ncfile,'/','last_calibration_date',char(lastPatternStr));
    ncwriteatt(ncfile,'/','calibration_link',char(calibration_linkStr));
    titleIndex = find(not(cellfun('isempty', strfind(networkFields, 'title'))));
    ncwriteatt(ncfile,'/','title',char(networkData{titleIndex}));
    summaryIndex = find(not(cellfun('isempty', strfind(networkFields, 'summary'))));
    ncwriteatt(ncfile,'/','summary',char(networkData{summaryIndex}));
    ncwriteatt(ncfile,'/','source',char('coastal structure'));
    ncwriteatt(ncfile,'/','source_platform_category_code',char('17'));
    ncwriteatt(ncfile,'/','institution',char(institution_nameStr));
    ncwriteatt(ncfile,'/','institution_edmo_code',char(num2str(NT_EDMO_code)));
    ncwriteatt(ncfile,'/','data_assembly_center',char('European HFR Node'));
    ncwriteatt(ncfile,'/','id',char(dataID));
    
    % Geo-spatial-temporal
    ncwriteatt(ncfile,'/','data_type', char('HF radar total data'));
    ncwriteatt(ncfile,'/','featureType',char('surface'));
    geospatial_lat_minIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lat_min'))));
    ncwriteatt(ncfile,'/','geospatial_lat_min',char(num2str(networkData{geospatial_lat_minIndex})));
    geospatial_lat_maxIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lat_max'))));
    ncwriteatt(ncfile,'/','geospatial_lat_max',char(num2str(networkData{geospatial_lat_maxIndex})));
    geospatial_lon_minIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lon_min'))));
    ncwriteatt(ncfile,'/','geospatial_lon_min',char(num2str(networkData{geospatial_lon_minIndex})));
    geospatial_lon_maxIndex = find(not(cellfun('isempty', strfind(networkFields, 'geospatial_lon_max'))));
    ncwriteatt(ncfile,'/','geospatial_lon_max',char(num2str(networkData{geospatial_lon_maxIndex})));
    ncwriteatt(ncfile,'/','geospatial_vertical_min', char('0'));
    ncwriteatt(ncfile,'/','geospatial_vertical_max', char(num2str(vertMax)));
    ncwriteatt(ncfile, '/','time_coverage_start',char(timeCoverageStart));
    ncwriteatt(ncfile, '/','time_coverage_end',char(timeCoverageEnd));
    % Conventions used
    ncwriteatt(ncfile,'/','format_version',char('v2.1.2'));
    ncwriteatt(ncfile,'/','Conventions',char('CF-1.6, OceanSITES Manual 1.2, SeaDataNet_1.0, INSPIRE'));
    % Publication information
    ncwriteatt(ncfile,'/','update_interval',char('void'));
    ncwriteatt(ncfile,'/','citation',char(citation_str));
    ncwriteatt(ncfile,'/','distribution_statement',char(distribution_str));
    ncwriteatt(ncfile,'/','publisher_name',char('European HFR Node'));
    ncwriteatt(ncfile,'/','publisher_url',char('http://eurogoos.eu/'));
    ncwriteatt(ncfile,'/','publisher_email',char('euhfrnode@azti.es'));
    licenseIndex = find(not(cellfun('isempty', strfind(networkFields, 'license'))));
    ncwriteatt(ncfile,'/','license',char(networkData{licenseIndex}));
    acknowledgmentIndex = find(not(cellfun('isempty', strfind(networkFields, 'acknowledgment'))));
    ncwriteatt(ncfile,'/','acknowledgment',char(networkData{acknowledgmentIndex}));
    % Provenance
    ncwriteatt(ncfile,'/','date_created',char(dateCreated));
    ncwriteatt(ncfile,'/','history',char(['Data collected from ' time_coll{1} ' to ' time_coll{2} '. ' dateCreated ' netCDF file created by the European HFR Node']));
    ncwriteatt(ncfile,'/','date_modified',char(dateCreated));
    ncwriteatt(ncfile,'/','date_update',char(dateCreated));
    ncwriteatt(ncfile,'/','processing_level',char('3B'));
    
    contributor_nameIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_name'))));
    ncwriteatt(ncfile,'/','contributor_name',char(networkData{contributor_nameIndex}));
    contributor_roleIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_role'))));
    ncwriteatt(ncfile,'/','contributor_role',char(networkData{contributor_roleIndex}));
    contributor_emailIndex = find(not(cellfun('isempty', strfind(networkFields, 'contributor_email'))));
    ncwriteatt(ncfile,'/','contributor_email',char(networkData{contributor_emailIndex}));
    
    % RECOMMENDED ATTRIBUTES
    % Discovery and Identification
    projectIndex = find(not(cellfun('isempty', strfind(networkFields, 'project'))));
    ncwriteatt(ncfile,'/','project',char(networkData{projectIndex}));
    ncwriteatt(ncfile,'/','naming_authority',char(naming_authorityStr));
    ncwriteatt(ncfile,'/','keywords',char('OCEAN CURRENTS, SURFACE WATER, RADAR, SCR-HF'));
    ncwriteatt(ncfile,'/','keywords_vocabulary',char('GCMD Science Keywords'));
    commentIndex = find(not(cellfun('isempty', strfind(networkFields, 'comment'))));
    ncwriteatt(ncfile,'/','comment',char(networkData{commentIndex}));
    ncwriteatt(ncfile,'/','data_language',char('eng'));
    ncwriteatt(ncfile,'/','data_character_set',char('utf8'));
    ncwriteatt(ncfile,'/','metadata_language',char('eng'));
    ncwriteatt(ncfile,'/','metadata_character_set',char('utf8'));
    ncwriteatt(ncfile,'/','topic_category',char('oceans'));
    network_nameIndex = find(not(cellfun('isempty', strfind(networkFields, 'network_name'))));
    ncwriteatt(ncfile,'/','network',char(networkData{network_nameIndex}));
    % Geo-spatial-temporal
    areaIndex = find(not(cellfun('isempty', strfind(networkFields, 'area'))));
    ncwriteatt(ncfile,'/','area',char(networkData{areaIndex}));
    ncwriteatt(ncfile,'/','geospatial_lat_units',char('degrees_north'));
    ncwriteatt(ncfile,'/','geospatial_lon_units',char('degrees_east'));
    ncwriteatt(ncfile,'/','geospatial_lat_resolution',char(num2str(latRes)));
    ncwriteatt(ncfile,'/','geospatial_lon_resolution',char(num2str(lonRes)));
    ncwriteatt(ncfile,'/','geospatial_vertical_resolution', char(num2str(vertMax)));
    ncwriteatt(ncfile,'/','geospatial_vertical_units', char('m'));
    ncwriteatt(ncfile,'/','geospatial_vertical_positive', char('down'));
    ncwriteatt(ncfile, '/','time_coverage_duration',char(timeCoverageDuration));
    ncwriteatt(ncfile, '/','time_coverage_resolution',char(timeCoverageResolution));
    ncwriteatt(ncfile,'/','reference_system',char('EPSG:4806'));
    grid_resolutionIndex = find(not(cellfun('isempty', strfind(networkFields, 'grid_resolution'))));
    ncwriteatt(ncfile,'/','grid_resolution',char(num2str(networkData{grid_resolutionIndex})));
    ncwriteatt(ncfile,'/','cdm_data_type',char('Grid'));
    % Conventions used
    ncwriteatt(ncfile,'/','netcdf_version',char(netcdf.inqLibVers));
    ncwriteatt(ncfile,'/','netcdf_format',char(ncfmt));
    
    % OTHER ATTRIBUTES
    ncwriteatt(ncfile,'/','metadata_contact',char('lorenzo.corgnati@sp.ismar.cnr.it'));
    ncwriteatt(ncfile,'/','metadata_date_stamp',char(dateCreated));
    ncwriteatt(ncfile,'/','standard_name_vocabulary',char('NetCDF Climate and Forecast (CF) Metadata Convention Standard Name Table Version 1.6'));
    ncwriteatt(ncfile,'/','sensor',char(sensorATT));
    ncwriteatt(ncfile,'/','institution_reference',char(institution_websiteStr));
    ncwriteatt(ncfile,'/','date_issued',char(dateCreated));
    ncwriteatt(ncfile,'/','software_name',char('HFR_Combiner'));
    ncwriteatt(ncfile,'/','software_version',char('v3.3'));
    ncwriteatt(ncfile,'/','references',char('High Frequency Radar European common data and metadata model Reference Card: all you need to know about High Frequency Radar (HFR) data harmonization at a glance. http://www.marineinsitu.eu/wp-content/uploads/2018/02/HFR_Data_Model_Reference_Card_v1.pdf'));
    
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

%% Retrieve information about the nc file
try
    ncfileInfo = dir(ncfile);
    ncFilesize = ncfileInfo.bytes/1024;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    tnA_err = 1;
end

%%

if(tnA_err==0)
    disp(['[' datestr(now) '] - - ' 'total_netCDF_aggregation_v212.m successfully executed.']);
end

return

