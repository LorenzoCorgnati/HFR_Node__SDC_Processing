%% SDC_radials.m
% This application builds the historical radial datasets to be distributed via the
% SeaDataNet infrastructure by reading hourly data from the EU HFR NODE
% THREDDS catalog via OpenDAP and aggregating them according to the
% European standard data model.
% The application also builds the CDI entry for each historical dataset.
% The information for assembling metadata are read from the HFR database.

% Author: Lorenzo Corgnati
% Date: July 18, 2019

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

warning('off', 'all');

SDCR_err = 0;

disp(['[' datestr(now) '] - - ' 'SDC_radials started.']);

%%

%% Connect to database

try
    conn = database(sqlConfig.database,sqlConfig.user,sqlConfig.password,'Vendor','MySQL','Server',sqlConfig.host);
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully established.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    SDCR_err = 1;
end

%%

%% Query the database for retrieving network data

% Set and exectute the query
try
    network_selectquery = 'SELECT * FROM network_tb WHERE SDC_distribution_flag=1';
    network_curs = exec(conn,network_selectquery);
    disp(['[' datestr(now) '] - - ' 'Query to network_tb table for retrieving network data successfully executed.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    SDCR_err = 1;
end

% Fetch data
try
    network_curs = fetch(network_curs);
    network_data = network_curs.Data;
    disp(['[' datestr(now) '] - - ' 'Network data successfully fetched from network_tb table.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    SDCR_err = 1;
end

% Retrieve column names
try
    network_columnNames = columnnames(network_curs,true);
    disp(['[' datestr(now) '] - - ' 'Column names from network_tb table successfully retrieved.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    SDCR_err = 1;
end

% Retrieve the number of networks
try
    numNetworks = rows(network_curs);
    disp(['[' datestr(now) '] - - ' 'Number of networks successfully retrieved from network_tb table.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    SDCR_err = 1;
end

% Close cursor
try
    close(network_curs);
    disp(['[' datestr(now) '] - - ' 'Cursor to network_tb table successfully closed.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    SDCR_err = 1;
end

%%

%% Scan the networks for generating aggregated datasets and related CDI entries

try
    % Find the index of the network_id field
    network_idIndexC = strfind(network_columnNames, 'network_id');
    network_idIndex = find(not(cellfun('isempty', network_idIndexC)));
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    SDCR_err = 1;
end

% Scan the networks
for network_idx=1:numNetworks
    SDCR_err = 0;
    % Retrieve information on the stations belonging to the current network
    try
        station_selectquery = ['SELECT * FROM station_tb WHERE network_id = ' '''' network_data{network_idx,network_idIndex} ''' AND SDC_distribution_flag=1'];
        station_curs = exec(conn,station_selectquery);
        disp(['[' datestr(now) '] - - ' 'Query to station_tb table for retrieving the stations of the ' network_data{network_idx,network_idIndex} ' network successfully executed.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        SDCR_err = 1;
    end
    
    % Fetch data
    try
        station_curs = fetch(station_curs);
        station_data = station_curs.Data;
        disp(['[' datestr(now) '] - - ' 'Data of the stations of the ' network_data{network_idx,network_idIndex} ' network successfully fetched from station_tb table.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        SDCR_err = 1;
    end
    
    % Retrieve column names
    try
        station_columnNames = columnnames(station_curs,true);
        disp(['[' datestr(now) '] - - ' 'Column names from station_tb table successfully retrieved.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        SDCR_err = 1;
    end
    
    % Retrieve the number of stations belonging to the current network
    try
        numStations = rows(station_curs);
        disp(['[' datestr(now) '] - - ' 'Number of stations belonging to the ' network_data{network_idx,network_idIndex} ' network successfully retrieved from station_tb table.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        SDCR_err = 1;
    end
    
    % Close cursor to station_tb table
    try
        close(station_curs);
        disp(['[' datestr(now) '] - - ' 'Cursor to station_tb table successfully closed.']);
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        SDCR_err = 1;
    end
    
    try
        % Find the index of the station_id field
        station_idIndexC = strfind(station_columnNames, 'station_id');
        station_idIndex = find(not(cellfun('isempty', station_idIndexC)));
    catch err
        disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        SDCR_err = 1;
    end
    
    for station_idx=1:numStations
        SDCR = 0;
        % Create the aggregated netCDF dataset
        % v2.1.2
        [rnA_err, datasetName, datasetSize, startDate, endDate, sdnLocalCDIid] = radial_netCDF_aggregation_v212(network_data(network_idx,:),network_columnNames,station_data(station_idx,:),station_columnNames,timeSpan);
        if(rnA_err==0)
            disp(['[' datestr(now) '] - - ' datasetName ' radial v2.1.2 dataset successfully created and stored.']);
        end
        
        % Insert generated radial dataset info in radial_SDCnetCDF_tb table
        try
            if((SDCR_err==0) && (rnA_err==0) && (exist('datasetName','var') ~= 0))
                % Delete the eventually present entry with the same name from the database
                radial_deletequery = ['DELETE FROM radial_SDCnetCDF_tb WHERE filename = ' '''' datasetName ''''];
                radial_deletecurs = exec(conn,radial_deletequery);
                
                % Define a cell array containing the column names to be added
                addColnames = {'filename' 'network_id' 'station_id' 'start_date' 'end_date' 'creation_date' 'filesize' 'sent_flag'};
                
                % Define a cell array that contains the data for insertion
                addData = {datasetName,network_data{network_idx,network_idIndex},station_data{station_idx,station_idIndex},(datestr(startDate,'yyyy-mm-dd HH:MM:SS')),(datestr(endDate-1/24,'yyyy-mm-dd HH:MM:SS')),(datestr(now,'yyyy-mm-dd HH:MM:SS')),datasetSize,0};
                
                % Append the product data into the total_HFRnetCDF_tb table on the database.
                tablename = 'radial_SDCnetCDF_tb';
                datainsert(conn,tablename,addColnames,addData);
                disp(['[' datestr(now) '] - - ' datasetName ' radial dataset information successfully inserted into total_SDCnetCDF_tb table.']);
            end
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            SDCR_err = 1;
        end
        
        % Insert metadata into DB for the generation of CDI entry
        try
            if((SDCR_err==0) &&(rnA_err==0))
                SDC_CDI_radials_metadata_to_DB;
            end
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            SDCR_err = 1;
        end
        
        % Create the CDI entry
        try
            if((SDCR_err==0)&& (rnA_err==0))
                % Set Mikado xml conf filename and output dir
                confFile = ['Mikado_conf_files' filesep ncreadatt(ncfile,'/','platform_code') '.xml'];
                outputDir = [station_data{station_idx,find(not(cellfun('isempty', strfind(station_columnNames, 'SDC_folder_path'))))} filesep station_data{station_idx,station_idIndex} filesep 'CDI'];
                if (exist(outputDir, 'dir') ~= 7)
                    mkdir(outputDir);
                end
                % Launch the Mikado batch command
                mikadoCommand = ['java -cp ' mikadoHome '/dist/*:dist/lib/* mikado.Mikado ' ...
                    'mikado-home=' mikadoHome '/ batch-type=XmlFiles  batch-mode=CDI19139  ' ...
                    'conf-file=' confFile ' output-dir=' outputDir ' continue-when-error=false ' ...
                    'log-file=/var/log/EU_HFR_Node_SDC_dataset_builder/CDI.log trace=false max-files-in-zip=1000 UpdateCenter=on'];
                [status,cmdout] = system(mikadoCommand);
                if((status==0) && (exist([outputDir filesep ncreadatt(ncfile,'/','id') '.xml'],'file')==2))
                    disp(['[' datestr(now) '] - - ' ncreadatt(ncfile,'/','id') '.xml CDI file successfully created and stored.']);
                else
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> Something went wrong with the CDI generation.']);
                    SDCR_err = 1;
                end
            end
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            SDCR_err = 1;
        end
        
        % Update radial_SDCnetCDF_tb table by adding SDN_LOCAL_CDI_ID
        try
            if((SDCR_err==0) && (rnA_err==0))
                if(strcmp(sdnLocalCDIid,ncreadatt(ncfile,'/','id')))
                    % Define a cell array containing the column names to be updated
                    updateColnames = {'SDN_LOCAL_CDI_ID'};
                    
                    % Define a cell array that contains the data for insertion
                    updateData = {sdnLocalCDIid};
                    
                    % Update the total_SDCnetCDF_tb table on the database
                    tablename = 'radial_SDCnetCDF_tb';
                    whereclause = ['WHERE filename = ' '''' datasetName ''''];
                    update(conn,tablename,updateColnames,updateData,whereclause);
                    disp(['[' datestr(now) '] - - ' 'radial_SDCnetCDF_tb table successfully updated with SDN_LOCAL_CDI_ID for dataset ' datasetName '.nc.']);
                else
                    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> Something went wrong with the SDN_LOCAL_CDI_ID.']);
                end
            end
        catch err
            disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            SDCR_err = 1;
        end
        
        
    end
    
end

%%

%% Close connection

try
    close(conn);
    disp(['[' datestr(now) '] - - ' 'Connection to database successfully closed.']);
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    SDCR_err = 1;
end

%%

if(SDCR_err==0)
    disp(['[' datestr(now) '] - - ' 'SDC_radials successfully executed.']);
end