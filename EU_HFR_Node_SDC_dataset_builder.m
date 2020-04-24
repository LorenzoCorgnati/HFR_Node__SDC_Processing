%% EU_HFR_Node_SDC_dataset_builder.m
% This wrapper launches the scripts for building the historical total and radial 
% datasets and the CDI entries for distribution on SeaDataNet infrastructure.

% Author: Lorenzo Corgnati
% Date: June 19, 2019

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

warning('off', 'all');

clear all
close all
clc

% Setup netCDF toolbox
setup_nctoolbox;

% Setup JBDC driver for MySQL
javaaddpath('/home/lorenz/Toolboxes/Matlab_HFR_AddOn/mysql-connector-java-5.1.17.jar');

% Setup map colormap
set(0,'DefaultFigureColormap',feval('jet'));

EHNSDB_err = 0;

disp(['[' datestr(now) '] - - ' 'EU_HFR_Node_SDC_dataset_builder started.']);

%%

%% Set Mikado home folder

% mikadoHome = '/opt/mikado_V3.5.2';
mikadoHome = '/opt/mikado_V3.5.3';
% mikadoHome = '/opt/mikado_V3.6';

%%

%% Set database parameters

sqlConfig.user = 'HFR_lorenzo';
sqlConfig.password = 'xWeLXHFQfvpBmDYO';
sqlConfig.host = '150.145.136.8';
sqlConfig.database = 'HFR_node_db';

%%

%% Processing

% Set the temporal interval of the datasets to be built
try
    % Time span expressed in number of months
    timeSpan = 1;
catch err
    disp(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    EHNSDB_err = 1;
end

%%

% TOTAL FILE PROCESSING
SDC_totals;

% RADIAL FILE PROCESSING
SDC_radials;

disp(['[' datestr(now) '] - - ' 'EU_HFR_Node_SDC_dataset_builder successfully executed.']);

%%