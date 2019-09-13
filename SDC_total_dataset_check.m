%% SDC_total_dataset_check.m
% This script compares the SDC total aggregated dataset with some hourly files

% Author: Lorenzo Corgnati
% Date: July 18, 2019

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

warning('off', 'all');

clear all
close all
clc

% Set non physical dimensions
maxSite_dim = 50;
string15_dim = 15;

% Setup netCDF toolbox
setup_nctoolbox;


%% CMEMS-INSTAC data

INSTACfile = '/home/lorenz/Downloads/HFR-TirLig-Total_2019_06_21_1700.nc';

% Read time and convert it to Matlab time
instac.time = ncread_cf_time(INSTACfile,'TIME');

% Read variables
% Coordinate variables
instac.latitude = ncread(INSTACfile,'LATITUDE');
instac.longitude = ncread(INSTACfile,'LONGITUDE');
instac.depth = ncread(INSTACfile,'DEPH');

% Data variables
instac.ewct = ncread(INSTACfile,'EWCT');
instac.nsct = ncread(INSTACfile,'NSCT');
instac.ewcs = ncread(INSTACfile,'EWCS');
instac.nscs = ncread(INSTACfile,'NSCS');
instac.ccov = ncread(INSTACfile,'CCOV');
instac.gdop = ncread(INSTACfile,'GDOP');
instac.narx = ncread(INSTACfile,'NARX');
instac.natx = ncread(INSTACfile,'NATX');
instac.sltr = ncread(INSTACfile,'SLTR');
instac.slnr = ncread(INSTACfile,'SLNR');
instac.sltt = ncread(INSTACfile,'SLTT');
instac.slnt = ncread(INSTACfile,'SLNT');
instac.scdr = ncread(INSTACfile,'SCDR');
instac.scdt = ncread(INSTACfile,'SCDT');

% QC variables
instac.time_seadatanet_qc = ncread(INSTACfile,'TIME_SEADATANET_QC');
instac.position_seadatanet_qc = ncread(INSTACfile,'POSITION_SEADATANET_QC');
instac.depth_seadatanet_qc = ncread(INSTACfile,'DEPTH_SEADATANET_QC');
instac.qcflag = ncread(INSTACfile,'QCflag');
instac.vart_qc = ncread(INSTACfile,'VART_QC');
instac.gdop_qc = ncread(INSTACfile,'GDOP_QC');
instac.ddns_qc = ncread(INSTACfile,'DDNS_QC');
instac.cspd_qc = ncread(INSTACfile,'CSPD_QC');

%%

%% SDC dataset

SDCfile = '/mnt/data/CNR/RADAR/DATI/Dati_HFR_TirLig/SDC/Totals/TV_HF_HFR-TirLig_TEST_201906.nc';

% Read time and convert it to Matlab time
sdc.time = ncread_cf_time(SDCfile,'TIME');

% Find the specific time index
iTime = find(sdc.time==datenum('21-Jun-2019 17:00:00'));

% Read variables
% Coordinate variables
sdc.latitude = ncread(SDCfile,'LATITUDE');
sdc.longitude = ncread(SDCfile,'LONGITUDE');
sdc.depth = ncread(SDCfile,'DEPTH');

% Data variables
sdc.ewct = ncread(SDCfile,'EWCT',[1,1,1,iTime],[length(sdc.longitude),length(sdc.latitude),length(sdc.depth),1]);
sdc.nsct = ncread(SDCfile,'NSCT',[1,1,1,iTime],[length(sdc.longitude),length(sdc.latitude),length(sdc.depth),1]);
sdc.ewcs = ncread(SDCfile,'EWCS',[1,1,1,iTime],[length(sdc.longitude),length(sdc.latitude),length(sdc.depth),1]);
sdc.nscs = ncread(SDCfile,'NSCS',[1,1,1,iTime],[length(sdc.longitude),length(sdc.latitude),length(sdc.depth),1]);
sdc.ccov = ncread(SDCfile,'CCOV',[1,1,1,iTime],[length(sdc.longitude),length(sdc.latitude),length(sdc.depth),1]);
sdc.gdop = ncread(SDCfile,'GDOP',[1,1,1,iTime],[length(sdc.longitude),length(sdc.latitude),length(sdc.depth),1]);
sdc.narx = ncread(SDCfile,'NARX',iTime,1);
sdc.natx = ncread(SDCfile,'NATX',iTime,1);
sdc.sltr = ncread(SDCfile,'SLTR',[1,iTime],[maxSite_dim,1]);
sdc.slnr = ncread(SDCfile,'SLNR',[1,iTime],[maxSite_dim,1]);
sdc.sltt = ncread(SDCfile,'SLTT',[1,iTime],[maxSite_dim,1]);
sdc.slnt = ncread(SDCfile,'SLNT',[1,iTime],[maxSite_dim,1]);
sdc.scdr = ncread(SDCfile,'SCDR',[1,1,iTime],[string15_dim,maxSite_dim,1]);
sdc.scdt = ncread(SDCfile,'SCDT',[1,1,iTime],[string15_dim,maxSite_dim,1]);

% QC variables
sdc.time_seadatanet_qc = ncread(SDCfile,'TIME_SEADATANET_QC',iTime,1);
sdc.position_seadatanet_qc = ncread(SDCfile,'POSITION_SEADATANET_QC',[1,1,1,iTime],[length(sdc.longitude),length(sdc.latitude),length(sdc.depth),1]);
sdc.depth_seadatanet_qc = ncread(SDCfile,'DEPTH_SEADATANET_QC',iTime,1);
sdc.qcflag = ncread(SDCfile,'QCflag',[1,1,1,iTime],[length(sdc.longitude),length(sdc.latitude),length(sdc.depth),1]);
sdc.vart_qc = ncread(SDCfile,'VART_QC',[1,1,1,iTime],[length(sdc.longitude),length(sdc.latitude),length(sdc.depth),1]);
sdc.gdop_qc = ncread(SDCfile,'GDOP_QC',[1,1,1,iTime],[length(sdc.longitude),length(sdc.latitude),length(sdc.depth),1]);
sdc.ddns_qc = ncread(SDCfile,'DDNS_QC',[1,1,1,iTime],[length(sdc.longitude),length(sdc.latitude),length(sdc.depth),1]);
sdc.cspd_qc = ncread(SDCfile,'CSPD_QC',[1,1,1,iTime],[length(sdc.longitude),length(sdc.latitude),length(sdc.depth),1]);


% Map QC variables to the SDC schema
sdc.time_seadatanet_qc(sdc.time_seadatanet_qc==49) = int8(1);
sdc.position_seadatanet_qc(sdc.position_seadatanet_qc==49) = int8(1);
sdc.depth_seadatanet_qc(sdc.depth_seadatanet_qc==49) = int8(1);

sdc.qcflag(sdc.qcflag==48) = int8(0);
sdc.qcflag(sdc.qcflag==49) = int8(1);
sdc.qcflag(sdc.qcflag==50) = int8(2);
sdc.qcflag(sdc.qcflag==51) = int8(3);
sdc.qcflag(sdc.qcflag==52) = int8(4);
sdc.qcflag(sdc.qcflag==56) = int8(8);

sdc.vart_qc(sdc.vart_qc==48) = int8(0);
sdc.vart_qc(sdc.vart_qc==49) = int8(1);
sdc.vart_qc(sdc.vart_qc==50) = int8(2);
sdc.vart_qc(sdc.vart_qc==51) = int8(3);
sdc.vart_qc(sdc.vart_qc==52) = int8(4);
sdc.vart_qc(sdc.vart_qc==56) = int8(8);

sdc.gdop_qc(sdc.gdop_qc==48) = int8(0);
sdc.gdop_qc(sdc.gdop_qc==49) = int8(1);
sdc.gdop_qc(sdc.gdop_qc==50) = int8(2);
sdc.gdop_qc(sdc.gdop_qc==51) = int8(3);
sdc.gdop_qc(sdc.gdop_qc==52) = int8(4);
sdc.gdop_qc(sdc.gdop_qc==56) = int8(8);

sdc.ddns_qc(sdc.ddns_qc==48) = int8(0);
sdc.ddns_qc(sdc.ddns_qc==49) = int8(1);
sdc.ddns_qc(sdc.ddns_qc==50) = int8(2);
sdc.ddns_qc(sdc.ddns_qc==51) = int8(3);
sdc.ddns_qc(sdc.ddns_qc==52) = int8(4);
sdc.ddns_qc(sdc.ddns_qc==56) = int8(8);

sdc.cspd_qc(sdc.cspd_qc==48) = int8(0);
sdc.cspd_qc(sdc.cspd_qc==49) = int8(1);
sdc.cspd_qc(sdc.cspd_qc==50) = int8(2);
sdc.cspd_qc(sdc.cspd_qc==51) = int8(3);
sdc.cspd_qc(sdc.cspd_qc==52) = int8(4);
sdc.cspd_qc(sdc.cspd_qc==56) = int8(8);

%%


%% Compare

timeDiff = sdc.time(iTime) - instac.time
latitudeDiff = sum(sdc.latitude - instac.latitude)
longitudeDiff = sum(sdc.longitude - instac.longitude)

ewctDiff = sum(sum(~isnan(sdc.ewct) - ~isnan(instac.ewct)))
nsctDiff = sum(sum(~isnan(sdc.nsct) - ~isnan(instac.nsct)))
ewcsDiff = sum(sum(~isnan(sdc.ewcs) - ~isnan(instac.ewcs)))
nscsDiff = sum(sum(~isnan(sdc.nscs) - ~isnan(instac.nscs)))
ccovDiff = sum(sum(~isnan(sdc.ccov) - ~isnan(instac.ccov)))
gdopDiff = sum(sum(~isnan(sdc.gdop) - ~isnan(instac.gdop)))

narxDiff = sdc.narx - instac.narx
natxDiff = sdc.natx - instac.natx

sltrDiff = sum(~isnan(sdc.sltr) - ~isnan(instac.sltr))
slnrDiff = sum(~isnan(sdc.slnr) - ~isnan(instac.slnr))
slttDiff = sum(~isnan(sdc.sltt) - ~isnan(instac.sltt))
slntDiff = sum(~isnan(sdc.slnt) - ~isnan(instac.slnt))

qcflagDiff = sum(sum(~isnan(sdc.qcflag) - ~isnan(instac.qcflag)))
vart_qcDiff = sum(sum(~isnan(sdc.vart_qc) - ~isnan(instac.vart_qc)))
gdop_qcDiff = sum(sum(~isnan(sdc.gdop_qc) - ~isnan(instac.gdop_qc)))
ddns_qcDiff = sum(sum(~isnan(sdc.ddns_qc) - ~isnan(instac.ddns_qc)))
cspd_qcDiff = sum(sum(~isnan(sdc.cspd_qc) - ~isnan(instac.cspd_qc)))

sdc.scdr
instac.scdr

sdc.scdt
instac.scdt

disp 'rock it';
