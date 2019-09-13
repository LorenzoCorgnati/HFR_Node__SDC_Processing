%% SDC_radial_dataset_check.m
% This script compares the SDC radial aggregated dataset with some hourly files

% Author: Lorenzo Corgnati
% Date: July 24, 2019

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

warning('off', 'all');

clear all
close all
clc

% Set non physical dimensions
maxSite_dim = 1;
string15_dim = 15;

% Setup netCDF toolbox
setup_nctoolbox;


%% CMEMS-INSTAC data

INSTACfile = '/home/lorenz/Downloads/HFR-TirLig-PCOR_2019_07_16_1900.nc';

% Read time and convert it to Matlab time
instac.time = ncread_cf_time(INSTACfile,'TIME');

% Read variables
% Coordinate variables
instac.bear = ncread(INSTACfile,'BEAR');
instac.rnge = ncread(INSTACfile,'RNGE');
instac.depth = ncread(INSTACfile,'DEPH');

instac.latitude = ncread(INSTACfile,'LATITUDE');
instac.longitude = ncread(INSTACfile,'LONGITUDE');

% Data variables
instac.rdva = ncread(INSTACfile,'RDVA');
instac.drva = ncread(INSTACfile,'DRVA');
instac.ewct = ncread(INSTACfile,'EWCT');
instac.nsct = ncread(INSTACfile,'NSCT');
instac.espc = ncread(INSTACfile,'ESPC');
instac.etmp = ncread(INSTACfile,'ETMP');
instac.maxv = ncread(INSTACfile,'MAXV');
instac.minv = ncread(INSTACfile,'MINV');
instac.ersc = ncread(INSTACfile,'ERSC');
instac.ertc = ncread(INSTACfile,'ERTC');
instac.xdst = ncread(INSTACfile,'XDST');
instac.ydst = ncread(INSTACfile,'YDST');
instac.sprc = ncread(INSTACfile,'SPRC');
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
instac.owtr_qc = ncread(INSTACfile,'OWTR_QC');
instac.mdfl_qc = ncread(INSTACfile,'MDFL_QC');
instac.vart_qc = ncread(INSTACfile,'VART_QC');
instac.cspd_qc = ncread(INSTACfile,'CSPD_QC');
instac.avrb_qc = ncread(INSTACfile,'AVRB_QC');
instac.rdct_qc = ncread(INSTACfile,'RDCT_QC');

%%

%% SDC dataset

SDCfile = '/mnt/data/CNR/RADAR/DATI/Dati_HFR_TirLig/SDC/Radials/PCOR_TEST/RV_HF_HFR-TirLig_TEST-PCOR_TEST_201907.nc';

% Read time and convert it to Matlab time
sdc.time = ncread_cf_time(SDCfile,'TIME');

% Find the specific time index
iTime = find(sdc.time==datenum('28-Jul-2019 10:00:00'));

% Read variables
% Coordinate variables
sdc.bear = ncread(SDCfile,'BEAR');
sdc.rnge = ncread(SDCfile,'RNGE');
sdc.depth = ncread(SDCfile,'DEPTH');

sdc.latitude = ncread(SDCfile,'LATITUDE');
sdc.longitude = ncread(SDCfile,'LONGITUDE');

% Data variables
sdc.rdva = ncread(SDCfile,'RDVA',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.drva = ncread(SDCfile,'DRVA',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.ewct = ncread(SDCfile,'EWCT',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.nsct = ncread(SDCfile,'NSCT',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.espc = ncread(SDCfile,'ESPC',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.etmp = ncread(SDCfile,'ETMP',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.maxv = ncread(SDCfile,'MAXV',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.minv = ncread(SDCfile,'MINV',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.ersc = ncread(SDCfile,'ERSC',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.ertc = ncread(SDCfile,'ERTC',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.xdst = ncread(SDCfile,'XDST');
sdc.ydst = ncread(SDCfile,'YDST');
sdc.sprc = ncread(SDCfile,'SPRC',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
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
sdc.position_seadatanet_qc = ncread(SDCfile,'POSITION_SEADATANET_QC',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.depth_seadatanet_qc = ncread(SDCfile,'DEPTH_SEADATANET_QC',iTime,1);

sdc.qcflag = ncread(SDCfile,'QCflag',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.owtr_qc = ncread(SDCfile,'OWTR_QC',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.mdfl_qc = ncread(SDCfile,'MDFL_QC',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.vart_qc = ncread(SDCfile,'VART_QC',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.cspd_qc = ncread(SDCfile,'CSPD_QC',[1,1,1,iTime],[length(sdc.rnge),length(sdc.bear),length(sdc.depth),1]);
sdc.avrb_qc = ncread(SDCfile,'AVRB_QC',iTime,1);
sdc.rdct_qc = ncread(SDCfile,'RDCT_QC',iTime,1);

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

sdc.owtr_qc(sdc.owtr_qc==48) = int8(0);
sdc.owtr_qc(sdc.owtr_qc==49) = int8(1);
sdc.owtr_qc(sdc.owtr_qc==50) = int8(2);
sdc.owtr_qc(sdc.owtr_qc==51) = int8(3);
sdc.owtr_qc(sdc.owtr_qc==52) = int8(4);
sdc.owtr_qc(sdc.owtr_qc==56) = int8(8);

sdc.mdfl_qc(sdc.mdfl_qc==48) = int8(0);
sdc.mdfl_qc(sdc.mdfl_qc==49) = int8(1);
sdc.mdfl_qc(sdc.mdfl_qc==50) = int8(2);
sdc.mdfl_qc(sdc.mdfl_qc==51) = int8(3);
sdc.mdfl_qc(sdc.mdfl_qc==52) = int8(4);
sdc.mdfl_qc(sdc.mdfl_qc==56) = int8(8);

sdc.cspd_qc(sdc.cspd_qc==48) = int8(0);
sdc.cspd_qc(sdc.cspd_qc==49) = int8(1);
sdc.cspd_qc(sdc.cspd_qc==50) = int8(2);
sdc.cspd_qc(sdc.cspd_qc==51) = int8(3);
sdc.cspd_qc(sdc.cspd_qc==52) = int8(4);
sdc.cspd_qc(sdc.cspd_qc==56) = int8(8);

sdc.avrb_qc(sdc.avrb_qc==48) = int8(0);
sdc.avrb_qc(sdc.avrb_qc==49) = int8(1);
sdc.avrb_qc(sdc.avrb_qc==50) = int8(2);
sdc.avrb_qc(sdc.avrb_qc==51) = int8(3);
sdc.avrb_qc(sdc.avrb_qc==52) = int8(4);
sdc.avrb_qc(sdc.avrb_qc==56) = int8(8);

sdc.rdct_qc(sdc.rdct_qc==48) = int8(0);
sdc.rdct_qc(sdc.rdct_qc==49) = int8(1);
sdc.rdct_qc(sdc.rdct_qc==50) = int8(2);
sdc.rdct_qc(sdc.rdct_qc==51) = int8(3);
sdc.rdct_qc(sdc.rdct_qc==52) = int8(4);
sdc.rdct_qc(sdc.rdct_qc==56) = int8(8);

%%


%% Compare

timeDiff = sdc.time(iTime) - instac.time
bearDiff = sum(sdc.bear - instac.bear)
rngeDiff = sum(sdc.rnge - instac.rnge)

latitudeDiff = sum(sum(sdc.latitude - instac.latitude))
longitudeDiff = sum(sum(sdc.longitude - instac.longitude))

rdvaDiff = sum(sum(~isnan(sdc.rdva) - ~isnan(instac.rdva)))
drvaDiff = sum(sum(~isnan(sdc.drva) - ~isnan(instac.drva)))
ewctDiff = sum(sum(~isnan(sdc.ewct) - ~isnan(instac.ewct)))
nsctDiff = sum(sum(~isnan(sdc.nsct) - ~isnan(instac.nsct)))
espcDiff = sum(sum(~isnan(sdc.espc) - ~isnan(instac.espc)))
etmpDiff = sum(sum(~isnan(sdc.etmp) - ~isnan(instac.etmp)))
maxvDiff = sum(sum(~isnan(sdc.maxv) - ~isnan(instac.maxv)))
minvDiff = sum(sum(~isnan(sdc.minv) - ~isnan(instac.minv)))
erscDiff = sum(sum(~isnan(sdc.ersc) - ~isnan(instac.ersc)))
ertcDiff = sum(sum(~isnan(sdc.ertc) - ~isnan(instac.ertc)))
xdstDiff = sum(sum(~isnan(sdc.xdst) - ~isnan(instac.xdst)))
ydstDiff = sum(sum(~isnan(sdc.ydst) - ~isnan(instac.ydst)))
sprcDiff = sum(sum(~isnan(sdc.sprc) - ~isnan(instac.sprc)))

narxDiff = sdc.narx - instac.narx
natxDiff = sdc.natx - instac.natx

sltrDiff = sum(~isnan(sdc.sltr) - ~isnan(instac.sltr))
slnrDiff = sum(~isnan(sdc.slnr) - ~isnan(instac.slnr))
slttDiff = sum(~isnan(sdc.sltt) - ~isnan(instac.sltt))
slntDiff = sum(~isnan(sdc.slnt) - ~isnan(instac.slnt))

qcflagDiff = sum(sum(~isnan(sdc.qcflag) - ~isnan(instac.qcflag)))
vart_qcDiff = sum(sum(~isnan(sdc.vart_qc) - ~isnan(instac.vart_qc)))
owtr_qcDiff = sum(sum(~isnan(sdc.owtr_qc) - ~isnan(instac.owtr_qc)))
mdfl_qcDiff = sum(sum(~isnan(sdc.mdfl_qc) - ~isnan(instac.mdfl_qc)))
cspd_qcDiff = sum(sum(~isnan(sdc.cspd_qc) - ~isnan(instac.cspd_qc)))
avrb_qcDiff = sum(sum(~isnan(sdc.avrb_qc) - ~isnan(instac.avrb_qc)))
rdct_qcDiff = sum(sum(~isnan(sdc.rdct_qc) - ~isnan(instac.rdct_qc)))

sdc.scdr
instac.scdr

sdc.scdt
instac.scdt

disp 'rock it';
