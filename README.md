# HFR_Node__SDC_Processing
Matlab scripts for the operational workflow for generating SeaDataCloud aggregated datasets and CDIs at the European HFR Node. Tools for the centralized processing

These applications are written in Matlab language and they are based on 'nctoolbox' toolbox and on the Mikado software (https://www.seadatanet.org/Software/MIKADO). The architecture of the workflow is based on a MySQL database containing information about data and metadata. The applications are designed for High Frequency Radar (HFR) data management according to the European HFR node processing workflow, thus generating radial and total velocity aggregated datasets in netCDF format according to the European standard data and metadata model for HFR current data and the related Common Data Index (CDI) files according to SeaDataNet erquirements.

The database is composed by the following tables:

	* network_tb: it contains the general information about the HFR network producing the radial and total files. These information will be used for the metadata content of the netCDF files.
	* station_tb: it contains the general information about the radar sites belonging to each HFR network producing the radial and total files. These information will be used for the metadata content of the netCDF files.
	* radial_CDIconf_tb: it contains all the information to be queried by the Mikado configuration files for automatic generation of the CDIs for the aggregated radial datasets.
	* radial_SDCnetCDF_tb: it contains information about the generated aggregated radial datasets.
	* total_CDIconf_tb: it contains all the information to be queried by the Mikado configuration files for automatic generation of the CDIs for the aggregated total datasets.
	* total_SDCnetCDF_tb: it contains information about the generated aggregated total datasets.

The applications are intended to:

	* load radial site and network information from the database tables network_tb and station_tb;
	* connect to the EU HFR Node THREDDS catalog and aggregate data via OpenDAP;
	* generate the netCDF aggregated datasets for radials and totals in netCDF format according to the European standard data and metadata model for HFR current data;
	* generate the CDIs for the aggregated datasets for radials and totals according to SeaDataNet erquirements.

General information for the tables network_tb and station_tb are loaded onto the database via a webform to be filled by the data providers. The webform is available at http://150.145.136.36/index.php

All generated radial and total netCDF datasets are quality controlled according the the QC tests defined as standard for the European HFR node and for the data distribution on CMEMS-INSTAC and SeaDataNet platforms.

The whole workflow is intended to run automatically to periodically generate aggregated datasets and the related CDIs. The wrapper CP_EU_HFR_Node_SDC_dataset_builer.m sets the aggregateion time span and lauches the processing applications.

The applications SDC_totals.m and SDC_radials.m call the aggregation functions, write metadata information to the database and call the Mikado application for generating the CDIs for total and radial data respectively.

The applications total_netCDF_aggregation_v211.m and radial_netCDF_aggregation_v211.m read aggregated data from the THREDDS catalog via OpenDAP and perform the temporal aggregation for total and radial data respectively.

The applications SDC_CDI_totals_metadata_to_DB.m and SDC_CDI_radials_metadata_to_DB.m write to the database all the information to be queried by the Mikado configuration files for automatic generation of the CDIs for the aggregated radial datasets for total and radial data respectively.

The folder Mikado_conf_file contains the xml configuration files to be used by Mikado for the automatic generation of the CDIs for the aggregated radial and total datasets. Please refer to the Mikado manual (https://www.seadatanet.org/content/download/651/file/sdn_Mikado_UserManual_V3.5.3.pdf) for the description of the configuration files and the automatic usage of Mikado.

The required toolboxes are:

    HFR_Progs-2.1.2 (https://github.com/rowg/hfrprogs);
    M_Map (https://www.eoas.ubc.ca/~rich/map.html);
    GSHHS (http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/);
    Nctoolbox-1.1.3 (https://github.com/nctoolbox/nctoolbox);
    mysql-connector-java-5.1.17 driver (https://mvnrepository.com/artifact/mysql/mysql-connector-java/5.1.17);
    Rdir (http://www.mathworks.com/matlabcentral/fileexchange/19550-recursive-directory-listing);
    uniqueStrCell (https://www.mathworks.com/matlabcentral/fileexchange/50476-unique-for-cell-array-of-string).
    
Cite as:
Lorenzo Corgnati. (2019). EU_HFR_NODE_SDC_Processing. Zenodo. https://doi.org/10.5281/zenodo.3855468


Author: Lorenzo Corgnati

Date: September 13, 2019

E-mail: lorenzo.corgnati@sp.ismar.cnr.it
