# GenZ_ADCP

## Project Information

ADCP data used to track the zooplankton dynamics in Lake Geneva (GenZ project).

Remote repository: https://github.com/tdoda/genz-adcp.git

Old remote repository: https://gitlab.renkulab.io/tomy.doda/genz-adcp.git

Renku session: https://renkulab.io/p/tomy.doda/genz-adcp/sessions/01JW8QSAD95XS2TJPJ91WEHV9G/start 

## Instruments

### ADCP

The ADCP is a Nortek Signature 1000 kHz: see [manual](notes/Signature_manuals/Signature_manual_2022.pdf) for more information. 

Typical setting for most of the GenZ campaigns: concurrent burst (echosounder + velocity) and averaging (velocity) plans. Details of each deployment can be found in the file `*.deploy` in Level 0 data.

The ADCP was either fixed on a floating platform attached to the Zodiac (test campaigns) or directly fixed to the LéM boat.

> **Data:**
> The ADCP data is stored in [data/Signature](data/Signature) (see [data section below](#folder-data)).

### GPS

GPS data was only collected for the most recent field campaigns. The GPS was a Garmin eTrex attached to the LéM boat.

> **Data:**
> The GPS data is stored in [data/GPS](data/GPS) (see [data section below](#folder-data)).

### RBR-CTD

<span style="color:red">To be added later.</span>

## Installation

You need to have [git](https://git-scm.com/downloads) and [git-lfs](https://git-lfs.github.com/) installed in order to successfully clone the repository.

- Clone the repository to your local machine using the command: 

 `git clone https://github.com/tdoda/genz-adcp.git `
 
 Note that the repository will be copied to your current working directory.

- Use Python 3 and install the requirements with:

 `pip install -r requirements.txt`

 The python version can be checked by running the command `python --version`. In case python is not installed or only an older version of it, it is recommend to install python through the anaconda distribution which can be downloaded [here](https://www.anaconda.com/products/individual).

## Usage

### Process new data

The raw and processed data are not uploaded to the remote git repository due to size limitations. To add it to your local repository:
1. Copy the raw ADCP and GPS data from the NAS sever to `data/*/*/Level0` (see the data structure in [data section below](#folder-data)).
2. Process the ADCP data by running the Python script `analysis/1-Export_Signature/main_export_Signature.py` (select the field campaigns of interest at the beginning of the script). The data is exported to netCDF files in Level1 and Level2 folders. 
3. Process the GPS data by running the Python script `analysis/2-Export_GPS/main_export_GPS.py` (select the field campaigns of interest at the beginning of the script). The data is exported to netCDF files in Level1 folder. 
4. Scripts performing further analysis (e.g., combining ADCP and GPS data) can be run in the `analysis` folder.

### Adapt/Extend data processing pipeline

New scripts can be added to the `analysis` and `figures` folders. Create a new folder for each new analysis step/new figure and follow the same structure.


## Folder `data`

The data is stored in the folder `data`. It is not uploaded to the remote git repository due to size limitations. The raw data can be added from the NAS server and processed with the scripts available in the `analysis` folder (see [scripts section below](#folder-analysis)).

The data is organized by data types and by dates of measurements for each data type. The data is structured in three levels:
- **Level 0**: Raw data downloaded from the instruments.
- **Level 1**: Raw data stored as netCDF files where attributes (such as sensors used, units, description of data, etc.) are added to the data. Column with quality flags are added: quality flag "1" indicates that the data point did not pass the 
quality checks and further investigation is needed, quality flag "0" indicates that no further investiagion is needed.
- **Level 2**: Processed data where quality flagged data has been removed and additional parameters have been computed.

### Subfolder `Signature` 

Contains the Signature ADCP data.

The `Level 0` folder contains the following files:
- `*.ad2cp`: all raw data (each ping), stored in binary format. 
- `*_avgd.ad2cp`: internally averaged data, stored in binary format. Only data where the correlation is above 50% will be included in the averaging, and the data will also include a percent-good value. Any bad data will be removed, or flagged with extreme values (e.g. -32.77 m/s). This averaging is only applicable to average mode; burst data will not be processed in this manner. "Percent Good" is defined as the percentage of data points above 50% correlation that go into an average ensemble. See [manual](notes/Signature_manuals/Signature_manual_2022.pdf), Sect. 5.5 "Recorder Data Retrieval".
- `*.meta`: metadata information following the json format.
- `*.deploy`: deployment file to be loaded/modified in the Nortek Deployment software to configure the Signature before the deployment.
- `*.cfg`: configuration file created by the Nortek Deployment software to save the Signature deployment settings.
- `*.ad2cp.index`: binary metadata file used by the Nortek Deployment software to structure the data.
- `Deployment_*.log`: instrument usage history provided by the Nortek Deployment software.
- `Instrument_*.log`: technical information about the Signature instrument provided by the Nortek Deployment software.
- `Sec_*.log`: communication events history between the computer and the Signature over both serial (RS‑232/422) and Ethernet connections, provided by the Nortek Deployment software.
- `System_*.log`: information about the firmware provided by the Nortek Deployment software.
- `NortekSupport_*.tar.gz`: compressed instrument support file containing debug information, provided by the Nortek Deployment software.

The processed data in `Level 2` have been quality-checked with both envass simple checks (numeric, bounds) and ADCP-specific checks (correlation, tilt). See the function `quality_flags` in `analysis/1-Export_Signature/adcp_signature.py` for more details.

### Subfolder `GPS` 

Contains the GPS data.

The `Level 0` folder contains the following files:
- `*.gpx`: GPS track.
- `*.meta`: metadata information following the json format.

The processed data in `Level 1` have not been quality-checked.

## Folder `analysis`
Contains the scripts processing the data and performing further analyses.

### Subfolder `1-Export_Signature`
ADCP Signature data export from raw data (Level 0) to processed data (Level 2).

- `main_export_Signature.py`: main script exporting new raw data (file *.ad2cp). Specify which campaign(s) to export at the beginning of the script. The script is designed to export the burst data that includes continuous echosounder measurements and repeated velocity bursts.
- `input_python.yaml`: name of the directories for the level-based data structure.
- `adcp_signature.py`: ADCP_Signature class with specific functions and parameters used to process the data.
- `functions.py`: additional functions used to process the data.
- `quality_checks_adcp.py`: ADCP-specific quality checks.
- `quality_assurance.json`: selected quality checks from the envass package.
- `log.txt`: log file produced during the data export process.

### Subfolder `2-Export_GPS`
GPS data export from raw data (Level 0) to processed data (Level 1).

- `main_export_GPS.py`: main script exporting new raw data (file *.gpx). Specify which campaign(s) to export at the beginning of the script. 
- `input_python.yaml`: name of the directories for the level-based data structure.
- `gps_device.py`: GPS class with specific functions and parameters used to process the data.

### Subfolder `3-Combine_Signature_GPS`
Combination of ADCP and GPS data.

- `combine_signature_GPS.py`: script combining the Level 2 Signature data with the Level 1 GPS data for a specific campaign.
- netCDF files produced by the script.

## Folder `functions`

Contains general Python functions used by various scripts in the `analysis`, `figures` and `notebooks` folder.

## Folder `notebooks`

Contains jupyter notebooks helping to visualize the processed data with simple figures.
- `vis_signature.ipynb`: plot the echosounder and burst velocity time series from the Level 2 Signature data.

## Folder `figures`

Contains scripts to make more advanced figures derived from the data analysis.

### Subfolder `1-Spatial_variability`

Spatial variability of backscattering and temperature.

- `plot_spatial_variability.py`: script plotting the spatial variability of backscattering and temperature.
- figures produced by the script as *.png files

## Collaborators 

- **Data collection**: Tomy Doda, Wyndel Sañoza
- **Scripts and analysis**: Tomy Doda, Wyndel Sañoza
- **Supervision**: Marie-Elodie Perga, Damien Bouffard

## Contact
- [Tomy Doda](mailto:tomy.doda@unil.ch)
- [Wyndel Sañoza](mailto:wyndel.sanoza@unil.ch)