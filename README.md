# GenZ_ADCP

## Project Information

ADCP data used to track the zooplankton dynamics in Lake Geneva (GenZ project).
Remote repository: https://gitlab.renkulab.io/tomy.doda/genz-adcp.git 

## Sensors

 <span style="color:red">Operation: combined burst-averaging plan</span>. 

## Installation

:warning You need to have [git](https://git-scm.com/downloads) and [git-lfs](https://git-lfs.github.com/) installed in order to successfully clone the repository.

- Clone the repository to your local machine using the command: 

 `git clone https://gitlab.renkulab.io/tomy.doda/genz-adcp.git `
 
 Note that the repository will be copied to your current working directory.

- Use Python 3 and install the requirements with:

 `pip install -r requirements.txt`

 The python version can be checked by running the command `python --version`. In case python is not installed or only an older version of it, it is recommend to install python through the anaconda distribution which can be downloaded [here](https://www.anaconda.com/products/individual).

 <span style="color:red">Explain how to install the develop branch of mhkit</span>. 

## Usage

### Process new data

### Adapt/Extend data processing pipeline


## Data

The data can be found in the folder `data`. The data is organized by data types and by dates of measurements for each data type. The data is structured in three levels:
- **Level 0**: Raw data downloaded from the instruments.
- **Level 1**: Raw data stored as netCDF files where attributes (such as sensors used, units, description of data, etc.) are added to the data. Column with quality flags are added: quality flag "1" indicates that the data point did not pass the 
quality checks and further investigation is needed, quality flag "0" indicates that no further investiagion is needed.
- **Level 2**: Processed data where quality flagged data has been removed and additional parameters have been computed.

### Signature

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


## Quality assurance

Quality checks include but are not limited to range validation, data type checking and flagging missing data.
<span style="color:red">Give more information on quality checks.</span>.

## Scripts
The scripts processing the data can be found in the folder `analysis`.
<span style="color:red">Give information on which scripts do what here! The script reads the *.ad2cp file (burst data = continuous echo measurements + 3 min long velocity burst every 10 min</span>.

## Data visualization
Jupyter notebooks helping to visualize the processed data with simple figures can be found in the folder `notebooks`.

## Figures
More complex figures derived from the data analysis can be found in the folder `figures` with the scripts used to create them.