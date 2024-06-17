EUMETSAT Meteosat Second Generation (MSG) Spinning Enhanced Visible Infra-Red Imager (SEVIRI) Level 1b Radiances
================================================================================================================

The data is available from EUMETSAT in a satellite-native file format. 
That native file format must be converted to NetCDF before applying the
obs converter. A Python-based native-to-netCDF conversion script is 
provided in ``observations/obs_converter/SEVIRI/work``.

The **convert_metesat_SEVIRI_L1b** program converts SEVIRI Level 1b 
Radiances in netCDF format to a DART observation sequence file with
``METEOSAT_8_IODC_SEVIRI_RADIANCE`` observations (there is a namelist 
option to select other MSG satellites, which will have the appropriate
observation type).


Specifying a vertical location
------------------------------

Jeff Steward added (PR 48) the capability to specify a vertical location
if desired. This allows for localization in the vertical. 

  It’s sometimes helpful, even though definitely wrong from a theoretical
  standpoint, to give a vertical location to satellite observations
  (which are integrated quantities). This has been an issue with
  observation-space localization for some time, and this is the standard
  workaround pioneered by Lili Lei and Jeff Whittaker.



A short description of the namelist options
-------------------------------------------

This table is meant to familiarize you with some of the options
available. Until we fully implement automatic documentation generation,
you would be well advised to familiarize yourself with the code. This is
not the full list of namelist variables …

+-------------------------+---------------------+---------------------+
| variable                | default             | meaning             |
+=========================+=====================+=====================+
| x_thin                  | 0                   | Skip this many per  |
|                         |                     | X scan.             |
+-------------------------+---------------------+---------------------+
| y_thin                  | 0                   | Skip this many per  |
|                         |                     | Y scan.             |
+-------------------------+---------------------+---------------------+
| meteosat_num            | 8                   | MSG Satellite      |
|                         |                     | number.             |
+-------------------------+---------------------+---------------------+
| reject_dqf_1            | .true.              | Bad scan rejection  |
|                         |                     | criteria. If .true. |
|                         |                     | and DQF /= 0, the   |
|                         |                     | scan is rejected.   |
|                         |                     | If .false. any DQF  |
|                         |                     | > 1 rejects the     |
|                         |                     | scan.               |
+-------------------------+---------------------+---------------------+
| verbose                 | .false.             | Run-time output     |
|                         |                     | verbosity           |
+-------------------------+---------------------+---------------------+
| obs_err                 | MISSING_R8          | The observation     |
|                         |                     | error standard      |
|                         |                     | deviation (std dev, |
|                         |                     | in radiance units)  |
|                         |                     | TODO: make this     |
|                         |                     | more sophisticated. |
|                         |                     | You must supply a   |
|                         |                     | value other than    |
|                         |                     | MISSING_R8. Be      |
|                         |                     | aware that the      |
|                         |                     | observation         |
|                         |                     | sequence files      |
|                         |                     | convert this to a   |
|                         |                     | variance.           |
+-------------------------+---------------------+---------------------+
| vloc_pres_hPa           | -1.0                | The vertical        |
|                         |                     | location of this    |
|                         |                     | observation (hPa).  |
|                         |                     | A negative means    |
|                         |                     | there is no         |
|                         |                     | vertical location   |
|                         |                     | (which is typical   |
|                         |                     | for a               |
|                         |                     | ve                  |
|                         |                     | rtically-integrated |
|                         |                     | quantity).          |
+-------------------------+---------------------+---------------------+
