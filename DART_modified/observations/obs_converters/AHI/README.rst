JAXA Himawari Advanced Himawari Imager (AHI) Level 1 Radiances
==============================================================

This README file is a shortened version of the README.rst for GOES observations.

The data are available from (requires registration)
`JAXA P-TREE <https://www.eorc.jaxa.jp/ptree/index.html>`__

The **convert_Himawari_AHI_L1** program converts AHI Level 1 Radiances in
netCDF format to a DART observation sequence file with
``HIMAWARI_9_AHI_RADIANCE`` observations (there is a namelist option to
select other GOES satellites, which will have the appropriate
observation type).

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
| goes_num                | 16                  | GOES Satellite      |
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
