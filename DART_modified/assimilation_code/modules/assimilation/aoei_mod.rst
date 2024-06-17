MODULE aoei_mod
===============

Overview
--------

This module implements the Adaptive Observation Error Inflation (AOEI) scheme of Minamide and Zhang (2017, Monthly Weather Review). 


How to use this module
-----------------------
By default, the module only applies AOEI on observations whose `QUANTITY` is `QTY_BRIGHTNESS_TEMPERATURE`.

To add other observation `QUANTITY`s, the steps are:

1) Open up `aoei_mod.f90`

2) Add the `QTY` codes to the `obs_kind_mod` use statement (line 42 of `aoei_mod.f90`)

3) Put those `QTY` codes into the `aoei_obs_qty_list` array in the subroutine `init_aoei`.




Note on AOEI implementation in DART
----------------------------------

This AOEI implementation differs from the original Minamide-Zhang implementation in the Pennsylvania State University EnKF (PSU-EnKF) system.
The PSU-EnKF system calls AOEI within the sequential assimilation loop. 
In other words, changing the order of the observations can change the behavior of AOEI.
To avoid this complication, DART's implementation of the AOEI happens outsideof the sequential assimilation loop.

