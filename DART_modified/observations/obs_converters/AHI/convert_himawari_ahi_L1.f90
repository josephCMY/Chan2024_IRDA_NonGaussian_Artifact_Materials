! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


program convert_advanced_himawari_imager_ir_obs_L1

! Program to convert Advanced Himawari Imager IR observations
! from JAXA P-tree NetCDF file (aka, ncfile) into DART observation
! format.
! Note that the JAXA P-tree data is L1. No parallex errors or geolocation
! errors are accounted for.
! TODO: Look for AHI data sources with parallel & geolocation errors handled?

! This program is created using the GOES observation converter program
! as a template.


use types_mod,        only : r8, deg2rad, PI, MISSING_R8
use obs_sequence_mod, only : obs_sequence_type, write_obs_seq, &
                             static_init_obs_sequence, destroy_obs_sequence, &
                             init_obs_sequence, print_obs_seq_summary
use himawari_ahi_mod, only : himawari_ahi_map_type, himawari_load_ahi_map,  &
                             make_obs_sequence
use    utilities_mod, only : initialize_utilities, register_module, &
                             error_handler, finalize_utilities, E_ERR, E_MSG, &
                             find_namelist_in_file, check_namelist_read, &
                             do_nml_file, do_nml_term, set_filename_list, &
                             nmlfileunit

implicit none

! ----------------------------------------------------------------------
! Declare local parameters
! ----------------------------------------------------------------------

integer                     :: filecount
type(himawari_ahi_map_type) :: map
type(obs_sequence_type)     :: seq

integer :: io, iunit, ifile, ichannel

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'convert_himawari_ahi_L1.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

! ----------------------------------------------------------------------
! Declare namelist parameters
! ----------------------------------------------------------------------

integer, parameter :: MAXFILES = 512

character(len=256) :: raw_ncfile       = ''   ! Name of JAXA P-TREE FILE
character(len=256) :: outputfile         = ''

real(r8) :: lon1 =   0.0_r8,  &   !  lower longitude bound
            lon2 = 360.0_r8,  &   !  upper longitude bound
            lat1 = -90.0_r8,  &   !  lower latitude bound
            lat2 =  90.0_r8       !  upper latitude bound

integer  :: x_thin = 0
integer  :: y_thin = 0
integer  :: himawari_satellite_number = 9
logical  :: reject_dqf_1 = .true.
logical  :: verbose = .false.
integer  :: ahi_channel_num = -999

real(r8) :: obs_err   = MISSING_R8     ! the fixed obs error (std dev, in radiance units)
                                       ! TODO: make this more sophisticated
real(r8) :: vloc_pres_hPa = -1.0_r8    ! the fixed location to "place" this observation
                                       ! in hPa. Negative means disable and use VERTISUNDEF

namelist /convert_himawari_AHI_L1_nml/  raw_ncfile, &
                                        outputfile, &
                                        lon1, lon2, lat1, lat2, &
                                        x_thin, y_thin, himawari_satellite_number, &
                                        reject_dqf_1, obs_err, verbose, &
                                        vloc_pres_hPa, ahi_channel_num

! ----------------------------------------------------------------------
! start of executable program code
! ----------------------------------------------------------------------

call initialize_utilities('convert_himawari_AHI_L1')
call register_module(source,revision,revdate)

! Initialize the obs_sequence module ...

call static_init_obs_sequence()


!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file('input.nml', 'convert_himawari_AHI_L1_nml', iunit)
read(iunit, nml = convert_himawari_AHI_L1_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_himawari_AHI_L1_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=convert_himawari_AHI_L1_nml)
if (do_nml_term()) write(    *      , nml=convert_himawari_AHI_L1_nml)

if (obs_err == MISSING_R8) then
  call error_handler(E_ERR, source, 'An obs_err value was not specified in the namelist, but is required.',source,revision,revdate)
end if

if (vloc_pres_hPa < 0.0_r8) then
  call error_handler(E_MSG, source, 'No vertical pressure level defined, using VERTISUNDEF.',source,revision,revdate)
end if

if ( (ahi_channel_num < 6) .or. (ahi_channel_num > 15) ) then
  call error_handler(E_ERR, source, 'Invalid AHI channel number specified',source,revision,revdate)
end if



! ! when this routine returns, the l1_files variable will have
! ! all the filenames, regardless of which way they were specified.
! filecount = set_filename_list(l1_files, l1_file_list, "convert_himawari_AHI_L1")


! ! for each input file
! do ifile=1, filecount

    ! read from netCDF file into a derived type that holds all the information
    call himawari_load_AHI_map(raw_ncfile, map, ahi_channel_num)
    ! call himawari_load_AHI_map(l1_files(ifile), map, ahi_channel_num)

    ! convert derived type information to DART sequence
    call make_obs_sequence(seq, map, lon1, lon2, lat1, lat2, &
                          x_thin, y_thin, himawari_satellite_number, reject_dqf_1, obs_err, &
                          vloc_pres_hPa)

    ! write the sequence to a disk file
    call write_obs_seq(seq, outputfile)

    if (verbose) call print_obs_seq_summary(seq,raw_ncfile)

    ! release the sequence memory
    call destroy_obs_sequence(seq)
! enddo

call error_handler(E_MSG, source, 'Finished successfully.',source,revision,revdate)
call finalize_utilities()

end program convert_advanced_himawari_imager_ir_obs_L1

