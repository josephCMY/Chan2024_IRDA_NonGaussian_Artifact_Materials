! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


! This mod is written by Man-Yau Chan by mimicking a mod written by Jeffrey
! Steward (see GOES observation converter)


module himawari_ahi_mod

use types_mod,     only : r4, r8, deg2rad, rad2deg, PI, digits12, &
                          MISSING_R8

use utilities_mod, only : error_handler, E_MSG, E_ERR, &
                          is_longitude_between, register_module

use time_manager_mod, only : time_type, get_date, set_date,            &
                             get_time, set_time, set_calendar_type,    &
                             GREGORIAN, print_date, print_time,        &
                             operator(+)

use obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                             set_obs_values, obs_sequence_type,              &
                             obs_type, set_copy_meta_data, set_qc_meta_data

use     location_mod, only : location_type, set_location, VERTISUNDEF, &
                             VERTISPRESSURE, get_location

use     obs_kind_mod,  only : get_index_for_type_of_obs

use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

use obs_def_rttov_mod, only : set_visir_metadata, &
                              get_rttov_option_logical

use netcdf_utilities_mod, only : nc_check

use netcdf

implicit none

private

public :: himawari_ahi_map_type, himawari_load_ahi_map, make_obs_sequence

!>@todo FIXME ... we should be using a lot more of the netcdf_utilities interfaces.
!>                they greatly simplify the code.

type himawari_ahi_map_type
   character(:), allocatable :: filename
   integer :: channel
   integer :: nx   ! Number of longitude points
   integer :: ny   ! Number of latitude points

   real(r8),    allocatable :: lat1d(:)     ! vector of observation latitudes 
   real(r8),    allocatable :: lon1d(:)     ! vector of observation longitudes
   

   ! Observation date (days since 1858-11-17 00:00:00)
   real(r8) :: obs_date


   ! Data quality flag
   integer(1),  allocatable :: dqf(:,:) ! data quality flag. 0 = good, 1 is conditionally okay, 2-4 is bad.
   
   ! Variables to read in
   real(r8),    allocatable :: bt(:,:)      ! Brightness temperature values
   real(r8),    allocatable :: sat_az(:,:)  ! Satellite azimuth angles
   real(r8),    allocatable :: sat_ze(:,:)  ! Satellite zenith  angles
   real(r8),    allocatable :: sun_az(:,:)  ! Solar azimuth angles
   real(r8),    allocatable :: sun_ze(:,:)  ! Solar zenith  angles
   real(r8),    allocatable :: lat(:,:)     ! Satellite pixel latitudes
   real(r8),    allocatable :: lon(:,:)     ! Satellite pixel longitudes

   ! Buffers used to interpret compressed float data in netCDF files
   integer,     allocatable :: raw_vals(:,:)
   real(r4)                 :: nc_scale
   real(r4)                 :: nc_offset

end type

character(len=512) :: msgstring

integer :: fid, varid

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'himawari_ahi_L1.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical, save :: module_initialized = .false.

contains

subroutine initialize_module

call register_module(source, revision, revdate)

call set_calendar_type(GREGORIAN)

module_initialized = .true.

end subroutine initialize_module

! ----------------------------------------------------------------------------

subroutine himawari_load_ahi_map(l1_file, map, ir_channel_num)

character(len=*),    intent(in)    :: l1_file
type(himawari_ahi_map_type), intent(inout) :: map
integer,             intent(in)    :: ir_channel_num

logical,          parameter :: little_endian = .true.

character(len=*), parameter :: routine = 'himawari_load_ahi_map'

character(len=512) :: string1
character(len=32) :: ir_channel_name

integer :: i, j
real(digits12) :: a, b, c, r_fac, xv, yv, r_s
real(digits12) :: sx, sy, sz
real(digits12) :: time_r8

real(digits12) :: obs_date_st
real(digits12) :: obs_date_ed

if ( .not. module_initialized ) call initialize_module

write(string1,*) 'Now loading L1 AHI data from file ',trim(l1_file)
call error_handler(E_MSG, routine, string1, source, revision, revdate)

allocate(character(len=len(l1_file)) :: map%filename)
map%filename = l1_file


! Open observation file
call nc_check( nf90_open(l1_file, nf90_nowrite, fid), 'file open', l1_file)


! Load observation longitude data vector from file
call nc_check( nf90_inq_dimid(fid, "longitude", varid),        'inq dimid longitude', l1_file)
call nc_check( nf90_inquire_dimension(fid, varid, len=map%nx), 'inq dim   longitude', l1_file)
allocate(map%lon1d( map%nx ))
call nc_check( nf90_get_var(fid, varid, map%lon1d),            'get var   longitude', l1_file)

! Load observation latitude data vector from file
call nc_check( nf90_inq_dimid(fid, "latitude", varid),         'inq dimid latitude', l1_file)
call nc_check( nf90_inquire_dimension(fid, varid, len=map%ny), 'inq dim   latitude', l1_file)
allocate(map%lat1d( map%ny ))
call nc_check( nf90_get_var(fid, varid, map%lat1d),            'get var   latitude', l1_file)



! Allocate 2d arrays to hold satellite data and decompression buffer
allocate( map%bt       ( map%nx, map%ny ) )
allocate( map%sat_az   ( map%nx, map%ny ) )
allocate( map%sat_ze   ( map%nx, map%ny ) )
allocate( map%sun_az   ( map%nx, map%ny ) )
allocate( map%sun_ze   ( map%nx, map%ny ) )
allocate( map%lat      ( map%nx, map%ny ) )
allocate( map%lon      ( map%nx, map%ny ) )
allocate( map%raw_vals ( map%nx, map%ny ) )


! Generate 2d array of longitude
do i=1,map%nx
   ! Longtiudes > 180 need special handling
   if ( map%lon1d(i) > 180. ) then
      map%lon(i,:) = -180. + (map%lon1d(i)-180)
   else
      map%lon(i,:) = map%lon1d(i)
   end if
end do

! Generate 2d array of latitude
do j=1,map%ny
    map%lat(:,j) = map%lat1d(j)
end do



! Determine ncfile variable name that corresponds to BT from desired channel
write(ir_channel_name, '(a,i0.2)') 'tbb_', ir_channel_num

! Store channel number into map 
map%channel = ir_channel_num


! Read brightness temperatures from ncfiles
map%raw_vals = 0
call nc_check( nf90_inq_varid(fid, ir_channel_name, varid) ,            'inq varid BT', l1_file)
call nc_check( nf90_get_var(fid, varid, map%raw_vals),                  'get var   BT', l1_file)
call nc_check( nf90_get_att(fid, varid, 'scale_factor', map%nc_scale) , 'get_att   BT scale', l1_file)
call nc_check( nf90_get_att(fid, varid, 'add_offset',   map%nc_offset), 'get_att   BT offset', l1_file)
map%bt = real(map%raw_vals)*map%nc_scale + map%nc_offset


! Read in the satellite azimuth angle
map%raw_vals = 0
call nc_check( nf90_inq_varid(fid, 'SAA', varid) ,                      'inq varid sat_az', l1_file)
call nc_check( nf90_get_var(fid, varid, map%raw_vals),                  'get var   sat_az', l1_file)
call nc_check( nf90_get_att(fid, varid, 'scale_factor', map%nc_scale) , 'get_att   sat_az scale', l1_file)
call nc_check( nf90_get_att(fid, varid, 'add_offset', map%nc_offset) ,  'get_att   sat_az offset', l1_file)
map%sat_az = real(map%raw_vals)*map%nc_scale + map%nc_offset

! Read in the satellite zenith angle
map%raw_vals = 0
call nc_check( nf90_inq_varid(fid, 'SAZ', varid) ,                      'inq varid sat_ze', l1_file)
call nc_check( nf90_get_var(fid, varid, map%raw_vals),                  'get var   sat_ze', l1_file)
call nc_check( nf90_get_att(fid, varid, 'scale_factor', map%nc_scale) , 'get_att   sat_ze scale', l1_file)
call nc_check( nf90_get_att(fid, varid, 'add_offset', map%nc_offset) ,  'get_att   sat_ze offset', l1_file)
map%sat_ze = real(map%raw_vals)*map%nc_scale + map%nc_offset

! Read in the solar azimuth angle
map%raw_vals = 0
call nc_check( nf90_inq_varid(fid, 'SOA', varid) ,                      'inq varid sun_az', l1_file)
call nc_check( nf90_get_var(fid, varid, map%raw_vals),                  'get var   sun_az', l1_file)
call nc_check( nf90_get_att(fid, varid, 'scale_factor', map%nc_scale) , 'get_att   sun_az scale', l1_file)
call nc_check( nf90_get_att(fid, varid, 'add_offset', map%nc_offset) ,  'get_att   sun_az offset', l1_file)
map%sun_az = real(map%raw_vals)*map%nc_scale + map%nc_offset


! Read in the solar zenith angle
map%raw_vals = 0
call nc_check( nf90_inq_varid(fid, 'SOZ', varid) ,                      'inq varid sun_ze', l1_file)
call nc_check( nf90_get_var(fid, varid, map%raw_vals),                  'get var   sun_ze', l1_file)
call nc_check( nf90_get_att(fid, varid, 'scale_factor', map%nc_scale) , 'get_att   sun_ze scale', l1_file)
call nc_check( nf90_get_att(fid, varid, 'add_offset', map%nc_offset) ,  'get_att   sun_ze offset', l1_file)
map%sun_ze = real(map%raw_vals)*map%nc_scale + map%nc_offset

! Read in observation date 
call nc_check( nf90_inq_varid(fid, 'start_time', varid) ,              'inq varid start_time', l1_file)
call nc_check( nf90_get_var(fid, varid, obs_date_st),                  'get var   start_time', l1_file)
call nc_check( nf90_inq_varid(fid, 'end_time', varid) ,                'inq varid end_time', l1_file)
call nc_check( nf90_get_var(fid, varid, obs_date_ed),                  'get var   end_time', l1_file)
map%obs_date = (obs_date_st+obs_date_ed)/2.



! close the file as we are now done reading it
call nc_check( nf90_close(fid) , 'close file', l1_file)


! Now handling missing values
allocate(map%dqf(map%nx,map%ny))
map%dqf(:,:) = 0.
where ( map%bt < 100 ) map%dqf = 4


end subroutine himawari_load_ahi_map











!------------------------------------------------------------------------------
!  extract the ABI channel observations from the map type
!  and convert to DART observation format.  allow caller to specify
!  a bounding box and only extract data within that region.

subroutine make_obs_sequence (seq, map, lon1, lon2, lat1, lat2, &
                              x_thin, y_thin, himawari_num, reject_dqf_1, &
                              obs_err_spec, vloc_pres_hPa)

type(obs_sequence_type),    intent(inout) :: seq
type(himawari_ahi_map_type),    intent(in)    :: map
real(r8),                   intent(in)    :: lon1, lon2, lat1, lat2
integer,                    intent(in)    :: x_thin, y_thin
integer,                    intent(in)    :: himawari_num
logical,                    intent(in)    :: reject_dqf_1
real(r8),                   intent(in)    :: obs_err_spec
real(r8),                   intent(in)    :: vloc_pres_hPa

type(obs_type)          :: obs, prev_obs

integer :: ix, iy
integer :: days, seconds
integer :: obs_num, key
integer :: which_vert 

real(r8) :: olon, olat, vloc
real(r8) :: obs_value, obs_err
real(r8) :: rqc
real(r8) :: latd, lond, beta

real(r8) :: lam1, lam2, phi1, phi2, r_fac

real(digits12) :: rdays, remainder

real(r8) :: sat_az, sat_ze, sun_az, sun_ze, specularity
integer :: platform_id, sat_id, sensor_id

type(time_type) :: obs_time, start_time

integer :: robstype
integer :: himawari_channel

integer :: num_copies, num_qc
! max possible obs from this one map. in practice if the
! real number of processed channels is very much smaller, make
! another parameter so we don't allocate all these unused obs
! (takes time & space) and then delete them at the end.
integer :: max_num

logical :: is_first_obs
type(time_type) :: pre_time

character(len=512) :: obs_type_name

character(len=*), parameter :: routine = 'make_obs_sequence'

if ( .not. module_initialized ) call initialize_module

! one observation data value and one quality control value
! per obs.  if you change these you have to set additional
! metadata for them below.
num_copies  = 1
num_qc      = 1

! Initialize an obs_sequence
max_num = map%nx*map%ny
call init_obs_sequence(seq, num_copies, num_qc, max_num)

! set meta data of obs_seq
call set_copy_meta_data(seq, 1, 'observation')
call set_qc_meta_data(seq, 1, 'Himawari QC')

! Initialize the obs variables
call init_obs(     obs, 1, 1)
call init_obs(prev_obs, 1, 1)

is_first_obs = .true.

! things known (and constant) about the input data and rttov
platform_id = 31   ! Himawari

select case(himawari_num)
    case (8)
        sat_id = himawari_num
        obs_type_name = 'HIMAWARI_8_AHI_RADIANCE'
    case (9)
        sat_id = himawari_num
        obs_type_name = 'HIMAWARI_9_AHI_RADIANCE'
    case default
        write(msgstring,*) 'Unknown Himawari number ',himawari_num,' should be either 8 or 9'
        call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
end select        

sensor_id   = 56  ! AHI

!------------------------------------------------------------------------------
!  loop over all observations within the file

obs_num = 1
if (vloc_pres_hPa < 0.0_r8) then
    which_vert = VERTISUNDEF
else
    which_vert = VERTISPRESSURE
end if

! assign each observation the correct observation type
robstype = get_index_for_type_of_obs(obs_type_name)
if (robstype < 1) then
   write(msgstring,*) 'unknown observation type ',trim(obs_type_name)
   call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
endif



xloop: do ix=1,map%nx
    ! if we're going to subset x, we will cycle here
    if (x_thin > 0) then
        if (modulo(ix, x_thin) /= 0) cycle xloop
    endif

   ! columns are across-track, varying faster than rows.
   yloop:  do iy=1,map%ny

      ! if we're going to subset y, ditto
      if (y_thin > 0) then
         if (modulo(iy, y_thin) /= 0) cycle yloop
      endif

      ! check channel quality control
      rqc = map%DQF(ix, iy)
      ! reject bad scans here, depending on whether conditional DQF=1 is okay
      if (reject_dqf_1) then
        if (rqc /= 0) cycle yloop 
      else
        if (rqc > 1) cycle yloop
      end if

      ! observation lat, lon:
      olat  = map%lat (ix,iy) ! valid range [ -90.00,  90.00]
      olon  = map%lon (ix,iy) ! valid range [-180.00, 180.00]

      ! verify the location is not outside valid limits.  AIRS  uses -180/180
      if((olon > 180.0_r8) .or. (olon < -180.0_r8) .or.  &
         (olat >  90.0_r8) .or. (olat <  -90.0_r8)) then
         write(*,*)'WARNING : invalid location.  x,y,lon,lat = ', ix,iy,olon,olat
         cycle yloop
      endif

      ! reject observations outside the bounding box (allowing wrapping)
      if(( olat < lat1) .or. ( olat > lat2 ) .or. &
         (.not. is_longitude_between(olon, lon1, lon2))) cycle yloop

      ! make sure lon is between 0 and 360
      if (olon < 0.0_r8) olon = olon + 360.0_r8

      ! set the zenith angle (aka earth incidence angle)
      ! see https://svn.ssec.wisc.edu/repos/cloud_team_cr/trunk/viewing_geometry_module.f90
      sat_ze = map%sat_ze(ix,iy) !* deg2rad

      ! calculate the bearing between the obs lat/lon and the SClat/lon
      sat_az = map%sat_az(ix,iy) !* deg2rad

      write(*,*) map%sat_ze(ix,iy), sat_ze, map%sat_az(ix,iy) , sat_az


      ! Convert obs_date to something DART understands
      ! Printout from ncdump of variable relating to obs_date
      ! double start_time(time) ;
      !     start_time:long_name = "observation start time" ;
      !     start_time:units = "days since 1858-11-17 0:0:0" ;
      !     start_time:standard_name = "time" ;
      rdays = floor(map%obs_date)
      remainder = (map%obs_date-rdays)*86400.0_digits12
      start_time = set_date(1858,11,17,0,0,0)
      obs_time = set_time(nint(remainder),nint(rdays)) + start_time
      call get_time(obs_time, seconds, days)

      ! create the radiance obs for this observation, add to sequence
      obs_value = map%bt(ix, iy)
      if (obs_value < 0.0_r8) cycle yloop

      ! TODO: specify this as a function of channel
      obs_err = obs_err_spec

      if (vloc_pres_hPa < 0.0_r8) then
          ! disable pressure location, use VERTISUNDEF, so no single vertical location
          vloc = 0.0_r8
      else
          ! assign the integrated column a vertical location
          ! can correspond to the height of the peak of the weighting function
          vloc = vloc_pres_hPa*100._r8 ! convert from hPa to Pa
      end if

      ! We don't yet have specularity data to add to the observations.
      if (get_rttov_option_logical('do_lambertian')) then
         write(msgstring,*) 'Advanced Himawari Imager IR obs do not yet support specularity or Lambertian'
         call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
      else if (get_rttov_option_logical('addsolar')) then
         write(msgstring,*) 'Advanced Himawari Imager IR obs do not yet support solar calculations'
         call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
      else
         sun_az = MISSING_R8
         sun_ze = MISSING_R8
         specularity = MISSING_R8
      end if

      ! the RTTOV ABI coefficients start from channel 7
      himawari_channel = map%channel !-6

      ! add additional metadata for this obs type.  returns key to use in create call
      call set_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, platform_id, sat_id, sensor_id, & 
         himawari_channel, specularity)

      call create_3d_obs(olat, olon, vloc, which_vert, obs_value, robstype, &
                             obs_err, days, seconds, rqc, obs, key)

      call add_obs_to_seq(seq, obs, obs_time, prev_obs, pre_time, is_first_obs)

      obs_num = obs_num + 1
   enddo yloop
enddo xloop

!! Print a little summary
!call print_obs_seq(seq, '')

write(msgstring,*) 'Finished loading ',obs_num-1,' of ',key,'total Advanced Himawari Imager IR observations for map ' // &
   trim(map%filename)
call error_handler(E_MSG, routine, msgstring, source, revision, revdate)

end subroutine make_obs_sequence

end module
