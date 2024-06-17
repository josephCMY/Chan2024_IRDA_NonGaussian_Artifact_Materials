! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module meteosat_seviri_L1b_mod

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
    
    public :: meteosat_seviri_map_type, meteosat_load_seviri_map, make_obs_sequence
    
    !>@todo FIXME ... we should be using a lot more of the netcdf_utilities interfaces.
    !>                they greatly simplify the code.
    
    type meteosat_seviri_map_type
       character(:), allocatable :: filename
       integer :: channel
       integer :: nx
       integer :: ny
    
    
       real(r8),    allocatable :: bt(:,:)  ! Brightness Temperature values
       integer(1),  allocatable :: dqf(:,:) ! data quality flag. 0 = good, 1 is conditionally okay, 2-4 is bad.
       real(r4)                 :: r_eq  ! radius of earth at equator
       real(r4)                 :: r_pol ! radius of earth at poles
       real(r4)                 :: H ! the total height of satellite 
       type(time_type)          :: scan_date  ! Mid-scan date
       real(r4)                 :: sc_lat ! space-craft latitude
       real(r4)                 :: sc_lon ! space-craft longitude
       real(r8),    allocatable :: lat(:,:) 
       real(r8),    allocatable :: lon(:,:)
    
    end type
    
    character(len=512) :: msgstring
    
    integer :: fid, varid, dimid
    
    ! version controlled file description for error handling, do not edit
    character(len=*), parameter :: source   = 'meteosat_seviri_L1b_mod.f90'
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
    
    subroutine meteosat_load_seviri_map(l1b_file, map, ir_channel_name, meteosat_num)
    
        character(len=*),    intent(in)    :: l1b_file
        type(meteosat_seviri_map_type), intent(inout) :: map
        character(len=*),    intent(in)    :: ir_channel_name
        integer,             intent(in)    :: meteosat_num
    
        logical,          parameter :: little_endian = .true.
    
        character(len=*), parameter :: routine = 'meteosat_load_seviri_map'
    
        character(len=512) :: string1
        
        character(len=19) :: scan_end_date_string
        integer :: scan_year, scan_month, scan_day, scan_hour, scan_minutes, scan_seconds
        
    
        integer :: i, j
        real(digits12) :: a, b, c, r_fac, xv, yv, r_s
        real(digits12) :: sx, sy, sz
        real(digits12) :: time_r8
    
        real(digits12), parameter :: r2d = 180.0_digits12/(atan(1.0_digits12)*4.0_digits12)
    
        if ( .not. module_initialized ) call initialize_module
    
        write(string1,*) 'Now loading L1B data from file ',trim(l1b_file)
        call error_handler(E_MSG, routine, string1, source, revision, revdate)
    
        allocate(character(len=len(l1b_file)) :: map%filename)
        map%filename = l1b_file


        ! Set useful constants (in meters)
        map%r_eq  = 6378137
        map%r_pol = 6356752
    
    
        ! Open file
        ! ---------
        call nc_check( nf90_open(l1b_file, nf90_nowrite, fid), 'file open', l1b_file)
    
    
    
        ! Read spatial dimensions
        ! -----------------------
        ! Load observation longitude dimension
        call nc_check( nf90_inq_dimid(fid, "lon", dimid),        'inq dimid longitude', l1b_file)
        call nc_check( nf90_inquire_dimension(fid, dimid, len=map%nx), 'inq dim   longitude', l1b_file)
    
        ! Load observation latitude dimension
        call nc_check( nf90_inq_dimid(fid, "lat", dimid),         'inq dimid latitude', l1b_file)
        call nc_check( nf90_inquire_dimension(fid, dimid, len=map%ny), 'inq dim   latitude', l1b_file)
    
        ! Load observation longitude data array
        allocate( map%lon( map%nx, map%ny ) )
        call nc_check( nf90_inq_varid(fid, 'longitude', varid) ,       'inq varid longitude', l1b_file)
        call nc_check( nf90_get_var(fid, varid, map%lon),              'get var   longitude', l1b_file)
    
        ! Load observation latitude data array
        allocate( map%lat( map%nx, map%ny ) )
        call nc_check( nf90_inq_varid(fid, 'latitude', varid) ,        'inq varid latitude', l1b_file)
        call nc_check( nf90_get_var(fid, varid, map%lat),              'get var   latitude', l1b_file)
    
    
    
        ! Read BT values
        ! ---------------
        allocate( map%bt( map%nx, map%ny ) )
        call nc_check( nf90_inq_varid(fid, ir_channel_name, varid) ,   'inq varid BT', l1b_file)
        call nc_check( nf90_get_var(fid, varid, map%bt),              'get var   BT', l1b_file)
    
    
    
        ! QC values for METEOSAT
        ! ----------------------
        allocate( map%dqf( map%nx, map%ny ))
        map%dqf(:,:) = 0
    
        ! Reject bad values of BT
        where ( map%bt > 500. ) map%dqf = 4
        where ( map%bt < 100. ) map%dqf = 4
    
        ! Reject bad locations
        where ( map%lat >   90 ) map%dqf = 4
        where ( map%lat <  -90 ) map%dqf = 4
        where ( map%lon >  180 ) map%dqf = 4
        where ( map%lon < -180 ) map%dqf = 4
    
    
    
        ! Determine date of scan
        ! -----------------------
        call nc_check( nf90_get_att(fid, NF90_GLOBAL, 'sensor_end_time'  , scan_end_date_string  ), 'inq scan end date'  , l1b_file )
        scan_end_date_string = replace_character_in_string( scan_end_date_string, (/'-','_',':'/), ' ' )
        read( scan_end_date_string , '(I4,X,I2,X,I2,X,I2,X,I2,X,I2)' ) scan_year, scan_month, scan_day, scan_hour, scan_minutes, scan_seconds
        map%scan_date = set_date( scan_year, scan_month, scan_day, scan_hour, scan_minutes, scan_seconds )
    
    
        ! close the file as we are now done reading it
        call nc_check( nf90_close(fid) , 'close file', l1b_file)


        ! Indicate channel number 
        ! -----------------------
        select case( ir_channel_name )
            case ("WV_062")
                map%channel = 5

            case default 
                write(msgstring,*) 'Unknown SEVIRI channel number ',ir_channel_name
                call error_handler(E_ERR,routine,msgstring,source,revision,revdate)

        end select


        ! Some satellite specifics
        select case(meteosat_num)
        case (9)
            map%sc_lat = 0.0
            map%sc_lon = 45.5
        case default
            write(msgstring,*) 'Unknown meteosat number ',meteosat_num
            call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
        end select   


    
    end subroutine meteosat_load_seviri_map
    
    
    
    
    
    
    
    
    
    
    !------------------------------------------------------------------------------
    !  extract the seviri channel observations from the map type
    !  and convert to DART observation format.  allow caller to specify
    !  a bounding box and only extract data within that region.
    
    subroutine make_obs_sequence (seq, map, lon1, lon2, lat1, lat2, &
                                  x_thin, y_thin, meteosat_num, reject_dqf_1, &
                                  obs_err_spec, vloc_pres_hPa)
    
        type(obs_sequence_type),    intent(inout) :: seq
        type(meteosat_seviri_map_type),    intent(in)    :: map
        real(r8),                   intent(in)    :: lon1, lon2, lat1, lat2
        integer,                    intent(in)    :: x_thin, y_thin
        integer,                    intent(in)    :: meteosat_num
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
        integer :: meteosat_channel
    
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
        call set_qc_meta_data(seq, 1, 'meteosat QC')
    
        ! Initialize the obs variables
        call init_obs(     obs, 1, 1)
        call init_obs(prev_obs, 1, 1)
    
        is_first_obs = .true.
    
        ! things known (and constant) about the input data and rttov   
        select case(meteosat_num)
            case (9)
                obs_type_name = 'MSG_2_SEVIRI_RADIANCE'
                platform_id = 12 ! RTTOV's platform ID for Meteosat Second Gen (MSG)
                sat_id = 2       ! Meteosat-9 corresponds to RTTOV's sat_id 2

            case default
                write(msgstring,*) 'Unknown meteosat number ',meteosat_num
                call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
        end select        
    
        sensor_id   = 21  ! seviri
    
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
    
        r_fac = (map%r_eq**2)/(map%r_pol**2)
    
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
    
            latd = (map%lat(ix,iy) - map%sc_lat)*deg2rad
            lond = (map%lon(ix,iy) - map%sc_lon)*deg2rad
            beta = acos(cos(latd)*cos(lond))
    
            sat_ze = map%H * sin(beta) / sqrt( r_fac + map%H**2 - 2.0_r8*map%r_eq*map%H*cos(beta) )
            sat_ze = max(-1.0_r8,min(1.0_r8, sat_ze))
            sat_ze = asin(sat_ze) / deg2rad
    
            lam1 = deg2rad*map%lon(ix,iy)*deg2rad
            lam2 = deg2rad*map%sc_lon*deg2rad
            phi1 = deg2rad*map%lat(ix,iy)*deg2rad
            phi2 = deg2rad*map%sc_lat*deg2rad
    
            ! calculate the bearing between the obs lat/lon and the SClat/lon
            sat_az = atan2(sin(lam2-lam1)*cos(phi2),&
                            cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(lam2-lam1))/deg2rad
    

            ! Converting scan time to seconds and days
            call get_time(map%scan_date, seconds, days)
    
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
                write(msgstring,*) 'meteosat observations do not yet support specularity or Lambertian'
                call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
            else if (get_rttov_option_logical('addsolar')) then
                write(msgstring,*) 'meteosat observations do not yet support solar calculations'
                call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
            else
                sun_az = MISSING_R8
                sun_ze = MISSING_R8
                specularity = MISSING_R8
            end if
    
            ! the RTTOV seviri channel
            meteosat_channel = map%channel
    
            ! add additional metadata for this obs type.  returns key to use in create call
            call set_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, platform_id, sat_id, sensor_id, & 
                meteosat_channel, specularity)
    
            call create_3d_obs(olat, olon, vloc, which_vert, obs_value, robstype, &
                                    obs_err, days, seconds, rqc, obs, key)
    
            call add_obs_to_seq(seq, obs, obs_time, prev_obs, pre_time, is_first_obs)
    
            obs_num = obs_num + 1
        enddo yloop
        enddo xloop
    
        !! Print a little summary
        !call print_obs_seq(seq, '')
    
        write(msgstring,*) 'Finished loading ',obs_num-1,' of ',key,'total meteosat seviri observations for map ' // &
        trim(map%filename)
        call error_handler(E_MSG, routine, msgstring, source, revision, revdate)
    
    end subroutine make_obs_sequence
    
    
    
    
    
    
    ! Simple character substitution function
    ! Copied from https://fortran-lang.github.io/fpm/proc/replace.html
    function replace_character_in_string(string, charset, target_char) result(res)
        character(*), intent(in) :: string
        character, intent(in) :: charset(:), target_char
        character(len(string)) :: res
        integer :: n
        res = string
        do n = 1, len(string)
            if (any(string(n:n) == charset)) then
                res(n:n) = target_char
            end if
        end do
    end function replace_character_in_string
    
    end module