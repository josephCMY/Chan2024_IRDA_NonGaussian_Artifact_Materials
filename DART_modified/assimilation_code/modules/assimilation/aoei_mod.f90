! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


! This module implements the Adaptive Observation Error Inflation (AOEI) 
! scheme of Minamide and Zhang (2017, Monthly Weather Review). 
!
! Implemented by Man-Yau (Joseph) Chan


! IMPORTANT NOTE:
! ---------------
! This AOEI implementation differs from the original Minamide-Zhang implementation
! in the Pennsylvania State University EnKF (PSU-EnKF) system.
! The PSU-EnKF system calls AOEI within the sequential assimilation loop. In other 
! words, changing the order of the observations can change the behavior of AOEI.
! To avoid this complication, DART's implementation of the AOEI happens outside
! of the sequential assimilation loop.



module adaptive_obs_error_inflation


    ! Section to use stuff from other DART modules
    ! --------------------------------------------

    use      types_mod,       only : r8, i8

    use         location_mod, only : location_type
    
    use  utilities_mod,       only : error_handler, E_MSG, find_namelist_in_file, check_namelist_read
    
    use obs_sequence_mod,     only : obs_sequence_type, obs_type, get_num_copies, get_num_qc, &
                                     init_obs, get_obs_from_key, get_obs_def, get_obs_values, &
                                     set_obs_def, set_obs
    
    use          obs_def_mod, only : obs_def_type, get_obs_def_error_variance, get_obs_def_type_of_obs,     &
                                     set_obs_def_error_variance
    
    use         obs_kind_mod, only : get_quantity_for_type_of_obs, QTY_BRIGHTNESS_TEMPERATURE
    
    use ensemble_manager_mod, only : ensemble_type, get_my_num_vars, get_my_vars,             &
                                     get_var_owner_index, map_pe_to_task
    
    use mpi_utilities_mod,    only : broadcast_send, broadcast_recv 

    use assim_model_mod,      only : get_state_meta_data
    
    



    implicit none
    private





    public :: apply_aoei, init_aoei, use_aoei




    ! Module-level variable definitions
    ! ---------------------------------
    character(len=*), parameter :: source = 'aoei_mod.f90'

    ! VERY IMPORTANT: QTYs to perform the AOEI for
    integer, parameter :: num_aoei_qty = 256
    integer, allocatable :: aoei_obs_qty_list(:)

    ! Parameter regulating strength of AOEI
    integer, parameter :: aoei_factor = 9  ! ensures d/sigma < 3

    ! Some convenient global storage items
    character(len=512)      :: msgstring


    ! Namelist quantity to turn on AOEI
    logical ::  use_aoei


    ! Namelist quantities
    namelist / aoei_nml / use_aoei



    ! Time to define subroutines!
    ! ---------------------------

    contains



    subroutine init_aoei()

        integer                 :: iunit, io

        ! Read in namelist
        call find_namelist_in_file("input.nml", "aoei_nml", iunit)
        read(iunit, nml = aoei_nml, iostat = io)
        call check_namelist_read(iunit, io, "aoei_nml")

        ! Init list of AOEI quantities
        allocate( aoei_obs_qty_list(num_aoei_qty) )
        aoei_obs_qty_list = -999

        ! Put in quantities that AOEI will be applied on.
        aoei_obs_qty_list(1) = QTY_BRIGHTNESS_TEMPERATURE
    
    end subroutine init_aoei

        




    subroutine apply_aoei( ens_handle, obs_ens_handle, obs_seq, keys, ens_size, obs_val_index )

        ! Define input variables
        type( ensemble_type ), intent(in)               :: ens_handle
        type( ensemble_type ), intent(in)               :: obs_ens_handle
        type( obs_sequence_type ), intent(inout)        :: obs_seq
        integer, intent(in)                             :: keys(:)
        integer, intent(in)                             :: ens_size
        integer, intent(in)                             :: obs_val_index


        integer             :: my_num_obs, i, j, owner, owners_index, my_num_state
        integer             :: base_obs_kind, base_obs_type, nth_obs

        integer(i8), allocatable    :: my_obs_indx(:)

        type(obs_type)       :: observation
        type(obs_def_type)   :: obs_def
        type(location_type)  :: dummyloc

        real(r8) :: obs(1), obs_err_var, obs_prior( ens_size )
        real(r8) :: aoei_yf_avg, aoei_innovation_squared, aoei_yf_var, aoei_obs_err

        logical :: skip_obs



        ! Allocation statements
        allocate( my_obs_indx(        obs_ens_handle%my_num_vars) )




        ! Get info on my number and indices for obs
        my_num_obs = get_my_num_vars(obs_ens_handle)
        call get_my_vars(obs_ens_handle, my_obs_indx)
        
        ! Construct an observation temporary
        call init_obs(observation, get_num_copies(obs_seq), get_num_qc(obs_seq))




        ! Loop through all the observations
        AOEI_LOOP_OVER_ALL_OBS : do i = 1, obs_ens_handle%num_vars


            ! Every pe has information about the global obs sequence
            call get_obs_from_key(obs_seq, keys(i), observation)
            call get_obs_def(observation, obs_def)
            obs_err_var = get_obs_def_error_variance(obs_def)
            base_obs_type = get_obs_def_type_of_obs(obs_def)
            if (base_obs_type > 0) then
                base_obs_kind = get_quantity_for_type_of_obs(base_obs_type)
            else
                call get_state_meta_data(-1 * int(base_obs_type,i8), dummyloc, base_obs_kind)  ! identity obs
            endif


            ! Check if current obs qty is targetted by AOEI
            skip_obs = .true. 
            CHECK_IF_AOEI_IS_USED : do j = 1, num_aoei_qty

                if ( base_obs_kind == aoei_obs_qty_list(j) ) skip_obs = .false.
            end do CHECK_IF_AOEI_IS_USED


            ! Skip current obs if it is not targetted by AOEI
            if ( skip_obs ) cycle AOEI_LOOP_OVER_ALL_OBS


            ! Get the value of the observation
            call get_obs_values(observation, obs, obs_val_index)
        

            ! Find out who has this observation and where it is
            call get_var_owner_index( ens_handle, int(i,i8), owner, owners_index )



            ! Prior observation ensemble is distributed across all processes. 
            ! Following block is done only by the owner of this observation
            !-----------------------------------------------------------------------
            if(ens_handle%my_pe == owner) then
            ! each task has its own subset of all obs.  if they were converted in the
            ! vertical up above, then we need to broadcast the new values to all the other

                ! Load prior obs values
                obs_prior = obs_ens_handle%copies(1:ens_size, owners_index)

                ! Broadcast prior obs info to all processes
                call broadcast_send(map_pe_to_task(ens_handle, owner), obs_prior )

            ! Next block is done by processes that do NOT own this observation
            !-----------------------------------------------------------------------
            else

                call broadcast_recv(map_pe_to_task(ens_handle, owner), obs_prior )

            endif
            !-----------------------------------------------------------------------



            ! Apply Adaptive Obs Err Inflation
            ! --------------------------------

            ! Compute innovation squared (OmB^2)
            aoei_yf_avg = sum(obs_prior(1:ens_size)) / ens_size
            aoei_innovation_squared = (obs(1) - aoei_yf_avg)**2

            ! Compute prior obs spread
            aoei_yf_var = sum( (obs_prior(1:ens_size)-aoei_yf_avg)**2 ) / ( ens_size - 1 )

            ! AOEI formula [Minamide and Zhang (2017; MWR), Eq 4]
            aoei_obs_err = max( obs_err_var, (aoei_innovation_squared - aoei_factor*aoei_yf_var)/aoei_factor )


            ! Overwrite obs_err_var if needed
            if ( aoei_obs_err > obs_err_var * 1.0001 ) then
            
                ! Print message to tell user that AOEI is evoked
                write(msgstring, '(2X, A20, I8)' ) 'AOEI called on obs ',  i
                call error_handler(E_MSG,'aoei_mod',msgstring)
                write(msgstring, '(4X, A25, F)') 'obs value = ',            obs(1)
                call error_handler(E_MSG,'aoei_mod',msgstring)
                write(msgstring, '(4X, A25, F)') 'avg obs prior = ',        aoei_yf_avg
                call error_handler(E_MSG,'aoei_mod',msgstring)
                write(msgstring, '(4X, A25, F)') 'obs_prior_var = ',        aoei_yf_var
                call error_handler(E_MSG,'aoei_mod',msgstring)
                write(msgstring, '(4X, A25, F)') 'original obs_err_var = ', obs_err_var
                call error_handler(E_MSG,'aoei_mod',msgstring)
                write(msgstring, '(4X, A25, F)') 'inflated obs_err_var = ', aoei_obs_err
                call error_handler(E_MSG,'aoei_mod',msgstring)
            
                ! Overwrite obs error
                call set_obs_def_error_variance( obs_def, aoei_obs_err )
                call set_obs_def( observation, obs_def )
                call set_obs( obs_seq, observation, keys(i))
            
            endif


        end do AOEI_LOOP_OVER_ALL_OBS


    

    end subroutine apply_aoei





end module adaptive_obs_error_inflation
