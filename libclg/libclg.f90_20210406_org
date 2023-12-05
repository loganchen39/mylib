! Description: useful self-defined functions and subroutines
! Author     : Ligang Chen, lchen2@umd.edu
! Modified   : 04/10/2012
! Created    : 04/10/2012


subroutine check(status)
  use netcdf
  integer, intent(in) :: status
  
  if (status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop 2
  end if
end subroutine check


! Return 1 for a leap year, and 0 for not a leap year
integer function isleapyear(iyr)
    logical y4 , y100 , y400

    y4   = mod(iyr, 4  ) .eq. 0
    y100 = mod(iyr, 100) .eq. 0
    y400 = mod(iyr, 400) .eq. 0

    if ((y4 .and. .not.y100) .or. y400) then
       isleapyear = 1
    else
       isleapyear = 0
    endif

    return
end function isleapyear


function days_in_month(year, month)
    implicit none
    integer, external :: isleapyear

    integer, intent(in) :: year, month
    integer             :: days_in_month

    integer :: dofm(12) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    if (isleapyear(year) /= 0) then 
        dofm(2) = 29
    else
        dofm(2) = 28
    end if

    days_in_month = dofm(month)

    return
end function days_in_month


function day_of_year(year, month, day)
    implicit none
    integer, external :: isleapyear

    integer, intent(in) :: year, month, day
    integer             :: day_of_year

    integer :: dofm(13) = (/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365/)


    day_of_year = -1

    if (year < 1000 .or. year > 9999 .or. month < 1 .or. month > 12 .or. day < 1 .or. day > 31 ) then
        write(*, *) "ERROR: in day_of_year, invalid parameters!"
        stop
    end if

    if (month <= 2) then
        day_of_year = dofm(month) + day
    else
        if (isleapyear(year) /= 0) then
            day_of_year = dofm(month) + 1 + day
        else
            day_of_year = dofm(month) + day
        end if
    end if

    return
end function day_of_year


function get_record_index(origin_year, origin_month, origin_day,  &
                          sample_year, sample_month, sample_day)
    implicit none
    integer, external :: isleapyear, day_of_year

    integer, intent(in) ::  &
        origin_year ,  &
        origin_month,  &
        origin_day  ,  &
        sample_year ,  &
        sample_month,  &
        sample_day

    integer :: get_record_index

    integer :: dofy_origin, dofy_sample
    integer :: i_yr, i_mon, i_day


    get_record_index = 0

    dofy_origin = day_of_year(origin_year, origin_month, origin_day)
    dofy_sample = day_of_year(sample_year, sample_month, sample_day)

    if (origin_year == sample_year) then
        get_record_index = dofy_sample - dofy_origin + 1
    else
        if (isleapyear(origin_year) /= 0) then
            get_record_index = get_record_index + 366 - dofy_origin + 1
        else
            get_record_index = get_record_index + 365 - dofy_origin + 1
        end if

        get_record_index = get_record_index + dofy_sample

        do i_yr = origin_year + 1, sample_year - 1
            if (isleapyear(i_yr) /= 0) then
                get_record_index = get_record_index + 366
            else
                get_record_index = get_record_index + 365
            end if
        end do
    end if

    return
end function get_record_index


function get_record_index_no_leap_year(origin_year, origin_month, origin_day,  &
                                       sample_year, sample_month, sample_day)
    implicit none
    integer, external :: isleapyear, day_of_year

    integer, intent(in) ::  &
        origin_year ,  &
        origin_month,  &
        origin_day  ,  &
        sample_year ,  &
        sample_month,  &
        sample_day

    integer :: get_record_index_no_leap_year

    integer :: dofy_origin, dofy_sample
    integer :: i_yr, i_mon, i_day


    get_record_index_no_leap_year = 0

    dofy_origin = day_of_year(origin_year, origin_month, origin_day)
    if (isleapyear(origin_year) .eq. 1 .and. origin_month .ge. 3) then
        dofy_origin = dofy_origin - 1
    end if

    dofy_sample = day_of_year(sample_year, sample_month, sample_day)
    if (isleapyear(sample_year) .eq. 1 .and. sample_month .ge. 3) then
        dofy_sample = dofy_sample - 1
    end if

    if (origin_year == sample_year) then
        get_record_index_no_leap_year = dofy_sample - dofy_origin + 1
    else
        if (isleapyear(origin_year) /= 0) then
            get_record_index_no_leap_year = get_record_index_no_leap_year + 365 - dofy_origin + 1
        else
            get_record_index_no_leap_year = get_record_index_no_leap_year + 365 - dofy_origin + 1
        end if

        get_record_index_no_leap_year = get_record_index_no_leap_year + dofy_sample

        do i_yr = origin_year + 1, sample_year - 1
            if (isleapyear(i_yr) /= 0) then
                get_record_index_no_leap_year = get_record_index_no_leap_year + 365
            else
                get_record_index_no_leap_year = get_record_index_no_leap_year + 365
            end if
        end do
    end if

    return
end function get_record_index_no_leap_year


subroutine pqsort(arr, p_arr, n_elem)
    implicit none 

    integer,                    intent(in   ) :: n_elem
    real   , dimension(n_elem), intent(inout) :: arr
    integer, dimension(n_elem), intent(out  ) :: p_arr  ! permutation array

    ! local vars
    integer :: i_elem, j_elem, idx_min, idx_tmp
    real    :: r_min, r_tmp

    do i_elem = 1, n_elem
        p_arr(i_elem) = i_elem
    end do

    ! sort and exchange
    do i_elem = 1, n_elem - 1
        r_min   = arr(i_elem)
        idx_min = i_elem

        do j_elem = i_elem + 1, n_elem
            if (arr(j_elem) .lt. r_min) then 
                r_min   = arr(j_elem)
                idx_min = j_elem
            end if
        end do
     
        if (idx_min .ne. i_elem) then 
            r_tmp        = arr(i_elem)
            arr(i_elem ) = r_min
            arr(idx_min) = r_tmp

            idx_tmp        = p_arr(i_elem)
            p_arr(i_elem ) = p_arr(idx_min)  ! lgchen: not just idx_min
            p_arr(idx_min) = idx_tmp
        end if
    end do

    return
end subroutine pqsort


subroutine qsort(arr, n_elem)
    implicit none 

    integer,                    intent(in ) :: n_elem
    real   , dimension(n_elem), intent(inout) :: arr

    ! local vars
    integer :: i_elem, j_elem, idx_min
    real    :: r_min, r_tmp

    ! sort and exchange
    do i_elem = 1, n_elem - 1
        r_min   = arr(i_elem)
        idx_min = i_elem

        do j_elem = i_elem + 1, n_elem
            if (arr(j_elem) .lt. r_min) then 
                r_min   = arr(j_elem)
                idx_min = j_elem
            end if
        end do
     
        if (idx_min .ne. i_elem) then 
            r_tmp        = arr(i_elem)
            arr(i_elem ) = r_min
            arr(idx_min) = r_tmp
        end if
    end do

    return
end subroutine qsort


! for NCL use, it should not be defined as a function or its returned value
! res_idx will be correct in fortran context like 789, but incorrect in NCL
! context and became 0.00288459. It should be a subroutine (see NCL manual). 
subroutine search_index_binary(sorted_arr, n_elem, key, res_idx)
    implicit none

    integer,                    intent(in) :: n_elem
    real   , dimension(n_elem), intent(in) :: sorted_arr
    real   ,                    intent(in) :: key
    integer,                    intent(out) :: res_idx


    ! local vars
    integer :: low, high, mid


    res_idx = -1

    low  = 1
    high = n_elem - 1

    do while (low <= high)
        mid = floor((low+high)/2.0)

        if (sorted_arr(mid) <= key .and. key < sorted_arr(mid+1)) then
            res_idx = mid
            return
        else if (key < sorted_arr(mid)) then
            high = mid - 1
        else
            low = mid + 1
        end if
    end do

    return
end subroutine search_index_binary


! pr_obs(NX, NY); pr_mod(NX, NY); N_ETS: 30 for seasonal mean
! return ets(N_ETS)
! if all of buffer_zone, ocean and us_land_only are false, then it's based on all_land
subroutine compute_pr_ets_us_cwrf(pr_obs, pr_mod, N_ETS  &
    , lcc, us_landmask                                   &
    , buffer_zone, ocean, us_land_only                   &
    , NX, NY, ets)
    implicit none

    integer,                    intent(in ) :: N_ETS, NX, NY
    real   , dimension(NX, NY), intent(in ) :: pr_obs, pr_mod, lcc, us_landmask
    logical,                    intent(in ) :: buffer_zone, ocean, us_land_only
    real   , dimension(N_ETS) , intent(out) :: ets


    ! local
    real, dimension(N_ETS) :: n_o, n_f, n_fo, n_nfo, n_fno, e_n_fo
    integer                :: n_all, i_y, i_x, i_ets


    n_o = 0.0
    n_all = 0

    n_f    = 0.0
    n_fo   = 0.0
    n_fno  = 0.0
    n_nfo  = 0.0
    e_n_fo = 0.0

    do i_y = 1, NY
        do i_x = 1, NX
            if (.not. buffer_zone) then
                if (i_y <= 14 .or. i_x <= 14 .or. NY - 13 <= i_y .or. NX - 13 <= i_x) then
                    cycle
                end if
            end if

            if (.not. ocean) then
                if (lcc(i_x, i_y) == 16.0) then
                    cycle
                end if
            end if

            if (us_land_only) then
                if (us_landmask(i_x, i_y) == 0.0) then
                    cycle
                end if
            end if

            if (pr_obs(i_x, i_y) > 1.0e5 .or. pr_mod(i_x, i_y) > 1.0e5) then
                cycle
            end if


            n_all = n_all + 1

            do i_ets = 1, N_ETS
                if (pr_obs(i_x, i_y) > 0.5*(i_ets)) then
                    n_o(i_ets) = n_o(i_ets) + 1
                end if

                if (pr_mod(i_x, i_y) > 0.5*(i_ets)) then
                    n_f(i_ets) = n_f(i_ets) + 1
                end if

                if (pr_obs(i_x, i_y) > 0.5*(i_ets) .and. pr_mod(i_x, i_y) > 0.5*(i_ets)) then
                    n_fo(i_ets) = n_fo(i_ets) + 1
                end if

                if (pr_obs(i_x, i_y) <= 0.5*(i_ets) .and. pr_mod(i_x, i_y) > 0.5*(i_ets)) then
                    n_fno(i_ets) = n_fno(i_ets) + 1
                end if

                if (pr_obs(i_x, i_y) > 0.5*(i_ets) .and. pr_mod(i_x, i_y) <= 0.5*(i_ets)) then
                    n_nfo(i_ets) = n_nfo(i_ets) + 1
                end if
            end do  ! end of "do i_ets = 1, N_ETS"

        end do  ! end of "do i_x = 1, NX"
    end do      ! end of "do i_y = 1, NY"

    ! write(*, *) "n_all = ", n_all, ", n_o = ", n_o, ", n_f = ", n_f, ", n_fo = ", n_fo, ", n_fno = ", n_fno, ", n_nfo = ", n_nfo

    do i_ets = 1, N_ETS
        e_n_fo(i_ets) = n_f(i_ets) * n_o(i_ets) / n_all

        if (n_fno(i_ets) + n_nfo(i_ets) + n_fo(i_ets) - e_n_fo(i_ets) > 0.0) then
            ets(i_ets) = (n_fo(i_ets) - e_n_fo(i_ets)) / (n_fno(i_ets) + n_nfo(i_ets) + n_fo(i_ets) - e_n_fo(i_ets))
        else
            ets(i_ets) = 0.0
        end if
    end do  ! end of "do i_ets = 1, N_ETS"

    return
end  subroutine compute_pr_ets_us_cwrf


! subroutine partition is for the quick sort subroutine below
subroutine Partition(A, arr_len, marker)
    ! real(kind = 4), intent(inout), dimension(:) :: A
    integer, intent(in   ) :: arr_len
    real   , intent(inout) :: A(arr_len)
    integer, intent(out  ) :: marker

    integer :: i, j
    real    :: temp
    real    :: x      ! pivot point

    x = A(1)
    i = 0
    j = arr_len + 1

    do
        j = j-1
        do
            if (A(j) <= x) exit
            j = j-1
        end do
        i = i+1
        do
            if (A(i) >= x) exit
            i = i+1
        end do
        if (i < j) then
            ! exchange A(i) and A(j)
            temp = A(i)
            A(i) = A(j)
            A(j) = temp
        elseif (i == j) then
            marker = i+1
            return
        else
            marker = i
            return
        endif
    end do
end subroutine Partition


! recursive quick sort algorithm
recursive subroutine QsortC(A, arr_len)
    implicit none

    integer, intent(in) :: arr_len
    real, intent(inout) :: A(arr_len)

    integer             :: iq

    if(size(A) > 1) then
        call Partition(A, arr_len, iq)
        call QsortC(A(1:iq-1), iq-1)
        call QsortC(A(iq:size(A)), size(A)-iq+1)
    endif
end subroutine QsortC


subroutine write_out_5dayAvg2monthlyAvg_SODA_MOM(curr_yr, curr_mon, avg_temp, avg_salt, avg_u  &
    , avg_v, avg_sea_level, avg_wt, avg_tau_x, avg_tau_y, XT_OCEAN, YT_OCEAN, ST_OCEAN)
    use netcdf
    implicit none

    integer,                                             intent(in) :: curr_yr, curr_mon, XT_OCEAN, YT_OCEAN, ST_OCEAN
    real   , dimension(XT_OCEAN, YT_OCEAN, ST_OCEAN, 1), intent(in) :: avg_temp, avg_salt, avg_u, avg_v, avg_wt
    real   , dimension(XT_OCEAN, YT_OCEAN,           1), intent(in) :: avg_sea_level, avg_tau_x, avg_tau_y

    ! local vars
    character(len = *), parameter :: DIR_ROOT_SODA_ANALYSIS = "\/glade/scratch/lgchen/project/MOM_run/SODA_Reanalysis_2011_2013/test3/MOM_SODA"  &
        , FN_NC_DIM = trim(DIR_ROOT_SODA_ANALYSIS) // "/20110101/history/20110101.ocean.nc"

    character(len=4) :: str_yr
    character(len=2) :: str_mon

    integer :: ncid_dim, ncid_monthly_avg, varid_temp, varid_salt, varid_u, varid_v, varid_sea_level, varid_wt, varid_tau_x  &
        , varid_tau_y, varid_time, varid_xt_ocean, varid_yt_ocean, varid_st_ocean, varid_st_edges_ocean, varid_xu_ocean, varid_yu_ocean  &
        , varid_sw_ocean, varid_sw_edges_ocean, varid_geolon_t, varid_geolat_t
    integer :: dimid_time, dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, dimid_st_edges_ocean, dimid_xu_ocean, dimid_yu_ocean  &
        , dimid_sw_ocean, dimid_sw_edges_ocean

    ! vars to receive dim data
    real(kind = 8) :: var_xt_ocean(XT_OCEAN), var_yt_ocean(YT_OCEAN), var_st_ocean(ST_OCEAN)  &
        , var_st_edges_ocean(51), var_xu_ocean(XT_OCEAN), var_yu_ocean(YT_OCEAN), var_sw_ocean(ST_OCEAN)  &
        , var_sw_edges_ocean(51), var_geolon_t(XT_OCEAN, YT_OCEAN), var_geolat_t(XT_OCEAN, YT_OCEAN)


    write(unit=str_yr , fmt="(I4.4)") curr_yr
    write(unit=str_mon, fmt="(I2.2)") curr_mon

    ! open and read the dim file data
    call check(nf90_open(trim(FN_NC_DIM), nf90_nowrite, ncid_dim))
    
    call check(nf90_inq_varid(ncid_dim, "xt_ocean"      , varid_xt_ocean      ))
    call check(nf90_inq_varid(ncid_dim, "yt_ocean"      , varid_yt_ocean      ))
    call check(nf90_inq_varid(ncid_dim, "st_ocean"      , varid_st_ocean      ))
    call check(nf90_inq_varid(ncid_dim, "st_edges_ocean", varid_st_edges_ocean))
    call check(nf90_inq_varid(ncid_dim, "xu_ocean"      , varid_xu_ocean      ))
    call check(nf90_inq_varid(ncid_dim, "yu_ocean"      , varid_yu_ocean      ))
    call check(nf90_inq_varid(ncid_dim, "sw_ocean"      , varid_sw_ocean      ))
    call check(nf90_inq_varid(ncid_dim, "sw_edges_ocean", varid_sw_edges_ocean))
    call check(nf90_inq_varid(ncid_dim, "geolon_t"      , varid_geolon_t      ))
    call check(nf90_inq_varid(ncid_dim, "geolat_t"      , varid_geolat_t      ))

    call check(nf90_get_var(ncid_dim, varid_xt_ocean      , var_xt_ocean      ))
    call check(nf90_get_var(ncid_dim, varid_yt_ocean      , var_yt_ocean      ))
    call check(nf90_get_var(ncid_dim, varid_st_ocean      , var_st_ocean      ))
    call check(nf90_get_var(ncid_dim, varid_st_edges_ocean, var_st_edges_ocean))
    call check(nf90_get_var(ncid_dim, varid_xu_ocean      , var_xu_ocean      ))
    call check(nf90_get_var(ncid_dim, varid_yu_ocean      , var_yu_ocean      ))
    call check(nf90_get_var(ncid_dim, varid_sw_ocean      , var_sw_ocean      ))
    call check(nf90_get_var(ncid_dim, varid_sw_edges_ocean, var_sw_edges_ocean))
    call check(nf90_get_var(ncid_dim, varid_geolon_t      , var_geolon_t      ))
    call check(nf90_get_var(ncid_dim, varid_geolat_t      , var_geolat_t      ))

!   call check(nf90_close(ncid_dim))


    ! write out the SODA-MOM monthly average result to a netCDF file
    ! first define related dimensions, variables and their associated attributes
    call check(nf90_create(DIR_ROOT_SODA_ANALYSIS // "/monthly_avg_ocean/" // str_yr // str_mon // ".ocean.nc", NF90_CLOBBER, ncid_monthly_avg))

    call check(nf90_def_dim(ncid_monthly_avg, "time"          , 1       , dimid_time          ))
    call check(nf90_def_dim(ncid_monthly_avg, "xt_ocean"      , XT_OCEAN, dimid_xt_ocean      ))
    call check(nf90_def_dim(ncid_monthly_avg, "yt_ocean"      , YT_OCEAN, dimid_yt_ocean      ))
    call check(nf90_def_dim(ncid_monthly_avg, "st_ocean"      , ST_OCEAN, dimid_st_ocean      ))
    call check(nf90_def_dim(ncid_monthly_avg, "st_edges_ocean", 51      , dimid_st_edges_ocean))
    call check(nf90_def_dim(ncid_monthly_avg, "xu_ocean"      , XT_OCEAN, dimid_xu_ocean      ))
    call check(nf90_def_dim(ncid_monthly_avg, "yu_ocean"      , YT_OCEAN, dimid_yu_ocean      ))
    call check(nf90_def_dim(ncid_monthly_avg, "sw_ocean"      , ST_OCEAN, dimid_sw_ocean      ))
    call check(nf90_def_dim(ncid_monthly_avg, "sw_edges_ocean", 51      , dimid_sw_edges_ocean))

!   dimid_4d = (/dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, dimid_time/)
!   dimid_3d = (/dimid_xt_ocean, dimid_yt_ocean,                 dimid_time/)
!   dimid_2d = (/dimid_xt_ocean, dimid_yt_ocean                            /)

    call check(nf90_def_var(ncid_monthly_avg, "time"          , NF90_DOUBLE, (/dimid_time          /), varid_time          ))
    call check(nf90_def_var(ncid_monthly_avg, "xt_ocean"      , NF90_DOUBLE, (/dimid_xt_ocean      /), varid_xt_ocean      ))
    call check(nf90_def_var(ncid_monthly_avg, "yt_ocean"      , NF90_DOUBLE, (/dimid_yt_ocean      /), varid_yt_ocean      ))
    call check(nf90_def_var(ncid_monthly_avg, "st_ocean"      , NF90_DOUBLE, (/dimid_st_ocean      /), varid_st_ocean      ))
    call check(nf90_def_var(ncid_monthly_avg, "st_edges_ocean", NF90_DOUBLE, (/dimid_st_edges_ocean/), varid_st_edges_ocean))
    call check(nf90_def_var(ncid_monthly_avg, "xu_ocean"      , NF90_DOUBLE, (/dimid_xu_ocean      /), varid_xu_ocean      ))
    call check(nf90_def_var(ncid_monthly_avg, "yu_ocean"      , NF90_DOUBLE, (/dimid_yu_ocean      /), varid_yu_ocean      ))
    call check(nf90_def_var(ncid_monthly_avg, "sw_ocean"      , NF90_DOUBLE, (/dimid_sw_ocean      /), varid_sw_ocean      ))
    call check(nf90_def_var(ncid_monthly_avg, "sw_edges_ocean", NF90_DOUBLE, (/dimid_sw_edges_ocean/), varid_sw_edges_ocean))

    call check(nf90_def_var(ncid_monthly_avg, "geolon_t", NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean/), varid_geolon_t))
    call check(nf90_def_var(ncid_monthly_avg, "geolat_t", NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean/), varid_geolat_t))

    call check(nf90_def_var(ncid_monthly_avg, "temp"     , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, 1/), varid_temp     ))
    call check(nf90_def_var(ncid_monthly_avg, "salt"     , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, 1/), varid_salt     ))
    call check(nf90_def_var(ncid_monthly_avg, "u"        , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, 1/), varid_u        ))
    call check(nf90_def_var(ncid_monthly_avg, "v"        , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, 1/), varid_v        ))
    call check(nf90_def_var(ncid_monthly_avg, "sea_level", NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_sea_level))
    call check(nf90_def_var(ncid_monthly_avg, "wt"       , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean, dimid_sw_ocean, 1/), varid_wt       ))
    call check(nf90_def_var(ncid_monthly_avg, "tau_x"    , NF90_FLOAT, (/dimid_xu_ocean, dimid_yu_ocean,                 1/), varid_tau_x    ))
    call check(nf90_def_var(ncid_monthly_avg, "tau_y"    , NF90_FLOAT, (/dimid_xu_ocean, dimid_yu_ocean,                 1/), varid_tau_y    ))


    call check(nf90_put_att(ncid_monthly_avg, NF90_GLOBAL, "description", "SODA-MOM5 monthly average analysis from 5-day average output."))    

    call check(nf90_put_att(ncid_monthly_avg, varid_xt_ocean, "long_name"     , "tcell longitude"))
    call check(nf90_put_att(ncid_monthly_avg, varid_xt_ocean, "units"         , "degrees_E"      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_xt_ocean, "cartesian_axis", "X"              ))

    call check(nf90_put_att(ncid_monthly_avg, varid_yt_ocean, "long_name"     , "tcell latitude"))
    call check(nf90_put_att(ncid_monthly_avg, varid_yt_ocean, "units"         , "degrees_N"     ))
    call check(nf90_put_att(ncid_monthly_avg, varid_yt_ocean, "cartesian_axis", "Y"             ))

    call check(nf90_put_att(ncid_monthly_avg, varid_time, "long_name"     , "time"       ))
    call check(nf90_put_att(ncid_monthly_avg, varid_time, "units"         , "month"      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_time, "cartesian_axis", "T"          ))
!   call check(nf90_put_att(ncid_monthly_avg, varid_time, "calendar_type" , "JULIAN"     ))
!   call check(nf90_put_att(ncid_monthly_avg, varid_time, "calendar"      , "JULIAN"     ))
!   call check(nf90_put_att(ncid_monthly_avg, varid_time, "bounds"        , "time_bounds"))

    call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "long_name"     , "tcell zstar depth"))
    call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "units"         , "meters"           ))
    call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "cartesian_axis", "Z"                ))
    call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "positive"      , "down"             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "edges"         , "st_edges_ocean"   ))

    call check(nf90_put_att(ncid_monthly_avg, varid_st_edges_ocean, "long_name"     , "tcell zstar depth edges"))
    call check(nf90_put_att(ncid_monthly_avg, varid_st_edges_ocean, "units"         , "meters"                 ))
    call check(nf90_put_att(ncid_monthly_avg, varid_st_edges_ocean, "cartesian_axis", "Z"                      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_st_edges_ocean, "positive"      , "down"                   ))

    call check(nf90_put_att(ncid_monthly_avg, varid_xu_ocean, "long_name"     , "ucell longitude"))
    call check(nf90_put_att(ncid_monthly_avg, varid_xu_ocean, "units"         , "degrees_E"      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_xu_ocean, "cartesian_axis", "X"              ))

    call check(nf90_put_att(ncid_monthly_avg, varid_yu_ocean, "long_name"     , "ucell latitude"))
    call check(nf90_put_att(ncid_monthly_avg, varid_yu_ocean, "units"         , "degrees_N"     ))
    call check(nf90_put_att(ncid_monthly_avg, varid_yu_ocean, "cartesian_axis", "Y"             ))

    call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "long_name"     , "ucell zstar depth"))
    call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "units"         , "meters"           ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "cartesian_axis", "Z"                ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "positive"      , "down"             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "edges"         , "sw_edges_ocean"   ))

    call check(nf90_put_att(ncid_monthly_avg, varid_sw_edges_ocean, "long_name"     , "ucell zstar depth edges"))
    call check(nf90_put_att(ncid_monthly_avg, varid_sw_edges_ocean, "units"         , "meters"                 ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sw_edges_ocean, "cartesian_axis", "Z"                      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sw_edges_ocean, "positive"      , "down"                   ))

    call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "long_name"    , "tracer longitude" ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "units"        , "degrees_E"        ))
!   call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "valid_range"  , "-281.f, 361.f"    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "valid_range"  , (/-281., 361./)    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "missing_value", 1.e+20             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "_FillValue"   , 1.e+20             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "cell_methods" , "time: point"      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "coordinates"  , "geolon_t geolat_t"))

    call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "long_name"    , "tracer latitude"  ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "units"        , "degrees_N"        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "valid_range"  , (/-91., 91./)      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "missing_value", 1.e+20             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "_FillValue"   , 1.e+20             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "cell_methods" , "time: point"      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "coordinates"  , "geolon_t geolat_t"))

    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "long_name"    , "Potential temperature"          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "units"        , "degrees C"                      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "valid_range"  , (/-10., 500./)                   ))
    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "missing_value", -1.e+20                          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "_FillValue"   , -1.e+20                          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "cell_methods" , "time: mean"                     ))
    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "coordinates"  , "geolon_t geolat_t"              ))
    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "standard_name", "sea_water_potential_temperature"))

    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "long_name"    , "Practical Salinity"))
    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "units"        , "psu"               ))
!   call check(nf90_put_att(ncid_monthly_avg, varid_salt, "valid_range"  , "-10.f, 100.f"      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "valid_range"  , (/-10., 100./)      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "missing_value", -1.e+20             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "_FillValue"   , -1.e+20             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "cell_methods" , "time: mean"        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "coordinates"  , "geolon_t geolat_t" ))
    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "standard_name", "sea_water_salinity"))

    call check(nf90_put_att(ncid_monthly_avg, varid_u, "long_name"    , "i-current"           ))
    call check(nf90_put_att(ncid_monthly_avg, varid_u, "units"        , "m/sec"               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_u, "valid_range"  , (/-10., 10./)         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_u, "missing_value", -1.e+20               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_u, "_FillValue"   , -1.e+20               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_u, "cell_methods" , "time: mean"          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_u, "coordinates"  , "geolon_t geolat_t"   ))
    call check(nf90_put_att(ncid_monthly_avg, varid_u, "standard_name", "sea_water_x_velocity"))

    call check(nf90_put_att(ncid_monthly_avg, varid_v, "long_name"    , "j-current"           ))
    call check(nf90_put_att(ncid_monthly_avg, varid_v, "units"        , "m/sec"               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_v, "valid_range"  , (/-10., 10./)         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_v, "missing_value", -1.e+20               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_v, "_FillValue"   , -1.e+20               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_v, "cell_methods" , "time: mean"          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_v, "coordinates"  , "geolon_t geolat_t"   ))
    call check(nf90_put_att(ncid_monthly_avg, varid_v, "standard_name", "sea_water_y_velocity"))

    call check(nf90_put_att(ncid_monthly_avg, varid_sea_level, "long_name"    , "effective sea level (eta_t + patm/(rho0*g)) on T cells"))
    call check(nf90_put_att(ncid_monthly_avg, varid_sea_level, "units"        , "meter"                         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sea_level, "valid_range"  , (/-1000., 1000./)               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sea_level, "missing_value", -1.e+20                         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sea_level, "_FillValue"   , -1.e+20                         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sea_level, "cell_methods" , "time: mean"                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sea_level, "coordinates"  , "geolon_t geolat_t"             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sea_level, "standard_name", "sea_surface_height_above_geoid"))

    call check(nf90_put_att(ncid_monthly_avg, varid_wt, "long_name"    , "dia-surface velocity T-points"))
    call check(nf90_put_att(ncid_monthly_avg, varid_wt, "units"        , "m/sec"                        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_wt, "valid_range"  , (/-100000., 100000./)          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_wt, "missing_value", -1.e+20                        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_wt, "_FillValue"   , -1.e+20                        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_wt, "cell_methods" , "time: mean"                   ))
    call check(nf90_put_att(ncid_monthly_avg, varid_wt, "coordinates"  , "geolon_t geolat_t"            ))

    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "long_name"    , "i-directed wind stress forcing u-velocity"))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "units"        , "N/m^2"                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "valid_range"  , (/-10., 10./)              ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "missing_value", -1.e+20                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "_FillValue"   , -1.e+20                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "cell_methods" , "time: mean"               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "coordinates"  , "geolon_t geolat_t"        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "standard_name", "surface_downward_x_stress"))

    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "long_name"    , "j-directed wind stress forcing v-velocity"))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "units"        , "N/m^2"                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "valid_range"  , (/-10., 10./)              ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "missing_value", -1.e+20                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "_FillValue"   , -1.e+20                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "cell_methods" , "time: mean"               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "coordinates"  , "geolon_t geolat_t"        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "standard_name", "surface_downward_y_stress"))

    call check(nf90_enddef(ncid_monthly_avg))


    ! write out all vars' data
    call check(nf90_put_var(ncid_monthly_avg, varid_time          , 0                 ))

    call check(nf90_put_var(ncid_monthly_avg, varid_xt_ocean      , var_xt_ocean      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_yt_ocean      , var_yt_ocean      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_st_ocean      , var_st_ocean      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_st_edges_ocean, var_st_edges_ocean))
    call check(nf90_put_var(ncid_monthly_avg, varid_xu_ocean      , var_xu_ocean      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_yu_ocean      , var_yu_ocean      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_sw_ocean      , var_sw_ocean      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_sw_edges_ocean, var_sw_edges_ocean))

    call check(nf90_put_var(ncid_monthly_avg, varid_geolon_t, var_geolon_t))
    call check(nf90_put_var(ncid_monthly_avg, varid_geolat_t, var_geolat_t))

    call check(nf90_put_var(ncid_monthly_avg, varid_temp     , avg_temp     ))
    call check(nf90_put_var(ncid_monthly_avg, varid_salt     , avg_salt     ))
    call check(nf90_put_var(ncid_monthly_avg, varid_u        , avg_u        ))
    call check(nf90_put_var(ncid_monthly_avg, varid_v        , avg_v        ))
    call check(nf90_put_var(ncid_monthly_avg, varid_sea_level, avg_sea_level))
    call check(nf90_put_var(ncid_monthly_avg, varid_wt       , avg_wt       ))
    call check(nf90_put_var(ncid_monthly_avg, varid_tau_x    , avg_tau_x    ))
    call check(nf90_put_var(ncid_monthly_avg, varid_tau_y    , avg_tau_y    ))

    call check(nf90_close(ncid_monthly_avg))

    call check(nf90_close(ncid_dim))

    write(*, *) "*** SUCCESS writing monthly average results for " // str_yr // str_mon // " ***"

    return
end subroutine write_out_5dayAvg2monthlyAvg_SODA_MOM


subroutine write_out_5dayAvg2monthlyAvg_SODA_MOM_ice(curr_yr, curr_mon, avg_CN, avg_HI, xt, yt, ct)
    use netcdf
    implicit none

    integer,                           intent(in) :: curr_yr, curr_mon, xt, yt, ct
    real   , dimension(xt, yt, ct, 1), intent(in) :: avg_CN
    real   , dimension(xt, yt,     1), intent(in) :: avg_HI

    ! local vars
    character(len = *), parameter :: DIR_ROOT_SODA_ANALYSIS = "/glade/scratch/lgchen/project/MOM_run/SODA_MOM5_sea_ice_2011_2013/test2/MOM_SODA"  &
        , FN_NC_DIM = trim(DIR_ROOT_SODA_ANALYSIS) // "/20110101/history/20110101.ice_month.nc"

    character(len=4) :: str_yr
    character(len=2) :: str_mon

    integer :: ncid_dim, ncid_monthly_avg, varid_CN, varid_HI  &
        , varid_time, varid_xt, varid_yt, varid_ct, varid_xb, varid_yb
    integer :: dimid_time, dimid_xt, dimid_yt, dimid_ct, dimid_xb, dimid_yb

    ! vars to receive dim data
    real(kind = 8) :: var_xt(xt), var_yt(yt), var_ct(ct), var_xb(1441), var_yb(1071)


    write(unit=str_yr , fmt="(I4.4)") curr_yr
    write(unit=str_mon, fmt="(I2.2)") curr_mon

    ! open and read the dim file data
    call check(nf90_open(trim(FN_NC_DIM), nf90_nowrite, ncid_dim))
    
    call check(nf90_inq_varid(ncid_dim, "xt", varid_xt))
    call check(nf90_inq_varid(ncid_dim, "yt", varid_yt))
    call check(nf90_inq_varid(ncid_dim, "ct", varid_ct))
    call check(nf90_inq_varid(ncid_dim, "xb", varid_xb))
    call check(nf90_inq_varid(ncid_dim, "yb", varid_yb))

    call check(nf90_get_var(ncid_dim, varid_xt, var_xt))
    call check(nf90_get_var(ncid_dim, varid_yt, var_yt))
    call check(nf90_get_var(ncid_dim, varid_ct, var_ct))
    call check(nf90_get_var(ncid_dim, varid_xb, var_xb))
    call check(nf90_get_var(ncid_dim, varid_yb, var_yb))

!   call check(nf90_close(ncid_dim))


    ! write out the SODA-MOM monthly average result to a netCDF file
    ! first define related dimensions, variables and their associated attributes
    call check(nf90_create(DIR_ROOT_SODA_ANALYSIS // "/monthly_avg/" // str_yr // str_mon // ".ice_month.nc", NF90_CLOBBER, ncid_monthly_avg))

    call check(nf90_def_dim(ncid_monthly_avg, "time", 1   , dimid_time))
    call check(nf90_def_dim(ncid_monthly_avg, "xt"  , xt  , dimid_xt  ))
    call check(nf90_def_dim(ncid_monthly_avg, "yt"  , yt  , dimid_yt  ))
    call check(nf90_def_dim(ncid_monthly_avg, "ct"  , ct  , dimid_ct  ))
    call check(nf90_def_dim(ncid_monthly_avg, "xb"  , 1441, dimid_xb  ))
    call check(nf90_def_dim(ncid_monthly_avg, "yb"  , 1071, dimid_yb  ))

!   dimid_4d = (/dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, dimid_time/)
!   dimid_3d = (/dimid_xt_ocean, dimid_yt_ocean,                 dimid_time/)
!   dimid_2d = (/dimid_xt_ocean, dimid_yt_ocean                            /)

    call check(nf90_def_var(ncid_monthly_avg, "time", NF90_DOUBLE, (/dimid_time/), varid_time))
    call check(nf90_def_var(ncid_monthly_avg, "xt"  , NF90_DOUBLE, (/dimid_xt  /), varid_xt  ))
    call check(nf90_def_var(ncid_monthly_avg, "yt"  , NF90_DOUBLE, (/dimid_yt  /), varid_yt  ))
    call check(nf90_def_var(ncid_monthly_avg, "ct"  , NF90_DOUBLE, (/dimid_ct  /), varid_ct  ))
    call check(nf90_def_var(ncid_monthly_avg, "xb"  , NF90_DOUBLE, (/dimid_xb  /), varid_xb  ))
    call check(nf90_def_var(ncid_monthly_avg, "yb"  , NF90_DOUBLE, (/dimid_yb  /), varid_yb  ))

    call check(nf90_def_var(ncid_monthly_avg, "CN"     , NF90_FLOAT, (/dimid_xt, dimid_yt, dimid_ct, 1/), varid_CN))
    call check(nf90_def_var(ncid_monthly_avg, "HI"     , NF90_FLOAT, (/dimid_xt, dimid_yt,           1/), varid_HI))


    call check(nf90_put_att(ncid_monthly_avg, NF90_GLOBAL, "description", "SODA-MOM5 ice monthly average from 5-day average output."))    

    call check(nf90_put_att(ncid_monthly_avg, varid_xt, "long_name"     , "longitude"))
    call check(nf90_put_att(ncid_monthly_avg, varid_xt, "units"         , "degrees_E"))
    call check(nf90_put_att(ncid_monthly_avg, varid_xt, "cartesian_axis", "X"        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_xt, "edges"         , "xb"       ))

    call check(nf90_put_att(ncid_monthly_avg, varid_yt, "long_name"     , "latitude" ))
    call check(nf90_put_att(ncid_monthly_avg, varid_yt, "units"         , "degrees_N"))
    call check(nf90_put_att(ncid_monthly_avg, varid_yt, "cartesian_axis", "Y"        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_yt, "edges"         , "yb"       ))

    call check(nf90_put_att(ncid_monthly_avg, varid_time, "long_name"     , "time"       ))
    call check(nf90_put_att(ncid_monthly_avg, varid_time, "units"         , "month"      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_time, "cartesian_axis", "T"          ))
!   call check(nf90_put_att(ncid_monthly_avg, varid_time, "calendar_type" , "JULIAN"     ))
!   call check(nf90_put_att(ncid_monthly_avg, varid_time, "calendar"      , "JULIAN"     ))
!   call check(nf90_put_att(ncid_monthly_avg, varid_time, "bounds"        , "time_bounds"))

    call check(nf90_put_att(ncid_monthly_avg, varid_ct, "long_name"     , "thickness"))
    call check(nf90_put_att(ncid_monthly_avg, varid_ct, "units"         , "meters"   ))
    call check(nf90_put_att(ncid_monthly_avg, varid_ct, "cartesian_axis", "Z"        ))

    call check(nf90_put_att(ncid_monthly_avg, varid_xb, "long_name"     , "longitude"))
    call check(nf90_put_att(ncid_monthly_avg, varid_xb, "units"         , "degrees_E"))
    call check(nf90_put_att(ncid_monthly_avg, varid_xb, "cartesian_axis", "X"        ))

    call check(nf90_put_att(ncid_monthly_avg, varid_yb, "long_name"     , "latitude" ))
    call check(nf90_put_att(ncid_monthly_avg, varid_yb, "units"         , "degrees_N"))
    call check(nf90_put_att(ncid_monthly_avg, varid_yb, "cartesian_axis", "Y"        ))



    call check(nf90_put_att(ncid_monthly_avg, varid_CN, "long_name"    , "ice concentration"))
    call check(nf90_put_att(ncid_monthly_avg, varid_CN, "units"        , "0-1"              ))
    call check(nf90_put_att(ncid_monthly_avg, varid_CN, "missing_value", -1.e+34            ))
    call check(nf90_put_att(ncid_monthly_avg, varid_CN, "_FillValue"   , -1.e+34            ))
    call check(nf90_put_att(ncid_monthly_avg, varid_CN, "cell_methods" , "time: mean"       ))

    call check(nf90_put_att(ncid_monthly_avg, varid_HI, "long_name"    , "ice thickness"))
    call check(nf90_put_att(ncid_monthly_avg, varid_HI, "units"        , "m-ice"        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_HI, "missing_value", -1.e+34        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_HI, "_FillValue"   , -1.e+34        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_HI, "cell_methods" , "time: mean"   ))

    call check(nf90_enddef(ncid_monthly_avg))


    ! write out all vars' data
    call check(nf90_put_var(ncid_monthly_avg, varid_time, 0     ))
    call check(nf90_put_var(ncid_monthly_avg, varid_xt  , var_xt))
    call check(nf90_put_var(ncid_monthly_avg, varid_yt  , var_yt))
    call check(nf90_put_var(ncid_monthly_avg, varid_ct  , var_ct))
    call check(nf90_put_var(ncid_monthly_avg, varid_xb  , var_xb))
    call check(nf90_put_var(ncid_monthly_avg, varid_yb  , var_yb))

    call check(nf90_put_var(ncid_monthly_avg, varid_CN, avg_CN))
    call check(nf90_put_var(ncid_monthly_avg, varid_HI, avg_HI))


    call check(nf90_close(ncid_monthly_avg))

    call check(nf90_close(ncid_dim))

    write(*, *) "*** SUCCESS writing ice monthly average results for " // str_yr // str_mon // " ***"

    return
end subroutine write_out_5dayAvg2monthlyAvg_SODA_MOM_ice




subroutine write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981(curr_yr, curr_mon, avg_temp, avg_salt, avg_u  &
    , avg_v, avg_ssh, avg_mlt, avg_mlp, avg_mls, avg_pbot,  avg_wt, avg_prho, avg_tau_x, avg_tau_y, XT_OCEAN, YT_OCEAN, ST_OCEAN)
    use netcdf
    implicit none

    integer,                                             intent(in) :: curr_yr, curr_mon, XT_OCEAN, YT_OCEAN, ST_OCEAN
    real   , dimension(XT_OCEAN, YT_OCEAN, ST_OCEAN, 1), intent(in) :: avg_temp, avg_salt, avg_u, avg_v, avg_wt, avg_prho
    real   , dimension(XT_OCEAN, YT_OCEAN,           1), intent(in) :: avg_ssh, avg_mlt, avg_mlp, avg_mls, avg_pbot, avg_tau_x, avg_tau_y

    ! local vars
    character(len = *), parameter :: DIR_ROOT_SODA_ANALYSIS = "/glade/scratch/lgchen/project/MOM_run/SODA_1979_2014/test1/MOM_SODA"  &
        , FN_NC_DIM = trim(DIR_ROOT_SODA_ANALYSIS) // "/19790101/history/19790101.5days.nc"

    character(len=4) :: str_yr
    character(len=2) :: str_mon

    integer :: ncid_dim, ncid_monthly_avg, varid_temp, varid_salt, varid_u, varid_v, varid_ssh, varid_mlt, varid_mlp  &
        , varid_mls, varid_pbot, varid_wt, varid_prho, varid_tau_x  &
        , varid_tau_y, varid_time, varid_xt_ocean, varid_yt_ocean, varid_st_ocean, varid_st_edges_ocean, varid_xu_ocean, varid_yu_ocean  &
        , varid_sw_ocean, varid_sw_edges_ocean, varid_geolon_t, varid_geolat_t
    integer :: dimid_time, dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, dimid_st_edges_ocean, dimid_xu_ocean, dimid_yu_ocean  &
        , dimid_sw_ocean, dimid_sw_edges_ocean

    ! vars to receive dim data
    real(kind = 8) :: var_xt_ocean(XT_OCEAN), var_yt_ocean(YT_OCEAN), var_st_ocean(ST_OCEAN)  &
        , var_st_edges_ocean(51), var_xu_ocean(XT_OCEAN), var_yu_ocean(YT_OCEAN), var_sw_ocean(ST_OCEAN)  &
        , var_sw_edges_ocean(51), var_geolon_t(XT_OCEAN, YT_OCEAN), var_geolat_t(XT_OCEAN, YT_OCEAN)


    write(unit=str_yr , fmt="(I4.4)") curr_yr
    write(unit=str_mon, fmt="(I2.2)") curr_mon

    ! open and read the dim file data
    call check(nf90_open(trim(FN_NC_DIM), nf90_nowrite, ncid_dim))
    
    call check(nf90_inq_varid(ncid_dim, "xt_ocean"      , varid_xt_ocean      ))
    call check(nf90_inq_varid(ncid_dim, "yt_ocean"      , varid_yt_ocean      ))
    call check(nf90_inq_varid(ncid_dim, "st_ocean"      , varid_st_ocean      ))
    call check(nf90_inq_varid(ncid_dim, "st_edges_ocean", varid_st_edges_ocean))
    call check(nf90_inq_varid(ncid_dim, "xu_ocean"      , varid_xu_ocean      ))
    call check(nf90_inq_varid(ncid_dim, "yu_ocean"      , varid_yu_ocean      ))
    call check(nf90_inq_varid(ncid_dim, "sw_ocean"      , varid_sw_ocean      ))
    call check(nf90_inq_varid(ncid_dim, "sw_edges_ocean", varid_sw_edges_ocean))
    call check(nf90_inq_varid(ncid_dim, "geolon_t"      , varid_geolon_t      ))
    call check(nf90_inq_varid(ncid_dim, "geolat_t"      , varid_geolat_t      ))

    call check(nf90_get_var(ncid_dim, varid_xt_ocean      , var_xt_ocean      ))
    call check(nf90_get_var(ncid_dim, varid_yt_ocean      , var_yt_ocean      ))
    call check(nf90_get_var(ncid_dim, varid_st_ocean      , var_st_ocean      ))
    call check(nf90_get_var(ncid_dim, varid_st_edges_ocean, var_st_edges_ocean))
    call check(nf90_get_var(ncid_dim, varid_xu_ocean      , var_xu_ocean      ))
    call check(nf90_get_var(ncid_dim, varid_yu_ocean      , var_yu_ocean      ))
    call check(nf90_get_var(ncid_dim, varid_sw_ocean      , var_sw_ocean      ))
    call check(nf90_get_var(ncid_dim, varid_sw_edges_ocean, var_sw_edges_ocean))
    call check(nf90_get_var(ncid_dim, varid_geolon_t      , var_geolon_t      ))
    call check(nf90_get_var(ncid_dim, varid_geolat_t      , var_geolat_t      ))

!   call check(nf90_close(ncid_dim))


    ! write out the SODA-MOM monthly average result to a netCDF file
    ! first define related dimensions, variables and their associated attributes
    call check(nf90_create(DIR_ROOT_SODA_ANALYSIS // "/monthly_avg/" // str_yr // str_mon // ".5days.nc", NF90_CLOBBER, ncid_monthly_avg))

    call check(nf90_def_dim(ncid_monthly_avg, "time"          , 1       , dimid_time          ))
    call check(nf90_def_dim(ncid_monthly_avg, "xt_ocean"      , XT_OCEAN, dimid_xt_ocean      ))
    call check(nf90_def_dim(ncid_monthly_avg, "yt_ocean"      , YT_OCEAN, dimid_yt_ocean      ))
    call check(nf90_def_dim(ncid_monthly_avg, "st_ocean"      , ST_OCEAN, dimid_st_ocean      ))
    call check(nf90_def_dim(ncid_monthly_avg, "st_edges_ocean", 51      , dimid_st_edges_ocean))
    call check(nf90_def_dim(ncid_monthly_avg, "xu_ocean"      , XT_OCEAN, dimid_xu_ocean      ))
    call check(nf90_def_dim(ncid_monthly_avg, "yu_ocean"      , YT_OCEAN, dimid_yu_ocean      ))
    call check(nf90_def_dim(ncid_monthly_avg, "sw_ocean"      , ST_OCEAN, dimid_sw_ocean      ))
    call check(nf90_def_dim(ncid_monthly_avg, "sw_edges_ocean", 51      , dimid_sw_edges_ocean))

!   dimid_4d = (/dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, dimid_time/)
!   dimid_3d = (/dimid_xt_ocean, dimid_yt_ocean,                 dimid_time/)
!   dimid_2d = (/dimid_xt_ocean, dimid_yt_ocean                            /)

    call check(nf90_def_var(ncid_monthly_avg, "time"          , NF90_DOUBLE, (/dimid_time          /), varid_time          ))
    call check(nf90_def_var(ncid_monthly_avg, "xt_ocean"      , NF90_DOUBLE, (/dimid_xt_ocean      /), varid_xt_ocean      ))
    call check(nf90_def_var(ncid_monthly_avg, "yt_ocean"      , NF90_DOUBLE, (/dimid_yt_ocean      /), varid_yt_ocean      ))
    call check(nf90_def_var(ncid_monthly_avg, "st_ocean"      , NF90_DOUBLE, (/dimid_st_ocean      /), varid_st_ocean      ))
    call check(nf90_def_var(ncid_monthly_avg, "st_edges_ocean", NF90_DOUBLE, (/dimid_st_edges_ocean/), varid_st_edges_ocean))
    call check(nf90_def_var(ncid_monthly_avg, "xu_ocean"      , NF90_DOUBLE, (/dimid_xu_ocean      /), varid_xu_ocean      ))
    call check(nf90_def_var(ncid_monthly_avg, "yu_ocean"      , NF90_DOUBLE, (/dimid_yu_ocean      /), varid_yu_ocean      ))
    call check(nf90_def_var(ncid_monthly_avg, "sw_ocean"      , NF90_DOUBLE, (/dimid_sw_ocean      /), varid_sw_ocean      ))
    call check(nf90_def_var(ncid_monthly_avg, "sw_edges_ocean", NF90_DOUBLE, (/dimid_sw_edges_ocean/), varid_sw_edges_ocean))

    call check(nf90_def_var(ncid_monthly_avg, "geolon_t", NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean/), varid_geolon_t))
    call check(nf90_def_var(ncid_monthly_avg, "geolat_t", NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean/), varid_geolat_t))

    call check(nf90_def_var(ncid_monthly_avg, "temp"     , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, 1/), varid_temp     ))
    call check(nf90_def_var(ncid_monthly_avg, "salt"     , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, 1/), varid_salt     ))
    call check(nf90_def_var(ncid_monthly_avg, "u"        , NF90_FLOAT, (/dimid_xu_ocean, dimid_yu_ocean, dimid_st_ocean, 1/), varid_u        ))
    call check(nf90_def_var(ncid_monthly_avg, "v"        , NF90_FLOAT, (/dimid_xu_ocean, dimid_yu_ocean, dimid_st_ocean, 1/), varid_v        ))
    call check(nf90_def_var(ncid_monthly_avg, "ssh"      , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_ssh      ))
    call check(nf90_def_var(ncid_monthly_avg, "mlt"      , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_mlt      ))
    call check(nf90_def_var(ncid_monthly_avg, "mlp"      , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_mlp      ))
    call check(nf90_def_var(ncid_monthly_avg, "mls"      , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_mls      ))
    call check(nf90_def_var(ncid_monthly_avg, "pbot"     , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_pbot     ))
    call check(nf90_def_var(ncid_monthly_avg, "wt"       , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean, dimid_sw_ocean, 1/), varid_wt       ))
    call check(nf90_def_var(ncid_monthly_avg, "prho"     , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, 1/), varid_prho     ))
    call check(nf90_def_var(ncid_monthly_avg, "tau_x"    , NF90_FLOAT, (/dimid_xu_ocean, dimid_yu_ocean,                 1/), varid_tau_x    ))
    call check(nf90_def_var(ncid_monthly_avg, "tau_y"    , NF90_FLOAT, (/dimid_xu_ocean, dimid_yu_ocean,                 1/), varid_tau_y    ))


    call check(nf90_put_att(ncid_monthly_avg, NF90_GLOBAL, "description", "SODA-MOM5 monthly average analysis from 5-day average output."))    

    call check(nf90_put_att(ncid_monthly_avg, varid_xt_ocean, "long_name"     , "tcell longitude"))
    call check(nf90_put_att(ncid_monthly_avg, varid_xt_ocean, "units"         , "degrees_E"      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_xt_ocean, "cartesian_axis", "X"              ))

    call check(nf90_put_att(ncid_monthly_avg, varid_yt_ocean, "long_name"     , "tcell latitude"))
    call check(nf90_put_att(ncid_monthly_avg, varid_yt_ocean, "units"         , "degrees_N"     ))
    call check(nf90_put_att(ncid_monthly_avg, varid_yt_ocean, "cartesian_axis", "Y"             ))

    call check(nf90_put_att(ncid_monthly_avg, varid_time, "long_name"     , "time"       ))
    call check(nf90_put_att(ncid_monthly_avg, varid_time, "units"         , "month"      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_time, "cartesian_axis", "T"          ))
!   call check(nf90_put_att(ncid_monthly_avg, varid_time, "calendar_type" , "JULIAN"     ))
!   call check(nf90_put_att(ncid_monthly_avg, varid_time, "calendar"      , "JULIAN"     ))
!   call check(nf90_put_att(ncid_monthly_avg, varid_time, "bounds"        , "time_bounds"))

    call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "long_name"     , "tcell zstar depth"))
    call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "units"         , "meters"           ))
    call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "cartesian_axis", "Z"                ))
    call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "positive"      , "down"             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "edges"         , "st_edges_ocean"   ))

    call check(nf90_put_att(ncid_monthly_avg, varid_st_edges_ocean, "long_name"     , "tcell zstar depth edges"))
    call check(nf90_put_att(ncid_monthly_avg, varid_st_edges_ocean, "units"         , "meters"                 ))
    call check(nf90_put_att(ncid_monthly_avg, varid_st_edges_ocean, "cartesian_axis", "Z"                      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_st_edges_ocean, "positive"      , "down"                   ))

    call check(nf90_put_att(ncid_monthly_avg, varid_xu_ocean, "long_name"     , "ucell longitude"))
    call check(nf90_put_att(ncid_monthly_avg, varid_xu_ocean, "units"         , "degrees_E"      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_xu_ocean, "cartesian_axis", "X"              ))

    call check(nf90_put_att(ncid_monthly_avg, varid_yu_ocean, "long_name"     , "ucell latitude"))
    call check(nf90_put_att(ncid_monthly_avg, varid_yu_ocean, "units"         , "degrees_N"     ))
    call check(nf90_put_att(ncid_monthly_avg, varid_yu_ocean, "cartesian_axis", "Y"             ))

    call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "long_name"     , "ucell zstar depth"))
    call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "units"         , "meters"           ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "cartesian_axis", "Z"                ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "positive"      , "down"             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "edges"         , "sw_edges_ocean"   ))

    call check(nf90_put_att(ncid_monthly_avg, varid_sw_edges_ocean, "long_name"     , "ucell zstar depth edges"))
    call check(nf90_put_att(ncid_monthly_avg, varid_sw_edges_ocean, "units"         , "meters"                 ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sw_edges_ocean, "cartesian_axis", "Z"                      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_sw_edges_ocean, "positive"      , "down"                   ))

    call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "long_name"    , "tracer longitude" ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "units"        , "degrees_E"        ))
!   call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "valid_range"  , "-281.f, 361.f"    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "valid_range"  , (/-281., 361./)    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "missing_value", 1.e+20             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "_FillValue"   , 1.e+20             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "cell_methods" , "time: point"      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "coordinates"  , "geolon_t geolat_t"))

    call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "long_name"    , "tracer latitude"  ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "units"        , "degrees_N"        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "valid_range"  , (/-91., 91./)      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "missing_value", 1.e+20             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "_FillValue"   , 1.e+20             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "cell_methods" , "time: point"      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "coordinates"  , "geolon_t geolat_t"))

    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "long_name"    , "Potential temperature"          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "units"        , "degrees C"                      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "valid_range"  , (/-10., 500./)                   ))
    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "missing_value", -1.e+20                          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "_FillValue"   , -1.e+20                          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "cell_methods" , "time: mean"                     ))
    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "coordinates"  , "geolon_t geolat_t"              ))
    call check(nf90_put_att(ncid_monthly_avg, varid_temp, "standard_name", "sea_water_potential_temperature"))

    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "long_name"    , "Practical Salinity"))
    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "units"        , "psu"               ))
!   call check(nf90_put_att(ncid_monthly_avg, varid_salt, "valid_range"  , "-10.f, 100.f"      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "valid_range"  , (/-10., 100./)      ))
    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "missing_value", -1.e+20             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "_FillValue"   , -1.e+20             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "cell_methods" , "time: mean"        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "coordinates"  , "geolon_t geolat_t" ))
    call check(nf90_put_att(ncid_monthly_avg, varid_salt, "standard_name", "sea_water_salinity"))

    call check(nf90_put_att(ncid_monthly_avg, varid_u, "long_name"    , "i-current"           ))
    call check(nf90_put_att(ncid_monthly_avg, varid_u, "units"        , "m/sec"               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_u, "valid_range"  , (/-10., 10./)         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_u, "missing_value", -1.e+20               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_u, "_FillValue"   , -1.e+20               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_u, "cell_methods" , "time: mean"          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_u, "coordinates"  , "geolon_t geolat_t"   ))
    call check(nf90_put_att(ncid_monthly_avg, varid_u, "standard_name", "sea_water_x_velocity"))

    call check(nf90_put_att(ncid_monthly_avg, varid_v, "long_name"    , "j-current"           ))
    call check(nf90_put_att(ncid_monthly_avg, varid_v, "units"        , "m/sec"               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_v, "valid_range"  , (/-10., 10./)         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_v, "missing_value", -1.e+20               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_v, "_FillValue"   , -1.e+20               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_v, "cell_methods" , "time: mean"          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_v, "coordinates"  , "geolon_t geolat_t"   ))
    call check(nf90_put_att(ncid_monthly_avg, varid_v, "standard_name", "sea_water_y_velocity"))

    call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "long_name"    , "effective sea level (eta_t + patm/(rho0*g)) on T cells"))
    call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "units"        , "meter"                         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "valid_range"  , (/-1000., 1000./)               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "missing_value", -1.e+20                         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "_FillValue"   , -1.e+20                         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "cell_methods" , "time: mean"                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "coordinates"  , "geolon_t geolat_t"             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "standard_name", "sea_surface_height_above_geoid"))

    call check(nf90_put_att(ncid_monthly_avg, varid_mlt, "long_name"    , "mixed layer depth determined by temperature criteria"))
    call check(nf90_put_att(ncid_monthly_avg, varid_mlt, "units"        , "m"                             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mlt, "valid_range"  , (/0., 1000000./)                ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mlt, "missing_value", -1.e+20                         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mlt, "_FillValue"   , -1.e+20                         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mlt, "cell_methods" , "time: mean"                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mlt, "coordinates"  , "geolon_t geolat_t"             ))

    call check(nf90_put_att(ncid_monthly_avg, varid_mlp, "long_name"    , "Depth of potential density mixed layer"))
    call check(nf90_put_att(ncid_monthly_avg, varid_mlp, "units"        , "m"                             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mlp, "valid_range"  , (/-1000000., 1000000./)         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mlp, "missing_value", -1.e+20                         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mlp, "_FillValue"   , -1.e+20                         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mlp, "cell_methods" , "time: mean"                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mlp, "coordinates"  , "geolon_t geolat_t"             ))

    call check(nf90_put_att(ncid_monthly_avg, varid_mls, "long_name"    , "mixed layer depth determined by salinity criteria"))
    call check(nf90_put_att(ncid_monthly_avg, varid_mls, "units"        , "m"                             ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mls, "valid_range"  , (/0., 1000000./)                ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mls, "missing_value", -1.e+20                         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mls, "_FillValue"   , -1.e+20                         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mls, "cell_methods" , "time: mean"                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_mls, "coordinates"  , "geolon_t geolat_t"             ))

    call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "long_name"    , "bottom pressure on T cells [Boussinesq (volume conserving) model]"))
    call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "units"        , "dbar"                           ))
    call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "valid_range"  , (/-1000000., 1000000./)          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "missing_value", -1.e+20                          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "_FillValue"   , -1.e+20                          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "cell_methods" , "time: mean"                     ))
    call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "coordinates"  , "geolon_t geolat_t"              ))
    call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "standard_name", "sea_water_pressure_at_sea_floor"))

    call check(nf90_put_att(ncid_monthly_avg, varid_wt, "long_name"    , "dia-surface velocity T-points"))
    call check(nf90_put_att(ncid_monthly_avg, varid_wt, "units"        , "m/sec"                        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_wt, "valid_range"  , (/-100000., 100000./)          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_wt, "missing_value", -1.e+20                        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_wt, "_FillValue"   , -1.e+20                        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_wt, "cell_methods" , "time: mean"                   ))
    call check(nf90_put_att(ncid_monthly_avg, varid_wt, "coordinates"  , "geolon_t geolat_t"            ))

    call check(nf90_put_att(ncid_monthly_avg, varid_prho, "long_name"    , "potential density referenced to 0 dbar"))
    call check(nf90_put_att(ncid_monthly_avg, varid_prho, "units"        , "kg/m^3"                         ))
    call check(nf90_put_att(ncid_monthly_avg, varid_prho, "valid_range"  , (/-10., 100000./)                ))
    call check(nf90_put_att(ncid_monthly_avg, varid_prho, "missing_value", -1.e+20                          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_prho, "_FillValue"   , -1.e+20                          ))
    call check(nf90_put_att(ncid_monthly_avg, varid_prho, "cell_methods" , "time: mean"                     ))
    call check(nf90_put_att(ncid_monthly_avg, varid_prho, "coordinates"  , "geolon_t geolat_t"              ))
    call check(nf90_put_att(ncid_monthly_avg, varid_prho, "standard_name", "sea_water_potential_density"    ))

    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "long_name"    , "i-directed wind stress forcing u-velocity"))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "units"        , "N/m^2"                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "valid_range"  , (/-10., 10./)              ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "missing_value", -1.e+20                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "_FillValue"   , -1.e+20                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "cell_methods" , "time: mean"               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "coordinates"  , "geolon_t geolat_t"        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_x, "standard_name", "surface_downward_x_stress"))

    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "long_name"    , "j-directed wind stress forcing v-velocity"))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "units"        , "N/m^2"                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "valid_range"  , (/-10., 10./)              ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "missing_value", -1.e+20                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "_FillValue"   , -1.e+20                    ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "cell_methods" , "time: mean"               ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "coordinates"  , "geolon_t geolat_t"        ))
    call check(nf90_put_att(ncid_monthly_avg, varid_tau_y, "standard_name", "surface_downward_y_stress"))

    call check(nf90_enddef(ncid_monthly_avg))


    ! write out all vars' data
    call check(nf90_put_var(ncid_monthly_avg, varid_time          , 0                 ))

    call check(nf90_put_var(ncid_monthly_avg, varid_xt_ocean      , var_xt_ocean      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_yt_ocean      , var_yt_ocean      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_st_ocean      , var_st_ocean      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_st_edges_ocean, var_st_edges_ocean))
    call check(nf90_put_var(ncid_monthly_avg, varid_xu_ocean      , var_xu_ocean      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_yu_ocean      , var_yu_ocean      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_sw_ocean      , var_sw_ocean      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_sw_edges_ocean, var_sw_edges_ocean))

    call check(nf90_put_var(ncid_monthly_avg, varid_geolon_t, var_geolon_t))
    call check(nf90_put_var(ncid_monthly_avg, varid_geolat_t, var_geolat_t))

    call check(nf90_put_var(ncid_monthly_avg, varid_temp     , avg_temp     ))
    call check(nf90_put_var(ncid_monthly_avg, varid_salt     , avg_salt     ))
    call check(nf90_put_var(ncid_monthly_avg, varid_u        , avg_u        ))
    call check(nf90_put_var(ncid_monthly_avg, varid_v        , avg_v        ))
    call check(nf90_put_var(ncid_monthly_avg, varid_ssh      , avg_ssh      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_mlt      , avg_mlt      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_mlp      , avg_mlp      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_mls      , avg_mls      ))
    call check(nf90_put_var(ncid_monthly_avg, varid_pbot     , avg_pbot     ))
    call check(nf90_put_var(ncid_monthly_avg, varid_wt       , avg_wt       ))
    call check(nf90_put_var(ncid_monthly_avg, varid_prho     , avg_prho     ))
    call check(nf90_put_var(ncid_monthly_avg, varid_tau_x    , avg_tau_x    ))
    call check(nf90_put_var(ncid_monthly_avg, varid_tau_y    , avg_tau_y    ))

    call check(nf90_close(ncid_monthly_avg))

    call check(nf90_close(ncid_dim))

    write(*, *) "*** SUCCESS writing monthly average results for " // str_yr // str_mon // " ***"

    return
end subroutine write_out_5dayAvg2monthlyAvg_SODA_MOM_1979_1981


! subroutine write_out_5dayAvg2monthlyAvg_SODA_MOM_339(curr_yr, curr_mon, avg_temp, avg_salt, avg_u  &
!     , avg_v, avg_ssh, avg_mlt, avg_mlp, avg_mls, avg_anompb, avg_sbd, avg_sbd_mld, avg_sbd_dhdt    &
!     , avg_sbd_horz, avg_sbd_vert, avg_wt, avg_prho, avg_taux, avg_tauy, avg_hflux_total, avg_net_heating  &
!     , avg_salt_flux_rstr, avg_salt_flux_total, XT_OCEAN, YT_OCEAN, ST_OCEAN)
! 
!     use netcdf
!     implicit none
! 
!     integer,                                             intent(in) :: curr_yr, curr_mon, XT_OCEAN, YT_OCEAN, ST_OCEAN
!     real   , dimension(XT_OCEAN, YT_OCEAN, ST_OCEAN, 1), intent(in) :: avg_temp, avg_salt, avg_u, avg_v, avg_wt, avg_prho
!     real   , dimension(XT_OCEAN, YT_OCEAN,           1), intent(in) :: avg_ssh, avg_mlt, avg_mlp, avg_mls, avg_anompb  &
!         , avg_sbd, avg_sbd_mld, avg_sbd_dhdt, avg_sbd_horz, avg_sbd_vert,  avg_taux, avg_tauy, avg_hflux_total         &
!         , avg_net_heating, avg_salt_flux_rstr, avg_salt_flux_total
! 
!     ! local vars
!     character(len = *), parameter :: DIR_ROOT_SODA_ANALYSIS = "/glade/scratch/lgchen/project/SODA_3.3.9"  &
!         , FN_NC_DIM = "/glade/p/umcp0006/SODA_338/ORIGINAL/ocean_5dy_1980_01_03.nc"
! 
!     character(len=4) :: str_yr
!     character(len=2) :: str_mon
! 
!     integer :: ncid_dim, ncid_monthly_avg, varid_temp, varid_salt, varid_u, varid_v, varid_ssh, varid_mlt, varid_mlp  &
!         , varid_mls, varid_anompb, avg_sbd, avg_sbd_mld, avg_sbd_dhdt, avg_sbd_horz, avg_sbd_vert, varid_wt, varid_prho, varid_taux  &
!         , varid_tauy, varid_hflux_total, varid_net_heating, varid_salt_flux_rstr, varid_salt_flux_total, varid_time, varid_xt_ocean  &
!         , varid_yt_ocean, varid_st_ocean, varid_st_edges_ocean, varid_xu_ocean, varid_yu_ocean  &
!         , varid_sw_ocean, varid_sw_edges_ocean
!     integer :: dimid_time, dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, dimid_st_edges_ocean, dimid_xu_ocean, dimid_yu_ocean  &
!         , dimid_sw_ocean, dimid_sw_edges_ocean
! 
!     ! vars to receive dim data
!     real(kind = 8) :: var_xt_ocean(XT_OCEAN), var_yt_ocean(YT_OCEAN), var_st_ocean(ST_OCEAN)  &
!         , var_st_edges_ocean(51), var_xu_ocean(XT_OCEAN), var_yu_ocean(YT_OCEAN), var_sw_ocean(ST_OCEAN)  &
!         , var_sw_edges_ocean(51), var_geolon_t(XT_OCEAN, YT_OCEAN), var_geolat_t(XT_OCEAN, YT_OCEAN)
! 
! 
!     write(unit=str_yr , fmt="(I4.4)") curr_yr
!     write(unit=str_mon, fmt="(I2.2)") curr_mon
! 
!     ! open and read the dim file data
!     call check(nf90_open(trim(FN_NC_DIM), nf90_nowrite, ncid_dim))
!     
!     call check(nf90_inq_varid(ncid_dim, "xt_ocean"      , varid_xt_ocean      ))
!     call check(nf90_inq_varid(ncid_dim, "yt_ocean"      , varid_yt_ocean      ))
!     call check(nf90_inq_varid(ncid_dim, "st_ocean"      , varid_st_ocean      ))
!     call check(nf90_inq_varid(ncid_dim, "st_edges_ocean", varid_st_edges_ocean))
!     call check(nf90_inq_varid(ncid_dim, "xu_ocean"      , varid_xu_ocean      ))
!     call check(nf90_inq_varid(ncid_dim, "yu_ocean"      , varid_yu_ocean      ))
!     call check(nf90_inq_varid(ncid_dim, "sw_ocean"      , varid_sw_ocean      ))
!     call check(nf90_inq_varid(ncid_dim, "sw_edges_ocean", varid_sw_edges_ocean))
! !   call check(nf90_inq_varid(ncid_dim, "geolon_t"      , varid_geolon_t      ))
! !   call check(nf90_inq_varid(ncid_dim, "geolat_t"      , varid_geolat_t      ))
! 
!     call check(nf90_get_var(ncid_dim, varid_xt_ocean      , var_xt_ocean      ))
!     call check(nf90_get_var(ncid_dim, varid_yt_ocean      , var_yt_ocean      ))
!     call check(nf90_get_var(ncid_dim, varid_st_ocean      , var_st_ocean      ))
!     call check(nf90_get_var(ncid_dim, varid_st_edges_ocean, var_st_edges_ocean))
!     call check(nf90_get_var(ncid_dim, varid_xu_ocean      , var_xu_ocean      ))
!     call check(nf90_get_var(ncid_dim, varid_yu_ocean      , var_yu_ocean      ))
!     call check(nf90_get_var(ncid_dim, varid_sw_ocean      , var_sw_ocean      ))
!     call check(nf90_get_var(ncid_dim, varid_sw_edges_ocean, var_sw_edges_ocean))
! !   call check(nf90_get_var(ncid_dim, varid_geolon_t      , var_geolon_t      ))
! !   call check(nf90_get_var(ncid_dim, varid_geolat_t      , var_geolat_t      ))
! 
! !   call check(nf90_close(ncid_dim))
! 
! 
!     ! write out the SODA-MOM monthly average result to a netCDF file
!     ! first define related dimensions, variables and their associated attributes
!     call check(nf90_create(DIR_ROOT_SODA_ANALYSIS // "/monthly_avg/ocean_" // str_yr // str_mon // ".nc", NF90_CLOBBER, ncid_monthly_avg))
! 
!     call check(nf90_def_dim(ncid_monthly_avg, "time"          , 1       , dimid_time          ))
!     call check(nf90_def_dim(ncid_monthly_avg, "xt_ocean"      , XT_OCEAN, dimid_xt_ocean      ))
!     call check(nf90_def_dim(ncid_monthly_avg, "yt_ocean"      , YT_OCEAN, dimid_yt_ocean      ))
!     call check(nf90_def_dim(ncid_monthly_avg, "st_ocean"      , ST_OCEAN, dimid_st_ocean      ))
!     call check(nf90_def_dim(ncid_monthly_avg, "st_edges_ocean", 51      , dimid_st_edges_ocean))
!     call check(nf90_def_dim(ncid_monthly_avg, "xu_ocean"      , XT_OCEAN, dimid_xu_ocean      ))
!     call check(nf90_def_dim(ncid_monthly_avg, "yu_ocean"      , YT_OCEAN, dimid_yu_ocean      ))
!     call check(nf90_def_dim(ncid_monthly_avg, "sw_ocean"      , ST_OCEAN, dimid_sw_ocean      ))
!     call check(nf90_def_dim(ncid_monthly_avg, "sw_edges_ocean", 51      , dimid_sw_edges_ocean))
! 
! !   dimid_4d = (/dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, dimid_time/)
! !   dimid_3d = (/dimid_xt_ocean, dimid_yt_ocean,                 dimid_time/)
! !   dimid_2d = (/dimid_xt_ocean, dimid_yt_ocean                            /)
! 
!     call check(nf90_def_var(ncid_monthly_avg, "time"          , NF90_DOUBLE, (/dimid_time          /), varid_time          ))
!     call check(nf90_def_var(ncid_monthly_avg, "xt_ocean"      , NF90_DOUBLE, (/dimid_xt_ocean      /), varid_xt_ocean      ))
!     call check(nf90_def_var(ncid_monthly_avg, "yt_ocean"      , NF90_DOUBLE, (/dimid_yt_ocean      /), varid_yt_ocean      ))
!     call check(nf90_def_var(ncid_monthly_avg, "st_ocean"      , NF90_DOUBLE, (/dimid_st_ocean      /), varid_st_ocean      ))
!     call check(nf90_def_var(ncid_monthly_avg, "st_edges_ocean", NF90_DOUBLE, (/dimid_st_edges_ocean/), varid_st_edges_ocean))
!     call check(nf90_def_var(ncid_monthly_avg, "xu_ocean"      , NF90_DOUBLE, (/dimid_xu_ocean      /), varid_xu_ocean      ))
!     call check(nf90_def_var(ncid_monthly_avg, "yu_ocean"      , NF90_DOUBLE, (/dimid_yu_ocean      /), varid_yu_ocean      ))
!     call check(nf90_def_var(ncid_monthly_avg, "sw_ocean"      , NF90_DOUBLE, (/dimid_sw_ocean      /), varid_sw_ocean      ))
!     call check(nf90_def_var(ncid_monthly_avg, "sw_edges_ocean", NF90_DOUBLE, (/dimid_sw_edges_ocean/), varid_sw_edges_ocean))
! 
! !   call check(nf90_def_var(ncid_monthly_avg, "geolon_t", NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean/), varid_geolon_t))
! !   call check(nf90_def_var(ncid_monthly_avg, "geolat_t", NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean/), varid_geolat_t))
! 
!     call check(nf90_def_var(ncid_monthly_avg, "temp"           , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, 1/), varid_temp           ))
!     call check(nf90_def_var(ncid_monthly_avg, "salt"           , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, 1/), varid_salt           ))
!     call check(nf90_def_var(ncid_monthly_avg, "u"              , NF90_FLOAT, (/dimid_xu_ocean, dimid_yu_ocean, dimid_st_ocean, 1/), varid_u              ))
!     call check(nf90_def_var(ncid_monthly_avg, "v"              , NF90_FLOAT, (/dimid_xu_ocean, dimid_yu_ocean, dimid_st_ocean, 1/), varid_v              ))
!     call check(nf90_def_var(ncid_monthly_avg, "ssh"            , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_ssh            ))
!     call check(nf90_def_var(ncid_monthly_avg, "mlt"            , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_mlt            ))
!     call check(nf90_def_var(ncid_monthly_avg, "mlp"            , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_mlp            ))
!     call check(nf90_def_var(ncid_monthly_avg, "mls"            , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_mls            ))
! 
!     call check(nf90_def_var(ncid_monthly_avg, "anompb"         , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_anompb         ))
!     call check(nf90_def_var(ncid_monthly_avg, "sbd"            , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_sbd            ))
!     call check(nf90_def_var(ncid_monthly_avg, "sbd_mld"        , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_sbd_mld        ))
!     call check(nf90_def_var(ncid_monthly_avg, "sbd_dhdt"       , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_sbd_dhdt       ))
!     call check(nf90_def_var(ncid_monthly_avg, "sbd_horz"       , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_sbd_horz       ))
!     call check(nf90_def_var(ncid_monthly_avg, "sbd_vert"       , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_sbd_vert       ))
! 
!     call check(nf90_def_var(ncid_monthly_avg, "wt"             , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean, dimid_sw_ocean, 1/), varid_wt             ))
!     call check(nf90_def_var(ncid_monthly_avg, "prho"           , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean, dimid_st_ocean, 1/), varid_prho           ))
!     call check(nf90_def_var(ncid_monthly_avg, "taux"           , NF90_FLOAT, (/dimid_xu_ocean, dimid_yu_ocean,                 1/), varid_taux           ))
!     call check(nf90_def_var(ncid_monthly_avg, "tauy"           , NF90_FLOAT, (/dimid_xu_ocean, dimid_yu_ocean,                 1/), varid_tauy           ))
! 
!     call check(nf90_def_var(ncid_monthly_avg, "hflux_total"    , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_hflux_total    ))
!     call check(nf90_def_var(ncid_monthly_avg, "net_heating"    , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_net_heating    ))
!     call check(nf90_def_var(ncid_monthly_avg, "salt_flux_rstr" , NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_salt_flux_rstr ))
!     call check(nf90_def_var(ncid_monthly_avg, "salt_flux_total", NF90_FLOAT, (/dimid_xt_ocean, dimid_yt_ocean,                 1/), varid_salt_flux_total))
! 
! 
! 
!     call check(nf90_put_att(ncid_monthly_avg, NF90_GLOBAL, "description", "SODA-MOM5 monthly average analysis from 5-day average output."))    
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_xt_ocean, "long_name"     , "tcell longitude"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_xt_ocean, "units"         , "degrees_E"      ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_xt_ocean, "cartesian_axis", "X"              ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_yt_ocean, "long_name"     , "tcell latitude"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_yt_ocean, "units"         , "degrees_N"     ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_yt_ocean, "cartesian_axis", "Y"             ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_time, "long_name"     , "time"       ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_time, "units"         , "month"      ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_time, "cartesian_axis", "T"          ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_time, "calendar_type" , "JULIAN"     ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_time, "calendar"      , "JULIAN"     ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_time, "bounds"        , "time_bounds"))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "long_name"     , "tcell zstar depth"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "units"         , "meters"           ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "cartesian_axis", "Z"                ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "positive"      , "down"             ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_st_ocean, "edges"         , "st_edges_ocean"   ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_st_edges_ocean, "long_name"     , "tcell zstar depth edges"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_st_edges_ocean, "units"         , "meters"                 ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_st_edges_ocean, "cartesian_axis", "Z"                      ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_st_edges_ocean, "positive"      , "down"                   ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_xu_ocean, "long_name"     , "ucell longitude"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_xu_ocean, "units"         , "degrees_E"      ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_xu_ocean, "cartesian_axis", "X"              ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_yu_ocean, "long_name"     , "ucell latitude"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_yu_ocean, "units"         , "degrees_N"     ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_yu_ocean, "cartesian_axis", "Y"             ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "long_name"     , "ucell zstar depth"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "units"         , "meters"           ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "cartesian_axis", "Z"                ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "positive"      , "down"             ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sw_ocean, "edges"         , "sw_edges_ocean"   ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_sw_edges_ocean, "long_name"     , "ucell zstar depth edges"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sw_edges_ocean, "units"         , "meters"                 ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sw_edges_ocean, "cartesian_axis", "Z"                      ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sw_edges_ocean, "positive"      , "down"                   ))
! 
! !   call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "long_name"    , "tracer longitude" ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "units"        , "degrees_E"        ))
! ! ! call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "valid_range"  , "-281.f, 361.f"    ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "valid_range"  , (/-281., 361./)    ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "missing_value", 1.e+20             ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "_FillValue"   , 1.e+20             ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "cell_methods" , "time: point"      ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_geolon_t, "coordinates"  , "geolon_t geolat_t"))
! 
! !   call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "long_name"    , "tracer latitude"  ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "units"        , "degrees_N"        ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "valid_range"  , (/-91., 91./)      ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "missing_value", 1.e+20             ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "_FillValue"   , 1.e+20             ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "cell_methods" , "time: point"      ))
! !   call check(nf90_put_att(ncid_monthly_avg, varid_geolat_t, "coordinates"  , "geolon_t geolat_t"))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_temp, "long_name"    , "Potential temperature"          ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_temp, "units"        , "degrees C"                      ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_temp, "valid_range"  , (/-10., 500./)                   ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_temp, "missing_value", -1.e+20                          ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_temp, "_FillValue"   , -1.e+20                          ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_temp, "cell_methods" , "time: mean"                     ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_temp, "coordinates"  , "geolon_t geolat_t"              ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_temp, "standard_name", "sea_water_potential_temperature"))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt, "long_name"    , "Practical Salinity"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt, "units"        , "psu"               ))
!   ! call check(nf90_put_att(ncid_monthly_avg, varid_salt, "valid_range"  , "-10.f, 100.f"      ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt, "valid_range"  , (/-10., 100./)      ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt, "missing_value", -1.e+20             ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt, "_FillValue"   , -1.e+20             ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt, "cell_methods" , "time: mean"        ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt, "coordinates"  , "geolon_t geolat_t" ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt, "standard_name", "sea_water_salinity"))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_u, "long_name"    , "i-current"           ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_u, "units"        , "m/sec"               ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_u, "valid_range"  , (/-10., 10./)         ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_u, "missing_value", -1.e+20               ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_u, "_FillValue"   , -1.e+20               ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_u, "cell_methods" , "time: mean"          ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_u, "coordinates"  , "geolon_t geolat_t"   ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_u, "standard_name", "sea_water_x_velocity"))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_v, "long_name"    , "j-current"           ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_v, "units"        , "m/sec"               ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_v, "valid_range"  , (/-10., 10./)         ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_v, "missing_value", -1.e+20               ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_v, "_FillValue"   , -1.e+20               ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_v, "cell_methods" , "time: mean"          ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_v, "coordinates"  , "geolon_t geolat_t"   ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_v, "standard_name", "sea_water_y_velocity"))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "long_name"    , "effective sea level (eta_t + patm/(rho0*g)) on T cells"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "units"        , "meter"                         ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "valid_range"  , (/-1000., 1000./)               ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "missing_value", -1.e+20                         ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "_FillValue"   , -1.e+20                         ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "cell_methods" , "time: mean"                    ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "coordinates"  , "geolon_t geolat_t"             ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_ssh, "standard_name", "sea_surface_height_above_geoid"))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_mlt, "long_name"    , "mixed layer depth determined by temperature criteria"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mlt, "units"        , "m"                             ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mlt, "valid_range"  , (/0., 1000000./)                ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mlt, "missing_value", -1.e+20                         ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mlt, "_FillValue"   , -1.e+20                         ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mlt, "cell_methods" , "time: mean"                    ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mlt, "coordinates"  , "geolon_t geolat_t"             ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_mlp, "long_name"    , "Depth of potential density mixed layer"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mlp, "units"        , "m"                             ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mlp, "valid_range"  , (/-1000000., 1000000./)         ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mlp, "missing_value", -1.e+20                         ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mlp, "_FillValue"   , -1.e+20                         ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mlp, "cell_methods" , "time: mean"                    ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mlp, "coordinates"  , "geolon_t geolat_t"             ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_mls, "long_name"    , "mixed layer depth determined by salinity criteria"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mls, "units"        , "m"                                                ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mls, "valid_range"  , (/0., 1000000./)                                   ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mls, "missing_value", -1.e+20                                            ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mls, "_FillValue"   , -1.e+20                                            ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mls, "cell_methods" , "time: mean"                                       ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_mls, "coordinates"  , "geolon_t geolat_t"                                ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_anompb, "long_name"    , "T-cell bottom pressure - rho0*grav*ht"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_anompb, "units"        , "dbar"                                 ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_anompb, "valid_range"  , (/-1000000., 1000000./)                ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_anompb, "missing_value", -1.e+20                                ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_anompb, "_FillValue"   , -1.e+20                                ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_anompb, "cell_methods" , "time: mean"                           ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_anompb, "coordinates"  , "geolon_t geolat_t"                    ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd, "long_name"    , "rate of mass transferred below the mixed layer base"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd, "units"        , "kg/sec"                                             ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd, "valid_range"  , (/-1.e+20, 1.e+20/)                                  ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd, "missing_value", -1.e+20                                              ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd, "_FillValue"   , -1.e+20                                              ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd, "cell_methods" , "time: mean"                                         ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd, "coordinates"  , "geolon_t geolat_t"                                  ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_mld, "long_name"    , "mixed layer depth used for subduction diagnostics"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_mld, "units"        , "m"                                                ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_mld, "valid_range"  , (/-100., 1.e+20/)                                  ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_mld, "missing_value", -1.e+20                                            ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_mld, "_FillValue"   , -1.e+20                                            ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_mld, "cell_methods" , "time: mean"                                       ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_mld, "coordinates"  , "geolon_t geolat_t"                                ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_dhdt, "long_name"    , "rate of mass transferred below the mixed layer base due to dh/dt"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_dhdt, "units"        , "kg/sec"                                                          ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_dhdt, "valid_range"  , (/-1.e+20, 1.e+20/)                                               ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_dhdt, "missing_value", -1.e+20                                                           ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_dhdt, "_FillValue"   , -1.e+20                                                           ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_dhdt, "cell_methods" , "time: mean"                                                      ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_dhdt, "coordinates"  , "geolon_t geolat_t"                                               ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_horz, "long_name"    , "rate of mass transferred below the mixed layer base due to horz advect"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_horz, "units"        , "kg/sec"                                                                ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_horz, "valid_range"  , (/-1.e+20, 1.e+20/)                                                     ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_horz, "missing_value", -1.e+20                                                                 ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_horz, "_FillValue"   , -1.e+20                                                                 ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_horz, "cell_methods" , "time: mean"                                                            ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_horz, "coordinates"  , "geolon_t geolat_t"                                                     ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_vert, "long_name"    , "rate of mass transferred below the mixed layer base due to vert advect"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_vert, "units"        , "kg/sec"                                                                ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_vert, "valid_range"  , (/-1.e+20, 1.e+20/)                                                     ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_vert, "missing_value", -1.e+20                                                                 ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_vert, "_FillValue"   , -1.e+20                                                                 ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_vert, "cell_methods" , "time: mean"                                                            ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_sbd_vert, "coordinates"  , "geolon_t geolat_t"                                                     ))
! 
! 
!   ! call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "long_name"    , "bottom pressure on T cells [Boussinesq (volume conserving) model]"))
!   ! call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "units"        , "dbar"                           ))
!   ! call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "valid_range"  , (/-1000000., 1000000./)          ))
!   ! call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "missing_value", -1.e+20                          ))
!   ! call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "_FillValue"   , -1.e+20                          ))
!   ! call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "cell_methods" , "time: mean"                     ))
!   ! call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "coordinates"  , "geolon_t geolat_t"              ))
!   ! call check(nf90_put_att(ncid_monthly_avg, varid_pbot, "standard_name", "sea_water_pressure_at_sea_floor"))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_wt, "long_name"    , "dia-surface velocity T-points"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_wt, "units"        , "m/sec"                        ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_wt, "valid_range"  , (/-100000., 100000./)          ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_wt, "missing_value", -1.e+20                        ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_wt, "_FillValue"   , -1.e+20                        ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_wt, "cell_methods" , "time: mean"                   ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_wt, "coordinates"  , "geolon_t geolat_t"            ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_prho, "long_name"    , "potential density referenced to 0 dbar"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_prho, "units"        , "kg/m^3"                                ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_prho, "valid_range"  , (/-10., 100000./)                       ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_prho, "missing_value", -1.e+20                                 ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_prho, "_FillValue"   , -1.e+20                                 ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_prho, "cell_methods" , "time: mean"                            ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_prho, "coordinates"  , "geolon_t geolat_t"                     ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_prho, "standard_name", "sea_water_potential_density"           ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_taux, "long_name"    , "i-directed wind stress forcing u-velocity"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_taux, "units"        , "N/m^2"                                    ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_taux, "valid_range"  , (/-10., 10./)                              ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_taux, "missing_value", -1.e+20                                    ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_taux, "_FillValue"   , -1.e+20                                    ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_taux, "cell_methods" , "time: mean"                               ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_taux, "coordinates"  , "geolon_t geolat_t"                        ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_taux, "standard_name", "surface_downward_x_stress"                ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_tauy, "long_name"    , "j-directed wind stress forcing v-velocity"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_tauy, "units"        , "N/m^2"                                    ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_tauy, "valid_range"  , (/-10., 10./)                              ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_tauy, "missing_value", -1.e+20                                    ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_tauy, "_FillValue"   , -1.e+20                                    ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_tauy, "cell_methods" , "time: mean"                               ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_tauy, "coordinates"  , "geolon_t geolat_t"                        ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_tauy, "standard_name", "surface_downward_y_stress"                ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_hflux_total, "long_name"    , "surface heat flux from coupler plus restore (omits mass transfer heating"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_hflux_total, "units"        , "Watts/m^2"                                                               ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_hflux_total, "valid_range"  , (/-10000., 10000./)                                                       ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_hflux_total, "missing_value", -1.e+20                                                                   ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_hflux_total, "_FillValue"   , -1.e+20                                                                   ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_hflux_total, "cell_methods" , "time: mean"                                                              ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_hflux_total, "coordinates"  , "geolon_t geolat_t"                                                       ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_net_heating, "long_name"    , "surface ocean heat flux coming through coupler and mass transfer"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_net_heating, "units"        , "Watts/m^2"                                                       ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_net_heating, "valid_range"  , (/-10000., 10000./)                                               ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_net_heating, "missing_value", -1.e+20                                                           ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_net_heating, "_FillValue"   , -1.e+20                                                           ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_net_heating, "cell_methods" , "time: mean"                                                      ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_net_heating, "coordinates"  , "geolon_t geolat_t"                                               ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt_flux_rstr, "long_name"    , "sfc_salt_flux_restore: flux from restoring term"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt_flux_rstr, "units"        , "kg/(m^2*sec"                                    ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt_flux_rstr, "valid_range"  , (/-10000., 10000./)                              ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt_flux_rstr, "missing_value", -1.e+20                                          ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt_flux_rstr, "_FillValue"   , -1.e+20                                          ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt_flux_rstr, "cell_methods" , "time: mean"                                     ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt_flux_rstr, "coordinates"  , "geolon_t geolat_t"                              ))
! 
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt_flux_total, "long_name"    , "sfc_salt_flux_total"))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt_flux_total, "units"        , "kg/(m^2*sec"                                    ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt_flux_total, "valid_range"  , (/-10000., 10000./)                              ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt_flux_total, "missing_value", -1.e+20                                          ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt_flux_total, "_FillValue"   , -1.e+20                                          ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt_flux_total, "cell_methods" , "time: mean"                                     ))
!     call check(nf90_put_att(ncid_monthly_avg, varid_salt_flux_total, "coordinates"  , "geolon_t geolat_t"                              ))
! 
! 
!     call check(nf90_enddef(ncid_monthly_avg))
! 
! 
!     ! write out all vars' data
!     call check(nf90_put_var(ncid_monthly_avg, varid_time          , 0                 ))
! 
!     call check(nf90_put_var(ncid_monthly_avg, varid_xt_ocean      , var_xt_ocean      ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_yt_ocean      , var_yt_ocean      ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_st_ocean      , var_st_ocean      ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_st_edges_ocean, var_st_edges_ocean))
!     call check(nf90_put_var(ncid_monthly_avg, varid_xu_ocean      , var_xu_ocean      ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_yu_ocean      , var_yu_ocean      ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_sw_ocean      , var_sw_ocean      ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_sw_edges_ocean, var_sw_edges_ocean))
! 
!   ! call check(nf90_put_var(ncid_monthly_avg, varid_geolon_t, var_geolon_t))
!   ! call check(nf90_put_var(ncid_monthly_avg, varid_geolat_t, var_geolat_t))
! 
!     call check(nf90_put_var(ncid_monthly_avg, varid_temp     , avg_temp     ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_salt     , avg_salt     ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_u        , avg_u        ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_v        , avg_v        ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_ssh      , avg_ssh      ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_mlt      , avg_mlt      ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_mlp      , avg_mlp      ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_mls      , avg_mls      ))
! 
!     call check(nf90_put_var(ncid_monthly_avg, varid_anompb   , avg_anompb   ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_sbd      , avg_sbd      ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_sbd_mld  , avg_sbd_mld  ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_sbd_dhdt , avg_sbd_dhdt ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_sbd_horz , avg_sbd_horz ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_sbd_vert , avg_sbd_vert ))
! 
!     call check(nf90_put_var(ncid_monthly_avg, varid_wt       , avg_wt       ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_prho     , avg_prho     ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_taux     , avg_taux     ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_tauy     , avg_tauy     ))
! 
!     call check(nf90_put_var(ncid_monthly_avg, varid_hflux_total    , avg_hflux_total    ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_net_heating    , avg_net_heating    ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_salt_flux_rstr , avg_salt_flux_rstr ))
!     call check(nf90_put_var(ncid_monthly_avg, varid_salt_flux_total, avg_salt_flux_total))
! 
! 
!     call check(nf90_close(ncid_monthly_avg))
! 
!     call check(nf90_close(ncid_dim))
! 
!     write(*, *) "*** SUCCESS writing monthly average results for " // str_yr // str_mon // " ***"
! 
!     return
! end subroutine write_out_5dayAvg2monthlyAvg_SODA_MOM_339


! compute julian day from SODA code
function jday(mon, iday, iyr)
    dimension dpm(12)
    data dpm /31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0/
    dpm(2) = 28.0
    if(mod(real(iyr),4.) .eq. 0.) dpm(2) = 29.0

    iyrs     = iyr-1970
    days     = 587.0+real(iyrs)*365.0
    num_leap = floor(real(iyrs+1)/4.)
    days     = days + real(num_leap)

    ! now sum up for days this year
    sum = 0.
    if(mon .gt. 1) then
        do l = 1, mon - 1
            sum = sum + dpm(l)
        end do
        days = days + sum
    end if
    
    jday = int(days) + iday
        
    return
end function jday


subroutine cal_temp_fcs_obs_diff_dist(yr, mon, day, jul_day, temp_tmp)
    implicit none

    integer,                    intent(in ) :: yr, mon, day
    integer,                    intent(out) :: jul_day
    integer, dimension(71, 34), intent(out) :: temp_tmp

    ! function
    integer :: jday

    ! local vars
    integer :: i, i_v_lev, i_int
    integer :: io_temp_var_wrk = 10
    character(len=* ), parameter :: DIR_SODA = "/glade/scratch/lgchen/project/MOM_run/SODA_1979_2014/test1/SODA/SODA_1979_1988"
    character(len=21) :: fn_temp
    character(len=4 ) :: str_yr
    character(len=2 ) :: str_mon
    character(len=2 ) :: str_day

    ! to read one diff record with all vertical levels
    integer :: jdatg, ilngg, jlatg, ncnt
    real    :: xlatg, xlngg, tempt(34)
    character(len=3) :: typeg
  
    write(unit=str_yr , fmt="(I4.4)") yr
    write(unit=str_mon, fmt="(I2.2)") mon
    write(unit=str_day, fmt="(I2.2)") day

    jul_day = jday(mon, day, yr)

    fn_temp = str_yr // str_mon // str_day // ".temp_var.wrk"

    open(unit=io_temp_var_wrk, file=trim(DIR_SODA) // "/" // trim(fn_temp), status="old")

    temp_tmp = 0
    do while (1)  ! not end of the file
        read(io_temp_var_wrk, 200, end=111)  &
            jdatg, xlatg, xlngg, ilngg, jlatg, typeg, ncnt, (tempt(i), i=1, ncnt)
        
        if (jul_day-2 .le. jdatg .and. jdatg .le. jul_day+2) then
          ! write(*, *) fn_temp, " jdatg = ", jdatg 

            do i_v_lev = 1, ncnt
                do i_int = 1, 71
                    if (-7.1+0.2*(i_int-1) .le. tempt(i_v_lev) .and. tempt(i_v_lev) .le. -6.9+0.2*(i_int-1)) then
                        temp_tmp(i_int, i_v_lev) = temp_tmp(i_int, i_v_lev) + 1

                        exit 
                    end if
                end do
            end do
        end if

        if (jul_day + 3 .le. jdatg) exit
    end do
111 continue

200 format(i6, 2f8.2, 2i4, a6, i4, 5f8.3/,3(10f8.3,/))

    return
end subroutine cal_temp_fcs_obs_diff_dist


subroutine cal_salt_fcs_obs_diff_dist(yr, mon, day, jul_day, salt_tmp)
    implicit none

    integer,                    intent(in ) :: yr, mon, day
    integer,                    intent(out) :: jul_day
    integer, dimension(21, 34), intent(out) :: salt_tmp

    ! function
    integer :: jday

    ! local vars
    integer :: i, i_v_lev, i_int
    integer :: io_salt_var_wrk = 10
    character(len=* ), parameter :: DIR_SODA = "/glade/scratch/lgchen/project/MOM_run/SODA_1979_2014/test1/SODA/SODA_1979_1988"
    character(len=21) :: fn_salt
    character(len=4 ) :: str_yr
    character(len=2 ) :: str_mon
    character(len=2 ) :: str_day

    ! to read one diff record with all vertical levels
    integer :: jdatg, ilngg, jlatg, ncnt
    real    :: xlatg, xlngg, saltt(34)
    character(len=3) :: typeg
  
    write(unit=str_yr , fmt="(I4.4)") yr
    write(unit=str_mon, fmt="(I2.2)") mon
    write(unit=str_day, fmt="(I2.2)") day

    jul_day = jday(mon, day, yr)

    fn_salt = str_yr // str_mon // str_day // ".salt_var.wrk"

    open(unit=io_salt_var_wrk, file=trim(DIR_SODA) // "/" // trim(fn_salt), status="old")

    salt_tmp = 0
    do while (1)  ! not end of the file
        read(io_salt_var_wrk, 200, end=111)  &
            jdatg, xlatg, xlngg, ilngg, jlatg, typeg, ncnt, (saltt(i), i=1, ncnt)
        
        if (jul_day-2 .le. jdatg .and. jdatg .le. jul_day+2) then
          ! write(*, *) fn_salt, " jdatg = ", jdatg 

            do i_v_lev = 1, ncnt
                do i_int = 1, 21
                    if (-2.1+0.2*(i_int-1) .le. saltt(i_v_lev) .and. saltt(i_v_lev) .le. -1.9+0.2*(i_int-1)) then
                        salt_tmp(i_int, i_v_lev) = salt_tmp(i_int, i_v_lev) + 1

                        exit 
                    end if
                end do
            end do
        end if

        if (jul_day + 3 .le. jdatg) exit
    end do
111 continue

200 format(i6, 2f8.2, 2i4, a6, i4, 5f8.3/,3(10f8.3,/))

    return
end subroutine cal_salt_fcs_obs_diff_dist


subroutine sst_noaa2soda_1hr(ni, nj, time, lon_noaa, lat_noaa, sst_noaa, sses_bias, l2p_flags  &
    , quality_level, sst_soda, n_sst_soda)
    implicit none

    integer        ,                      intent(in ) :: ni, nj, time
    real   (kind=4), dimension(ni , nj ), intent(in ) :: lon_noaa, lat_noaa
    real   (kind=4), dimension(ni , nj ), intent(in ) :: sst_noaa, sses_bias
    integer(kind=2), dimension(ni , nj ), intent(in ) :: l2p_flags
    integer(kind=1), dimension(ni , nj ), intent(in ) :: quality_level
  ! integer        , dimension(ni , nj ), intent(in ) :: l2p_flags
  ! integer        , dimension(ni , nj ), intent(in ) :: quality_level
    real   (kind=4), dimension(360, 180), intent(out) :: sst_soda
    integer(kind=4), dimension(360, 180), intent(out) :: n_sst_soda


    ! local vars
    real(kind=4) :: lon_soda(360), lat_soda(180)
    real(kind=4) :: lon_noaa2soda

    integer :: i, j, i_idx, j_idx


    ! set up lon & lat for SODA
  ! write(*, *) "Entering subroutine sst_noaa2soda_1hr, setting up lon and lat for SODA"
    do i = 1, 360
        lon_soda(i) = float(i)
    end do

    do i = 1, 180
        lat_soda(i) = -89.5 + (i - 1)*1.0
    end do


    ! extract valid obs sst data to soda grid
    sst_soda   = 0.0
    n_sst_soda = 0

    do j = 1, nj; do i = 1, ni
      ! write(*, *) "j = ", j, ", i = ", i
        if (quality_level(i, j) .eq. 5 .and. (.not. btest(l2p_flags(i, j), 10))  &
            .and. 200.0 <= sst_noaa(i, j) .and. sst_noaa(i, j) <= 350.0) then
          ! write(*, *) "valid, j = ", j, ", i = ", i

            j_idx = floor(lat_noaa(i, j)) + 91
            if (lon_noaa(i, j) < 0.0) then
                i_idx = floor(lon_noaa(i, j) + 360.0 + 0.5)
            else
                i_idx = floor(lon_noaa(i, j) + 0.5)
            end if
            if (i_idx < 1) then
                i_idx = 1
            end if

            sst_soda  (i_idx, j_idx) = sst_soda(i_idx, j_idx) + sst_noaa(i, j) - sses_bias(i, j)
            n_sst_soda(i_idx, j_idx) = n_sst_soda(i_idx, j_idx) + 1
        end if
    end do; end do

  ! write(*, *) "returning from subroutine sst_noaa2soda_1hr"

    return
end subroutine sst_noaa2soda_1hr


subroutine sst_noaa2soda_10min(lon, lat, time, lon_noaa, lat_noaa, sst_noaa, sses_bias, l2p_flags  &
    , quality_level, sst_soda, n_sst_soda)
    implicit none

    integer        ,                      intent(in ) :: lon, lat, time
    real   (kind=4), dimension(lon     ), intent(in ) :: lon_noaa
    real   (kind=4), dimension(     lat), intent(in ) :: lat_noaa
    real   (kind=4), dimension(lon, lat), intent(in ) :: sst_noaa, sses_bias
    integer(kind=2), dimension(lon, lat), intent(in ) :: l2p_flags
    integer(kind=1), dimension(lon, lat), intent(in ) :: quality_level
    real   (kind=4), dimension(360, 180), intent(out) :: sst_soda
    integer(kind=4), dimension(360, 180), intent(out) :: n_sst_soda


    ! local vars
    real(kind=4) :: lon_soda(360), lat_soda(180)
    real(kind=4) :: lon_noaa2soda

    integer :: i, j, i_idx, j_idx


    ! set up lon & lat for SODA
  ! write(*, *) "Entering subroutine sst_noaa2soda_1hr, setting up lon and lat for SODA"
    do i = 1, 360
        lon_soda(i) = float(i)
    end do

    do i = 1, 180
        lat_soda(i) = -89.5 + (i - 1)*1.0
    end do


    ! extract valid obs sst data to soda grid
    sst_soda   = 0.0
    n_sst_soda = 0

    do j = 1, lat; do i = 1, lon
        if (quality_level(i, j) .eq. 5 .and. (.not. btest(l2p_flags(i, j), 10))  &
            .and. 200.0 <= sst_noaa(i, j) .and. sst_noaa(i, j) <= 350.0) then
          ! write(*, *) "valid, j = ", j, ", i = ", i

          ! j_idx = floor(lat_noaa(i, j)) + 91
            j_idx = floor(lat_noaa(j)) + 91
            if (lon_noaa(i) < 0.0) then
                i_idx = floor(lon_noaa(i) + 360.0 + 0.5)
            else
                i_idx = floor(lon_noaa(i) + 0.5)
            end if
            if (i_idx < 1) then
                i_idx = 1
            end if

            sst_soda  (i_idx, j_idx) = sst_soda(i_idx, j_idx) + sst_noaa(i, j) - sses_bias(i, j)
            n_sst_soda(i_idx, j_idx) = n_sst_soda(i_idx, j_idx) + 1
        end if
    end do; end do

  ! write(*, *) "returning from subroutine sst_noaa2soda_1hr"

    return
end subroutine sst_noaa2soda_10min


subroutine write_sst_noaa2soda_5day_avg(ndat, sst_soda)
    implicit none

    integer     ,                      intent(in) :: ndat
    real(kind=4), dimension(360, 180), intent(in) :: sst_soda


    ! local vars
    real(kind=4) :: lon_soda(360), lat_soda(180)
    character(len=256) :: fn_sst
    integer :: i, j


    ! set up lon & lat for SODA
    write(*, *) "Entering subroutine write_sst_noaa2soda_5day_avg ..."

    do i = 1, 360
        lon_soda(i) = float(i)
    end do
    do i = 1, 180
        lat_soda(i) = -89.5 + (i - 1)*1.0
    end do

    ! open 
    fn_sst = "/gpfs/fs1/p/univ/umcp0009/OBS4SODA/nst_sst_tmp.bin"
    open(30, file = fn_sst, status = 'unknown')

    ! write out sst 5day_avg on soda grid
    write(*, *) "Writing out sst 5day_avg on soda grid ..."
    do i = 1, 360; do j = 1, 180
        if (273.15 - 2.0 <= sst_soda(i, j) .and. sst_soda(i, j) <= 273.15 + 32.0) then
            write(30, 111) ndat, lat_soda(j), lon_soda(i), sst_soda(i, j) - 273.15
        end if
    end do; end do

111 format(1x, i5, 1x, f5.1, 1x, f5.1, 2x, f5.2)

    return
end subroutine write_sst_noaa2soda_5day_avg




subroutine write_out_monthlyAvg_nc_SODA(curr_yr, curr_mon, i_var, str_var, var)
    use netcdf
    implicit none

    integer, parameter :: N_MONTH = 12, NY = 330, NX = 720, NZ = 50

    integer,                        intent(in) :: curr_yr, curr_mon, i_var
  ! character(len=128),             intent(in) :: str_var
    character*(*),             intent(in) :: str_var
    real   , dimension(NX, NY, NZ), intent(in) :: var

    ! local vars
  ! character(len = *), parameter :: DIR_DATA_NC = "/glade/p/umcp0006/lgchen/project/MOM_run/SODA_3.3_1/monthly/nc_ligang"
  ! character(len = *), parameter :: DIR_DATA_NC = "/glade/p/umcp0006/lgchen/project/MOM_run/SODA_3.5_1/monthly/nc_ligang"
    character(len = *), parameter :: DIR_DATA_NC = "/glade/scratch/lgchen/data/SODA_3.4_monthly_regrided_nc"

    character(len=4  ) :: str_yr
    character(len=2  ) :: str_mon
    character(len=256) :: fn_nc  

    integer :: ncid_dim, varid


    write(unit=str_yr , fmt="(I4.4)") curr_yr
    write(unit=str_mon, fmt="(I2.2)") curr_mon

  ! fn_nc = "soda_3.3_mn_ocean_" // str_yr // ".nc"
  ! fn_nc = "soda_3.5_mn_ocean_" // str_yr // ".nc"
    fn_nc = "soda_3.4_mn_ocean_" // str_yr // ".nc"


    ! open and read the dim file data
    call check(nf90_open(trim(DIR_DATA_NC) // "/" // trim(fn_nc), nf90_write, ncid_dim))
    
    call check(nf90_inq_varid(ncid_dim, trim(str_var), varid))

    ! write out all vars' data
    if (i_var <= 11) then
        call check(nf90_put_var(ncid_dim, varid, var(:, :, 1), start = (/1, 1, curr_mon/)  &
            , count = (/NX, NY, 1/)))
    else if (12 <= i_var) then
        call check(nf90_put_var(ncid_dim, varid, var, start = (/1, 1, 1, curr_mon/)  &
            , count = (/NX, NY, NZ, 1/)))
    end if

    call check(nf90_close(ncid_dim))

    write(*, *) "*** SUCCESS writing monthly average results for " // str_yr // str_mon // " " // trim(str_var) // " ***"

    return
end subroutine write_out_monthlyAvg_nc_SODA




! subroutine write_out_monthlyAvg_nc_SODA_339(curr_yr, curr_mon, i_var, str_var, var)
!     use netcdf
!     implicit none
! 
!     integer, parameter :: N_MONTH = 12, NY = 330, NX = 720, NZ = 50
! 
!     integer,                        intent(in) :: curr_yr, curr_mon, i_var
!   ! character(len=128),             intent(in) :: str_var
!     character*(*),             intent(in) :: str_var
!     real   , dimension(NX, NY, NZ), intent(in) :: var
! 
!     ! local vars
!   ! character(len = *), parameter :: DIR_DATA_NC = "/glade/p/umcp0006/lgchen/project/MOM_run/SODA_3.3_1/monthly/nc_ligang"
!   ! character(len = *), parameter :: DIR_DATA_NC = "/glade/p/umcp0006/lgchen/project/MOM_run/SODA_3.5_1/monthly/nc_ligang"
!   ! character(len = *), parameter :: DIR_DATA_NC = "/glade/scratch/lgchen/data/SODA_3.4_monthly_regrided_nc"
!     character(len = *), parameter :: DIR_DATA_NC = "/glade/scratch/lgchen/project/SODA_3.3.9/month_avg/netcdf"
! 
!     character(len=4  ) :: str_yr
!     character(len=2  ) :: str_mon
!     character(len=256) :: fn_nc  
! 
!     integer :: ncid_dim, varid
! 
! 
!     write(unit=str_yr , fmt="(I4.4)") curr_yr
!     write(unit=str_mon, fmt="(I2.2)") curr_mon
! 
!   ! fn_nc = "soda_3.3_mn_ocean_" // str_yr // ".nc"
!   ! fn_nc = "soda_3.5_mn_ocean_" // str_yr // ".nc"
!   ! fn_nc = "soda_3.4_mn_ocean_" // str_yr // ".nc"
!     fn_nc = "SODA_3.3.0_mn_ocean_" // str_yr // ".nc"
! 
! 
!     ! open and read the dim file data
!     call check(nf90_open(trim(DIR_DATA_NC) // "/" // trim(fn_nc), nf90_write, ncid_dim))
!     
!     call check(nf90_inq_varid(ncid_dim, trim(str_var), varid))
! 
!     ! write out all vars' data
!     if (i_var <= 8) then
!         call check(nf90_put_var(ncid_dim, varid, var(:, :, 1), start = (/1, 1, curr_mon/)  &
!             , count = (/NX, NY, 1/)))
!     else if (9 <= i_var) then
!         call check(nf90_put_var(ncid_dim, varid, var, start = (/1, 1, 1, curr_mon/)  &
!             , count = (/NX, NY, NZ, 1/)))
!     end if
! 
!     call check(nf90_close(ncid_dim))
! 
!     write(*, *) "*** SUCCESS writing monthly average results for " // str_yr // str_mon // " " // trim(str_var) // " ***"
! 
!     return
! end subroutine write_out_monthlyAvg_nc_SODA_339




subroutine write_out_monthlyAvg_nc_SODA_339(curr_yr, curr_mon, is_var_3d, str_var, var)
    use netcdf
    implicit none

    integer, parameter :: N_MONTH = 12, NY = 330, NX = 720, NZ = 50

    integer,                        intent(in) :: curr_yr, curr_mon, is_var_3d
  ! character(len=128),             intent(in) :: str_var
    character*(*),             intent(in) :: str_var
    real   , dimension(NX, NY, NZ), intent(in) :: var

    ! local vars
  ! character(len = *), parameter :: DIR_DATA_NC = "/glade/p/umcp0006/lgchen/project/MOM_run/SODA_3.3_1/monthly/nc_ligang"
  ! character(len = *), parameter :: DIR_DATA_NC = "/glade/p/umcp0006/lgchen/project/MOM_run/SODA_3.5_1/monthly/nc_ligang"
  ! character(len = *), parameter :: DIR_DATA_NC = "/glade/scratch/lgchen/data/SODA_3.4_monthly_regrided_nc"
  ! character(len = *), parameter :: DIR_DATA_NC = "/glade/scratch/lgchen/project/SODA_3.3.9/month_avg/netcdf"
    character(len = *), parameter :: DIR_DATA_NC = "/glade2/scratch2/lgchen/project/soda3.4.0_mn_bin2nc_fromAOSC/nc"

    character(len=4  ) :: str_yr
    character(len=2  ) :: str_mon
    character(len=256) :: fn_nc  

    integer :: ncid_dim, varid


    write(unit=str_yr , fmt="(I4.4)") curr_yr
    write(unit=str_mon, fmt="(I2.2)") curr_mon

  ! fn_nc = "soda_3.3_mn_ocean_" // str_yr // ".nc"
  ! fn_nc = "soda_3.5_mn_ocean_" // str_yr // ".nc"
  ! fn_nc = "soda_3.4_mn_ocean_" // str_yr // ".nc"
  ! fn_nc = "soda3.3.1_mn_ocean_" // str_yr // ".nc"
    fn_nc = "soda3.4.0_mn_ocean_reg_" // str_yr // ".nc"


    ! open and read the dim file data
    call check(nf90_open(trim(DIR_DATA_NC) // "/" // trim(fn_nc), nf90_write, ncid_dim))
    
    call check(nf90_inq_varid(ncid_dim, trim(str_var), varid))

    ! write out all vars' data
    if (is_var_3d == 0) then  ! 2d
        call check(nf90_put_var(ncid_dim, varid, var(:, :, 1), start = (/1, 1, curr_mon/), count = (/NX, NY, 1/)))
    else  ! 3d
        call check(nf90_put_var(ncid_dim, varid, var, start = (/1, 1, 1, curr_mon/), count = (/NX, NY, NZ, 1/)))
    end if

    call check(nf90_close(ncid_dim))

    write(*, *) "*** SUCCESS writing monthly average results for " // str_yr // str_mon // " " // trim(str_var) // " ***"

    return
end subroutine write_out_monthlyAvg_nc_SODA_339




subroutine write_out_monthlyAvg_isopycn_nc_SODA_339(curr_yr, curr_mon, str_var, var)
    use netcdf
    implicit none

    integer, parameter :: N_MONTH = 12, NY = 330, NX = 720, NZ = 16

    integer,                        intent(in) :: curr_yr, curr_mon
  ! character(len=128),             intent(in) :: str_var
    character*(*),             intent(in) :: str_var
    real   , dimension(NX, NY, NZ), intent(in) :: var

    ! local vars
    character(len = *), parameter :: DIR_DATA_NC = "/glade/scratch/lgchen/project/SODA_3.3.9/month_avg/netcdf"
  ! character(len = *), parameter :: DIR_DATA_NC = "/glade/scratch/lgchen/project/SODA_3.3.0/month_avg/netcdf/isopycn"

    character(len=4  ) :: str_yr
    character(len=2  ) :: str_mon
    character(len=256) :: fn_nc  

    integer :: ncid_dim, varid


    write(unit=str_yr , fmt="(I4.4)") curr_yr
    write(unit=str_mon, fmt="(I2.2)") curr_mon

    fn_nc = "soda3.3.1_mn_isopycn_reg_" // str_yr // ".nc"
  ! fn_nc = "soda3.3.0_mn_isopycn_reg_" // str_yr // ".nc"


    ! open and read the dim file data
    call check(nf90_open(trim(DIR_DATA_NC) // "/" // trim(fn_nc), nf90_write, ncid_dim))
    
    call check(nf90_inq_varid(ncid_dim, trim(str_var), varid))

    ! write out all vars' data
    call check(nf90_put_var(ncid_dim, varid, var, start = (/1, 1, 1, curr_mon/), count = (/NX, NY, NZ, 1/)))

    call check(nf90_close(ncid_dim))

    write(*, *) "*** SUCCESS writing monthly average results for " // str_yr // str_mon // " " // trim(str_var) // " ***"

    return
end subroutine write_out_monthlyAvg_isopycn_nc_SODA_339




subroutine write_out_monthlyAvg_ice_nc_SODA_339(curr_yr, curr_mon, is_var_3d, str_var, var)
    use netcdf
    implicit none

    integer, parameter :: N_MONTH = 12, NY = 330, NX = 720, NZ = 5

    integer,                        intent(in) :: curr_yr, curr_mon, is_var_3d
  ! character(len=128),             intent(in) :: str_var
    character*(*),             intent(in) :: str_var
    real   , dimension(NX, NY, NZ), intent(in) :: var

    ! local vars
    character(len = *), parameter :: DIR_DATA_NC = "/glade/scratch/lgchen/project/SODA_3.3.9/month_avg/netcdf/ice"

    character(len=4  ) :: str_yr
    character(len=2  ) :: str_mon
    character(len=256) :: fn_nc  

    integer :: ncid_dim, varid


    write(unit=str_yr , fmt="(I4.4)") curr_yr
    write(unit=str_mon, fmt="(I2.2)") curr_mon

    fn_nc = "soda3.3.1_mn_ice_reg_" // str_yr // ".nc"


    ! open and read the dim file data
    call check(nf90_open(trim(DIR_DATA_NC) // "/" // trim(fn_nc), nf90_write, ncid_dim))
    
    call check(nf90_inq_varid(ncid_dim, trim(str_var), varid))

    ! write out all vars' data
    if (is_var_3d == 0) then  ! 2d
        call check(nf90_put_var(ncid_dim, varid, var(:, :, 1), start = (/1, 1, curr_mon/), count = (/NX, NY, 1/)))
    else  ! 3d
        call check(nf90_put_var(ncid_dim, varid, var, start = (/1, 1, 1, curr_mon/), count = (/NX, NY, NZ, 1/)))
    end if

    call check(nf90_close(ncid_dim))

    write(*, *) "*** SUCCESS writing monthly average results for " // str_yr // str_mon // " " // trim(str_var) // " ***"

    return
end subroutine write_out_monthlyAvg_ice_nc_SODA_339




subroutine write_out_monthlyAvg_nc_SODA_330(curr_yr, curr_mon, is_var_3d, str_var, var)
    use netcdf
    implicit none

    integer, parameter :: N_MONTH = 12, NY = 330, NX = 720, NZ = 50

    integer,                        intent(in) :: curr_yr, curr_mon, is_var_3d
  ! character(len=128),             intent(in) :: str_var
    character*(*),             intent(in) :: str_var
    real   , dimension(NX, NY, NZ), intent(in) :: var

    ! local vars
    character(len = *), parameter :: DIR_DATA_NC = "/glade/scratch/lgchen/project/SODA_3.3.0/month_avg/netcdf/ocean"

    character(len=4  ) :: str_yr
    character(len=2  ) :: str_mon
    character(len=256) :: fn_nc  

    integer :: ncid_dim, varid


    write(unit=str_yr , fmt="(I4.4)") curr_yr
    write(unit=str_mon, fmt="(I2.2)") curr_mon

    fn_nc = "soda3.3.0_mn_ocean_reg_" // str_yr // ".nc"


    ! open and read the dim file data
    call check(nf90_open(trim(DIR_DATA_NC) // "/" // trim(fn_nc), nf90_write, ncid_dim))
    
    call check(nf90_inq_varid(ncid_dim, trim(str_var), varid))

    ! write out all vars' data
    if (is_var_3d == 0) then  ! 2d
        call check(nf90_put_var(ncid_dim, varid, var(:, :, 1), start = (/1, 1, curr_mon/), count = (/NX, NY, 1/)))
    else  ! 3d
        call check(nf90_put_var(ncid_dim, varid, var, start = (/1, 1, 1, curr_mon/), count = (/NX, NY, NZ, 1/)))
    end if

    call check(nf90_close(ncid_dim))

    write(*, *) "*** SUCCESS writing monthly average results for " // str_yr // str_mon // " " // trim(str_var) // " ***"

    return
end subroutine write_out_monthlyAvg_nc_SODA_330




subroutine write_out_monthlyAvg_ice_nc_SODA_330(curr_yr, curr_mon, is_var_3d, str_var, var)
    use netcdf
    implicit none

    integer, parameter :: N_MONTH = 12, NY = 330, NX = 720, NZ = 5

    integer,                        intent(in) :: curr_yr, curr_mon, is_var_3d
  ! character(len=128),             intent(in) :: str_var
    character*(*),             intent(in) :: str_var
    real   , dimension(NX, NY, NZ), intent(in) :: var

    ! local vars
    character(len = *), parameter :: DIR_DATA_NC = "/glade/scratch/lgchen/project/SODA_3.3.0/month_avg/netcdf/ice"

    character(len=4  ) :: str_yr
    character(len=2  ) :: str_mon
    character(len=256) :: fn_nc  

    integer :: ncid_dim, varid


    write(unit=str_yr , fmt="(I4.4)") curr_yr
    write(unit=str_mon, fmt="(I2.2)") curr_mon

    fn_nc = "soda3.3.0_mn_ice_reg_" // str_yr // ".nc"


    ! open and read the dim file data
    call check(nf90_open(trim(DIR_DATA_NC) // "/" // trim(fn_nc), nf90_write, ncid_dim))
    
    call check(nf90_inq_varid(ncid_dim, trim(str_var), varid))

    ! write out all vars' data
    if (is_var_3d == 0) then  ! 2d
        call check(nf90_put_var(ncid_dim, varid, var(:, :, 1), start = (/1, 1, curr_mon/), count = (/NX, NY, 1/)))
    else  ! 3d
        call check(nf90_put_var(ncid_dim, varid, var, start = (/1, 1, 1, curr_mon/), count = (/NX, NY, NZ, 1/)))
    end if

    call check(nf90_close(ncid_dim))

    write(*, *) "*** SUCCESS writing monthly average results for " // str_yr // str_mon // " " // trim(str_var) // " ***"

    return
end subroutine write_out_monthlyAvg_ice_nc_SODA_330


subroutine write_out_monthlyAvg_nc_SODA_341(curr_yr, curr_mon, is_var_3d, str_var, var)
    use netcdf
    implicit none

    integer, parameter :: N_MONTH = 12, NY = 330, NX = 720, NZ = 50

    integer,                        intent(in) :: curr_yr, curr_mon, is_var_3d
  ! character(len=128),             intent(in) :: str_var
    character*(*),             intent(in) :: str_var
    real   , dimension(NX, NY, NZ), intent(in) :: var

    ! local vars
  ! character(len = *), parameter :: DIR_DATA_NC = "/glade/p/umcp0006/lgchen/project/MOM_run/SODA_3.3_1/monthly/nc_ligang"
  ! character(len = *), parameter :: DIR_DATA_NC = "/glade/p/umcp0006/lgchen/project/MOM_run/SODA_3.5_1/monthly/nc_ligang"
  ! character(len = *), parameter :: DIR_DATA_NC = "/glade/scratch/lgchen/data/SODA_3.4_monthly_regrided_nc"
  ! character(len = *), parameter :: DIR_DATA_NC = "/glade/scratch/lgchen/project/SODA_3.3.9/month_avg/netcdf"
    character(len = *), parameter :: DIR_DATA_NC = "/glade/scratch/lgchen/tmp/20170203_soda3.4.1_mn_bin/netcdf"

    character(len=4  ) :: str_yr
    character(len=2  ) :: str_mon
    character(len=256) :: fn_nc  

    integer :: ncid_dim, varid


    write(unit=str_yr , fmt="(I4.4)") curr_yr
    write(unit=str_mon, fmt="(I2.2)") curr_mon

  ! fn_nc = "soda_3.3_mn_ocean_" // str_yr // ".nc"
  ! fn_nc = "soda_3.5_mn_ocean_" // str_yr // ".nc"
  ! fn_nc = "soda_3.4_mn_ocean_" // str_yr // ".nc"
  ! fn_nc = "soda3.3.1_mn_ocean_" // str_yr // ".nc"
    fn_nc = "soda3.4.1_mn_ocean_reg_" // str_yr // ".nc"


    ! open and read the dim file data
    call check(nf90_open(trim(DIR_DATA_NC) // "/" // trim(fn_nc), nf90_write, ncid_dim))
    
    call check(nf90_inq_varid(ncid_dim, trim(str_var), varid))

    ! write out all vars' data
    if (is_var_3d == 0) then  ! 2d
        call check(nf90_put_var(ncid_dim, varid, var(:, :, 1), start = (/1, 1, curr_mon/), count = (/NX, NY, 1/)))
    else  ! 3d
        call check(nf90_put_var(ncid_dim, varid, var, start = (/1, 1, 1, curr_mon/), count = (/NX, NY, NZ, 1/)))
    end if

    call check(nf90_close(ncid_dim))

    write(*, *) "*** SUCCESS writing monthly average results for " // str_yr // str_mon // " " // trim(str_var) // " ***"

    return
end subroutine write_out_monthlyAvg_nc_SODA_341




