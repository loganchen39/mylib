subroutine ssh_1dto2d_1day_1sat(ssh_1d, ssh_2d, num_ssh_2d, idx_lat_1dto2d, idx_lon_1dto2d)
    implicit none

    ! can not have intent(in) or f2py compiling failed due to double intent(in)
    !   which could be due to dimension(:).
    real   (kind=4), dimension(:) :: ssh_1d  
    real   (kind=4), dimension(1440, 720) :: ssh_2d
    integer(kind=4), dimension(1440, 720) :: num_ssh_2d
  ! integer(kind=4), dimension(1440, 720) :: landmask_oi  
    integer(kind=4), dimension(:) :: idx_lat_1dto2d
    integer(kind=4), dimension(:) :: idx_lon_1dto2d

! f2py instructions, for f90 use !f2py not Cf2py.
!f2py intent(in) ssh_1d
!f2py intent(out) ssh_2d
!f2py intent(out) num_ssh_2d
!f2py intent(in) idx_lat_1dto2d
!f2py intent(in) idx_lon_1dto2d

    ! local vars
    integer :: i, idx_lo, idx_la, OBS_SIZE
  ! parameter(OBS_SIZE = size(sst_iq))
    OBS_SIZE = size(ssh_1d)
    write(*, *) "OBS_SIZE=", OBS_SIZE

    ssh_2d     = 0.0
    num_ssh_2d = 0

    do i = 1, OBS_SIZE
        idx_lo = idx_lon_1dto2d(i)
        idx_la = idx_lat_1dto2d(i)

        ssh_2d    (idx_lo, idx_la) = ssh_2d    (idx_lo, idx_la) + ssh_1d(i)
        num_ssh_2d(idx_lo, idx_la) = num_ssh_2d(idx_lo, idx_la) + 1
    end do

    return
end subroutine ssh_1dto2d_1day_1sat


subroutine convert_sobsc_1dto2d(fn, sobsc, landmask_oi)
    implicit none

    character(len=*),                      intent(in) :: fn
    real   (kind=4), dimension(1440, 720), intent(out) :: sobsc 
    integer(kind=4), dimension(1440, 720), intent(in) :: landmask_oi

! f2py instructions, for f90 not to use Cf2py.
!f2py intent(out) sobsc
!f2py intent(in) landmask_oi

    ! local vars
    integer, parameter :: MAX_OBS_BOX = 10000
    real, parameter :: LAT_START = -89.875, LON_START = 0.125, GRID_SIZE = 0.25

    real, dimension(MAX_OBS_BOX) :: wsst, wlat, wlon, wnum
    character(len=8) :: ctype
    integer :: jyr, jmo, jda
    integer :: i, iostat_val, ict = 0, iv = 10
    integer :: ila, ilo

    real, parameter, dimension(720 ) :: LAT = (/ (LAT_START+(i-1)*GRID_SIZE, i=1,720 ) /) 
    real, parameter, dimension(1440) :: LON = (/ (LON_START+(i-1)*GRID_SIZE, i=1,1440) /) 

    !\without this, sobsc will be initialized to 0.0 because of intent(out),
    ! even if it's been assigned to -9999 in the python script! so you can also
    ! use intent(inout) to keep the -9999 value. 
    sobsc = -9999  

    !\the OISST uses access='sequential' as well, so it would be easier to use
    ! it here too. if you use 'direct' then you need to read the 'head' and
    ! 'tail' for each record!
    open(unit=iv, file=fn, form='unformatted', status='old', access='sequential', iostat=iostat_val, convert='big_endian')
    if (iostat_val /= 0) then
        write(*, *) "In Fortran subroutine convert_sobsc_1dto2d, open file failed!"
        write(*, *) "Failed file name: " // fn
        return
    end if

    ict = 0

  ! write(*, *) "anaysis_year=", analysis_year, "analysis_mon=", analysis_mon, "analysis_day=", analysis_day, "ctype=", ctype

    read(iv) jyr, jmo, jda, ctype
    read(iv, iostat=iostat_val) ict  ! each 'read' corresponds to a written record
    if (iostat_val/=0 .or. ict==0) then
        write(*, *) "Read data failed, iostat_vali", iostat_val, ", ict", ict
        close(iv)
        return
    endif
  ! write(*, *) "jyr=", jyr, "jmo=", jmo, "jda=", jda, "ctype=", ctype, "ict=", ict

    do while (ict > 0)
        read(iv) wnum(1:ict) 
        read(iv) wsst(1:ict)
        read(iv) wlat(1:ict)
        read(iv) wlon(1:ict)
      ! write(*, *) "wsst=", wsst(1:10)

        do i = 1, ict
            ila = nint( (wlat(i) + 89.875)*4 + 1 )
            ilo = nint( (wlon(i) - 0.125 )*4 + 1 )

            if (ila < 1) then
                write(*, *) "ila < 1", "ila=", ila, "wlat=", wlat(i)
                ila = 1
            else if (ila > 720) then
                write(*, *) "ila > 720", "ila=", ila, "wlat=", wlat(i)
                ila = 720
            end if

            if (ilo < 1) then
                write(*, *) "ilo < 1", "ilo=", ilo, "wlon=", wlon(i)
                ilo = 1
            else if (ilo > 1440) then
                write(*, *) "ilo > 720", "ilo=", ilo, "wlon=", wlon(i)
                ilo = 1440
            end if

            if (landmask_oi(ilo, ila) == 1) then
                sobsc(ilo, ila) = wsst(i)
            end if
        end do


        ! Read next record
        read(iv, iostat=iostat_val) ict
        if (iostat_val/=0 .or. ict==0) then
            write(*, *) "Read data failed, inside the loop, iostat_val=", iostat_val, ", ict=", ict
            close(iv)
            return
        endif

      ! write(*, *) 'ict=', ict
    end do

    
    close(iv)

    return
end subroutine convert_sobsc_1dto2d




subroutine write_1d_sat_sst_super_obs(fn, sst_oi, num_sst_oi, landmask_oi, analysis_year, analysis_mon, analysis_day, ctype)
    implicit none

    character(len=*),                      intent(in) :: fn
    real   (kind=4), dimension(1440, 720), intent(in) :: sst_oi
    integer(kind=4), dimension(1440, 720), intent(in) :: num_sst_oi
    integer(kind=4), dimension(1440, 720), intent(in) :: landmask_oi
    integer(kind=4),                       intent(in) :: analysis_year, analysis_mon, analysis_day
    character(len=8),                      intent(in) :: ctype   

! f2py instructions, for f90 not to use Cf2py.
!f2py intent(in) sst_oi
!f2py intent(in) num_sst_oi
!f2py intent(in) landmask_oi
!f2py intent(in) year
!f2py intent(in) month
!f2py intent(in) day
!f2py intent(in) ctype

    ! local vars
    integer, parameter :: MAX_OBS_BOX = 10000
    real, dimension(MAX_OBS_BOX) :: wsst, wlat, wlon, wnum
    integer :: i, j, ii, icts, iostat_val, ict = 0, iwrt = 0, isup = 0, iv = 10
    real(kind=4) :: lat, lon

    open(unit=iv, file=fn, form='unformatted', status='replace', iostat=iostat_val, convert='big_endian')
    if (iostat_val /= 0) then
        write(*, *) "In Fortran subroutine write_1d_sat_sst_super_obs, open file failed!"
        write(*, *) "Failed file name: " // fn
        return
    end if

    ict = 0
    iwrt = 0
    isup = 0

  ! write(*, *) "anaysis_year=", analysis_year, "analysis_mon=", analysis_mon, "analysis_day=", analysis_day, "ctype=", ctype

    write(iv) analysis_year, analysis_mon, analysis_day, ctype

    do j = 1, 720 
        lat = -89.875 + (j-1)*0.25
        do i = 1, 1440
            if (landmask_oi(i, j) == 1 .and. num_sst_oi(i, j) > 0.9) then
                lon = 0.125 + (i-1)*0.25

                ict = ict + 1 
                wnum(ict) = float(num_sst_oi(i, j) )
                wsst(ict) = sst_oi(i, j)
                wlat(ict) = lat
                wlon(ict) = lon
            end if

            if (ict == MAX_OBS_BOX) then
                iwrt = iwrt + 1

                write(iv) ict
                write(iv) (wnum (ii), ii = 1, ict)
                write(iv) (wsst (ii), ii = 1, ict)
                write(iv) (wlat (ii), ii = 1, ict)
                write(iv) (wlon (ii), ii = 1, ict)

                isup = isup + ict
                ict = 0
            end if
        end do 
    end do

    if (ict > 0) then
        isup = isup + ict
        iwrt = iwrt + 1
        icts = ict

        write(iv) ict
        write(iv) (wnum (ii), ii = 1, ict)
        write(iv) (wsst (ii), ii = 1, ict)
        write(iv) (wlat (ii), ii = 1, ict)
        write(iv) (wlon (ii), ii = 1, ict)
    endif

  ! write(*, *) "isup = ", isup
    close(iv)

    return
end subroutine write_1d_sat_sst_super_obs


subroutine sst_acspo2oi_1rec(sst_acspo, sst_oi, num_sst_oi, landmask_oi  &
    , idx_lat_oi2acspo, idx_lon_oi2acspo)
    implicit none

    real   (kind=4), dimension(18000, 9000), intent(in ) :: sst_acspo
    real   (kind=4), dimension(1440, 720), intent(out) :: sst_oi
    integer(kind=4), dimension(1440, 720), intent(out) :: num_sst_oi
    integer(kind=4), dimension(1440, 720), intent(in) :: landmask_oi
    integer(kind=4), dimension(2, 720 ), intent(in) :: idx_lat_oi2acspo
    integer(kind=4), dimension(2, 1440), intent(in) :: idx_lon_oi2acspo

! f2py instructions, for f90 not to use Cf2py.
!f2py intent(in) sst_acspo
!f2py intent(out) sst_oi
!f2py intent(out) num_sst_oi
!f2py intent(in) landmask_oi
!f2py intent(in) idx_lat_oi2acspo
!f2py intent(in) idx_lon_oi2acspo

    ! local vars
    integer :: i, j, lon_min, lon_max, lat_min, lat_max

    sst_oi     = 0.0
    num_sst_oi = 0

    do j = 1, 720; do i = 1, 1440
        if (landmask_oi(i, j) == 1) then
            lon_min = idx_lon_oi2acspo(1, i); lon_max = idx_lon_oi2acspo(2, i)
            lat_min = idx_lat_oi2acspo(1, j); lat_max = idx_lat_oi2acspo(2, j)
            sst_oi(i, j)     = sum  ( sst_acspo(lon_min:lon_max, lat_min:lat_max) )
            num_sst_oi(i, j) = count( sst_acspo(lon_min:lon_max, lat_min:lat_max) > 1.0)
        end if
    end do; end do

    return
end subroutine sst_acspo2oi_1rec


subroutine sst_pf2oi_1rec(sst_pf, sst_oi, num_sst_oi, landmask_oi  &
    , idx_lat_oi2pf, idx_lon_oi2pf)
    implicit none

    real   (kind=4), dimension(8640, 4320), intent(in ) :: sst_pf
    real   (kind=4), dimension(1440, 720), intent(out) :: sst_oi
    integer(kind=4), dimension(1440, 720), intent(out) :: num_sst_oi
    integer(kind=4), dimension(1440, 720), intent(in) :: landmask_oi
    integer(kind=4), dimension(2, 720 ), intent(in) :: idx_lat_oi2pf
    integer(kind=4), dimension(2, 1440), intent(in) :: idx_lon_oi2pf

! f2py instructions, for f90 use !f2py not Cf2py.
!f2py intent(in) sst_pf
!f2py intent(out) sst_oi
!f2py intent(out) num_sst_oi
!f2py intent(in) landmask_oi
!f2py intent(in) idx_lat_oi2pf
!f2py intent(in) idx_lon_oi2pf

    ! local vars
    integer :: i, j, lon_min, lon_max, lat_min, lat_max

    sst_oi     = 0.0
    num_sst_oi = 0

    do j = 1, 720; do i = 1, 1440
        if (landmask_oi(i, j) == 1) then
            lon_min = idx_lon_oi2pf(1, i); lon_max = idx_lon_oi2pf(2, i)
            lat_min = idx_lat_oi2pf(1, j); lat_max = idx_lat_oi2pf(2, j)
            sst_oi(i, j)     = sum  ( sst_pf(lon_min:lon_max, lat_min:lat_max) )
            num_sst_oi(i, j) = count( sst_pf(lon_min:lon_max, lat_min:lat_max) > 1.0)
        end if
    end do; end do

    return
end subroutine sst_pf2oi_1rec


subroutine sst_iquam2oi_1mon(sst_iq, sst_oi, num_sst_oi, landmask_oi  &
    , day, idx_lat_iq2oi, idx_lon_iq2oi, is_daytime, quality_level)
    implicit none

    ! can not have intent(in) or f2py compiling failed due to double intent(in)
    !   which could be due to dimension(:).
    real   (kind=8), dimension(:) :: sst_iq  
    real   (kind=4), dimension(1440, 720, 3, 31) :: sst_oi
    integer(kind=4), dimension(1440, 720, 3, 31) :: num_sst_oi
    integer(kind=4), dimension(1440, 720) :: landmask_oi  
    integer(kind=1), dimension(:) :: day
    integer(kind=4), dimension(:) :: idx_lat_iq2oi
    integer(kind=4), dimension(:) :: idx_lon_iq2oi
    integer(kind=2), dimension(:) :: is_daytime  ! f90 no unsigned type
    integer(kind=1), dimension(:) :: quality_level

! f2py instructions, for f90 use !f2py not Cf2py.
!f2py intent(in) sst_iq
!f2py intent(out) sst_oi
!f2py intent(out) num_sst_oi
!f2py intent(in) landmask_oi
!f2py intent(in) day
!f2py intent(in) is_daytime

    ! local vars
    integer :: i, idx_lo, idx_la, dy, is_dy, ql, OBS_SIZE
  ! parameter(OBS_SIZE = size(sst_iq))
    OBS_SIZE = size(sst_iq)
    write(*, *) "OBS_SIZE=", OBS_SIZE

    sst_oi     = 0.0
    num_sst_oi = 0

    do i = 1, OBS_SIZE
        idx_lo = idx_lon_iq2oi(i)
        idx_la = idx_lat_iq2oi(i)
        dy     = day(i)
        is_dy  = is_daytime(i)
        ql     = quality_level(i)

        if (landmask_oi(idx_lo, idx_la)==1 .and. ql==5) then
            if (sst_iq(i)<268.15 .or. sst_iq(i)>323.15) then
          ! if (sst_iq(i)<271.15 .or. sst_iq(i)>313.15) then
                write(*, *) "bad sst from iquam, sst=", sst_iq(i)
              ! cycle
            end if

            sst_oi    (idx_lo, idx_la, 3, dy) = sst_oi    (idx_lo, idx_la, 3, dy) + sst_iq(i)
            num_sst_oi(idx_lo, idx_la, 3, dy) = num_sst_oi(idx_lo, idx_la, 3, dy) + 1
            if (is_dy .ne. 0) then
                sst_oi    (idx_lo, idx_la, 1, dy) = sst_oi    (idx_lo, idx_la, 1, dy) + sst_iq(i)
                num_sst_oi(idx_lo, idx_la, 1, dy) = num_sst_oi(idx_lo, idx_la, 1, dy) + 1
            else
                sst_oi    (idx_lo, idx_la, 2, dy) = sst_oi    (idx_lo, idx_la, 2, dy) + sst_iq(i)
                num_sst_oi(idx_lo, idx_la, 2, dy) = num_sst_oi(idx_lo, idx_la, 2, dy) + 1
            end if
        end if
    end do

    return
end subroutine sst_iquam2oi_1mon


subroutine sst_acspo_l3u2oi_1rec(sst_ac, lon_ac, lat_ac, sst_oi, sst_stat_oi, landmask_oi, is_daytime)
    implicit none

    ! can not have intent(in) or f2py compiling failed due to double intent(in)
    !   which could be due to dimension(:).
    real   (kind=4), dimension(:) :: sst_ac
    real   (kind=4), dimension(:) :: lon_ac 
    real   (kind=4), dimension(:) :: lat_ac 
    real   (kind=4), dimension(2000, 1440, 720, 2) :: sst_oi
    real   (kind=4), dimension(1440, 720, 6, 2) :: sst_stat_oi  ! stats: num; mean; median; lowest; highest; std;
    integer(kind=4), dimension(1440, 720) :: landmask_oi
    integer(kind=2), dimension(:) :: is_daytime

! f2py instructions, for f90 not to use Cf2py.
!f2py intent(in) sst_ac
!f2py intent(in) lat_ac 
!f2py intent(in) lon_ac 
!f2py intent(inout) sst_oi
!f2py intent(inout) sst_stat_oi
!f2py intent(in) landmask_oi
!f2py intent(in) is_daytime

    ! local vars
    integer :: i_sst, i, j, k, sst_size
    sst_size = size(sst_ac)
    write(*, *) "In Fortran sst_acspo_l3u2oi_1rec, sst_size = ", sst_size 

    do i_sst = 1, sst_size
        j = nint(4*(lat_ac(i_sst)+89.875)) + 1 
        i = nint(4*(lon_ac(i_sst)-0.125 )) + 1

      ! write(*, *) 'lat_ac=', lat_ac(i_sst), ', lon_ac=', lon_ac(i_sst), ', j=', j, ', i=', i
 
        if (j < 1) then
            write(*, *) "In Fortran sst_acspo_l3u2oi_1rec, j < 1, lat = ", lat_ac(i_sst)
            j = 1
        end if
        if (j > 720) then
            write(*, *) "In Fortran sst_acspo_l3u2oi_1rec, j > 720, lat = ", lat_ac(i_sst)
            j = 720
        end if
        if (i < 1) then
            write(*, *) "In Fortran sst_acspo_l3u2oi_1rec, i < 1, lon = ", lon_ac(i_sst)
            i = 1
        end if
        if (i > 1440) then
            write(*, *) "In Fortran sst_acspo_l3u2oi_1rec, i > 1440, lon = ", lon_ac(i_sst)
            i = 1440
        end if

        if (landmask_oi(i, j) == 1) then
            if (is_daytime(i_sst) .ne. 0) then  ! daytime
                sst_stat_oi(i, j, 1, 1) = sst_stat_oi(i, j, 1, 1) + 1
                k = nint(sst_stat_oi(i, j, 1, 1))
                if (k .le. 2000) then
                    sst_oi(k, i, j, 1) = sst_ac(i_sst)
                else
                    write(*, *) "k .gt. 2000, exit!"
                    call exit(1)
                end if
            else  ! nighttime
                sst_stat_oi(i, j, 1, 2) = sst_stat_oi(i, j, 1, 2) + 1
                k = nint(sst_stat_oi(i, j, 1, 2))
                if (k .le. 2000) then
                    sst_oi(k, i, j, 2) = sst_ac(i_sst)
                else
                    write(*, *) "k .gt. 2000, exit!"
                    call exit(1)
                end if           
            end if
        end if
    end do  ! end of "do i_sst = 1, sst_size"

  ! write(*, *) "In Fortran sst_acspo_l3u2oi_1rec, returning"
    return
end subroutine sst_acspo_l3u2oi_1rec


subroutine calc_sst_stat_acspo_l3u2oi_1day(sst_oi, sst_stat_oi, landmask_oi)
    use support_subs

    implicit none
    real   (kind=4), dimension(2000, 1440, 720, 2) :: sst_oi
  ! real   (kind=4), dimension(6 , 1440, 720, 2) :: sst_stat_oi  ! stats: num; mean; median; lowest; highest; std;
    real   (kind=4), dimension(1440, 720, 6, 2) :: sst_stat_oi  ! stats: num; mean; median; lowest; highest; std;
    integer(kind=4), dimension(1440, 720) :: landmask_oi

! f2py instructions, for f90 not to use Cf2py.
!f2py intent(inout) :: sst_oi
!f2py intent(inout) :: sst_stat_oi
!f2py intent(in)    :: landmask_oi

! f2py intent(callback, hide) qsort


    ! local vars
    integer :: i, j, k, n_sst, i_sst
    real    :: mean, sumsqr
  ! write(*, *) "In Fortran sst_acspo_l3u2oi_1rec, sst_size = ", sst_size 

    do k = 1, 2; do j = 1, 720; do i = 1, 1440
        n_sst = nint(sst_stat_oi(i, j, 1, k))
        if (landmask_oi(i, j) /= 1 .or. n_sst == 0) then
            cycle
        end if  ! has to have then and end if

        ! mean
        sst_stat_oi(i, j, 2, k) = sum(sst_oi(1:n_sst, i, j, k)) / (1.0*n_sst)
        mean = sst_stat_oi(i, j, 2, k)

        ! std
        sumsqr = 0.0
        do i_sst = 1, n_sst
            sumsqr = sumsqr + (sst_oi(i_sst, i, j, k) - mean)**2
        end do
        sst_stat_oi(i, j, 6, k) = sqrt(sumsqr/n_sst)



        call qsort(sst_oi(1:n_sst, i, j, k))

      ! ! for debug
      ! if (n_sst >= 10) then
      !     write(*, *) "k=", k, ", j=", j, ", i=", i, ", n_sst=", n_sst, ", sorted_sst:"
      !     write(*, *) sst_oi(1:n_sst, i, j, k)
      ! end if

        ! lowest, highest
        sst_stat_oi(i, j, 4, k) = sst_oi(1    , i, j, k)  ! lowest
        sst_stat_oi(i, j, 5, k) = sst_oi(n_sst, i, j, k)  ! highest
        
        ! median
        if (mod(n_sst, 2) == 0) then
           sst_stat_oi(i, j, 3, k) = (sst_oi(n_sst/2, i, j, k) + sst_oi(n_sst/2 + 1, i, j, k)) / 2.0
        else
           sst_stat_oi(i, j, 3, k) = sst_oi(n_sst/2 + 1, i, j, k)
        end if

      ! ! for debug
      ! if (n_sst >= 50) then
      !     write(*, *) "k=", k, ", j=", j, ", i=", i, ", n_sst=", n_sst, ", mean=", mean, ", std=", sst_stat_oi(6, i, j, k)
      !     write(*, *) ", lowest=", sst_stat_oi(4, i, j, k), ", highest=", sst_stat_oi(5, i, j, k)
      ! end if     
    end do; end do; end do  ! end of "do j = 1, 720; do i = 1, 1440"

    return
end subroutine calc_sst_stat_acspo_l3u2oi_1day


