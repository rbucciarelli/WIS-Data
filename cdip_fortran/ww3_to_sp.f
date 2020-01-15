c-------------------------------------------------------------------------------
c  WW3_TO_SP - Converts ww3 spectral output into standard CDIP spectral files.
c    A direct port to f90 of bor's 'ww3_2_dfc' f77 code.
c-------------------------------------------------------------------------------
      program ww3_to_sp

      use archive_info
      use dates
      use elemental_stats
      use file_ops
      use lookup_utils
      use parameter
      use spectral
      use strings
      use unit_conversions

      integer::  bands, dbins, Dp, err_code, in_unit=97, peak_band, sec_length=0

      real
     *  a0(SP_max_coeffs), a1(SP_max_coeffs), a2(SP_max_coeffs),
     *  b1(SP_max_coeffs), b2(SP_max_coeffs), bw(SP_max_coeffs),
     *  check(SP_max_coeffs),
     *  depth, dir(SP_max_coeffs), ds(SP_max_coeffs,24),
     *  energy(SP_max_coeffs),
     *  freq(SP_max_coeffs),
     *  Hs,
     *  M0, M1,
     *  rad_dir, rdir(24), 
     *  Ta, Tp, total_energy

      logical::  has_direction

      character*3    station
      character*5    nstation
      character*14   sp_timestr
      character*19   source
      character*100  comment, fname, fpath, sp_name, sp_path
      character*500  dline

      type(ai_time_frame)  frame
      type(date_block)     sp_time
      type(pm_record)      pm_rec
      type(sp_data_block)  spec
      type(sp_hdr_block)   sp_hdr


c-- Check for command-line arguments, get ww3 filename

        num_args = iargc()
        if (num_args .eq. 9) then
          write(6,'(/,a)') 'Input the source file name: '
          read(5,'(a19)') source
        else
          call getarg(1,source)
        end if


c-- Set paths and filenames, open data file

        fname = source
        fpath = '/project/data06/ww3/'

        call open_read(in_unit, fpath, fname, err_code, 6)
        if (err_code .ne. 0) then
          write(6,'(a)') 'ERROR, aborting: can not open ww3 file'
          call exit(1)
        end if

        do j = 1, SP_max_coeffs
          check(j) = -1
        end do


c--   Begin reading wavewatch file spectrum. The array rdir contains the true 
c--   compass head from which the waves are heading in radians

        read(in_unit,'(a)') dline
        bands = get_field_int(dline,' ', 4)
        dbins = get_field_int(dline,' ', 5)	!* assumed equal to 24
        if (dbins .ne. 24) then
          write(6,'(a)') 'ERROR, aborting: dir bins != 24'
          call exit(1)
        end if

        read(in_unit,*) (freq(j), j=1, bands)        
        read(in_unit,*) (rdir(j), j=1, dbins)        


c--   Calculate frequency bandwidths 

        xfr=1.1

        do j=1,bands
          bw(j)=0.5*(xfr-1./xfr)*freq(j)
        end do


c--   Loop over spectal files within the WW3 output

        icnt=0

        read(in_unit, '(a)', iostat=err_code) dline
        do while (err_code .eq. 0)
          icnt=icnt+1

          do i=1,bands
            a0(i)=0
            a1(i)=0
            b1(i)=0
            a2(i)=0
            b2(i)=0
          end do

          read(dline,'(i4,i2,i2,1x,i2)') iy, im, id, ih
          sp_time = init_date(iy, im, id, ih, 0, 0)
          sp_timestr = make_datestring(sp_time)


c--   Read in forecast site name,latitude, E longitude, depth,
c--   wind speed (m/s), wind direction (T compass from),current
c--   speed, current direction

          read(in_unit,'(a)') dline
          nstation = dline(2:6)
          if (nstation == 'MAJUR') then
            station = '863'
          else if (nstation(1:3) == '460') then
            station = '7'//dline(5:6)
          else 
            station = '8'//dline(5:6)
          end if
          flat = get_field_real(dline, ' ', 3)
          flon = get_field_real(dline, ' ', 4)
          flon = flon - 360.0
          depth = get_field_real(dline, ' ', 5)
          wspd = get_field_real(dline, ' ', 6)
          wdir = get_field_real(dline, ' ', 7)
          curs = get_field_real(dline, ' ', 8)
          curd = get_field_real(dline, ' ', 9)


c--   Load archive, check for open time frame

          sp_name = 'sp'//station//'01'//sp_timestr(1:12)
          write(6,'(/,2a)') 'Processing file ', TRIM(sp_name)
          sp_path = get_file_path(sp_name, err_info%code)
          sp_path = './'

          frame = get_file_frame(sp_name, err_code)

          if (err_code .ne. 0) then
            write(6,'(a)') 'ERROR, aborting: no archive time frame found'
            call exit(1)
          end if


c--   Check position, depth values in WW3 file header

          if (NINT(depth) .ne. NINT(frame%water_depth/100.0)) 
     *      write(6,'(a)') 'WARNING: ww3 depth does not agree with archive'
          if (flat .ne. frame%deploy_site%lat) 
     *      write(6,'(a)') 'WARNING: ww3 latitude does not agree with archive'
          if (flon .ne. frame%deploy_site%long) 
     *      write(6,'(a)') 'WARNING: ww3 longitude does not agree with archive'
 

c--   Read in directional spectrum (m^2/hz-rad)

          read(in_unit,*,iostat=err_code) ((ds(i,j),i=1,bands),j=1,dbins)


c--   When calculating moments, rotate directions by pi so that the resulting 
c--   coefficients are in true compass "arriving from" coordinates.

          pi = 3.1415927

          do i=1, bands
            do ll=1, dbins


c--   Little floating point exception check here to catch 
c--   very small double-precision energies from ww3 model

              if(ds(i,ll).lt.1.0e-15) ds(i,ll) = 0.

              ds(i,ll) = ds(i,ll)*0.261799
              a0(i) = a0(i)+ds(i,ll)
              a1(i) = a1(i)+ds(i,ll)*cos(rdir(ll)+pi)
              b1(i) = b1(i)+ds(i,ll)*sin(rdir(ll)+pi)
              a2(i) = a2(i)+ds(i,ll)*cos(2.*(rdir(ll)+pi))
              b2(i) = b2(i)+ds(i,ll)*sin(2.*(rdir(ll)+pi))
            end do
            

c--   Normalize fourier coeffiecients

            if(a0(i).gt.0) then
              a1(i) = a1(i)/a0(i)
              b1(i) = b1(i)/a0(i)
              a2(i) = a2(i)/a0(i)
              b2(i) = b2(i)/a0(i)
            endif

c... convert freq spectrum from m^2/hz-rad to m^2/hz by multiplying times 2*pi
c              a0(i) = a0(i)*6.28318531


c--   Calculate direction

              rad_dir = ATAN2(b1(i),a1(i))
              dir(i) = to_degrees(rad_dir)
              if (dir(i) .lt. 0) dir(i) = dir(i) + 360

          end do


c--   Calculate Hs, Tp, Dp, Ta

          total_energy = 0
          M0 = 0
          M1 = 0

          do i = 1, bands
            energy(i) = a0(i) * bw(i)
            total_energy = total_energy + energy(i)
            M0 = M0 + energy(i)
            M1 = M1 + energy(i) * freq(i)
          end do

          Hs = 4 * SQRT(total_energy)
          peak_band = maxloc_lh(a0, 1, bands, SP_max_coeffs)
          Tp = 1.0 / freq(peak_band)
          Dp = NINT(dir(peak_band))
          Ta = M0 / M1


c--   Update the monthly xml file - first create sp data and header blocks

           spec%bands = bands
           do i = 1, bands
             spec%freq(i) = freq(i)
             spec%band_width(i) = bw(i)
             spec%ener_dens(i) = a0(i)
             spec%dir(i) = dir(i)
             spec%a1(i) = a1(i)
             spec%b1(i) = b1(i)
             spec%a2(i) = a2(i)
             spec%b2(i) = b2(i)
             spec%check(i) = check(i)
           end do

           sp_hdr%filename = sp_name(1:19)
           sp_hdr%sample_rate = frame%sample_rate
           sp_hdr%sample_length = 2400
           sp_hdr%position = frame%deploy_site
           sp_hdr%Hs = Hs
           sp_hdr%Tp = Tp
           sp_hdr%Ta = Ta
           sp_hdr%Dp = Dp

           call update_monthly_xml(spec, sp_hdr, .false., err_info%code)


c--   Create spectral file

          has_direction = .true.
          comment = 'Wind speed, dir: '
c         call create_spfile(sp_name, source, frame, sec_length, Hs, 
c    *      Tp, Dp, Ta, frame%water_depth/100.0, -1, freq, bw, a0, bands, has_direction, 
c    *      dir, a1, b1, a2, b2, check, comment)


c--   Update the parameter file

           call init_pmrecord(pm_rec)
           pm_rec%date = sp_time
           pm_rec%hs = Hs
           pm_rec%tp = Tp
           pm_rec%ta = Ta
           write(pm_rec%dp(1:3),'(i3)') Dp
c          call modify_parameter(station, 'p1', .false., pm_rec, err_code, .false.)

          read(in_unit, '(a)', iostat=err_code) dline

        end do

        close(in_unit)

      end
