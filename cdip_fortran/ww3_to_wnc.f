c-------------------------------------------------------------------------------
c  WW3_TO_WNC - Converts ww3 spectral output into netcdf files for net_model 
c    input. A direct port to f90 of bor's 'ww3_2_dfc' f77 code, with a switch 
c    from sp-format output to wnc output.
c-------------------------------------------------------------------------------
      program ww3_to_wnc

      use archive_info
      use dates
      use elemental_stats
      use file_ops
      use metadata_utils
      use spectral
      use strings
      use unit_conversions
      use wavecdf5_data
      use wmo_utils

      integer,parameter::  ww3_max_recs = 1000

      integer::  angle, bands, deg_per_bin, depth_diff, dbins, dir_idx(72), 
     *  err_code, in_unit=97, ncid, ogroups, peak_band, sec_length=0, times(ww3_max_recs)

      real
     *  a0(SP_max_coeffs, ww3_max_recs), a1(SP_max_coeffs, ww3_max_recs), 
     *  a2(SP_max_coeffs, ww3_max_recs),
     *  b1(SP_max_coeffs, ww3_max_recs), b2(SP_max_coeffs, ww3_max_recs), 
     *  bw(SP_max_coeffs), bin_dirs(72),
     *  check(SP_max_coeffs),
     *  depth, dist_diff, dir(SP_max_coeffs, ww3_max_recs), Dp(ww3_max_recs),
     *  ds(SP_max_coeffs,72), ds_set(72, SP_max_coeffs, ww3_max_recs),
     *  energy(SP_max_coeffs),
     *  freq(SP_max_coeffs), 
     *  Hs(ww3_max_recs),
     *  M0, M1,
     *  rad_dir, rdir(72), 
     *  Ta(ww3_max_recs), Tp(ww3_max_recs), total_energy, tesths

      logical exists, historic, zero_energy(ww3_max_recs)

      character*6    cdip_id
      character*6    yearmo
      character*10   stnid, pstnid
      character*20   histflag, id_string, nc_timestr
      character*100  fname, fpath, source, wnc_file, wnc_path, wnc_full
      character*500  comment, dline, nc_history

      type(ai_time_frame)  frame
      type(date_block)     sp_time
      type(location)       ww3_site
      type(meta_attribute_list)   meta_atts
      type(meta_variable_list)    meta_vars
      type(wc5_dataset)    wnc_set, onc_set

        write(6,'(a)') 'WW3_TO_WNC: WaveWatchIII netCDF processing'

c-- Check for command-line arguments, get ww3 filename

        num_args = iargc()
        if (num_args .eq. 0) then
          write(6,'(/,a)') 'Input the source file name: '
          read(5,'(a100)') source
        else
          call getarg(1,source)
        end if

        historic = .false.
        if (num_args .gt. 1) then
          call getarg(2,histflag)
          if (histflag(1:1) .eq. 'h') historic = .true.
        end if



c-- Set paths and filenames, open data file

        fname = source
        if (historic) then
          yearmo = fname(LEN_TRIM(fname)-5:LEN_TRIM(fname)+1)
          fpath = '/project/wave_hc/CFSR_phase_2/WW3_DATA/'//yearmo//'/'
        else
          fpath = '/project/data06/ww3/'
        end if

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
        dbins = get_field_int(dline,' ', 5)
        if (dbins .gt. 72) then
          write(6,'(a)') 'ERROR, aborting: dir bins != 72'
          call exit(1)
        end if
        deg_per_bin = 360/dbins

        read(in_unit,*) (freq(j), j=1, bands)        
        read(in_unit,*) (rdir(j), j=1, dbins)        


c--   Assign WW3 directions to correct bins (direction coming from, 
c--   starting with 5 deg bin)

        do i = 1, dbins
          angle = NINT(to_degrees(rdir(i))) + 180
          if (angle .ge. 360) angle = angle - 360
          dir_idx(i) = angle/deg_per_bin + 1
          bin_dirs(dir_idx(i)) = REAL(angle)
        end do


c--   Calculate frequency bandwidths 

c       xfr=1.1
c       do j=1,bands
c         bw(j)=0.5*(xfr-1./xfr)*freq(j)
c       end do

        do j=2, bands-1
          bw(j) = (freq(j+1)-freq(j-1))/2.0
        end do
 
        bw(1) = freq(1) * (bw(2)/freq(2))
        bw(bands) = freq(bands) * (bw(2)/freq(2))

        icnt=0
        pstnid = 'NULL'
        a0 = 0
        a1 = 0
        b1 = 0
        a2 = 0
        b2 = 0


c--   Loop over spectal files within the WW3 output

        read(in_unit, '(a)', iostat=err_code) dline
        do while (err_code .eq. 0)
          icnt=icnt+1

          read(dline,'(i4,i2,i2,1x,i2)') iy, im, id, ih
          times(icnt) = date_to_timestamp(init_date(iy, im, id, ih, 0, 0))


c--   Read in forecast site name,latitude, E longitude, depth,
c--   wind speed (m/s), wind direction (T compass from),current
c--   speed, current direction

          read(in_unit,'(a)') dline
          stnid = dline(2:11)
          if (stnid .ne. pstnid) then
            if (pstnid .ne. 'NULL ') then
              write(6,'(a)') 'ERROR, aborting: inconsistent station numbers'
              call exit(1)
            end if
          end if

          pstnid = stnid


          lon_start = INDEX(dline(16:28), '-', .true.)
          if (lon_start .eq. 0) lon_start = INDEX(dline(16:28), ' ') + 1
          lon_start = lon_start + 15

          read(dline(15:lon_start-1),'(f10.10)') flat
          read(dline(lon_start:28),'(f10.10)') flon
          if (flon .gt. 180) flon = flon - 360.0

          read(dline(30:38),'(f10.10)') depth

c         wspd = get_field_real(dline, ' ', 6)
c         wdir = get_field_real(dline, ' ', 7)
c         curs = get_field_real(dline, ' ', 8)
c         curd = get_field_real(dline, ' ', 9)




c--   Read in directional spectrum (m^2/hz-rad)

          read(in_unit,*,iostat=err_code) ((ds(i,j),i=1,bands),j=1,dbins)

          do i = 1, bands
            do j = 1, dbins
              ds_set(dir_idx(j),i,icnt) = ds(i,j)
            end do
          end do


c--   When calculating moments, rotate directions by pi so that the resulting 
c--   coefficients are in true compass "arriving from" coordinates.

          pi = 3.1415927

          zero_energy(icnt) = .true.
          do i=1, bands
            do j=1, dbins


c--   Little floating point exception check here to catch 
c--   very small double-precision energies from ww3 model

              if(ds_set(j,i,icnt).lt.1.0e-15) ds_set(j,i,icnt) = 0.

c             ds_set(j,i,icnt) = ds_set(j,i,icnt)*0.261799
              ds_set(j,i,icnt) = ds_set(j,i,icnt)*(2*3.14159/dbins)
              a0(i,icnt) = a0(i,icnt)+ds_set(j,i,icnt)
              a1(i,icnt) = a1(i,icnt)+ds_set(j,i,icnt)*COS(to_radians(bin_dirs(j)))
              b1(i,icnt) = b1(i,icnt)+ds_set(j,i,icnt)*SIN(to_radians(bin_dirs(j)))
              a2(i,icnt) = a2(i,icnt)+ds_set(j,i,icnt)*COS(to_radians(2.*bin_dirs(j)))
              b2(i,icnt) = b2(i,icnt)+ds_set(j,i,icnt)*SIN(to_radians(2.*bin_dirs(j)))
              ds_set(j,i,icnt) = ds_set(j,i,icnt) / deg_per_bin
            end do
            

c--   Normalize fourier coeffiecients

            if(a0(i,icnt).gt.0) then
              zero_energy(icnt) = .false.
              a1(i,icnt) = a1(i,icnt)/a0(i,icnt)
              b1(i,icnt) = b1(i,icnt)/a0(i,icnt)
              a2(i,icnt) = a2(i,icnt)/a0(i,icnt)
              b2(i,icnt) = b2(i,icnt)/a0(i,icnt)
              dir(i,icnt) = to_degrees(ATAN2(b1(i,icnt),a1(i,icnt)))
              if (dir(i,icnt) .lt. 0) dir(i,icnt) = dir(i,icnt) + 360
            endif

c... convert freq spectrum from m^2/hz-rad to m^2/hz by multiplying times 2*pi
c              a0(i) = a0(i)*6.28318531

          end do


c--   Calculate Hs, Tp, Dp, Ta

          if (.not. zero_energy(icnt)) then
            total_energy = 0
            M0 = 0
            M1 = 0

            do i = 1, bands
              energy(i) = a0(i,icnt) * bw(i)
              total_energy = total_energy + energy(i)
              M0 = M0 + energy(i)
              M1 = M1 + energy(i) * freq(i)
            end do

            Hs(icnt) = 4 * SQRT(total_energy)
            peak_band = maxloc_lh(a0(:,icnt), 1, bands, SP_max_coeffs)
            Tp(icnt) = 1.0 / freq(peak_band)
            Dp(icnt) = dir(peak_band,icnt)
            Ta(icnt) = M0 / M1
  

            total_energy = 0.0
            do i=1, bands
              do j=1, dbins
                total_energy = total_energy + ds_set(j,i,icnt)*deg_per_bin*bw(i)
              end do
            end do
            tesths = 4*total_energy**0.5
          else
            write(6,'(a)') 'NO PARAMS: zero-energy record'
          end if
          read(in_unit, '(a)', iostat=err_code) dline
        end do

        close(in_unit)

        write(6,'(3a,i4)') '  Processing ', TRIM(stnid), ', records: ', icnt

 
c--   Create wc5_dataset

        call wc5_initialize_set(wnc_set, 22)	!* gauge type 22 = ww3 output
        wnc_set%unlimited = .true.
        wnc_set%wave%time_count = icnt
        wnc_set%wave%freq_count = bands
        wnc_set%wave%dir_count = dbins
        call wc5_allocate_set(wnc_set)
        wnc_set%wave%dirs = bin_dirs(1:dbins)
        wnc_set%wave%freqs = freq(1:bands)
        wnc_set%wave%bw = bw(1:bands)
        wnc_set%wave%fflags = 1
        wnc_set%wave%fflags2 = 0
        call meta_initialize_attributes(meta_atts)
        nc_history = 'Ww3_to_wnc version 0.1. Runtime arguments:'//TRIM(source)
        call meta_add_attribute(meta_atts, 'history', nc_history, META_char_type)
        nc_timestr = write_iso_8601_date(current_utc())
        call meta_add_attribute(meta_atts, 'date_created', nc_timestr, META_char_type)
        call meta_add_attribute(meta_atts, 'date_modified', nc_timestr, META_char_type)
        call meta_initialize_variables(meta_vars)
        call meta_assign_model_variables(meta_vars, stnid, flat, flon, depth)

        do i = 1, icnt
          wnc_set%wave%times(i) = times(i)
          if (zero_energy(i)) then
            wnc_set%wave%flags(i) = 4
            wnc_set%wave%flags2(i) = 2
          else
            wnc_set%wave%flags(i) = 1
            wnc_set%wave%flags2(i) = 0
            wnc_set%wave%hs(i) = Hs(i)
            wnc_set%wave%tp(i)= Tp(i)
            wnc_set%wave%dp(i) = Dp(i)
            wnc_set%wave%ta(i) = Ta(i)
            
            do j = 1, bands
              wnc_set%wave%mdir(j,i) = dir(j,i)
              wnc_set%wave%a0(j,i) = a0(j,i)
              wnc_set%wave%a1(j,i) = a1(j,i)
              wnc_set%wave%b1(j,i) = b1(j,i)
              wnc_set%wave%a2(j,i) = a2(j,i)
              wnc_set%wave%b2(j,i) = b2(j,i)

              do k = 1, dbins
                wnc_set%wave%dirspec(k,j,i) = ds_set(k,j,i)
              end do
            end do
          end if
        end do


c--   Load extant nc file, combine

        if (historic) then
          wnc_file = TRIM(stnid)//'_CFSR2_historic.nc'
        else
          wnc_file = TRIM(stnid)//'_WW3_realtime.nc'
        end if
c       wnc_path = '/project/wvutil/ww3_to_wnc/output/'
        wnc_path = '/project/WNC/WNC_DATA/WW3/'
        wnc_full = TRIM(wnc_path)//TRIM(wnc_file)
        inquire(file=wnc_full, exist=exists) 

        if (exists) then

          wnc_set%classic = .true.
          wnc_set%is_2d_model = .true.
          call wc5_update_ncfile(wnc_full, wnc_set, err_code)

          call nc_call_func(nf90_open(wnc_full, NF90_WRITE, ncid))
          call nc_call_func(nf90_redef(ncid))
          nc_timestr = write_iso_8601_date(current_utc())
          call nc_call_func(nf90_put_att(ncid, NF90_GLOBAL, 'date_modified', TRIM(nc_timestr)))
          call nc_call_func(nf90_close(ncid))

c         depth_diff = ABS(NINT(depth-meta_depth))
c         write(6,'(2a,i4)') 'WARNING: site depth inconsistency, diff(m) = ', depth_diff
c         dist_diff = get_distance(meta_site, ww3_site)
c         write(6,'(2a,i4)') 'WARNING: site location inconsistency, diff(m) = ', dist_diff

        else
          call wc5_create_ncfile(wnc_set, wnc_full, err_code, enddef=.false.)
          call wc5_def_metadata_variables(wnc_set%ncid, meta_vars)
          call wc5_add_global_attributes(wnc_set%ncid, meta_atts)

          call load_wmo_array(err_code)
          if (err_code .ne. 0) then
            write(6,'(a)') '  Error - loading CDIP/WMO id table, aborting'
            call exit(1)
          end if

          cdip_id = get_cdip_id(TRIM(stnid))
          if (cdip_id .eq. 'NULL') cdip_id = '6'//stnid(4:5)
          cdip_id = cdip_id(1:3)//'w3'

          id_string = 'CDIP_'//TRIM(cdip_id)//'_WW3'
          call nc_call_func(nf90_put_att(wnc_set%ncid, NF90_GLOBAL, 'id', TRIM(id_string)))
          comment = 'This dataset is intended for use only as input to CDIP wave models. '//
     *      'For WAVEWATCH III data access please visit http://polar.ncep.noaa.gov/waves/.'
          call nc_call_func(nf90_put_att(wnc_set%ncid, NF90_GLOBAL, 'comment', TRIM(comment)))

          call nc_call_func(nf90_enddef(wnc_set%ncid))

          call wc5_put_metadata_variables(wnc_set%ncid, meta_vars)
          call wc5_fill_ncfile(wnc_set, err_code)
        end if

        write(6,'(a)') '  Complete; netCDF file updated.'

      end
