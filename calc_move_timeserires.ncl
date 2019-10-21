;-------------------------------------------------------------------------------------------------------------------------------------
; Bob Tomas
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; Code reads in a POP cross-section at locations where a j-index has ~constant latitude over a range of i-indices.  This is done to 
; reduce the volume of data processed in subsequent steps, which is important when looking at at multiple years of monthly 0.1 degree
; oean data.
;
; For use in doing analyses in the Atlantic, the RAPID (26.5N) and MOVE (16 N) sections.
; This approach has evolved over the last few months having worked on these analyses. It pares the full 4d data to the subset
; for further processing and plotting:  data(time, depth, latitude, longitude) --> data(time, depth, subset of longitude)
; 
; 0) Assumes that one wants a cross-section in the Atlantic at a ~constant latitude, over a specified longitude range
;    given by lat_in, lon_w_in, lon_e_in.
; 1) Because the left edge of the native grid starts ~40W and I am interested in cross sections in the Atlantic,  
;    rotate the data so that the longitudes ~range from -180 to 180 and change the longitude coordinats from 0 to 360 to -180 to 180.
; 2) Determine the i's & j(s) for a longitude range and latitude.  
; 3) Read in data at j or average over j's if there is more than 1 j ( of the form data(time, depth, longitude) )
; 4) Extract the longitude range of interest from the full longitude range ( of the form data(time, depth, subset of longitude) )
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;          18 October   2019 - I did some development work and found that the primary reason VGEO_TD and VGEO_TD_boundaries differed
;                              is that the longitude indices were not consistent between the two.  I replace all the coordinate 
;                              subscripting with index subscipting and chose ind_lon_loc_ge to be one greater than used for averaging
;                              vvel and MOVE time series calculated both ways agree pretty closely.  VGEO_BU and VGEO_BU_boundaries
;                              also agree closely if the lowest level specified for analysis does not intersect topography.
; 
;          14 October   2019 - I did some development work and found that the primary reason VGEO_TD and VGEO_BU differ is the choice
;                              of lower BC in VGEO_BU.  This conclusion is based on analyses I did comparing the two, and finding
;                              that when I remove VGEO_TD(bottom) from every level in a column, the resulting MOVE time series
;                              compares pretty well with VGEO_BU.
;          10 September 2019 - need to check that time_year is correct and why using data available at all longitudes gives a different
;                              mean value than using the data at just the boundaries
;           8 September 2019 - I need to put in extra boundaries to demark the beginning and end of where transport calculated from 
;                             (1) model velocitities and where the (2) geostrophic approximation begins and ends.  I found the
;                             calculations are sensitive to where the western edge of the geostrophic approximation is located.
;           7 September 2019 - update to handle LR ocean data
;           3 September 2019 - update to handle monthly data to produce time series of transport.
;          27 August 2019 - update to handle t-grid/u-grids correctly.  Prevously I averaged TEMP & SALT and SSH onto the u-grid
;                           and calculated density with those values.  Now I retain the two latitudes flanking VVEL, calculate
;                           density with those, and average density. - update: it makes very little difference.
;          26 August 2019 - updated to handle HR ocean data
;           2 August 2019
;
; 6 Aug 19 9:52 a.m. Status:
; I am working with a single LR decadal average.
; I just finished calcuating density and looked at the results.  The values are in line in an average density profile
; published in Pond & Pickard (1983), p. 9, Table 2.1.
; 7 Aug 19 2:53 p.m.  Status
; I calculated vvel geostrophic using points along y-index j; I should also do y-index j+1 and average with those at j, if I understand
; the POP grid stencil correctly.  My ininital calculations look a lot like vvel, so I am at least close.  There is some vertical 
; striping in the geostrophic velocity with small amplitude.
; I think I'll calculate a time series using the section integrated vvel and vvel_geo and compare the two.
; 9 Aug 19 7:42 a.m.  Status
; The vgeo cross section has vertical stripes at grid level scale using monthly means.  I talked with Gokhan, and he 
; suggested a couple of things.  One was to look at the SSH contribution.  When I removed SSH from the calculation,
; the striping seemed to go away.  He suggested that I try doing the calculation using all 4 points that surround the
; vvel point --- like one is suppose to instead of the points just to the south.  So I'm going to try that today.
; 12 Aug 19 8:51 a.m. Status:
; I got the 4-point stencil coded up and the stripes are not a problem any more.  I'm going to clean up the code today before I
; move on to the next steps. The next steps are to output 10 years of VGEO and apply the same metric I'm applying to VVEL 
; and compare the two results.
;
;=====================================================================================================================================
begin
;----------
; Model resolution
   cres = "LR"
;
; User specified lat, lon range
;
   lat_loc_in    =  16.0
;
; longitudes range from -180 to 180 after shifting (done below) for a sub-section
; These longitudes are used to pull out an Atlatntic full basin section.  The transport section is
; defined using longitude values located in code below 
   lon_loc_w   = -91.0
   lon_loc_e   = -10.0

; longitude and depth range for transport calculation.  Latitude is at lat_loc_in, above and longitude range must
; fall within the sub-section range specified, above.
   if(cres .eq. "LR")
      lon_loc_trans_cw = -61.0   ; -180 to 180.  Location of western edge of where model currents are used to calculate transport
      lon_loc_trans_ce = -60.0   ; -180 to 180.  Location of eastern edge of where model currents are used to calculate transport
      lon_loc_trans_gw = -59.0   ; -180 to 180.  Location of western edge of where geostrophic calculation is used to calculate transport
      lon_loc_trans_ge = -51.5   ; -180 to 180.  Location of eastern edge of where geostrophic calculation is used to calculate transport
   end if
   if(cres .eq. "HR")
      lon_loc_trans_cw = -61.5   ; -180 to 180.  Location of western edge of where model currents are used to calculate transport
      lon_loc_trans_ce = -60.1   ; -180 to 180.  Location of eastern edge of where model currents are used to calculate transport
      lon_loc_trans_gw = -52.0   ; -180 to 180.  Location of western edge of where geostrophic calculation is used to calculate transport
      lon_loc_trans_ge = -51.0   ; -180 to 180.  Location of eastern edge of where geostrophic calculation is used to calculate transport
   end if
;
   z_loc_lower     = 500000. ; cm
   z_loc_upper     = 120000. ; cm
;
; gravitational acceleration 
   grav = 9.80616  ; m s^-2
   rot  = 7.29212351699037d-05 ; rad s-1
   f    = 2.*rot*sin(lat_loc_in * 3.141592765d0/180. )
   print("====================================================================================================")
   print("f at " + lat_loc_in + "  = " + f)
   print("====================================================================================================")

; multiply factor for decibars to Pascals
   pa2dbar = 10000. ;  Pa to dbar

; multiply factor cm^3 to Sv
   cm3pstosv = 1./(100.*100*100*1.e6) ; cm->m  &  1.e6 m^3/s->Sv
;=====================================================================================================================================
; Read in coordinate data

   if (cres .ne. "LR" .and. cres .ne. "HR" ) then
      print("cres must be set to LR or HR.  It is set to: " + cres)
      exit
   end if
;
   if(cres .eq. "LR") then   
;      file_in = "/glade/work/tomas/RapidMoc/data/g20a10.GIAF_JRA.gx1v7.C03/g20a10.GIAF_JRA.gx1v7.C03.pop.h.004301_avg_005212.nc"
      file_in = "/glade/p/cgd/oce/projects/JRA55/IAF/g20a10.GIAF_JRA.gx1v7.C03/ocn/hist/g20a10.GIAF_JRA.gx1v7.C03.pop.h.0026-02.nc"
   end if
;
   if(cres .eq. "HR") then   
;      file_in = "/glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/hist/g.e20.G.TL319_t13.control.001_hfreq.pop.h.0026-01-10.nc"
      file_in = "/glade/scratch/altuntas/archive/g.e20.G.TL319_t13.control.001_hfreq/ocn/hist/mavg/g.e20.G.TL319_t13.control.001_hfreq.pop.h.0026-01.nc"
   end if
;
   f_in    = addfile(file_in, "r")
   ulat = f_in->ULAT
   ulon = f_in->ULONG
   tlat = f_in->TLAT
   tlon = f_in->TLONG
   dxu  = f_in->DXU
   dxt  = f_in->DXT
   dz   = f_in->dz
   z_t  = f_in->z_t
   depth = z_t
   depth = depth/100.
   depth@units = "meters"
   depth@valid_min = 5
   depth@valid_max = 5375
   depth!0 = "depth"
   depth&depth = depth
;=====================================================================================================================================
;----------
; initialize some indicies for shift the data
;----------
   if(cres .eq. "LR") then
      imaxgrid = 320
      ipivot_u = 194 ; for ULONG
      ipivot_t = 195 ; for TLONG
   end if
;
   if(cres .eq. "HR") then
      imaxgrid = 3200
      ipivot_u = 2899 ; for ULONG
      ipivot_t = 2900 ; for TLONG
   end if
;----------
; shift longitude arrays and change coordinate range from 0 to 360 to -180 to 180
;----------
;
   ulon_shift = ulon
   ulon_shift(:,imaxgrid-ipivot_u-1:imaxgrid-1) = ulon(:,0:ipivot_u)
   ulon_shift(:,0:imaxgrid-ipivot_u-2)      = ulon(:,ipivot_u+1:imaxgrid-1)
   ulon_shift = where(ulon_shift .gt. 180.d0, ulon_shift - 360.d0, ulon_shift)
;---
   tlon_shift = tlon
   tlon_shift(:,imaxgrid-ipivot_t-1:imaxgrid-1) = tlon(:,0:ipivot_t)
   tlon_shift(:,0:imaxgrid-ipivot_t-2)      = tlon(:,ipivot_t+1:imaxgrid-1)
   tlon_shift = where(tlon_shift .gt. 180.d0, tlon_shift - 360.d0, tlon_shift)
;----------
; shift latitude arrays
;----------
;
   ulat_shift = ulat
   ulat_shift(:,imaxgrid-ipivot_u-1:imaxgrid-1) = ulat(:,0:ipivot_u)
   ulat_shift(:,0:imaxgrid-ipivot_u-2)      = ulat(:,ipivot_u+1:imaxgrid-1)
;---
   tlat_shift = tlat
   tlat_shift(:,imaxgrid-ipivot_t-1:imaxgrid-1) = tlat(:,0:ipivot_t)
   tlat_shift(:,0:imaxgrid-ipivot_t-2)      = tlat(:,ipivot_t+1:imaxgrid-1)
;
;----------
; shift dx arrays
;----------
   dxu_shift = dxu
   dxu_shift(:,imaxgrid-ipivot_u-1:imaxgrid-1) = dxu(:,0:ipivot_u)
   dxu_shift(:,0:imaxgrid-ipivot_u-2)      = dxu(:,ipivot_u+1:imaxgrid-1)
;---
   dxt_shift = dxt
   dxt_shift(:,imaxgrid-ipivot_u-1:imaxgrid-1) = dxt(:,0:ipivot_u)
   dxt_shift(:,0:imaxgrid-ipivot_u-2)      = dxt(:,ipivot_u+1:imaxgrid-1)
;----------
; take a look...
;----------
   if(0 .eq. 1) then
      system("rm " + dir_out + "latlon.nc")
      f_out = addfile(dir_out + "/latlon.nc", "c")
      filedimdef(f_out,"time",-1,True)
      f_out->ulon       = ulon
      f_out->ulon_shift = ulon_shift
      f_out->ulat       = ulat
      f_out->ulat_shift = ulat_shift

      f_out->tlon       = tlon
      f_out->tlon_shift = tlon_shift
      f_out->tlat       = tlat
      f_out->tlat_shift = tlat_shift
   end if
;=====================================================================================================================================
; Find nearest i & j indices using ulon, ulat (not tlon, tlat).  
; These indices correspond to a larger cross-section than used to calculate the transport time series, for example, 
; a cross-section of the full Atlantic basin.  
;----------
; western most point
;----------
   lon_loc = lon_loc_w
;-
   dist_array = sqrt ( (ulat_shift - lat_loc_in)^2 + (ulon_shift - lon_loc)^2 )
   dims = dimsizes(dist_array)
;
   oned_dist_array = ndtooned(dist_array)
   oned_min_loc = minind(oned_dist_array)
   indices_min_loc = ind_resolve(oned_min_loc,dims)
   delete(dims)
   i_loc_ulon   = indices_min_loc(0,1)
   j_loc_ulat   = indices_min_loc(0,0)
   lon_loc_ulon = ulon_shift(indices_min_loc(0,0), indices_min_loc(0,1))
   lat_loc_ulat = ulat_shift(indices_min_loc(0,0), indices_min_loc(0,1))
;
; The i, j & lon, lat found in the ulon, ulat array that are closest to the lon, lat desired
;
   i_loc_ulon_w   = i_loc_ulon
   j_loc_ulat_w   = j_loc_ulat
   lon_loc_ulon_w = lon_loc_ulon
   lat_loc_ulat_w = lat_loc_ulat
;
;----------
; eastern most point
;----------
   lon_loc = lon_loc_e
;-
   dist_array = sqrt ( (ulat_shift - lat_loc_in)^2 + (ulon_shift - lon_loc)^2 )
   dims = dimsizes(dist_array)
;
   oned_dist_array = ndtooned(dist_array)
   oned_min_loc = minind(oned_dist_array)
   indices_min_loc = ind_resolve(oned_min_loc,dims)
   delete(dims)
   i_loc_ulon   = indices_min_loc(0,1)
   j_loc_ulat   = indices_min_loc(0,0)
   lon_loc_ulon = ulon_shift(indices_min_loc(0,0), indices_min_loc(0,1))
   lat_loc_ulat = ulat_shift(indices_min_loc(0,0), indices_min_loc(0,1))
;
; The i, j & lon, lat found in the ulon, ulat array that are closest to the lon, lat desired
;
   i_loc_ulon_e   = i_loc_ulon
   j_loc_ulat_e   = j_loc_ulat
   lon_loc_ulon_e = lon_loc_ulon
   lat_loc_ulat_e = lat_loc_ulat
;----------
; See if j_loc_ulat_e = j_loc_ulat_w.  For now, j_loc_ulat_e = j_loc_ulat_w.  It would be better to evaluate the sum of the 
; distances along each section and choose the j_loc which has the smallest distance
   njs_ulat   = abs(j_loc_ulat_w - j_loc_ulat_e) + 1
   if ( njs_ulat .ne. 1) then
      print("")
      print("for j's found using ulat, j at the west most point is different than j at the east most point")
      print("j_loc_ulat_w: " + j_loc_ulat_w + " lat_loc_ulat_w: " + lat_loc_ulat_w)
      print("j_loc_ulat_e: " + j_loc_ulat_e + " lat_loc_ulat_e: " + lat_loc_ulat_e)
;      
      print("using j at the west most point for both locations,")
      print("it is probably better to determine which j results in the smallest deviation from the desired latitude across the domain")
; use west lat location at east lat location (force j_loc_ulat_e to be the same as j_loc_ulat_w)
      j_loc_ulat_e = j_loc_ulat_w
      lat_loc_ulat_e =  ulat_shift(j_loc_ulat_w, i_loc_ulon_e)
   end if
;----------
   print("----------------------------------------------------------------------------------------------------")
   print("A basin-scale data section will be made using the coordinate values, below.")
   print("The coordinates were found using ulon and ulat, not tlong and tlat.")
   print("Verify that the lat/lon values and indices are to your liking.")
   print("These j/i, lat/lon values are from an array modified to put the GM near the middle of the domain, not the native ulat, ulon grids.")
   print("")
   print("You asked for:  lat_loc_in: " + lat_loc_in + "   lon_loc_w: " + lon_loc_w + "   lon_loc_e: " + lon_loc_e)
   print("Values from ulat, ulon:")
   print("j from ulat, west: " + j_loc_ulat_w + ", i from ulon, west: " + i_loc_ulon_w)
   print("lat from ulat, west: " + lat_loc_ulat_w + ", lon from ulon, west: " + lon_loc_ulon_w)
   print("---")
   print("j from ulat, east: " + j_loc_ulat_e + ", i from ulon, east: " + i_loc_ulon_e)
   print("lat from ulat, east: " + lat_loc_ulat_e + ", lon from ulon, east: " + lon_loc_ulon_e)
   print("----------------------------------------------------------------------------------------------------")
;=====================================================================================================================================
; (1) Read in data at desired latitude given by the j-index determined above, all longitudes.
; Note:  Although there is a j_w and  a j_e values, they should be forced to be the same, as the code is written (search for "force j_loc_ulat_e to be the same as j_loc_ulat_w", above ~20 lines).  
; At some point, you need to read in the j+1 value to properly calculate all
; the quantities needed to determine d(rho)/dx at the same location as the V-velocity.
; (2) Rotate as was done above for coordinate variables, so that longitudes range ~-180 to 180
;     This yields the "_all_lons" variables. 
; (3) Take the subset of desired longitudes given by the i-indices determined above
;     This yields the "_section" variables
; 
;----------------------------------------------------------------------
; tlon -  shifted form of variable already exists.  Extract sub-range of lat/lon to form a _section varaible, compatible with other _section variables
   dummy = dim_avg_n_Wrap(tlon_shift(j_loc_ulat_w:j_loc_ulat_e,i_loc_ulon_w:i_loc_ulon_e),0) ; dxt is 2d, lat,lon array read in early on in the code. dim_avg_n_Wrap removes first dimension which is 1.
   delete(tlon_shift)
   tlon_section = 1.d0*(dummy) ; float to double
   copy_VarAtts(dummy,tlon_section)
   delete(dummy)
;   print(tlon_section)
;----------
; ulon -  shifted form of variable already exists.  Extract sub-range of lat/lon to form a _section varaible, compatible with other _section variables
   dummy = dim_avg_n_Wrap(ulon_shift(j_loc_ulat_w:j_loc_ulat_e,i_loc_ulon_w:i_loc_ulon_e),0) ; dxt is 2d, lat,lon array read in early on in the code. dim_avg_n_Wrap removes first dimension which is 1.
   delete(ulon_shift)
   ulon_section = 1.d0*(dummy) ; float to double
   copy_VarAtts(dummy,ulon_section)
   delete(dummy)
;   printVarSummary(ulon_section)
;=====================================================================================================================================
; A diversion.  Calculate indicies for transport region.
;----------
; 14 October 2019 Some code to determine and print the values of indices used based on on the lon_loc's specified at the top of the code
;     lon_loc_trans_cw ; -180 to 180.  Location of western edge of where model currents are used to calculate transport
;     lon_loc_trans_ce ; -180 to 180.  Location of eastern edge of where model currents are used to calculate transport
;     lon_loc_trans_gw ; -180 to 180.  Location of western edge of where geostrophic calculation is used to calculate transport
;     lon_loc_trans_ge ; -180 to 180.  Location of eastern edge of where geostrophic calculation is used to calculate transport
;
   print("------------------------------------------------------------------------------------------------------------------------")
   print("Make sure the western wedge east most index and gestrophic west most index do not overlap")
   ind_lon_loc_trans_cw = ind_nearest_coord(lon_loc_trans_cw, ulon_section,0)
   print("ind_lon_loc_trans_cw: " + ind_lon_loc_trans_cw + "  lon from ulon(ind_lon_loc_trans_cw): " + ulon_section(ind_lon_loc_trans_cw))
   ind_lon_loc_trans_ce = ind_nearest_coord(lon_loc_trans_ce, ulon_section,0)
   print("ind_lon_loc_trans_ce: " + ind_lon_loc_trans_ce + "  lon from ulon(ind_lon_loc_trans_ce): " + ulon_section(ind_lon_loc_trans_ce))
   print("")
;
   ind_lon_loc_trans_gw = ind_nearest_coord(lon_loc_trans_gw, ulon_section,0)
   print("ind_lon_loc_trans_gw: " + ind_lon_loc_trans_gw + "  lon from ulon(ind_lon_loc_trans_gw): " + ulon_section(ind_lon_loc_trans_gw))
   ind_lon_loc_trans_ge = ind_nearest_coord(lon_loc_trans_ge, ulon_section,0)
   print("ind_lon_loc_trans_ge: " + ind_lon_loc_trans_ge + "  lon from ulon(ind_lon_loc_trans_ge): " + ulon_section(ind_lon_loc_trans_ge))
   print("")
;
   ind_z_loc_upper = ind_nearest_coord(z_loc_upper, z_t,0)
   print("ind_z_loc_upper: " + ind_z_loc_upper + "  z_t(ind_z_loc_upper): " + z_t(ind_z_loc_upper) + "cm")
;
   ind_z_loc_lower = ind_nearest_coord(z_loc_lower, z_t,0)
   print("ind_z_loc_lower: " + ind_z_loc_lower + "  z_t(ind_z_loc_lower): " + z_t(ind_z_loc_lower) + "cm")
;
   print("------------------------------------------------------------------------------------------------------------------------")
;
;=====================================================================================================================================
;----------
; dxt - shifted form of variable already exists.  Extract sub-range of lat/lon to form a _section varaible, compatible with other _section variables
   dummy = dim_avg_n_Wrap(dxt_shift(j_loc_ulat_w:j_loc_ulat_e,i_loc_ulon_w:i_loc_ulon_e),0) ; dxt is 2d, lat,lon array read in early on in the code.  dim_avg_n_Wrap removes first dimension which is 1.
   delete(dxt_shift)
   dxt_section = 1.d0*(dummy) ; float to double
   delete(dummy)
   dxt_section!0 = "tlon"
   dxt_section&tlon = tlon_section
;   printVarSummary(dxt_section)
;----------
; dxu - shifted form of variable already exists.  Extract sub-range of lat/lon to form a _section varaible, compatible with other _section variables
   dummy = dim_avg_n_Wrap(dxu_shift(j_loc_ulat_w:j_loc_ulat_e,i_loc_ulon_w:i_loc_ulon_e),0) ; dxu is 2d, lat,lon array read in early on in the code.  dim_avg_n_Wrap removes first dimension which is 1.
   delete(dxu_shift)
   dxu_section = 1.d0*(dummy) ; float to double
   copy_VarAtts(dummy,dxu_section)
   delete(dummy)
   dxu_section!0 = "ulon"
   dxu_section&ulon = ulon_section
;   printVarSummary(dxu_section)
;----------------------------------------------------------------------------------------------------------------------------------
; vvel - on u-grid.  Read in 1-latitude.
;----------
;   print("j_loc_ulat_w: " + j_loc_ulat_w + "   j_loc_ulat_e: " + j_loc_ulat_e)
;-----
   if(cres .eq. "LR") then
;      file_in = "/glade/work/tomas/RapidMoc/data/g20a10.GIAF_JRA.gx1v7.C03/g20a10.GIAF_JRA.gx1v7.C03.pop.h.004301_avg_005212.nc"
;      file_in = "/glade/p/cgd/oce/projects/JRA55/IAF/g20a10.GIAF_JRA.gx1v7.C03/ocn/hist/g20a10.GIAF_JRA.gx1v7.C03.pop.h.0026-02.nc"
      file_in = "/glade/work/tomas/RapidMoc/data/g20a10.GIAF_JRA.gx1v7.C03/g20a10.GIAF_JRA.gx1v7.C03.pop.h.VVEL.025001-026612.nc"
;-----
      f_in = addfile(file_in, "r")
      vvel_all_lons  = f_in->VVEL(:,:,j_loc_ulat_w:j_loc_ulat_e,:)
      vvel_all_lons = where( fabs(vvel_all_lons ) .gt. 9.e36, vvel_all_lons@_FillValue, vvel_all_lons)
;      printVarSummary(vvel_all_lons)
;-----
   time           = f_in->time
;-----
      dir_out = "/glade/work/tomas/RapidMoc/data/g20a10.GIAF_JRA.gx1v7.C03"
   end if ; LR case
;----------
   if(cres .eq. "HR") then   
;      file_in = "/glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/hist/g.e20.G.TL319_t13.control.001_hfreq.pop.h.0026-01-10.nc"
;      file_in = "/glade/scratch/altuntas/archive/g.e20.G.TL319_t13.control.001_hfreq/ocn/hist/mavg/g.e20.G.TL319_t13.control.001_hfreq.pop.h.0026-01.nc"
;
;      files_in = systemfunc("ls /glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/proc/tseries/monthly/g.e20.G.TL319_t13.control.001_hfreq.pop.h.VVEL.00{43,44}01-00??12.nc")
      files_in = systemfunc("ls /glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/proc/tseries/monthly/g.e20.G.TL319_t13.control.001_hfreq.pop.h.VVEL.00{43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59}01-00??12.nc")
;
      fs_in = addfiles(files_in, "r")
      vvel_all_lons  = fs_in[:]->VVEL(:,:,j_loc_ulat_w:j_loc_ulat_e,:) ; don't dim_avg over latitudes.  Retain latitude.
      vvel_all_lons = where( fabs(vvel_all_lons ) .gt. 9.e36, vvel_all_lons@_FillValue, vvel_all_lons)
;      printVarSummary(vvel_all_lons)
;-----
      time           = fs_in[:]->time
      dir_out = "/glade/work/tomas/RapidMoc/data/g.e20.G.TL319_t13.control.001_hfreq"
;-----
   end if ; HR case
;----------
;
   print("CHECK that the time_year is correct!!!")
   print("I'm not sure that it is correct!!!")
   time_year = cd_calendar(time,4)
   if(cres .eq. "LR") then
      time_year = time_year + 1750 ; year 2000 forcing corresponds to cycle 5, year 251 
   end if
   if(cres .eq. "HR") then
      time_year = time_year + 1957 ; year 2000 forcing corresponds to cycle 1, year 43
   end if
   time_year!0 = "time"
   time_year&time = time
; 
; [time | n] x [z_t | 62] x [nlat | 1] x [nlon | 3600] HR
;
; rotate the data and extract the subsection of interest
   vvel_all_lons_shift = vvel_all_lons
   vvel_all_lons_shift(:,:,:,imaxgrid-ipivot_u-1:imaxgrid-1) = vvel_all_lons(:,:,:,0:ipivot_u)
   vvel_all_lons_shift(:,:,:,0:imaxgrid-ipivot_u-2)      = vvel_all_lons(:,:,:,ipivot_u+1:imaxgrid-1)
   vvel_section = vvel_all_lons_shift(:,:,:,i_loc_ulon_w:i_loc_ulon_e)
   delete(vvel_all_lons_shift)
 ;  printVarSummary(vvel_section)
;----------
; a small diversion.  create coord_template_section_u
;----------
   coord_template_section_u = vvel_section
   coord_template_section_u = vvel_section@_FillValue
   coord_template_section_u!0 = "time"
   coord_template_section_u&time = time
   coord_template_section_u!1 = "z_t"
   coord_template_section_u&z_t = z_t
   coord_template_section_u!2 =  "ulat"
   coord_template_section_u&ulat = ulat_shift(j_loc_ulat_w:j_loc_ulat_e,i_loc_ulon_w)
   coord_template_section_u!3 =  "ulon"
   coord_template_section_u&ulon = ulon_section
;----------
; use coord_template_sction to assign coordinate variables
   copy_VarCoords(coord_template_section_u, vvel_section)
;   printVarSummary(vvel_section)
;----------------------------------------------------------------------------------------------------------------------------------
;----------
; temp - on t-grid. Read in 2-latitudes.
;----------
   if(cres .eq. "LR") then
;      file_in = "/glade/work/tomas/RapidMoc/data/g20a10.GIAF_JRA.gx1v7.C03/g20a10.GIAF_JRA.gx1v7.C03.pop.h.004301_avg_005212.nc"
;      file_in = "/glade/p/cgd/oce/projects/JRA55/IAF/g20a10.GIAF_JRA.gx1v7.C03/ocn/hist/g20a10.GIAF_JRA.gx1v7.C03.pop.h.0026-02.nc"
      file_in = "/glade/work/tomas/RapidMoc/data/g20a10.GIAF_JRA.gx1v7.C03/g20a10.GIAF_JRA.gx1v7.C03.pop.h.TEMP.025001-026612.nc"
;-----
      f_in = addfile(file_in, "r")
      temp_all_lons  = f_in->TEMP(:,:,j_loc_ulat_w:j_loc_ulat_e+1,:) ; don't dim_avg over latitudes.  Retain both latitudes to calculate the density and average those values.
      temp_all_lons = where( fabs(temp_all_lons ) .gt. 9.e36, temp_all_lons@_FillValue, temp_all_lons)
;      printVarSummary(temp_all_lons)
;-----
   end if ; LR case
;----------
   if(cres .eq. "HR") then   
;      file_in = "/glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/hist/g.e20.G.TL319_t13.control.001_hfreq.pop.h.0026-01-10.nc"
;      file_in = "/glade/scratch/altuntas/archive/g.e20.G.TL319_t13.control.001_hfreq/ocn/hist/mavg/g.e20.G.TL319_t13.control.001_hfreq.pop.h.0026-01.nc"
;
;      files_in = systemfunc("ls /glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/proc/tseries/monthly/g.e20.G.TL319_t13.control.001_hfreq.pop.h.TEMP.00{43,44,45,46,47}01-00??12.nc")

;      files_in = systemfunc("ls /glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/proc/tseries/monthly/g.e20.G.TL319_t13.control.001_hfreq.pop.h.TEMP.00{43,44}01-00??12.nc")
      files_in = systemfunc("ls /glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/proc/tseries/monthly/g.e20.G.TL319_t13.control.001_hfreq.pop.h.TEMP.00{43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59}01-00??12.nc")

      fs_in = addfiles(files_in, "r")
      temp_all_lons  = fs_in[:]->TEMP(:,:,j_loc_ulat_w:j_loc_ulat_e+1,:) ; don't dim_avg over latitudes.  Retain 2-latitudes to average onto u-grid.
      temp_all_lons = where( fabs(temp_all_lons ) .gt. 9.e36, temp_all_lons@_FillValue, temp_all_lons)
;-----
   end if ; HR case
;----------
; rotate the data and extract the subsection of interest
   temp_all_lons_shift = temp_all_lons
   temp_all_lons_shift(:,:,:,imaxgrid-ipivot_u-1:imaxgrid-1) = temp_all_lons(:,:,:,0:ipivot_u)
   temp_all_lons_shift(:,:,:,0:imaxgrid-ipivot_u-2)      = temp_all_lons(:,:,:,ipivot_u+1:imaxgrid-1)
   dummy = temp_all_lons_shift(:,:,:,i_loc_ulon_w:i_loc_ulon_e)
   delete(temp_all_lons_shift)
   temp_section = 1.d0*dummy ; float to double
   copy_VarAtts(dummy,temp_section)
   delete(dummy)
;----------
; a small diversion.  create coord_template_section_t with two latitudes
;----------
   coord_template_section_t = temp_section
   coord_template_section_t = temp_section@_FillValue
   coord_template_section_t!0 = "time"
   coord_template_section_t&time = time
   coord_template_section_t!1 = "z_t"
   coord_template_section_t&z_t = z_t
   coord_template_section_t!2 =  "tlat"
;   printVarSummary(coord_template_section_t)
;   printVarSummary( tlat_shift(j_loc_ulat_w:j_loc_ulat_e+1,i_loc_ulon_w) )
   coord_template_section_t&tlat = tlat_shift(j_loc_ulat_w:j_loc_ulat_e+1,i_loc_ulon_w)
   coord_template_section_t!3 =  "tlon"
   coord_template_section_t&tlon = tlon_section
;----------
; use coord_template_sction to assign coordinate variables
   copy_VarCoords(coord_template_section_t, temp_section)
;   printVarSummary(temp_section)
;Dimensions and sizes:    [time | 1] x [z_t | 62] x [tlat | 2] x [tlon | 811]
;Coordinates:
;            time: [9500..9500]
;            z_t: [500..587499.1]
;            tlat: [25.92276271463276..26.0127011319162]
;            tlon: [-90.94999999999999..-9.949999999999989]
;----------------------------------------------------------------------------------------------------------------------------------
; Extract the cross-sections of the data.
;----------------------------------------------------------------------------------------------------------------------------------
; salt - on t-grid. Read in 2-latitudes.
;----------
   if(cres .eq. "LR") then
;      file_in = "/glade/work/tomas/RapidMoc/data/g20a10.GIAF_JRA.gx1v7.C03/g20a10.GIAF_JRA.gx1v7.C03.pop.h.004301_avg_005212.nc"
;      file_in = "/glade/p/cgd/oce/projects/JRA55/IAF/g20a10.GIAF_JRA.gx1v7.C03/ocn/hist/g20a10.GIAF_JRA.gx1v7.C03.pop.h.0026-02.nc"
      file_in = "/glade/work/tomas/RapidMoc/data/g20a10.GIAF_JRA.gx1v7.C03/g20a10.GIAF_JRA.gx1v7.C03.pop.h.SALT.025001-026612.nc"
;-----
      f_in = addfile(file_in, "r")
      salt_all_lons  = f_in->SALT(:,:,j_loc_ulat_w:j_loc_ulat_e+1,:) ; don't dim_avg over latitudes.  Retain both latitudes to calculate the density and average those values.
      salt_all_lons = where( fabs(salt_all_lons ) .gt. 9.e36, salt_all_lons@_FillValue, salt_all_lons)
;      printVarSummary(salt_all_lons)
;-----
   end if ; LR case
;----------
   if(cres .eq. "HR") then
;      file_in = "/glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/hist/g.e20.G.TL319_t13.control.001_hfreq.pop.h.0026-01-10.nc"
;      file_in = "/glade/scratch/altuntas/archive/g.e20.G.TL319_t13.control.001_hfreq/ocn/hist/mavg/g.e20.G.TL319_t13.control.001_hfreq.pop.h.0026-01.nc"
;
;      files_in = systemfunc("ls /glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/proc/tseries/monthly/g.e20.G.TL319_t13.control.001_hfreq.pop.h.SALT.00{43,44,45,46,47}01-00??12.nc")

;      files_in = systemfunc("ls /glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/proc/tseries/monthly/g.e20.G.TL319_t13.control.001_hfreq.pop.h.SALT.00{43,44}01-00??12.nc")
      files_in = systemfunc("ls /glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/proc/tseries/monthly/g.e20.G.TL319_t13.control.001_hfreq.pop.h.SALT.00{43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59}01-00??12.nc")

      fs_in = addfiles(files_in, "r")
      salt_all_lons  = fs_in[:]->SALT(:,:,j_loc_ulat_w:j_loc_ulat_e+1,:) ; don't dim_avg over latitudes.  Retain 2-latitudes to average onto u-grid.
      salt_all_lons = where( fabs(salt_all_lons ) .gt. 9.e36, salt_all_lons@_FillValue, salt_all_lons)
;-----
   end if ; HR case
;----------
; check for scale factor on SALT and apply and delete it if it is present
   if( isatt(salt_all_lons, "scale_factor") ) then
      salt_all_lons = salt_all_lons*salt_all_lons@scale_factor
      delete(salt_all_lons@scale_factor)
   end if
;----------
; rotate the data and extract the subsection of interest
   salt_all_lons_shift = salt_all_lons
   salt_all_lons_shift(:,:,:,imaxgrid-ipivot_u-1:imaxgrid-1) = salt_all_lons(:,:,:,0:ipivot_u)
   salt_all_lons_shift(:,:,:,0:imaxgrid-ipivot_u-2)      = salt_all_lons(:,:,:,ipivot_u+1:imaxgrid-1)
   dummy = salt_all_lons_shift(:,:,:,i_loc_ulon_w:i_loc_ulon_e)
   delete(salt_all_lons_shift)
   salt_section = 1.d0*(dummy) ; float to double
   copy_VarAtts(dummy,salt_section)
   delete(dummy)
;----------
; use coord_template_sction to assign coordinate variables
   copy_VarCoords(coord_template_section_t, salt_section)
;   printVarSummary(salt_section)
;----------
; ssh  - on t-grid. Read in 2-latitudes.
;----------
   if(cres .eq. "LR") then
;      file_in = "/glade/work/tomas/RapidMoc/data/g20a10.GIAF_JRA.gx1v7.C03/g20a10.GIAF_JRA.gx1v7.C03.pop.h.004301_avg_005212.nc"
;      file_in = "/glade/p/cgd/oce/projects/JRA55/IAF/g20a10.GIAF_JRA.gx1v7.C03/ocn/hist/g20a10.GIAF_JRA.gx1v7.C03.pop.h.0026-02.nc"
      file_in = "/glade/work/tomas/RapidMoc/data/g20a10.GIAF_JRA.gx1v7.C03/g20a10.GIAF_JRA.gx1v7.C03.pop.h.SSH.025001-026612.nc"
;-----
      f_in = addfile(file_in, "r")
      ssh_all_lons  = f_in->SSH(:,j_loc_ulat_w:j_loc_ulat_e+1,:) ; don't dim_avg over latitudes.  Retain both latitudes to calculate the density and average those values.
      ssh_all_lons = where( fabs(ssh_all_lons ) .gt. 9.e36, ssh_all_lons@_FillValue, ssh_all_lons)
;      printVarSummary(ssh_all_lons)
;-----
   end if ; LR case
;----------
   if(cres .eq. "HR") then
;      file_in = "/glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/hist/g.e20.G.TL319_t13.control.001_hfreq.pop.h.0026-01-10.nc"
;      file_in = "/glade/scratch/altuntas/archive/g.e20.G.TL319_t13.control.001_hfreq/ocn/hist/mavg/g.e20.G.TL319_t13.control.001_hfreq.pop.h.0026-01.nc"
;
;      files_in = systemfunc("ls /glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/proc/tseries/monthly/g.e20.G.TL319_t13.control.001_hfreq.pop.h.SSH.00{43,44,45,46,47}01-00??12.nc")

;      files_in = systemfunc("ls /glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/proc/tseries/monthly/g.e20.G.TL319_t13.control.001_hfreq.pop.h.SSH.00{43,44}01-00??12.nc")
      files_in = systemfunc("ls /glade/scratch/tomas/g.e20.G.TL319_t13.control.001_hfreq/ocn/proc/tseries/monthly/g.e20.G.TL319_t13.control.001_hfreq.pop.h.SSH.00{43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59}01-00??12.nc")

      fs_in = addfiles(files_in, "r")
      ssh_all_lons  = fs_in[:]->SSH(:,j_loc_ulat_w:j_loc_ulat_e+1,:) ; don't dim_avg over latitudes.  Retain 2-latitudes to average onto u-grid.
      ssh_all_lons = where( fabs(ssh_all_lons ) .gt. 9.e36, ssh_all_lons@_FillValue, ssh_all_lons)
;-----
   end if ; HR case
;----------
; rotate the data and extract the subsection of interest
   ssh_all_lons_shift = ssh_all_lons
;   printVarSummary(ssh_all_lons)

   ssh_all_lons_shift(:,:,imaxgrid-ipivot_u-1:imaxgrid-1) = ssh_all_lons(:,:,0:ipivot_u)
   ssh_all_lons_shift(:,:,0:imaxgrid-ipivot_u-2)      = ssh_all_lons(:,:,ipivot_u+1:imaxgrid-1)


   dummy = ssh_all_lons_shift(:,:,i_loc_ulon_w:i_loc_ulon_e)
   delete(ssh_all_lons_shift)
   ssh_section = 1.d0*(dummy) ; float to double
   copy_VarAtts(dummy,ssh_section)
   delete(dummy)
;----------
; use coord_template_sction to assign coordinate variables
   copy_VarCoords(coord_template_section_t(:,0,:,:), ssh_section) ; don't assign depth coordinate
;   printVarSummary(ssh_section)
;----------
   if(1 .eq. 1) then
      system("rm " + dir_out + "/state.nc")
      f_out = addfile(dir_out + "/state.nc", "c")
      filedimdef(f_out,"time",-1,True)
      f_out->vvel_section=vvel_section
      f_out->temp_section=temp_section
      f_out->salt_section=salt_section
      f_out->ssh_section=ssh_section
      f_out->dxu_section=dxu_section
      f_out->dxt_section=dxu_section
      f_out->dz=dz
   end if
;-------------------------------------------------------------------------------------------------------------------------------------
; calculate density 
   dims  = dimsizes(temp_section)
   ntime  = dims(0)
   ndepth = dims(1)
   nlat  = dims(2)
   nlon  = dims(3)
   delete(dims)
;
   rho_section = temp_section
   rho_section = rho_section@_FillValue
   rho_section@units = "kg/m3"
   rho_section@long_name = "Density"
;   printVarSummary(rho_section)
; For some reason, this function has t & s as a fuction of depth, longitude but depth is a scalar
; I guess because it was written to calculate potential density (depth = 0) and not the actual density.
   do nt = 0, ntime - 1
      do nj = 0, nlat - 1
         do nz = 0, ndepth - 1
            dummy = rho_mwjf(temp_section(nt, nz:nz, nj, :), salt_section(nt, nz:nz, nj,:), depth(nz)) ; just calculate density at the level of interest
            rho_section(nt, nz, nj, :) = (/ 1000.*dummy(0,:) /); multiply by 1000 to put into kg/m3, although the ncl documentation 
                                                               ;  says something different:  pd = 1000.*(pd-1.)
         end do
      end do
   end do
;   printVarSummary(rho_section)
;   Dimensions and sizes:    [time | 1] x [z_t | 62] x [tlat | 2] x [tlon | 811]
;----------
; use coord_template_sction to assign coordinate variables
   copy_VarCoords(coord_template_section_t, rho_section)
;   printVarSummary(rho_section)
;----------
;
   if(1 .eq. 1) then
      system("rm " + dir_out + "/density.nc")
      f_out = addfile(dir_out + "/density.nc", "c")
      filedimdef(f_out,"time",-1,True)
      f_out->rho_section=rho_section
   end if
;-----------
; integrate w.r.t. dz.  From the top to bottom for now.
   pres_section = rho_section
   pres_section@long_name = "Pressure"
   pres_section@units     = "decibar"
   pres_section = pres_section@_FillValue
;Dimensions and sizes:    [time | 1] x [z_t | 62] x [tlat | 2] x [tlon | 811]
;   printVarSummary(pres_section)
;-----------
; Pressure; include the SSH contribution by adding SSH to the top layer depth
   thickness_top = ssh_section/100. + dz(0)/100. ; convert cm to m and add to the top layer thickness
;
;commmented out so no kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge 
;   print("OMITTING SSH!!!!!")
;   print("OMITTING SSH!!!!!")
;   print("OMITTING SSH!!!!!")
;   print("OMITTING SSH!!!!!")
;   print("OMITTING SSH!!!!!")
;   thickness_top = dz(0)/100. ; convert cm to m and add to the top layer thickness
;commmented out so no kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge kludge 
;
   pres_section(:,0,:,:) = rho_section(:,0,:,:)*thickness_top*grav/pa2dbar
;
; 
;-----------
; 3oct 19 Set top level pressure to 0.
;   pres_section(:,0,:,:) = 0.
;----------
; Try and omit rho y-variations (set to rho=rho_average in y )
;-----------
; 3oct19
;   print("?????????? experimenting ??????????")
;   print("Setting 1000m  levels pressure contribution to 0")
;   print("?????????? experimenting ??????????")
;   printVarSummary(rho_section)
;
; make a copy of rho_section to doctor up and use in the pressure calculation, immediately below
   rho_section_copy = rho_section
;-----------
; Now integrate downward
;-----------
   do nz = 1, ndepth - 1
      pres_section(:, nz, :, :) = pres_section(:,nz-1, :, :) + rho_section_copy(:, nz, :, :)*(dz(nz)/100.)*grav/pa2dbar
   end do 
;----------
; use coord_template_sction to assign coordinate variables
   copy_VarCoords(coord_template_section_t, pres_section)
;   printVarSummary(pres_section)
;----------
;
   if(1 .eq. 1) then
      system("rm " + dir_out + "/pressure.nc")
      f_out = addfile(dir_out + "/pressure.nc", "c")
      filedimdef(f_out,"time",-1,True)
      f_out->pres_section=pres_section
   end if
;-------------------------------------------------------------------------------------------------------------------------------------
; Calculate d/dx of the latitude average pressure (there should be two t-grid latitudes that lie on each side of the u-grid latitude and average to the u-grid value).  
   p1 = pres_section&tlat
   p2 = vvel_section&ulat
   print("")
   print("------------------------------------------------------------------------------------------------------------------------")
   print("The average of the t-grid latitudes should equal the u-grid latitude to ~6 digits.  t-grid lat avg: " + dim_avg(p1) + "  u-grid lat: " + p2)
   print("------------------------------------------------------------------------------------------------------------------------")
   print("")
   delete(p1)
   delete(p2)
;
; Average pres_section over both latitudes effectively putting it onto the u-grid and calculate the derivative using dxu
;
   pres_section_avg = 1.d0*vvel_section
   pres_section_avg(:,:,0,:) = dim_avg_n(pres_section, 2)
   copy_VarCoords(coord_template_section_u, pres_section_avg)
;   printVarSummary(pres_section_avg)
;
; 3oct19 Take a look at rho_section_avg  and set to a constant - doesn't make much difference in the index value
   rho_section_avg = 1.d0*vvel_section
   rho_section_avg(:,:,0,:) = dim_avg_n(rho_section, 2)
   copy_VarCoords(coord_template_section_u, rho_section_avg)
;   file_out = "rho_section_avg.nc"
;   system("rm " + file_out)
;   f_out = addfile(file_out, "c")
;   f_out->rho_section_avg = rho_section_avg
;
;   16oct19 set rho_section_avg = 1030.
;   print("setting rho_section_avg = 1030. everywhere")
;   rho_section_avg = 1030.
;
   dpdx_section = pres_section_avg
   copy_VarCoords(coord_template_section_u, rho_section_avg)
   do ni = 0, nlon - 2
      dpdx_section(:,:,:,ni) = pa2dbar*( pres_section_avg(:,:,:,ni+1) - pres_section_avg(:,:,:,ni) ) / ((dxu_section(ni+1)+dxu_section(ni))/2./100. ) ; dbar to Pa in the numerator and cm to m in the denominator 
   end do
; assign value at nlon-2 to nlon-1  for now.  Include an extra i-value in temp, salt, ssh, then rho and calculate the correct value in the future.
   ni = nlon - 1
   dpdx_section(:,:,:,ni) = dpdx_section(:,:,:, ni-1)
;------------
; calculate the top-down meridional geostrophic velocity
;------------
   vgeo_td_section = dpdx_section
   vgeo_td_section = (/  100.*(1./rho_section_avg)*(1./f) * dpdx_section /) ; in units of cm s-1 to compare with vvel
   copy_VarCoords(coord_template_section_u, vgeo_td_section)
;   printVarSummary(vgeo_td_section)
;------------
; 3oct19 try matching the BU calculation by subtracting the bottom velocity from all levels (level of no motion is the bottom)
;------------
;   printVarSummary(vgeo_td_section)
;----------
; for the cross section, find the lowest level where there is a valid value (not missing value)
; make an array of size nlon that has the k-index of the deepest valid ocean point at that lon
; This is time invariant, so use time index = 0
   max_depth_ind_section = new( nlon, "integer")
;
   do ni = 0, nlon - 1
      l_valid_vals =  .not.ismissing(vgeo_td_section(0,:,0,ni))
;
      if(.not.all(l_valid_vals)) then
         inds = ind(l_valid_vals) 
         maximum_ind = max(inds)
         max_depth_ind_section(ni) = maximum_ind
;         print("maximum_ind: " + maximum_ind)
         delete(inds)
      end if
   end do
;   print("max_depth_ind_section: " + max_depth_ind_section)
;----------
; Remove the vertical average or the bottom value from each vertical profile, all times
;   do ni = 0, nlon - 1
;      if(.not.ismissing(max_depth_ind_section(ni))) then
;         printVarSummary(vgeo_td_section(:,:,:,ni))
;         printVarSummary(vgeo_td_section(:,max_depth_ind_section(ni),:,ni))
;
; 14 October 2019 try subtracting the weighted vertical mean (a.k.a. barotropic component); I'm not sure why, but this results in something very different than omitting the SSH from the VGEO_TD 
; calculation
;            vgeo_td_section(:,:,:,ni) = vgeo_td_section(:,:,:,ni) - conform( vgeo_td_section(:,:,:,ni), dim_avg_wgt_n(vgeo_td_section(:,:,:,ni),dz,1,1), (/0,2/))
;
;  try subtracting the bottom velocity value at all levels
;             vgeo_td_section(:,:,:,ni) = vgeo_td_section(:,:,:,ni) - conform( vgeo_td_section(:,:,:,ni), vgeo_td_section(:,max_depth_ind_section(ni),:,ni), (/0,2/))
;         
;      end if
;   end do
;
;   file_out = "no_motion.nc"
;   system("rm " + file_out)
;   f_out = addfile(file_out, "c")
;   f_out->vgeo_td_section=vgeo_td_section
;-----------------------------------------------------------------------------------------------------------------------------
; Calculate d/dx of latitude average density.  (there should be two t-grid latitudes that lie on each side of the u-grid latitude and average to the u-grid value).
   p1 = rho_section&tlat
   p2 = vvel_section&ulat
   print("")
   print("------------------------------------------------------------------------------------------------------------------------")
   print("The average of the t-grid latitudes should equal the u-grid latitude to ~6 digits.  t-grid lat avg: " + dim_avg(p1) + "  u-grid lat: " + p2)
   print("------------------------------------------------------------------------------------------------------------------------")
   print("")
   delete(p1)
   delete(p2)
;
; Average rho_section over both latitudes effectively putting it onto the u-grid and calculate the derivative using dxu
;
   rho_section_avg = 1.d0*vvel_section
   rho_section_avg(:,:,0,:) = dim_avg_n(rho_section, 2)
   copy_VarCoords(coord_template_section_u, rho_section_avg)
;   printVarSummary(rho_section_avg)
;
   drhodx_section = rho_section_avg
   copy_VarCoords(coord_template_section_u, rho_section_avg)
   do ni = 0, nlon - 2
      drhodx_section(:,:,:,ni) = ( rho_section_avg(:,:,:,ni+1) - rho_section_avg(:,:,:,ni) ) / ((dxu_section(ni+1)+dxu_section(ni))/2. )
   end do
; assign value at nlon-2 to nlon-1  for now.  Include an extra i-value in temp, salt, ssh, then rho and calculate the correct value in the future.
   ni = nlon - 1
   drhodx_section(:,:,:,ni) = drhodx_section(:,:,:, ni-1)
;------------
; calculate the vertical gradient of the meridional geostrophic velocity
;------------
   dvgeodz_section = drhodx_section
   dvgeodz_section = (/  (-1.*grav/rho_section_avg)*(1./f) * drhodx_section /) ; factor of 100 in numerator to convert from m to cm
   copy_VarCoords(coord_template_section_u, dvgeodz_section)
   copy_VarCoords(coord_template_section_u, drhodx_section)
;-----------
; Now bottom-up taking vgeo = 0 at the bottom
;-----------
   vgeo_bu_section = dvgeodz_section
; put zeros where there are not missing values
   vgeo_bu_section = where(ismissing(vgeo_bu_section),vgeo_bu_section@_FillValue,0.)
;   printVarSummary(vgeo_bu_section)
;
; Integrate from the bottom up.  Integrate only when the layer below is not a missing value. Use 1/2 dvdz from below and 1/2 dvdz at the same level.
   do nz = ndepth - 2, 0, 1
      vgeo_bu_section(:,nz,:,:) = where( ismissing(vgeo_bu_section(:,nz+1,:,:)), vgeo_bu_section(:,nz,:,:), vgeo_bu_section(:,nz+1,:,:) + (0.5*dvgeodz_section(:,nz+1,:,:)*dz(nz+1) + 0.5*dvgeodz_section(:,nz,:,:)*dz(nz)))
   end do
   vgeo_bu_section = (/ 100.*vgeo_bu_section /)
;----------
;
   if(1 .eq. 1) then
      system("rm " + dir_out + "/dpdx_vgeo_td_bu.nc")
      f_out = addfile(dir_out + "/dpdx_vgeo_td_bu.nc", "c")
      filedimdef(f_out,"time",-1,True)
      f_out->dpdx_section=dpdx_section
      f_out->vgeo_td_section=vgeo_td_section
      f_out->vgeo_bu_section=vgeo_bu_section
      f_out->dxu_section=dxu_section
      f_out->dz=dz
   end if
;-------------------------------------------------------------------------------------------------------------------------------------
;
; This is a quick and dirty calculation.  I have not double checked to make sure the indexing is correct; I think it
; is, but it could be off by a little in the x and/or y directtion. That is, vvel and dxu do not line up properly
   print("Indexing of vvel, dxu, dz is at least close; need to verify that they are correct")
   flux_vvel_section = vvel_section
   print("about to calculate flux")
;   printVarSummary(vvel_section)
;   printVarSummary(dxu_section)
;   printVarSummary(dz)
   flux_vvel_section = (/ doubletofloat(vvel_section * conform(vvel_section, dxu_section, 3) * conform(vvel_section, dz, 1)) /)
   print("calculated flux_vvel_section")
   flux_vvel_section@long_name = "flux"
   flux_vvel_section@units     = "cm^3/s"
;
   print("about to calculate vgeo_td flux")
   flux_vgeo_td_section = flux_vvel_section
   flux_vgeo_td_section = (/ doubletofloat(vgeo_td_section * conform(vgeo_td_section, dxu_section, 3) * conform(vgeo_td_section, dz, 1)) /)
   print("calculated flux_vgeo_td_section")
;   printVarSummary(flux_vvel_section)
;   Dimensions and sizes:    [time | 120] x [z_t | 60] x [ulat | 1] x [ulon | 73]
;
   print("about to calculate vgeo_bu flux")
   flux_vgeo_bu_section = flux_vvel_section
   flux_vgeo_bu_section = (/ doubletofloat(vgeo_bu_section * conform(vgeo_bu_section, dxu_section, 3) * conform(vgeo_bu_section, dz, 1)) /)
   print("calculated flux_vgeo_bu_section")
;   printVarSummary(flux_vvel_section)
;   Dimensions and sizes:    [time | 120] x [z_t | 60] x [ulat | 1] x [ulon | 73]
;-------------------------------------------------------------------------------------------------------------------------------------
; Calculate transport (flux) over the geostrophic domain using model velocities and  geostosphic velocities (td and bu)  everywhere 
   ind_lon_loc_trans_w = ind_lon_loc_trans_gw
   ind_lon_loc_trans_e = ind_lon_loc_trans_ge
;
   flux_move_vvel = dim_sum_Wrap(dim_sum_Wrap(dim_sum_Wrap(flux_vvel_section(:,ind_z_loc_lower:ind_z_loc_upper,:,ind_lon_loc_trans_w:ind_lon_loc_trans_e))))
   flux_move_vvel_sv = flux_move_vvel
   flux_move_vvel_sv = (/ cm3pstosv*flux_move_vvel /)
   flux_move_vvel_sv@units = "Sv"
;
   flux_move_vgeo_td = dim_sum_Wrap(dim_sum_Wrap(dim_sum_Wrap(flux_vgeo_td_section(:,ind_z_loc_lower:ind_z_loc_upper,:,ind_lon_loc_trans_w:ind_lon_loc_trans_e))))
   flux_move_vgeo_td_sv = flux_move_vgeo_td
   flux_move_vgeo_td_sv = (/ cm3pstosv*flux_move_vgeo_td /)
   flux_move_vgeo_td_sv@units = "Sv"
   flux_move_vgeo_td_sv@units = "Sv"
;
   flux_move_vgeo_bu = dim_sum_Wrap(dim_sum_Wrap(dim_sum_Wrap(flux_vgeo_bu_section(:,ind_z_loc_lower:ind_z_loc_upper,:,ind_lon_loc_trans_w:ind_lon_loc_trans_e))))
   flux_move_vgeo_bu_sv = flux_move_vgeo_bu
   flux_move_vgeo_bu_sv = (/ cm3pstosv*flux_move_vgeo_bu /)
   flux_move_vgeo_bu_sv@units = "Sv"
   flux_move_vgeo_bu_sv@units = "Sv"
;
;   printVarSummary(flux_move_vgeo_td_sv)
;----------
; Calculate transport (flux) over the domain using model velocities and  geostosphic velocities (td and bu)  just along the western part of the domain,
; to the west of where the geostrophic section calculation will be done
; Calculate transport (flux) over the domain using model velocities and  geostosphic velocities (td and bu)  everywhere 
;
   ind_lon_loc_trans_w = ind_lon_loc_trans_cw
   ind_lon_loc_trans_e = ind_lon_loc_trans_ce
;
;   print("ind_lon_loc_trans_cw: " + ind_lon_loc_trans_w)
;   print("ind_lon_loc_trans_ce: " + ind_lon_loc_trans_e)
;   printVarSummary( flux_vvel_section(:,ind_z_loc_lower:ind_z_loc_upper,:,ind_lon_loc_trans_w:ind_lon_loc_trans_e)  )
;
   flux_westwedge_vvel = dim_sum_Wrap(dim_sum_Wrap(dim_sum_Wrap(flux_vvel_section(:,ind_z_loc_lower:ind_z_loc_upper,:,ind_lon_loc_trans_w:ind_lon_loc_trans_e))))
   flux_westwedge_vvel_sv = flux_westwedge_vvel
   flux_westwedge_vvel_sv = (/ cm3pstosv*flux_westwedge_vvel /)
   flux_westwedge_vvel_sv@units = "Sv"
;
;   print("flux_westwedge_vvel_sv: " + flux_westwedge_vvel_sv)
;   exit
;
   flux_westwedge_vgeo_td = dim_sum_Wrap(dim_sum_Wrap(dim_sum_Wrap(flux_vgeo_td_section(:,ind_z_loc_lower:ind_z_loc_upper,:,ind_lon_loc_trans_w:ind_lon_loc_trans_e))))
   flux_westwedge_vgeo_td_sv = flux_westwedge_vgeo_td
   flux_westwedge_vgeo_td_sv = (/ cm3pstosv*flux_westwedge_vgeo_td /)
   flux_westwedge_vgeo_td_sv@units = "Sv"
   flux_westwedge_vgeo_td_sv@units = "Sv"
;
   flux_westwedge_vgeo_bu = dim_sum_Wrap(dim_sum_Wrap(dim_sum_Wrap(flux_vgeo_bu_section(:,ind_z_loc_lower:ind_z_loc_upper,:,ind_lon_loc_trans_w:ind_lon_loc_trans_e))))
   flux_westwedge_vgeo_bu_sv = flux_westwedge_vgeo_bu
   flux_westwedge_vgeo_bu_sv = (/ cm3pstosv*flux_westwedge_vgeo_bu /)
   flux_westwedge_vgeo_bu_sv@units = "Sv"
   flux_westwedge_vgeo_bu_sv@units = "Sv"
;
;   printVarSummary(flux_westwedge_vgeo_td)
;----------------------------------------------------------------------------------------------------------------------
; calculate the transport (flux) using top down geostrophic velocity calculated using pressure at the west and east coordinates 
; (lon_loc_trans_gw, lon_loc_trans_ge) for the LR case. 
;----------
; 
; 18 Oct 19 Use one more point to the east for the geostrophic calculations
   ind_lon_loc_trans_w = ind_lon_loc_trans_gw
   ind_lon_loc_trans_e = ind_lon_loc_trans_ge + 1
;
;   printVarSummary(pres_section_avg(:,:,:,ind_lon_loc_trans_w))
;   printVarSummary(pres_section_avg(:,:,:,ind_lon_loc_trans_e))
   dpdx_move_bounds = pa2dbar*(pres_section_avg(:,:,:,ind_lon_loc_trans_e) - pres_section_avg(:,:,:,ind_lon_loc_trans_w))/(dim_sum(dxu_section(ind_lon_loc_trans_w:ind_lon_loc_trans_e) )/100. ) ; (pres_east - pres_west) / sum of dx's
;                                                                                                                                                                                      ; dbar to Pa in the numerator and cm to m in the denominator 
   rho_move_bounds  = (rho_section_avg(:,:,:,ind_lon_loc_trans_e) + rho_section_avg(:,:,:,ind_lon_loc_trans_w))/2.                                                            ; (rho_east + rho_west) / 2.
   vgeo_td_move_bounds = (/  100.*(1./rho_move_bounds)*(1./f) * dpdx_move_bounds /)                                                                                       ; meridional geostrophic velocity in units of cm s-1 to compare with vvel
;
;   printVarSummary(pres_section_avg)
;   printVarSummary(vgeo_td_move_bounds)
;
   copy_VarCoords_1(pres_section_avg, vgeo_td_move_bounds)                                                                                                                ; copy over some coordinate names and variables

;
; 18 Oct 19 to compare, calculate the average vvel, vgeo_td, vgeo_bu over the same move region, but use one lest point
; It doesn't seem to make much difference
   ind_lon_loc_trans_w = ind_lon_loc_trans_gw
   ind_lon_loc_trans_e = ind_lon_loc_trans_ge
;
   vvel_move_bounds_avg =  dim_avg_Wrap(vvel_section(:,:,:,ind_lon_loc_trans_w:ind_lon_loc_trans_e))
   vgeo_td_move_bounds_avg =  dim_avg_Wrap(vgeo_td_section(:,:,:,ind_lon_loc_trans_w:ind_lon_loc_trans_e))
   vgeo_bu_move_bounds_avg =  dim_avg_Wrap(vgeo_bu_section(:,:,:,ind_lon_loc_trans_w:ind_lon_loc_trans_e))
;
;   printVarSummary(vgeo_td_move_bounds)
;   printVarSummary(vvel_move_bounds_avg)
;  vgeo_td_move_bounds:  Dimensions and sizes:    [time | 120] x [z_t | 60] x [ulat | 1]
;----------
; 18 Oct 19 Use one more point to the east for the geostrophic calculations
   ind_lon_loc_trans_w = ind_lon_loc_trans_gw
   ind_lon_loc_trans_e = ind_lon_loc_trans_ge + 1
;
   flux_move_vgeo_td_bounds_sv = cm3pstosv*dim_sum(dim_sum(vgeo_td_move_bounds(:,ind_z_loc_upper:ind_z_loc_lower,:) * dim_sum(dxu_section(ind_lon_loc_trans_w:ind_lon_loc_trans_e)) * conform(vgeo_td_move_bounds(:,ind_z_loc_upper:ind_z_loc_lower,:),dz(ind_z_loc_upper:ind_z_loc_lower),1) ) )
   flux_move_vgeo_td_bounds_sv!0 = "time"
   flux_move_vgeo_td_bounds_sv&time = time
   flux_move_vgeo_td_bounds_sv@long_name = flux_move_vvel_sv@long_name
   flux_move_vgeo_td_bounds_sv@units     = flux_move_vvel_sv@units
;   printVarSummary(flux_move_vgeo_td_bounds_sv)
;--------------------
; calculate the transport (flux) using bottom up geostrophic velocity calculated using pressure at the west and east coordinates 
; (lon_loc_trans_gw, lon_loc_trans_ge) for the LR case. 
;----------
   drhodx_move_bounds = (rho_section_avg(:,:,:,ind_lon_loc_trans_e) - rho_section_avg(:,:,:,ind_lon_loc_trans_w))/(dim_sum(dxu_section(ind_lon_loc_trans_w:ind_lon_loc_trans_e) ) ) ; (rho_east - rho_west) / sum of dx's
;------------
; calculate the vertical gradient of the meridional geostrophic velocity
;------------
   rho_move_bounds  = (rho_section_avg(:,:,:,ind_lon_loc_trans_e) + rho_section_avg(:,:,:,ind_lon_loc_trans_w))/2.                                                         ; (rho_east + rho_west) / 2.
   dvgeodz_move_bounds = drhodx_move_bounds
   dvgeodz_move_bounds = (/  (-100.*grav/rho_move_bounds)*(1./f) * drhodx_move_bounds /)                                                                               ; in units of cm s-1 to compare with vvel
   copy_VarCoords(coord_template_section_u, dvgeodz_move_bounds)
   copy_VarCoords(coord_template_section_u, drhodx_move_bounds)
;   printVarSummary(dvgeodz_move_bounds)
; Dimensions and sizes:    [time | 120] x [z_t | 60] x [ulat | 1]
;-----------
; Now bottom to top taking vgeo = 0 at the bottom
;-----------
  vgeo_bu_move_bounds = dvgeodz_move_bounds
; put zeros where there are not missing values
   vgeo_bu_move_bounds = where(ismissing(vgeo_bu_move_bounds),vgeo_bu_move_bounds@_FillValue,0.)
; replace "section" with "move_bounds" below to "stopped" and remove a dimension
; Integrate from the bottom up.  Integrate only when the layer below is not a missing value. Use 1/2 dvdz from below and 1/2 dvdz at the same level.
   do nz = ndepth - 2, 0, 1
      vgeo_bu_move_bounds(:,nz,:) = where( ismissing(vgeo_bu_move_bounds(:,nz+1,:)), vgeo_bu_move_bounds(:,nz,:), vgeo_bu_move_bounds(:,nz+1,:) + (0.5*dvgeodz_move_bounds(:,nz+1,:)*dz(nz+1) + 0.5*dvgeodz_move_bounds(:,nz,:)*dz(nz)))
   end do
;   printVarSummary(vgeo_bu_move_bounds)
;----------
   flux_move_vgeo_bu_bounds_sv = cm3pstosv*dim_sum(dim_sum(vgeo_bu_move_bounds(:,ind_z_loc_upper:ind_z_loc_lower,:) * dim_sum(dxu_section(ind_lon_loc_trans_w:ind_lon_loc_trans_e)) * conform(vgeo_bu_move_bounds(:,ind_z_loc_upper:ind_z_loc_lower,:),dz(ind_z_loc_upper:ind_z_loc_lower),1) ) )
   flux_move_vgeo_bu_bounds_sv!0 = "time"
   flux_move_vgeo_bu_bounds_sv&time = time
   flux_move_vgeo_bu_bounds_sv@long_name = flux_move_vvel_sv@long_name
   flux_move_vgeo_bu_bounds_sv@units     = flux_move_vvel_sv@units
;   printVarSummary(flux_move_vgeo_bu_bounds_sv)
;==================
; HERE
;   file_out = "vvel_comparison.nc"
;   system("rm " + file_out)
;   f_out = addfile(file_out, "c")
;   ind_lon_loc_w = ind_nearest_coord(lon_loc_trans_gw, ulon_section,0)
;   ind_lon_loc_e = ind_nearest_coord(lon_loc_trans_ge, ulon_section,0)
;   print("ind_lon_loc_w: " + ind_lon_loc_w + "   ind_lon_loc_e: " + ind_lon_loc_e)
;
;   f_out->vgeo_bu_subsection = vgeo_bu_section(:,:,:,ind_lon_loc_trans_gw:ind_lon_loc_trans_ge)
;   f_out->vgeo_td_subsection = vgeo_td_section(:,:,:,ind_lon_loc_trans_gw:ind_lon_loc_trans_ge)
;   f_out->vgeo_bu_subsection = vgeo_bu_section(:,:,:,ind_lon_loc_w:ind_lon_loc_e)
;   f_out->vgeo_td_subsection = vgeo_td_section(:,:,:,ind_lon_loc_w:ind_lon_loc_e)
;
;   f_out->vgeo_bu_move_bounds = vgeo_bu_move_bounds
;   f_out->vgeo_td_move_bounds = vgeo_td_move_bounds
;
;   delete(ind_lon_loc_w)
;   delete(ind_lon_loc_e)
;
;   print(vgeo_bu_section&ulon-ulon_section)
; END HERE
;==================
;
;----------
; take a look...
;----------
   if(1 .eq. 1) then
      system("rm " + dir_out + "/vgeo_tdbu_move_bounds.nc")
      f_out = addfile(dir_out + "/vgeo_tdbu_move_bounds.nc", "c")
      filedimdef(f_out,"time",-1,True)
;
      flux_move_vgeo_td_m_flux_move_vgeo_bu = flux_move_vgeo_td
      flux_move_vgeo_td_m_flux_move_vgeo_bu = (/ flux_move_vgeo_td - flux_move_vgeo_bu /)
;
      f_out->dpdx_move_bounds = dpdx_move_bounds
      f_out->drhodx_move_bounds = drhodx_move_bounds
      f_out->vgeo_td_move_bounds = vgeo_td_move_bounds
      f_out->vgeo_bu_move_bounds = vgeo_bu_move_bounds
   end if
;
   if(1 .eq. 1) then
      system("rm " + dir_out + "/fluxes.nc")
      f_out = addfile(dir_out + "/fluxes.nc", "c")
      filedimdef(f_out,"time",-1,True)
;
      f_out->flux_westwedge_vvel_sv=flux_westwedge_vvel_sv
      f_out->flux_westwedge_vgeo_bu_sv=flux_westwedge_vvel_sv
      f_out->flux_westwedge_vvel_sv=flux_westwedge_vvel_sv
      f_out->flux_move_vgeo_td_bounds_sv = flux_move_vgeo_td_bounds_sv
      f_out->flux_move_vgeo_bu_bounds_sv = flux_move_vgeo_bu_bounds_sv
      f_out->flux_move_vvel_sv = flux_move_vvel_sv
      f_out->flux_move_vgeo_td_sv = flux_move_vgeo_td_sv
      f_out->flux_move_vgeo_bu_sv = flux_move_vgeo_bu_sv
      f_out->flux_move_vgeo_td_m_flux_move_vgeo_bu = flux_move_vgeo_td_m_flux_move_vgeo_bu
   end if
;
   if(1 .eq. 1) then
      vgeo_td_m_bu_move_bounds_avg = vgeo_td_move_bounds_avg
      vgeo_td_m_bu_move_bounds_avg = (/ vgeo_td_move_bounds_avg - vgeo_bu_move_bounds_avg /)
;
      vvel_m_td_move_bounds_avg = vgeo_td_move_bounds_avg
      vvel_m_td_move_bounds_avg = (/ vvel_move_bounds_avg - vgeo_td_move_bounds_avg /)
;
      vvel_m_bu_move_bounds_avg = vgeo_bu_move_bounds_avg
      vvel_m_bu_move_bounds_avg = (/ vvel_move_bounds_avg - vgeo_bu_move_bounds_avg /)
;      
      system("rm " + dir_out + "/move_bounds_averages.nc")
      f_out = addfile(dir_out + "/move_bounds_averages.nc", "c")
      filedimdef(f_out,"time",-1,True)
;
; vgeo at all longitudes and averaged over the boundary regions to compare with the ones without _avg
      f_out->vvel_move_bounds_avg = vvel_move_bounds_avg
      f_out->vgeo_td_move_bounds_avg = vgeo_td_move_bounds_avg
      f_out->vgeo_bu_move_bounds_avg = vgeo_bu_move_bounds_avg
;
; vgeo calculated using the boundary values only
      f_out->vgeo_td_move_bounds = vgeo_td_move_bounds
      f_out->vgeo_bu_move_bounds = vgeo_bu_move_bounds
;
      f_out->vgeo_td_m_bu_move_bounds_avg = vgeo_td_m_bu_move_bounds_avg
      f_out->vvel_m_td_move_bounds_avg = vvel_m_td_move_bounds_avg
      f_out->vvel_m_bu_move_bounds_avg = vvel_m_bu_move_bounds_avg
   end if
;
;-------------------------------------------------------------------------------------------------------------------------------------
   script_name = get_script_prefix_name()
   time_date = systemfunc("date")
   file_out = script_name + "_" + cres + ".nc"
   system("rm " + file_out)
   f_out = addfile(file_out,"c")

   fatt = True
   fatt@creation_date = systemfunc("date")
   script_name = get_script_name()
   pwd = systemfunc("echo $PWD")
   fatt@creation_script = pwd + "/" + script_name
   fatt@author = "Bob Tomas"
;
   fatt@lat_loc_in = lat_loc_in
   fatt@lon_loc_trans_cw = lon_loc_trans_cw
   fatt@lon_loc_trans_ce = lon_loc_trans_ce
   fatt@lon_loc_trans_gw = lon_loc_trans_gw
   fatt@lon_loc_trans_ge = lon_loc_trans_ge
   fatt@z_loc_lower = z_loc_lower
   fatt@z_loc_upper = z_loc_upper
   fileattdef( f_out, fatt )
;
   f_out->time      = time
   f_out->time_year = time_year
   f_out->flux_move_vvel_sv = flux_move_vvel_sv
   f_out->flux_move_vgeo_td_sv = flux_move_vgeo_td_sv
   f_out->flux_move_vgeo_td_bounds_sv = flux_move_vgeo_td_bounds_sv
   f_out->flux_move_vgeo_bu_sv = flux_move_vgeo_bu_sv
   f_out->flux_move_vgeo_bu_bounds_sv = flux_move_vgeo_bu_bounds_sv
   f_out->flux_westwedge_vvel_sv = flux_westwedge_vvel_sv
;
   lat_loc_in!0      = "ncl_scalar"
   lon_loc_trans_cw!0 = "ncl_scalar"
   lon_loc_trans_ce!0 = "ncl_scalar"
   lon_loc_trans_gw!0 = "ncl_scalar"
   lon_loc_trans_ge!0 = "ncl_scalar"
   z_loc_lower!0     = "ncl_scalar"
   z_loc_upper!0     = "ncl_scalar"
;
   f_out->lat_loc_in = lat_loc_in
   f_out->lon_loc_trans_cw = lon_loc_trans_cw
   f_out->lon_loc_trans_ce = lon_loc_trans_ce
   f_out->lon_loc_trans_gw = lon_loc_trans_gw
   f_out->lon_loc_trans_ge = lon_loc_trans_ge
   f_out->z_loc_lower = z_loc_lower
   f_out->z_loc_upper = z_loc_upper
;======================================================================================================================
   print("Plotting")
   script_name = get_script_prefix_name()
   time_date = systemfunc("date")
   file_out = script_name + "_" + cres + ".png"
   wks_type = "png"
   wks_type@wkWidth = 1024   ; 1024 is default
   wks_type@wkHeight = 1024  ; 1024 is default
   wks = gsn_open_wks(wks_type,file_out)
;
   res=True
   res@gsnDraw = False
  res@gsnFrame = False
;
   res@trYReverse = False
   res@vpWidthF   = 0.6
   res@vpHeightF  = 0.2
;
   res@tiXAxisString = "Year"
   res@tiYAxisString = "Transport(Sv)"
   res@tiMainString = "POP " + cres + " Ocean Volume Transport Across MOVE Section (16~S~o~N~)"
;
; add a legend
   res@pmLegendDisplayMode    = "Always"              ; turn on legend
;
   res@pmLegendSide           = "Bottom"                 ; Change location of
   res@pmLegendParallelPosF   =  0.5                  ; move units right
   res@pmLegendOrthogonalPosF = -0.00                 ; more neg = down
;
   res@pmLegendWidthF         = 0.25                  ; Change width and
   res@pmLegendHeightF        = 0.15                  ; height of legend.
   res@lgLabelFontHeightF     = .012                  ; change font height
   res@lgPerimOn              = False                 ; no box around
;
   nplots = 3
   plot   = new(1, graphic)
   plot_data = new( (/nplots,dimsizes(flux_move_vvel_sv)/), typeof(flux_move_vvel_sv) )
;
   plot_data!1 = "time"
   plot_data&time = time_year
   plot_data(0,:) = flux_move_vvel_sv
;   plot_data(1,:) = (/ dble2flt(flux_move_vgeo_td_sv) /)
;   plot_data(2,:) = (/ dble2flt(flux_move_vgeo_bu_sv) /)
;
   plot_data(1,:) = (/ dble2flt(flux_move_vgeo_td_sv + flux_westwedge_vvel_sv) /)
   plot_data(2,:) = (/ dble2flt(flux_move_vgeo_bu_sv + flux_westwedge_vvel_sv) /)
;
   res@xyLineThicknesses := (/ 5.0, 5.0, 5.0 /) ; (/ 2.0, 5.0, 5.0, 5.0 /)
   res@xyDashPatterns := (/ 0,0,0 /) ;    2, 12
   res@xyLineColors   := (/  "black", "red","blue" /)
   res@trXMinF = 2000.
   res@trXMaxF = 2017.
;   res@gsnXRefLine = 0.
;
;------------------------------------------------------------------------------------------

   sd_plot_data0 = sprintf("%3.1f", dim_stddev(plot_data(0,:)))
   mn_plot_data0 = sprintf("%3.1f", dim_avg(plot_data(0,:)))
   sd_plot_data1 = sprintf("%3.1f", dim_stddev(plot_data(1,:)))
   mn_plot_data1 = sprintf("%3.1f", dim_avg(plot_data(1,:)))
   sd_plot_data2 = sprintf("%3.1f", dim_stddev(plot_data(2,:)))
   mn_plot_data2 = sprintf("%3.1f", dim_avg(plot_data(2,:)))

;
;
   res@xyExplicitLegendLabels := (/ "VVEL~B~Domain Avg~N~   Mean: " + mn_plot_data0 + " SD: " + sd_plot_data0,   \
                                    "VGEO~B~TD~N~   Mean: " + mn_plot_data1 + " SD: " + sd_plot_data1 ,          \
                                    "VGEO~B~BU~N~   Mean: " + mn_plot_data2 + " SD: " + sd_plot_data2 /)
;
   plot  = gsn_csm_xy(wks, time_year, plot_data, res)
;----------
;
   pres = True
   gsn_panel(wks,plot,(/1,1/),pres)             ;
;-------------------------------------------------------------------------------------------------------------------------------------
end
