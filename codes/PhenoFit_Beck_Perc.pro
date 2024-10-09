PRO PhenoFit_Beck_Perc_V1
  ;**********************************************************************************************;
  ; This script is to extract leaf phenology from PlanetScope-Derived normalized differential 
  ; greenness index (NDGI). The key process include in this code include:
  ; 1. clean observations that are within a defined winter season during the surface is fully
  ;    coverd my snow
  ; 2. replace pixels with ngdi values below 0.05 using background value 0.05
  ; 3. extend the temporal range of time series data to the day 1 and day 365 of the year using 
  ;    background value 0.05 and a white noise. this extension is to minimize the effects of the
  ;    curve shape during shoulder seasons on phenology extraction
  ; 4. conduct Savitzky–Golay filter to smooth the time-series NDGI data
  ; 5. apply Beck et al. (2006) method to smooth the time-series NDGI
  ; 6. extract phenoloyg (upturn day (UD), stop day (SD), downturn day (DD), and return day (RD)), 
  ;    as well as corresponding fitted NDGI at the four days, and annual maximum and minimum 
  ;    NDGI values

  ; Critial notes:
  ; 1. the code uses output NDGI raster file processed using "Generate_Time_Series_VIs_Multiple.R"
  ;    or "Generate_Time_Series_VIs_Single.R", as well corresponding day of year of each band layer
  ;    stored in a csv file. 
  ; 2. this code require MPFIT function and resistent_main function. please download and place the
  ;    corresponding code in ENVI/IDL lib
  
  ;           Latest updated on 2024/10/01 by Daryl Yang <<yangd@ornl.gov>>
  ;**********************************************************************************************;

  ;********************************* SET OUPUT DIRECTORY ****************************************;
  ; please define output directory for storing leaf phenological maps and fitted NDVI curve data
  out_dir = ''
  file_mkdir, out_dir
  if out_dir then begin
    print, '--- Creating --- "', out_dir, '" --- Successful ---'
  endif else begin
    print, '--- Creating --- "', out_dir, '" --- Failed ---'
  endelse
  ;**********************************************************************************************;

  ;********************************* SET USER PARAMETER *****************************************;
  ; read in csv file that contains the day of year of each image layer that in the layerstack of
  ; NDGI index. This file is an output from the R script "Generate_Time_Series_VIs_Single.R" or
  ; "Generate_Time_Series_VIs_Multiple.R"
  date_stamp_orig = read_csv('', $
    header = header)
  doy_stamp = date_stamp_orig.field2
  
  ; define the range of day of year (DOY) outside of which the vegetation index will be exclude. 
  ; these paremeters to remove observations that are apparently known as snow coverage periods
  ; especially in cold region like the Arctic. if not sure about the snow cover date, please set
  ; expBEG to 1 and expEND to 365
  expBEG = 100   
  expEND = 300  
  ;**********************************************************************************************;

  ;************************************** LOAD DATA *********************************************;
  ; load in ndgi time series data 
  ndgi_series_dir = ''
  ndgi_raster = envi.OpenRaster(ndgi_series_dir)
  ndgi_img_data = ndgi_raster.GetData()
  ndgi_img_data = transpose(ndgi_img_data, [2, 1, 0]) ; this line may not be necessary depending on
                                                      ; the format of orginal data. bands should be 
                                                      ; the last 

  ; remove areas that with very low NDGI values, corresponding to non-vegetation surfaces
  mask = max(ndgi_img_data, dimension = 3)
  mask[where(mask lt 0.3)] = !VALUES.F_NAN
  mask[where(mask ge 0.3)] = 1
  ;**********************************************************************************************;


  ;********************************* Process Phenology ******************************************;
  ; the day that the time-series NDGI will be extended to using backgound value. please do not
  ; change this parameter if not necessary
  doy_start = 1
  ; the last that the time-series NDGI will be extended to using backgound value. please do not
  ; change this parameter if not necessary
  doy_end = 365
  ; define the interval for extending the time series. please do not change this parameter if 
  ; not necessary
  intv = 3

  ; get the size the orginal raster NDGI file
  dims_ndgi = size(ndgi_img_data, /dimensions)
  ; create empty matrix to store the output phenological map and fitted daily NDGI time series 
  pheno_pars = fltarr(dims_ndgi[0], dims_ndgi[1], 10)
  pheno_ndgi = fltarr(dims_ndgi[0], dims_ndgi[1], 365)
  ; process pixel-wise phenology
  for lon_col=0L, dims_ndgi[0]-1 do begin ; dims_ndgi[0]-1
    print, lon_col
    for lon_row = 0L, dims_ndgi[1]-1 do begin ;dims_ndgi[1]-1
      ; extract pixel-wise NDGI time series
      pixel_ndgi_in = reform(ndgi_img_data[lon_col, lon_row, *])
      ; fit phenology for the pixel
      pixl_pheno = phenofit(pixel_ndgi_in, doy_stamp, expBEG, expEND, doy_start, doy_end, intv, $
        pheno_par=pheno_par, fitted_ndgi=fitted_ndgi)
      pheno_pars[lon_col,lon_row,*] = pheno_par
      pheno_ndgi[lon_col,lon_row,*] = fitted_ndgi
      ; remore current data cycle to save space
      delvar, pixel_ndgi_in
    endfor
  endfor

  ; mask out areas define by low-NDGI mask and run a simple resistant_mean to remove outliers in the 
  ; phenological maps
  for i=0, 9 do begin
    band = pheno_pars[*,*,i]
    band_smoothed = simple_smooth(band, 5)
    band_smoothed = band_smoothed*mask
    pheno_pars[*,*,i] = band_smoothed
  endfor

  ; transform map back to the orignal structure same to input time-series NDGI. this line may not be
  ; necessary depending on the format of orginal data. 
  pheno_pars_out = transpose(pheno_pars, [1, 0, 2])
  pheno_ndgi_out = transpose(pheno_ndgi, [1, 0, 2])

  print, '...... saving Double Logistic results'
  envi_open_file, ndgi_series_dir, r_fid=fid
  map_info = envi_get_map_info(fid = fid)
  file_name = file_basename(ndgi_series_dir, '.tif')
  out_name = out_dir + '\' + 'Double_Logistic_Parms_Full'
  envi_write_envi_file, pheno_pars_out, out_name = out_name, map_info = map_info, $
    bnames = ['ud', 'sd', 'dd', 'rd', 'ndgi_ud', 'ndgi_sd', 'ndgi_dd', 'ndgi_rd', 'ndgimin', 'ndgimax']

  out_name = out_dir + '\' + 'Double_Logistic_Fitted_NDGI_Full'
  envi_write_envi_file, pheno_ndgi_out, out_name = out_name, map_info = map_info, $
    bnames = [1:365:1]

END

Function phenofit, pixel_vi, date_stamp, expBEG, expEND, doy_start, doy_end, intv, pheno_par=pheno_par, fitted_vi=fitted_vi
  subsnow = where(date_stamp ge expBEG and date_stamp le expEND)
  pixel_vi_subsnow = pixel_vi[subsnow]
  date_stamp_subsnow = date_stamp[subsnow]
  ; remove outliers
  pixel_vi_subsnow[where(pixel_vi_subsnow lt 0.0 or pixel_vi_subsnow gt 1)] = 0.05
  pixel_vi_subsnow[where(~finite(pixel_vi_subsnow))] = 0
  
  date_stamp = date_stamp_subsnow
  pixel_vi = pixel_vi_subsnow
  ; filter vi with sg
  num = n_elements(pixel_vi)
  if num gt 20 then begin
    ; initialize phenofit parameters
    front_extend_dates = [doy_start:(expBEG-1):intv]
    end_extend_dates = [doy_end:(expEND+1):-intv]
    end_extend_dates = end_extend_dates[sort(end_extend_dates)]
    ext_doy = [front_extend_dates, date_stamp, end_extend_dates]
    parms = [0.01, 1.0, 0.05, 0.02, 140, 260]
    yerr = randomn(seed, n_elements(ext_doy))/100
    pixel_t_vi = extend(date_stamp, pixel_vi, ext_doy, expBEG, expEND, ext_vi = ext_vi)
    pixel_vi_sgfilterred =  sgfilter(ext_vi)
    result = MPFITFUN('DoubleLogBeck', ext_doy, pixel_vi_sgfilterred, yerr, parms, yfit = yfit, $
      STATUS=ST, ERRMSG=ERR, /QUIET) ;
    ;doublog_parms[lon_col, lon_row, *] = result

    ; generate smooth fitted vi time series
    t = [1:365:1]
    if n_elements(result) eq 6 then begin
      mn = result[0]
      mx = result[1]
      rsp = result[2]
      rau = result[3]
      sos = result[4]
      eos = result[5]
      fitted_vi = mn + (mx-mn)*(1/(1+exp(-rsp*(t-sos))) + 1/(1+exp(rau*(t-eos))))
      ; determine pheno dates based on percentage
      vi_green = fitted_vi[100:213] ;213
      vi_brown = fitted_vi[213:300]

      vimax = max(fitted_vi[expBEG:expEND])
      vimin = min(fitted_vi[expBEG:expEND])
      vi25 = 0.15*(vimax-vimin) + vimin
      vi75 = 0.85*(vimax-vimin) + vimin

      viUD= vi25
      UD = where(abs(vi_green - viUD) eq min(abs(vi_green - viUD))) + 100
      if n_elements(UD) gt 1 then begin
        UD = !VALUES.F_NAN
      endif
      viSD = vi75
      SD = where(abs(vi_green - viSD) eq min(abs(vi_green - viSD))) + 100
      if n_elements(SD) gt 1 then begin
        SD = !VALUES.F_NAN
      endif
      viDD = vi75
      DD = where(abs(vi_brown - viDD) eq min(abs(vi_brown - viDD))) + 213
      if n_elements(DD) gt 1 then begin
        DD = !VALUES.F_NAN
      endif
      viRD = vi25
      RD = where(abs(vi_brown - viRD) eq min(abs(vi_brown - viRD))) + 213
      if n_elements(RD) gt 1 then begin
        RD = !VALUES.F_NAN
      endif
      ; store pheno results
      pheno_par = [UD, SD, DD, RD, viUD, viSD, viDD, viRD, vimin, vimax]
      fitted_vi = fitted_vi
    endif
    if n_elements(result) ne 6 then begin
      pheno_par = replicate(0, 10)
      fitted_vi = replicate(0, 365)
    endif
  endif
  if num le 20 then begin
    ; if the number of valid vi is less than tha half of input then skip
    pheno_par = replicate(0, 10)
    fitted_vi = replicate(0, 365)
  endif

  delvar, pixel_vi, date_stamp
End
  
Function double_logistic, t, parms
  pred_y = parms[0] + parms[1]*(1/(1+exp(parms[2]-parms[3]*t))-1/(1+exp(parms[4]-parms[5]*t)))
  return, pred_y
End

Function DoubleLogBeck, t, parms
  mn = parms[0]
  mx = parms[1]
  rsp = parms[2]
  rau = parms[3]
  sos = parms[4]
  eos = parms[5]
  pred_y = mn + (mx-mn)*(1/(1+exp(-rsp*(t-sos))) + 1/(1+exp(rau*(t-eos))))
  return, pred_y
End

Function simple_smooth, pheno_pars, windsize
  dims = size(pheno_pars, /dimensions)
  ns = dims[0]
  nl = dims[1]
  num_point_sample= ns mod windsize
  num_point_line= ns mod windsize
  if num_point_sample ne 0 then begin
    num_count_sample=fix(ns/windsize)+1
  endif else begin
    num_count_sample=fix(ns/windsize)
  endelse

  if num_point_line ne 0 then begin
    num_count_line=fix(nl/windsize)+1
  endif else begin
    num_count_line=fix(nl/windsize)
  endelse

  pheno_par_corrected=fltarr(dims)

  for i=0L,num_count_sample-1 do begin
    for j=0L,num_count_line-1 do begin
      end_point_sample=(i+1)*windsize-1
      end_point_line=(j+1)*windsize-1

      if end_point_sample le ns-1 and end_point_line le nl-1 then begin
        cut_img=pheno_pars[[i*windsize]:[(i+1)*windsize-1],[j*windsize]:[(j+1)*windsize-1]]
        resistant_mean, cut_img, 1, meanv, meansig
        cut_img[where(abs(cut_img-meanv) gt 3*stddev(cut_img))] = meanv
        pheno_par_corrected[[i*windsize]:[(i+1)*windsize-1],[j*windsize]:[(j+1)*windsize-1]] = cut_img
      endif
      if end_point_sample le ns-1 and end_point_line gt nl-1 then begin
        cut_img=pheno_pars[[i*windsize]:[(i+1)*windsize-1],[j*windsize]:[nl-1]]
        resistant_mean, cut_img, 1, meanv, meansig
        cut_img[where(abs(cut_img-meanv) gt 3*stddev(cut_img))] = meanv
        pheno_par_corrected[[i*windsize]:[(i+1)*windsize-1],[j*windsize]:[nl-1]] = cut_img
      endif
      if end_point_sample gt ns-1 and end_point_line le nl-1 then begin
        cut_img=pheno_pars[[i*windsize]:[ns-1],[j*windsize]:[(j+1)*windsize-1]]
        resistant_mean, cut_img, 1, meanv, meansig
        cut_img[where(abs(cut_img-meanv) gt 3*stddev(cut_img))] = meanv
        pheno_par_corrected[[i*windsize]:[ns-1],[j*windsize]:[(j+1)*windsize-1]] = cut_img
      endif
      if end_point_sample gt ns-1 and end_point_line gt nl-1 then begin
        cut_img=pheno_pars[[i*windsize]:[ns-1],[j*windsize]:[nl-1]]
        resistant_mean, cut_img, 1, meanv, meansig
        cut_img[where(abs(cut_img-meanv) gt 3*stddev(cut_img))] = meanv
        pheno_par_corrected[[i*windsize]:[ns-1],[j*windsize]:[nl-1]] = cut_img
      endif
    endfor
  endfor
  return, pheno_par_corrected
End


Function fill_up, vector_in, maxNDVI               ;remove the cloudy values

  num_elements=n_elements(vector_in)

  ;remove  continuous 0
  in_pixel = 0
  if vector_in[in_pixel] EQ 0 then begin
    pixel_start = in_pixel
    while vector_in[in_pixel] EQ 0 && in_pixel LT num_elements -2 do begin
      in_pixel = in_pixel +1
    endwhile
    pixel_end = in_pixel
    vector_in[pixel_start:pixel_end-1] = vector_in[pixel_end]
  endif

  in_pixel = num_elements-1
  if vector_in[in_pixel] EQ 0 then begin
    pixel_end = in_pixel
    while vector_in[in_pixel] EQ 0 && in_pixel GT 2 do begin
      in_pixel = in_pixel -1
    endwhile
    pixel_start = in_pixel
    vector_in[pixel_start+1:pixel_end] = vector_in[pixel_start]
  endif
  ;-----------
  in_pixel = 1
  while in_pixel LT num_elements - 2 do begin
    if vector_in[in_pixel] EQ 0 then begin
      pixel_start = in_pixel
      num_pixel = 1
      while vector_in[in_pixel] EQ 0 && in_pixel LT num_elements -2 do begin
        in_pixel = in_pixel +1
        num_pixel = num_pixel +1
      endwhile
      pixel_end = in_pixel
      temp = (vector_in[pixel_end] - vector_in[pixel_start-1]) / num_pixel
      for pixel = 0, num_pixel-2 do begin
        vector_in[pixel_start + pixel] = vector_in[pixel_start -1] + (pixel+1)*temp
      endfor
    endif
    in_pixel = in_pixel + 1
  endwhile

  in_pixel = 1
  while  in_pixel LT num_elements -2 do begin
    if (vector_in[in_pixel] - vector_in[in_pixel -1]) GE 0.2*maxNDVI $
      && (vector_in[in_pixel] - vector_in[in_pixel +1]) GE 0.2*maxNDVI then begin
      vector_in[in_pixel] = (vector_in[in_pixel -1] + vector_in[in_pixel +1]) /2.0
      in_pixel = in_pixel +2
    endif else begin
      in_pixel = in_pixel +1
    endelse
  endwhile

  return, vector_in

End

Function sgfilter,vector_in                 ;S-G filter

  num_elements = n_elements(vector_in)
  ; The first Savitzky-Golay fitting
  vector_in=reform(vector_in,num_elements)                         ; num_elements is the number of values of time-series
  savgolFilter = SAVGOL(5,5,0,2)                          ;set the window width(4,4) and degree (2) for computing trend curve
  rst = CONVOL(vector_in, savgolFilter, /EDGE_TRUNCATE)

  ; Calculate the threshold for loop control, so that the fit is maximize
  gdis = 0.0
  fl = IntARR(num_elements)

  for i =0,(num_elements-1) do begin
    fl[i] = (vector_in[i] ge rst[i])
    gdis = gdis + fl[i]*abs(vector_in[i]-rst[i])*abs(vector_in[i]-rst[i])
  endfor

  ra4 = fltARR(num_elements)
  pre = fltARR(num_elements)

  ormax = gdis
  num   = 0

  loop_times = 0l
  while (gdis le ormax) && loop_times LT 10 do begin
    loop_times = loop_times +1
    for i =0,(num_elements-1) do begin
      ra4[i] = (vector_in[i] ge rst[i]) ? vector_in[i] : rst[i]
      pre[i] = rst[i]
    endfor

    ; The Savitzky-Golay fitting
    savgolFilter = SAVGOL(5, 5, 0, 6)        ;set the window width(4,4) and degree (6) for repetition
    rst = CONVOL(ra4, savgolFilter, /EDGE_TRUNCATE)
    ormax = gdis
    ; Calculate the fitting-effect index
    gdis = 0.0
    for i =0,(num_elements-1) do begin
      gdis = gdis + fl[i]*abs(vector_in[i]-rst[i])*abs(vector_in[i]-rst[i])
    endfor
  endwhile

  if loop_times GE 1000 then begin
    print, 'loop times is: ', loop_times
  endif

  return, pre

End ; of function sgfilter


Function getfilename, pathname
  ;**********************************************************************************************;
  ; This function is for extracting the file name from a file directory
  ;**********************************************************************************************;
  ; get the filename in a directory
  IF (N_PARAMS() NE 1) THEN RETURN, ''
  idx = STRPOS(pathname,'\', /REVERSE_SEARCH)
  filename = STRMID(pathname, idx + 1)
  RETURN, filename
End

Function extend, time_stamp, vi, ext_doy, expBEG, expEND, ext_vi = ext_vi
  ; define the first date extend to
  ;vi = vi[where(time_stamp ge expBEG and time_stamp le expEND)]
  ;time_stamp = time_stamp[where(time_stamp ge expBEG and time_stamp le expEND)]

  n_front_dates = n_elements(where(ext_doy lt time_stamp[0]))
  front_extend_vi = replicate(min(vi), n_front_dates) + randomu(seed, n_front_dates)*0.02 ;
  n_end_dates = n_elements(where(ext_doy gt time_stamp[n_elements(time_stamp)-1]))
  end_extend_vi = replicate(min(vi), n_end_dates) + randomu(seed, n_end_dates)*0.02 ;min(vi)

  ext_vi = [front_extend_vi, vi, end_extend_vi]
End

