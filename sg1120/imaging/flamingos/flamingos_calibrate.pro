function read_twomass, twomassfile

    data = djs_readlines(twomassfile,nhead=4,head=head)

    keep = where(strmatch(data,'#*') eq 0B)
    data = data[keep]
    ndata = n_elements(data)
    
    twomass = {star: '', ra: 0.0D, dec: 0.0D, mag: -999.0D, magerr: -999.0D}
    twomass = replicate(twomass,ndata)

    for ii = 0L, ndata-1L do begin
       twomass[ii].star   = strmid(data[ii],36,18)
       twomass[ii].ra     = strmid(data[ii],0,10)
       twomass[ii].dec    = strmid(data[ii],11,10)
       twomass[ii].mag    = strmid(data[ii],114,6)
       magerr = strmid(data[ii],120,6)
       if strmatch(magerr,'*---*') then twomass[ii].magerr = -999.0 else $
           twomass[ii].magerr = strmid(data[ii],120,6)
    endfor

    twomass = twomass[where((twomass.mag gt -900.0) and (twomass.magerr gt -900.0))]
    
return, twomass
end

pro flamingos_calibrate, sextractor=sextractor, findtwomass=findtwomass, postscript=postscript
; jm07jan27nyu - calibrate the FLAMINGOS mosaic 
;   11:20:22.11 -12:04:09.5, 16.8'x20.5'

    datapath = flamingos_path(/mosaics)
    sexpath = sg1120_path(/sex)
    
    sexcatfile = datapath+'sg1120_ks.cat'
    mosaic_file = datapath+'sg1120_ks.fits'
    mosaic_weightfile = datapath+'sg1120_ks.weight.fits'
       
    twomassfile = datapath+'twomass.dat'
    if keyword_set(findtwomass) then begin

       extast, headfits(mosaic_file), astr
       radec = strtrim(im_dec2hms(astr.crval[0]/15.0,/colon),2)+$
         strtrim(im_dec2hms(astr.crval[1],/colon),2)
       spawn, 'find2mass -r 20 -m 1000 -c "'+radec+'" > twomass.dat'

    endif

    twomass = read_twomass(twomassfile)
;   struct_print, twomass

; run SE on the mosaic

    if keyword_set(sextractor) then begin

       sexconfig = sexpath+'default.sex'
       sexparam = sexpath+'sg1120.sex.param'
       sexconv = sexpath+'default.conv'
       sexnnw = sexpath+'default.nnw'

       spawn, 'sex '+mosaic_file+' -c '+sexconfig+' -CATALOG_NAME '+sexcatfile+' -CATALOG_TYPE FITS_LDAC'+$
         ' -PARAMETERS_NAME '+sexparam+' -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+$
         ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+mosaic_weightfile+' -WEIGHT_THRESH 0'+$
         ' -BACK_TYPE MANUAL -BACK_VALUE 0.0'+$
         ' -MEMORY_BUFSIZE 8192 -VERBOSE_TYPE NORMAL -NTHREADS 4', /sh

    endif

; read the SE catalog and match

    hdr = headfits(mosaic_file)
    extast, hdr, astr
    
    cat = mrdfits(sexcatfile,2,/silent)
    
    xy2ad, cat.xwin_image, cat.ywin_image, astr, aa, dd
    spherematch, aa, dd, twomass.ra, twomass.dec,1.5/3600.0, $
      catmatch, twomassmatch, distance12, maxmatch=1
    sortindex = uniq(twomassmatch,sort(twomassmatch))
    twomassmatch = twomassmatch[sortindex]
    catmatch = catmatch[sortindex]
    ntwomassmatch = n_elements(twomassmatch)

    splog, 'Matched '+string(ntwomassmatch,format='(I0)')+' TWOMASS stars.'
;   struct_print, twomass[twomassmatch]

    result = {x: -1.0D, y: -1.0D, ra: -1.0D, dec: -1.0D, mag: 0.0D, magerr: 0.0D, $
      mag_auto: 0.0D, magerr_auto: 0.0D, flux_auto: 0.0D, fluxerr_auto: 0.0D, class_star: 0.0}
    result = replicate(result,ntwomassmatch)
    
    result.x           = cat[catmatch].xwin_image
    result.y           = cat[catmatch].ywin_image
    result.ra          = cat[catmatch].xwin_world
    result.dec         = cat[catmatch].ywin_world
    result.mag         = twomass[twomassmatch].mag
    result.magerr      = twomass[twomassmatch].magerr
    result.mag_auto    = cat[catmatch].mag_auto
    result.magerr_auto = cat[catmatch].magerr_auto
    result.flux_auto    = cat[catmatch].flux_auto
    result.fluxerr_auto = cat[catmatch].fluxerr_auto
    result.class_star   = cat[catmatch].class_star

    struct_print, result

; compute the zero-point

    keep = where(result.class_star gt 0.15 and result.mag_auto gt 11.5)
    djs_iterstat, result[keep].mag-result[keep].mag_auto, sigrej=3.0, median=zpt, sigma=zpt_err

    coeff = linfit(result[keep].mag_auto,result[keep].mag)

    splog, coeff
    splog, 'Zero-point = ', zpt, zpt_err

    struct_print, result, tarray=tarray
    
    openw, lun, repstr(mosaic_file,'.fits','_STARS.dat'), /get_lun
    printf, lun, '# 2MASS calibration standards.'
    printf, lun, '# ZPT    = '+strtrim(string(zpt,format='(F12.4)'),2)
    printf, lun, '# ZPTERR = '+strtrim(string(zpt_err,format='(F12.4)'),2)
    printf, lun, '# '+strmid(tarray[0],2)
    printf, lun, '# '+strmid(tarray[1],2)
    struct_print, result, lun=lun, /no_head
    free_lun, lun

; QA plot    
    
    magaxis = findgen((25.0-5.0)/0.1)*0.1+5.0
    plotsym, 0, 1, /fill

    if keyword_set(postscript) then begin
       postthick1 = 4.0
    endif else begin
       postthick1 = 2.0
       im_window, 0, xratio=0.45, /square
    endelse

    psname = 'sg1120_ks.calib'
    im_openclose, datapath+psname, postscript=postscript, xsize=6.9, ysize=7.5, encapsulated=encapsulated
    pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=5.0, height=[3.5,2.5], $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=6.5, ypage=7.5, $
      position=pos, /normal

    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=postthick1, ythick=postthick1, $
      charsize=1.8, charthick=postthick1, xr=[9.5,16.5], yr=[9.5,16.5], $
      xtitle='', ytitle='2MASS mag', position=pos[*,0], xtickname=replicate(' ',10)
    djs_oplot, magaxis, poly(magaxis,[0.0,1.0]), thick=postthick1, line=0
;   djs_oplot, magaxis, poly(magaxis,[zpt,1.0]), thick=postthick1, line=2
;   djs_oplot, magaxis, poly(magaxis,coeff), thick=postthick1, line=0
    djs_oplot, result.mag_auto, result.mag, ps=8
    djs_oplot, result[keep].mag_auto, result[keep].mag, ps=8, color='dark green'

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xthick=postthick1, ythick=postthick1, $
      charsize=1.8, charthick=postthick1, xr=[9.5,16.5], yr=[-1,1], $
      xtitle='SE mag auto', ytitle='Residuals', position=pos[*,1]
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick1
    djs_oplot, result.mag_auto, result.mag-result.mag_auto, ps=8
    djs_oplot, result[keep].mag_auto, result[keep].mag-result[keep].mag_auto, ps=8, color='dark green'

    im_openclose, postscript=postscript, /close

stop       

return
end
    
