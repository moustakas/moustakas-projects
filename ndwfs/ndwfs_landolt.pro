function zeropoint_func, mtrue, param, minst=minst, $
  airmass=airmass, color=color
; m0 = m + Z - k1*AM - k2*C - k3*C*AM
;
;   m0 = Standard magnitude
;   m  = instrumental magnitude
;   Z  = zero point
;   k1 = extinction coefficient [mag/airmass]
;   AM = airmass
;   k2 = color term coefficient
;   C  = Standard system color of the object
;   k3 = second order extinction coefficient

    return, mtrue - param[0] + param[1]*airmass - $
      param[2]*color - param[3]*color*airmass
    
end

pro ndwfs_landolt, makelinks=makelinks, sextractor=sextractor, $
  zeropoint=zeropoint
; jm10feb18ucsd - parse the Landolt observations of the NDWFS

    ndwfspath = getenv('RESEARCHPATH')+'/data/ndwfs/landolt/'
    datapath = ndwfspath+'data/'
    outpath = ndwfspath+'results/'

    sexpath = ndwfspath+'sex/'
    sexconfig = sexpath+'default.sex'
    sexparam = sexpath+'sex.param.standards'
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'

    pixscale = 0.258 ; MOSAIC pixel scale [arcsec/pixel]
    pixscl = '0.258'

    photaper_arcsec = findgen(20)+1.0          ; [diameter,arcsec]
    photaper_pixel = photaper_arcsec/pixscale ; [diameter,pixel]
    photaper = strjoin(strtrim(string(photaper_pixel,format='(F10.2)'),2),',')

    bandpass = ['Bw','R','I']
    extcoeff = [0.25,0.08,0.03] 
    extcoeff_err = [0.1,0.08,0.03]
    
; read the Stetson (2000) photometric catalog    
    stet1 = read_00stetson() 

; --------------------------------------------------
; build the soft links    
    if keyword_set(makelinks) then begin
       spawn, '/bin/rm -f '+datapath+'*'
       pushd, ndwfspath
       allfits = file_search('*/*.fits',count=nfits)
       imroot = repstr(file_basename(allfits),'.fits','')
       for ii = 0, nfits-1 do begin
          hdr = headfits(allfits[ii])
          filter = repstr(strtrim(sxpar(hdr,'FILTER'),2),' ','_')
          date = strmid(sxpar(hdr,'DATE-OBS'),0,10)
          outfile = imroot[ii]+'_'+filter+'_'+date+'.fits'
          print, outfile
;         obj = sxpar(hdr,'OBJECT')
;         outfile = imroot[ii]+'_'+repstr(repstr(strtrim(obj,2),'.',''),' ','_')+'.fits'
          spawn, 'ln -sfv ../'+allfits[ii]+' data/'+outfile, /sh
       endfor
       popd
    endif

; ---------------------------------------------------------------------------    
; generate SE catalogs
    if keyword_set(sextractor) then begin
       imagelist = file_search(datapath+'*.fits',count=nimage)
;      imagelist = imagelist[1] & nimage = 1
       catlist = outpath+file_basename(repstr(imagelist,'.fits','.cat'))
       radecreglist = repstr(catlist,'.cat','.reg')

;      info = im_headerforage(imagelist,ext=0) ; grab header info
       
; initialize the SE configuration parameters
       config = init_sex_config(nimage)
       configfile = outpath+'sex.config'
       
       config.catalog_name = catlist
       config.catalog_type = 'FITS_LDAC'
       config.parameters_name = sexparam
       config.filter_name = sexconv
       config.starnnw_name = sexnnw

       config.detect_thresh = 20.0
       config.analysis_thresh = 20.0
       config.detect_minarea = 10.0
       config.deblend_mincont = 0.1
       config.interp_type = 'NONE'

       config.seeing_fwhm = 1.0
;      config.mag_zeropoint = 0.0

       config.checkimage_type = 'NONE' ; SEGMENTATION
       config.phot_apertures = photaper
       
       im_mwrfits, config, configfile+'.fits', /clobber

; do it!       
       t0 = systime(1)
       im_sex, imagelist, config, silent=silent ; do not pass CONFIGFILE
       splog, 'Total time to generate SE catalog = ', $
         (systime(1)-t0)/60.0, ' minutes'

; write a regions file          
       for ii = 0, nimage-1 do begin
          splog, 'Witing DS9 region file '+file_basename(radecreglist[ii])
          cat = mrdfits(catlist[ii],2)
          write_ds9_regionfile, cat.xwin_world, cat.ywin_world, $
            filename=radecreglist[ii], color='red'
       endfor 
    endif

; ---------------------------------------------------------------------------    
; generate the photometry table and compute the zeropoints 
    if keyword_set(zeropoint) then begin
; loop on each bandpass
;      for ib = 2, 2 do begin
       for ib = 0, n_elements(bandpass)-1 do begin
          zpt = {band: '', nallstar: 0, nstar: 0, $
            zpt: 0.0, zpt_err: 0.0, extcoeff: 0.0, $
            extcoeff_err: 0.0, colorcoeff: 0.0, $
            colorcoeff_err: 0.0}
          zpt.band = bandpass[ib]
          delvarx, info
;         imagelist = file_search(datapath+'*'+bandpass[ib]+'*.fits',count=nimage)
;         catlist = outpath+file_basename(repstr(imagelist,'.fits','.cat'))
          catlist = file_search(outpath+'*'+bandpass[ib]+'*.cat')
          imagelist = datapath+file_basename(repstr(catlist,'.cat','.fits'))
          nimage = n_elements(imagelist)

          for ii = 0, nimage-1 do begin
; match against Stetson             
             cat = mrdfits(catlist[ii],2)
             spherematch, cat.xwin_world, cat.ywin_world, $
               stet1.ra, stet1.dec, 5.0/3600.0, m1, m2
             if (m1[0] ne -1L) then begin
                stet = stet1[m2]
                info1 = struct_trimtags(cat[m1],select=[$
                  'XWIN_IMAGE','YWIN_IMAGE','XWIN_WORLD','YWIN_WORLD',$
                  'MAG_AUTO','MAGERR_AUTO','MAG_APER','MAGERR_APER',$
                  'FWHM_WORLD','FLAGS'])
                moretags = replicate({file: file_basename(imagelist[ii]), $
                  date: '', exptime: 0.0, airmass: 0.0},n_elements(info1))
                info1 = struct_addtags(struct_addtags(moretags,info1),stet)
                info1.date = sxpar(headfits(imagelist[ii]),'DATE-OBS')
                info1.exptime = sxpar(headfits(imagelist[ii]),'EXPTIME')
                info1.airmass = sxpar(headfits(imagelist[ii]),'AIRMASS')
                if (n_elements(info) eq 0L) then info = info1 else info = [info,info1]
             endif 
          endfor 
          
; apply some cuts to the standard-star observations: non-zero SE
; flags; small magnitude errors; and those with the best seeing
;         good = where((info.flags eq 0) and (info.magerr_auto lt 0.1) and $
;           (info.fwhm_world lt weighted_quantile(info.fwhm_world,quant=0.90)),ngood)
;         if (ngood eq 0L) then message, 'No good stars!'
;         info = info[good]
;         struct_print, struct_trimtags(info,except=['*aper*','arcfile'])

; write out the FITS catalog             
          outfile = outpath+'ndwfs_landolt_'+bandpass[ib]+'.fits'
          im_mwrfits, info, outfile, /clobber

          nstar = n_elements(info)
          splog, 'Bandpass '+bandpass[ib]+': '+string(nstar,format='(I0)')+' stars'
          
; compute the zeropoint
          case bandpass[ib] of
             'Bw': begin
                these = where((info.mag_b lt 90.0) and (info.mag_v lt 90.0),nthese)
                mtrue = info[these].mag_b & color = info[these].mag_b-info[these].mag_v
             end
             'R': begin
                these = where((info.mag_r lt 90.0) and (info.mag_v lt 90.0),nthese)
                mtrue = info[these].mag_r & color = info[these].mag_v-info[these].mag_r
             end
             'I': begin
                these = where((info.mag_i lt 90.0) and (info.mag_r lt 90.0),nthese)
                mtrue = info[these].mag_i & color = info[these].mag_r-info[these].mag_i
             end
          endcase

; iteratively reject outliers
          minst = info[these].mag_auto + 2.5*alog10(info[these].exptime)

          djs_iterstat, mtrue-minst, median=zpt_guess, $
            mask=mask, sigrej=2.0, maxiter=100
          goodstars = where(mask,nstar)
          if (nstar eq 0L) then message, 'This is very bad'

          minst = minst[goodstars]
          mtrue = mtrue[goodstars]
          color = color[goodstars]
          airmass = info[these[goodstars]].airmass
          minst_err = info[these[goodstars]].magerr_auto

          parinfo = replicate({value: 0.0D, fixed: 0},4)         ; 4-parameter fit
          parinfo[0].value = zpt_guess & parinfo[0].fixed = 0    ; zeropoint
             
          parinfo[1].value = extcoeff[ib] & parinfo[1].fixed = 1
;         parinfo[1].value = 0.1D & parinfo[1].fixed = 0 ; linear extinction coeff

;         parinfo[2].value = colorcoeff[iq,ib] & parinfo[2].fixed = 1 ; color term

          parinfo[2].value = 0.1D & parinfo[2].fixed = 0 ; color term
          parinfo[3].value = 0.0D & parinfo[3].fixed = 1 ; quadratic extinction coeff

          functargs = {minst: minst, airmass: airmass, color: color}
          param = mpfitfun('zeropoint_func',mtrue,minst,minst_err,$
            functargs=functargs,perror=perror,quiet=0,parinfo=parinfo,$
            dof=dof,bestnorm=bestnorm,yfit=minst_best)
;         niceprint, minst_best, minst, mtrue-param[0]

          zpt.zpt = param[0]
          zpt.zpt_err = djsig(minst-minst_best,sigrej=2.0)
          zpt.extcoeff = param[1]
;         zpt.extcoeff_err = extcoeff_err[ib]
          zpt.colorcoeff = param[2] ; note!
          zpt.colorcoeff_err = 0.0

; generate a QA plot
          djs_plot, [0], [0], /nodata, xthick=2.0, ythick=2.0, $
            charsize=2.0, charthick=2.0, xsty=1, ysty=1, $
            xrange=[-0.3,2.0], xtitle='R-I', ytitle='I-i', $
            yrange=minmax(mtrue-minst)+[-0.1,0.1]
          plotsym, 8, 1.5, fill=1
          djs_oplot, color, mtrue-minst, ps=8
;         cc = get_kbrd(1)

; write out the FITS table for this bandpass
          outfile = outpath+'ndwfs_landolt_zpt_'+bandpass[ib]+'.fits'
          im_mwrfits, zpt, outfile, /clobber
          struct_print, zpt
       endfor ; bandpass loop

; generate a curve-of-growth plot for all bandpasses/quadrants
       qafile = outpath+'qa_curve_of_growth.ps'
       im_plotconfig, 0, pos, psfile=qafile
;      for ib = 2, 2 do begin
       for ib = 0, n_elements(bandpass)-1 do begin
          info = mrdfits(outpath+'ndwfs_landolt_'+bandpass[ib]+'.fits.gz',1)
          plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
            xrange=[0,max(photaper_arcsec)], yrange=[0.2,1.3], $
            xtitle='Aperture Diameter (arcsec)', $
            ytitle='Light Fraction (relative to MAG_AUTO)', $
            title=bandpass[ib]
          djs_oplot, !x.crange, [1,1], line=0, thick=2.0
          for kk = 0, n_elements(info)-1 do djs_oplot, $
            photaper_arcsec, 10^(-0.4*(info[kk].mag_aper-info[kk].mag_auto)), $
            psym=6
       endfor
       im_plotconfig, psfile=qafile, /psclose, /gzip
    endif 
    
stop    
    
; average BVR extinction coefficients, assumed constant for all
; quadrants; see, e.g., 
; http://www.eso.org/observing/dfo/quality/FORS2/qc/photcoeff/photcoeffs_fors2.html
    
    colorcoeff = [$
      [-0.074,-0.076,-0.076,-0.081],$ ; B
      [-0.021,-0.029,-0.020,-0.035],$ ; V
      [-0.110,-0.112,-0.068,-0.087] ] ; R
    colorcoeff_err = [$
      [0.005,0.005,0.006,0.006],$
      [0.005,0.006,0.005,0.004],$
      [0.007,0.007,0.008,0.010] ]

return
end
