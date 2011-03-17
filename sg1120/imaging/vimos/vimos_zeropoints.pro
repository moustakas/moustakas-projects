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

pro vimos_zeropoints, dec03=dec03, feb06=feb06, sextractor=sextractor, $
  zeropoint=zeropoint, update_headers=update_headers, qc_zeropoint=qc_zeropoint
; jm08jul18nyu - derive the photometric zeropoints for the VIMOS data 
    
    if ((not keyword_set(dec03)) and (not keyword_set(feb06))) or $
      (keyword_set(dec03) and keyword_set(feb06)) then begin
       splog, 'Must specify either DEC03 or FEB06 keyword!'
       return
    endif
    
    sexpath = sg1120_path(/sex)
    qcpath = vimos_path(/qc)
    qaplotspath = vimos_path(/qaplots)
    datapath = vimos_path(dec03=dec03,feb06=feb06)+'sg1120/'
    standardspath = vimos_path(dec03=dec03,feb06=feb06)+'standards/'

    if keyword_set(dec03) then suffix = '03dec'
    if keyword_set(feb06) then suffix = '06feb'

;   bandpass = ['R']
    bandpass = ['B','V','R']
    quadrant = ['Q1','Q2','Q3','Q4']
    gain = '1.73,1.86,1.95,1.80' ; [e/ADU]

; average BVR extinction coefficients, assumed constant for all
; quadrants; see, e.g., 
; http://www.eso.org/observing/dfo/quality/FORS2/qc/photcoeff/photcoeffs_fors2.html
    
    extcoeff = [0.25,0.15,0.10] 
    extcoeff_err = extcoeff*0.1

    colorcoeff = [$
      [-0.074,-0.076,-0.076,-0.081],$ ; B
      [-0.021,-0.029,-0.020,-0.035],$ ; V
      [-0.110,-0.112,-0.068,-0.087] ] ; R
    colorcoeff_err = [$
      [0.005,0.005,0.006,0.006],$
      [0.005,0.006,0.005,0.004],$
      [0.007,0.007,0.008,0.010] ]
    
    sexconfig = sexpath+'default.sex'
    sexparam = sexpath+'sg1120.sex.param.standards'
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'

    pixscale = 0.205 ; pixel scale [arcsec/pixel]
    photaper_arcsec = findgen(15)+1.0         ; [diameter,arcsec]
    photaper_pixel = photaper_arcsec/pixscale ; [diameter,pixel]
    photaper = strjoin(strtrim(string(photaper_pixel,format='(F10.2)'),2),',')

    stet1 = read_00stetson() ; read the Stetson (2000) photometric catalog

; ---------------------------------------------------------------------------    
; compare the zeropoints downloaded from the VIMOS/QC website to the
; zeropoints placed in the image headers; note that this does not tell
; us if our science observations were observed during photometric
; conditions, only whether our *standards* were observed during
; photometric conditions
; ---------------------------------------------------------------------------    
    
    if keyword_set(qc_zeropoint) then begin

       pagemaker, nx=2, ny=2, xpage=8.5, ypage=8.5, $
         xspace=0.0, yspace=0.0, xmargin=[1.2,0.2], $
         ymargin=[0.3,1.1], width=3.55*[1,1], $
         height=3.55*[1,1], /normal, position=pos
       dfpsplot, qaplotspath+'qc_zeropoint_'+suffix+'.ps', $
         xsize=8.5, ysize=8.5, /color
       
       for ib = 0L, n_elements(bandpass)-1L do begin

          imagelist = file_search(datapath+'ra.sg1120*_'+$
            bandpass[ib]+'.fits',count=nimage)
          
          for iq = 0L, n_elements(quadrant)-1L do begin

; read the zeropoint computed by the VIMOS pipeline
             obs = replicate({mjd: 0.0D, zpt: 0.0, zpt_err: 0.0},nimage)
             for ii = 0L, nimage-1L do begin
                hdr = headfits(imagelist[ii],ext=iq+1L)
                obs[ii].mjd = sxpar(hdr,'MJD-OBS')
                obs[ii].zpt = sxpar(hdr,'MAG0')
                obs[ii].zpt_err = sxpar(hdr,'EMAG0')
             endfor

; read the QC historical zerotpoints, which were presumbly only
; measured during photometric conditions; these were downloaded for
; each bandpass and quadrant from
; http://www.eso.org/observing/dfo/quality/VIMOS/qc/zeropoints.html   

             readcol, qcpath+'vimos_zpt_'+bandpass[ib]+'_'+quadrant[iq]+'.dat', $
               mjd, qczpt, filt, dat, qq, skipline=2, format='D,F,A,A,I', /silent

             mjdpivot = long(mean(obs.mjd))
             zptpivot = mean(obs.zpt)

             xrange = 49*[-1,1] ; days
;            yrange = zptpivot+0.5*[-1,1]
             if keyword_set(dec03) then begin
                case bandpass[ib] of 
                   'B': yrange = [26.9,27.9]
                   'V': yrange = [26.6,27.6]
                   'R': yrange = [26.6,27.6]
                endcase
             endif
             if keyword_set(feb06) then begin
                case bandpass[ib] of 
                   'B': yrange = [27.1,28.1]
                   'V': yrange = [26.7,27.7]
                   'R': yrange = [26.7,27.7]
                endcase
             endif

             if (iq eq 0L) or (iq eq 1L) then begin
                xtitle = ''
                xtickname = replicate(' ',10)
             endif else begin
                xtitle = 'MJD - '+string(mjdpivot,format='(I0)')+' (days)'
                delvarx, xtickname
             endelse

             if (iq eq 0L) or (iq eq 2L) then begin
                ytitle = 'Zeropoint (mag)'
                delvarx, ytickname
             endif else begin
                ytitle = ''
                ytickname = replicate(' ',10)
             endelse
             
             djs_plot, [0], [0], /nodata, noerase=(iq gt 0L), xsty=1, ysty=1, $
               xrange=xrange, yrange=yrange, charsize=1.7, charthick=3.0, $
               xthick=4.0, ythick=4.0, position=pos[*,iq], $
               xtitle=xtitle, ytitle=ytitle, xtickname=xtickname, $
               ytickname=ytickname
             djs_oplot, !x.crange, zptpivot*[1,1]+0.1, line=5, thick=3.0
             djs_oplot, !x.crange, zptpivot*[1,1]-0.1, line=5, thick=3.0
             legend, bandpass[ib]+'/'+quadrant[iq], /left, /top, box=0, $
               charsize=1.7, charthick=3.0, margin=0
             djs_oplot, mjd-mjdpivot, qczpt, ps=4, thick=4.0
             djs_oplot, obs.mjd-mjdpivot, obs.zpt, ps=6, sym=2.5, $
               color='red', thick=5.0

          endfor 

       endfor

       dfpsclose
       
    endif
                   
; ---------------------------------------------------------------------------    
; generate SE catalogs
; ---------------------------------------------------------------------------    
    
    if keyword_set(sextractor) then begin

       t0 = systime(1)
       for ib = 0L, n_elements(bandpass)-1L do begin
          for iq = 0L, n_elements(quadrant)-1L do begin

             imagelist = file_search(standardspath+'ra.[PG,SA]*_'+bandpass[ib]+'_'+quadrant[iq]+'.fits')
             weightlist = repstr(imagelist,'.fits','.weight.fits')
             rmslist = repstr(imagelist,'.fits','.rms.fits')
             flaglist = repstr(imagelist,'.fits','.flag.fits')

             catlist = repstr(imagelist,'.fits','.cat')
             xyreglist = repstr(catlist,'.cat','.xy.reg')
             radecreglist = repstr(catlist,'.cat','.radec.reg')
             radecreglist_stetson = repstr(catlist,'.cat','.radec.stetson.reg')

             for ii = 0L, n_elements(imagelist)-1L do begin
                print, 'SExtracting '+imagelist[ii]
                spawn, 'sex '+imagelist[ii]+' -c '+sexconfig+' -CATALOG_NAME '+$
                  catlist[ii]+' -CATALOG_TYPE FITS_LDAC'+$
                  ' -PARAMETERS_NAME '+sexparam+' -DETECT_THRESH 3.0 -ANALYSIS_THRESH 3.0 '+$
                  ' -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+$
                  ' -WEIGHT_TYPE MAP_RMS -WEIGHT_IMAGE '+rmslist[ii]+' -WEIGHT_GAIN N'+$
                  ' -FLAG_IMAGE '+flaglist[ii]+' -FLAG_TYPE OR -INTERP_TYPE NONE'+$
                  ' -PHOT_APERTURES '+photaper+' -VERBOSE_TYPE NORMAL -NTHREADS 4', /sh

; write out regions files of all the sources and of the Stetson
; standards; this is just for a quick-look comparison
                
                cat = mrdfits(catlist[ii],2,/silent)
                spherematch, cat.xwin_world, cat.ywin_world, $
                  stet1.ra, stet1.dec, 5.0/3600.0, m1, m2
;               splog, 'Witing DS9 region file '+file_basename(radecreglist[ii])
;               write_ds9_regionfile, cat.xwin_world, cat.ywin_world, $
;                 filename=radecreglist[ii], color='red'
                if (m1[0] eq -1L) then begin
                   splog, 'No standard stars found!'
                endif else begin
                   splog, 'Witing DS9 region file '+file_basename(radecreglist_stetson[ii])
                   write_ds9_regionfile, cat[m1].xwin_world, cat[m1].ywin_world, $
                     filename=radecreglist_stetson[ii], color='green', symbol='cross'
                endelse
             endfor
             
          endfor
       endfor
       splog, 'Total time to generate SE catalogs = ', (systime(1)-t0)/60.0, ' minutes.'
          
    endif

; ---------------------------------------------------------------------------    
; generate the photometry table and compute the zeropoints 
; ---------------------------------------------------------------------------    
    
    if keyword_set(zeropoint) then begin
; loop on each bandpass
       for ib = 0L, n_elements(bandpass)-1L do begin
          zpt = {band: '', quadrant: '', nallstar: 0, nstar: 0, $
            zpt: 0.0, zpt_err: 0.0, $
            extcoeff: 0.0, extcoeff_err: 0.0, $
            colorcoeff: 0.0, colorcoeff_err: 0.0}
          zpt = replicate(zpt,4) ; one per quadrant
          zpt.band = bandpass[ib]
          zpt.quadrant = quadrant
; loop on each quadrant
          delvarx, info
          for iq = 0L, n_elements(quadrant)-1L do begin
             imagelist = file_search(standardspath+'ra.[PG,SA]*_'+$
               bandpass[ib]+'_'+quadrant[iq]+'.fits')
             catlist = repstr(imagelist,'.fits','.cat')
             for ii = 0L, n_elements(imagelist)-1L do begin
; non-photometric standard-star observations
;               if $
;                 (strmatch(imagelist[ii],'*ra.SA98_2006-02-02T00*_B*') eq 0B) and $
;                 (strmatch(imagelist[ii],'*ra.SA98_2006-02-02T00*_V*') eq 0B) and $
;                 (strmatch(imagelist[ii],'*ra.SA98_2006-02-02T00*_R*') eq 0B) then begin
;                  splog, 'Reading '+catlist[ii]
                   cat = mrdfits(catlist[ii],2,/silent) ; match against Stetson
                   spherematch, cat.xwin_world, cat.ywin_world, $
                     stet1.ra, stet1.dec, 5.0/3600.0, m1, m2
                   if (m1[0] ne -1L) then begin
                      stet = stet1[m2]
                      info1 = struct_trimtags(cat[m1],select=[$
                        'XWIN_IMAGE','YWIN_IMAGE','XWIN_WORLD','YWIN_WORLD',$
                        'MAG_BEST','MAGERR_BEST','MAG_APER','MAGERR_APER',$
                        'FWHM_WORLD','FLAGS'])
                      moretags = replicate({file: file_basename(imagelist[ii]), $
                        arcfile: '', exptime: 0.0, airmass: 0.0},n_elements(info1))
                      info1 = struct_addtags(struct_addtags(moretags,info1),stet)
                      info1.arcfile = sxpar(headfits(imagelist[ii]),'ARCFILE')
                      info1.exptime = sxpar(headfits(imagelist[ii]),'EXPTIME')
                      info1.airmass = sxpar(headfits(imagelist[ii]),'AIRMASS')
                      if (n_elements(info) eq 0L) then info = info1 else info = [info,info1]
                   endif
;               endif
             endfor

             zpt[iq].nallstar = n_elements(info)
             
; apply some cuts to the standard-star observations: non-zero SE
; flags; small magnitude errors; and those with the best seeing

             good = where((info.flags eq 0) and (info.magerr_best lt 0.1) and $
               (info.fwhm_world lt weighted_quantile(info.fwhm_world,quant=0.90)),ngood)
             if (ngood eq 0L) then message, 'No good stars!'
             info = info[good]
;            struct_print, struct_trimtags(info,except=['*aper*','arcfile'])

; write out the FITS catalog             

             outfile = standardspath+'photo_standards_'+$
               bandpass[ib]+'_'+quadrant[iq]+'.fits'
             splog, 'Writing '+outfile
             mwrfits, info, outfile, /create

; compute the zeropoint
             case bandpass[ib] of
                'B': begin
                   these = where((info.mag_b lt 90.0) and (info.mag_v lt 90.0),nthese)
                   mtrue = info[these].mag_b & color = info[these].mag_b-info[these].mag_v
                end
                'V': begin
                   these = where((info.mag_b lt 90.0) and (info.mag_v lt 90.0),nthese)
                   mtrue = info[these].mag_v & color = info[these].mag_b-info[these].mag_v
                end
                'R': begin
                   these = where((info.mag_r lt 90.0) and (info.mag_v lt 90.0),nthese)
                   mtrue = info[these].mag_r & color = info[these].mag_v-info[these].mag_r
                end
             endcase

; iteratively reject outliers
             
             minst = info[these].mag_best + 2.5*alog10(info[these].exptime)

             djs_iterstat, mtrue-minst, median=zpt_guess, $
               mask=mask, sigrej=2.0, maxiter=100
             goodstars = where(mask,nstar)
             if (nstar eq 0L) then message, 'This is very bad'
             zpt[iq].nstar = nstar

             minst = minst[goodstars]
             mtrue = mtrue[goodstars]
             color = color[goodstars]
             airmass = info[these[goodstars]].airmass
             minst_err = info[these[goodstars]].magerr_best

             parinfo = replicate({value: 0.0D, fixed: 0},4) ; 4-parameter fit
             parinfo[0].value = zpt_guess & parinfo[0].fixed = 0 ; zeropoint
             
             parinfo[1].value = extcoeff[ib] & parinfo[1].fixed = 1 ; linear extinction coeff
;            parinfo[2].value = colorcoeff[iq,ib] & parinfo[2].fixed = 1 ; color term
;            parinfo[1].value = 0.0D  & parinfo[1].fixed = 1
             parinfo[2].value = 0.0D  & parinfo[2].fixed = 1
             parinfo[3].value = 0.0D  & parinfo[3].fixed = 1 ; quadratic extinction coeff

             functargs = {minst: minst, airmass: airmass, color: color}
             param = mpfitfun('zeropoint_func',mtrue,minst,minst_err,$
               functargs=functargs,perror=perror,quiet=1,parinfo=parinfo,$
               dof=dof,bestnorm=bestnorm,yfit=minst_best)
;            niceprint, minst_best, minst, mtrue-param[0]

             zpt[iq].zpt = param[0]
             zpt[iq].zpt_err = djsig(minst-minst_best,sigrej=2.0)
             zpt[iq].extcoeff = param[1]
             zpt[iq].extcoeff_err = extcoeff_err[ib]
             zpt[iq].colorcoeff = param[2] ; note!
             zpt[iq].colorcoeff_err = 0.0

; generate a QA plot

;            djs_plot, [0], [0], /nodata, xthick=2.0, ythick=2.0, $
;              charsize=2.0, charthick=2.0, xsty=1, ysty=1, $
;              xrange=[-0.3,2.0], xtitle='B-V', ytitle='B-b', $
;              yrange=minmax(mtrue-minst)+[-0.1,0.1]
;            plotsym, 8, 1.5, fill=1
;            djs_oplot, color, mtrue-minst, ps=8
;            zpt[iq].zpt+zpt[iq].zpt_err*10.0*[-1,1]
             
          endfor 

; write out the FITS table for this bandpass

          outfile = standardspath+'photo_zpt_table_'+bandpass[ib]+'.fits'
          splog, 'Writing '+outfile
          mwrfits, zpt, outfile, /create
          struct_print, zpt & print

       endfor 

; generate a curve-of-growth plot for all bandpasses/quadrants

       qafile = standardspath+'qa_curve_of_growth.ps'
       splog, 'Writing '+qafile
       dfpsplot, qafile, /square, /color
       for ib = 0L, n_elements(bandpass)-1L do begin
          for iq = 0L, n_elements(quadrant)-1L do begin
             info = mrdfits(standardspath+'photo_standards_'+$
               bandpass[ib]+'_'+quadrant[iq]+'.fits',1,/silent)
             plot, [0], [0], /nodata, xthick=4.0, ythick=4.0, $
               charsize=1.8, charthick=3.0, xsty=1, ysty=1, $
               xrange=[0,max(photaper_arcsec)], yrange=[0.2,1.3], $
               xtitle='Aperture Diameter (arcsec)', $
;              ytitle='F(<D)/F(MAG_BEST)', $
               ytitle='Light Fraction (relative to MAG_BEST)', $
               title=bandpass[ib]+'/'+quadrant[iq]
             djs_oplot, !x.crange, [1,1], line=0, thick=2.0
             for kk = 0L, n_elements(info)-1L do djs_oplot, $
               photaper_arcsec, 10^(-0.4*(info[kk].mag_aper-info[kk].mag_best)), $
               psym=6; line=0, thick=4.0

; generate a histogram of seeing values

;            seeing = info.fwhm_world
;            stats = im_stats(seeing,sigrej=3.0)
;            xstats = strtrim(string(stats.mean_rej,format='(F12.2)'),2)+'\pm'+$
;              strtrim(string(stats.sigma_rej,format='(F12.3)'),2)
;            
;            im_plothist, seeing, bin=0.1, xx, yy, /noplot
;            xrange = [0.4,max(xx)*1.05]
;            yrange = [0.0,max(yy)*1.05]
;
;            plot, [0], [0], /nodata, /noerase, xthick=4.0, ythick=4.0, $
;              charsize=1.4, charthick=3.0, xsty=1, ysty=1, $
;              xrange=xrange, yrange=yrange, position=[0.5,0.3,0.9,0.6], $
;              xtitle='Seeing (FWHM, arcsec)', ytitle='Number of Stars'
;            im_plothist, seeing, bin=0.1, thick=4.0, /overplot
;            legend, textoidl(xstats), /right, /top, box=0, $
;              charsize=1.4, charthick=3.0

          endfor
       endfor
       dfpsclose 
    endif 

; ---------------------------------------------------------------------------    
; update the FITS headers
; ---------------------------------------------------------------------------    
    
    if keyword_set(update_headers) then begin
       for ib = 0L, n_elements(bandpass)-1L do begin
          imagelist = file_search(datapath+'ra.sg1120*_'+bandpass[ib]+'.fits')
          weightlist = repstr(imagelist,'.fits','.weight.fits')
          rmslist = repstr(imagelist,'.fits','.rms.fits')
          flaglist = repstr(imagelist,'.fits','.flag.fits')
          splog, 'Reading '+standardspath+'photo_zpt_table_'+bandpass[ib]+'.fits'
          zpt = mrdfits(standardspath+'photo_zpt_table_'+bandpass[ib]+'.fits',1,/silent)
          for ii = 0L, n_elements(imagelist)-1L do begin
; add a scamp keyword indicating which observations were taken in
; 'photometric' conditions; these were determined by looking at the
; scamp output plots (see vimos_mosaics, /scamp) and are *only*
; appropriate when using both epochs simultaneously
             case file_basename(imagelist[ii]) of
                'ra.sg1120_2003-12-27T07:37:16.456_V.fits': photflag = 'T'
                'ra.sg1120_2003-12-27T06:57:46.975_B.fits': photflag = 'T'
                'ra.sg1120_2006-02-04T08:16:30.340_R.fits': photflag = 'T'
                else: photflag = 'F'
             endcase
             if (photflag eq 'T') then splog, '***Photometric field '+$
               file_basename(imagelist[ii])+'***'
             for iq = 1L, n_elements(quadrant) do begin
                thesefits = [imagelist[ii],weightlist[ii],rmslist[ii],flaglist[ii]]
                if (iq eq 1L) then print, '--------------------------------------------------'
                for jj = 0L, n_elements(thesefits)-1L do begin
                   if (iq eq 1L) then splog, 'Updating '+file_basename(thesefits[jj])
                   hdr = headfits(thesefits[jj],ext=iq)
                   sxaddpar, hdr, 'MAGZERO', float(zpt[iq-1L].zpt), $
                     ' magnitude zero-point (counts/s)', after='EMAG0'
                   sxaddpar, hdr, 'EMAGZERO', float(zpt[iq-1L].zpt_err), $
                     ' error in MAGZERO', after='MAGZERO'
                   sxaddpar, hdr, 'NZPTSTAR', fix(zpt[iq-1L].nstar), ' number '+$
                     'of stars used to compute MAGZERO', after='EMAGZERO'
                   sxaddpar, hdr, 'PHOT_K', float(zpt[iq-1L].extcoeff), ' atmospheric '+$
                     'extinction coefficient (mag/airmass)', after='NZPTSTAR'
                   sxaddpar, hdr, 'PHOTFLAG', photflag, ' photometric? (T/F)', after='PHOT_K'
;                  sxaddhist, "'Photometry parameters added "+hogg_iso_date()+"'", hdr
                   modfits, thesefits[jj], 0, hdr, ext=iq
                endfor
             endfor
          endfor 
       endfor
    endif

return
end
