pro ldss3_catalogs, make_flag_image=make_flag_image, sextractor=sextractor, $
  clean_catalogs=clean_catalogs, psf=psf, qaplots=qaplots, zmatch=zmatch, $
  dcr=dcr, band=band, reread=reread
; jm07aug20nyu - generate the LDSS3 catalogs using SE in double-image
;                mode with the chi2-band image as the detection image
;                (based on earlier code)
; jm09mar23nyu - major rewrite to the v2.0 catalogs, with calibration!

    common ldss3_cat, gcat, rcat, gcat_chi2, rcat_chi2
    
    version = ldss3_catalogs_version()
    
    match_radius_zcat = 1.0     ; [arcsec]
    if (n_elements(band) eq 0L) then band = ['gprime','rprime']
    nband = n_elements(band)
    
    sexpath = sg1120_path(/sex)
    catpath = ldss3_path(/catalogs)
    mosaicpath = ldss3_path(/mosaics)
    
    sexconfig = sexpath+'default.sex'
    sexparam = sexpath+'sg1120.sex.param.mosaic'
;   sexconv = sexpath+'gauss_5.0_9x9.conv' 
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'

    http = 'http://sdss.physics.nyu.edu/ioannis/research/sg1120/sex/'
    xslsex = http+'sex.xsl'

    pixscale = 0.188 ; pixel scale [arcsec/pixel]
    pixscl = '0.188'

    photaper_arcsec = findgen(5)+1.0          ; [diameter,arcsec]
    photaper_pixel = photaper_arcsec/pixscale ; [diameter,pixel]
    photaper = strjoin(strtrim(string(photaper_pixel,format='(F10.2)'),2),',')

; read the catalogs    
    
    if (not keyword_set(sextractor)) and (not keyword_set(clean_catalogs)) then begin
       catnames = catpath+'sg1120_'+['gprime','rprime']+'_'+version+'.cat'
       if (n_elements(gcat) eq 0L) or keyword_set(reread) then $
         gcat = rsex(catnames[0])
       if (n_elements(rcat) eq 0L) or keyword_set(reread) then $
         rcat = rsex(catnames[1])

       catnames_chi2 = catpath+'sg1120_'+['gprime','rprime']+'_'+version+'.chi2.cat'
       if (n_elements(gcat_chi2) eq 0L) or keyword_set(reread) then $
         gcat_chi2 = rsex(catnames_chi2[0])
       if (n_elements(rcat_chi2) eq 0L) or keyword_set(reread) then $
         rcat_chi2 = rsex(catnames_chi2[1])
    endif
    
; ---------------------------------------------------------------------------    
; generate a flag image for each of the mosaics; for now, just set a
; single bit (2^0) for pixels where the weight is zero
; ---------------------------------------------------------------------------    

    if keyword_set(make_flag_image) then begin
       for ib = 0L, nband-1L do begin
          weightname = mosaicpath+'sg1120_'+band[ib]+'.weight.fits'
          flagname = repstr(weightname,'weight','flag')

          splog, 'Reading '+weightname
          weight = mrdfits(weightname,0,hdr,/silent)
          flag = fix(weight eq 0)

          splog, 'Writing '+flagname
          mwrfits, flag, flagname, hdr, /create
       endfor
    endif
       
; ---------------------------------------------------------------------------
; for each [gr] bandpass generate two catalogs: one in which detection
; and measurements are done on the same image and another in which the
; chi2 image is used for detection; note that in the latter case all
; three catalogs will have the same number of sources, whereas that
; will not be true in the former case

    if keyword_set(sextractor) then begin

; Note:
;   * the mosaics are already in physical units, so we don't need
;     MAGZERO
;   * the seeing for each mosaic was determined iteratively (see PSF
;     keyword)
;   * note that if the CHI2 image is used for detection FWHM_WORLD and
;     CLASS_STAR are not very meaningful
;   * setting INTERP_TYPE to NONE prevents spurious sources near the
;     edges from being identified
;   * PHOT_FLUXFRAC is set to 0.2,0.5,0.9 so that we can measure the
;     SDSS-style 'concentration'

; initialize the configuration parameters       
       config = init_sex_config()

       config.parameters_name = sexparam
       config.filter_name = sexconv
       config.starnnw_name = sexnnw

       config.catalog_type = 'ASCII_HEAD'
       config.detect_thresh = 1.5
       config.analysis_thresh = 1.5
;      config.deblend_mincont = 0.0025

       config.weight_type = 'MAP_WEIGHT'
       config.weight_gain = 'Y'
       config.interp_type = 'NONE'
       config.nthreads = 8
       config.write_xml = 'N' ; 'Y' ; xml output lacks style sheet
       config.xsl_url = xslsex

       config.checkimage_type = 'SEGMENTATION,APERTURES' ; NONE

       config.seeing_fwhm = 1.0 ; update this!
       config.phot_apertures = photaper
       config.phot_fluxfrac = '0.2,0.5,0.9'
       config.pixel_scale = pixscale

; build the catalogs       
       for ib = 0L, nband-1L do begin

          imagename = mosaicpath+'sg1120_'+band[ib]+'.fits'
          weightname = repstr(imagename,'.fits','.weight.fits')
          flagname = repstr(imagename,'.fits','.flag.fits')

          detect_imagename = [imagename,mosaicpath+'sg1120_gr_chi2.fits']
          detect_weightname = repstr(detect_imagename,'.fits','.weight.fits')
          detect_suffix = ['','.chi2']

;         for idetect = 1, 1 do begin 
          for idetect = 0, 1 do begin 
             
             config.weight_image = detect_weightname[idetect]+','+weightname

             catname = catpath+'sg1120_'+band[ib]+'_'+$
               version+detect_suffix[idetect]+'.cat'
             
             config.catalog_name = catname
             config.xml_name = repstr(catname,'.cat','.xml')
             config.checkimage_name = strjoin(repstr(catname,'.cat','')+$
               ['.seg','.aper']+'.fits',',')
             config.flag_image = flagname
             
             configfile = repstr(catname,'.cat','.config')
             mwrfits, config, configfile+'.fits', /create

             im_sex, imagename, config, detect_imagelist=detect_imagename[idetect], $
               configfile=configfile, silent=silent

; write out a region file for all sources
             radecregname = repstr(catname,'.cat','.reg')
             splog, 'Witing DS9 region file'+radecregname
             cat = rsex(catname)
             write_ds9_regionfile, cat.xwin_world, cat.ywin_world, $
               filename=radecregname, color='red', symbol='box'
          endfor
       endfor
    endif

; ---------------------------------------------------------------------------
; clean up the catalogs: identify multiple (poorly deblended) sources
; and also objects with specified non-zero SE flags

    if keyword_set(clean_catalogs) then begin
       if (n_elements(dcr) eq 0L) then dcr = 0.5D ; [arcsec]

       detect_suffix = ['','.chi2']
       for idetect = 0L, 1L do begin

          catnames = catpath+'sg1120_'+band+'_'+version+$
            detect_suffix[idetect]+'.cat'
          radecregnames = repstr(catnames,'.cat','.reg')
          ncat = n_elements(catnames)

          for icat = 0L, ncat-1L do begin

             splog, 'Reading '+catnames[icat]
             cat = rsex(catnames[icat])
             ngalaxy = n_elements(cat)

             if tag_exist(cat,'FLAG_SPLIT') then begin
                cat.flag_split = 0L ; reset
             endif else begin
                moretags = replicate({flag_split: 0L},ngalaxy)
                cat = struct_addtags(temporary(cat),moretags)
             endelse
             
             ingroup = spheregroup(cat.xwin_world,cat.ywin_world,dcr/3600D0,$
               multgroup=mult,firstgroup=first,nextgroup=next)
             ugroup = ingroup[uniq(ingroup,sort(ingroup))] ; unique groups
             finalindx = lindgen(ngalaxy)

             for igroup = 0L, n_elements(ugroup)-1L do begin
                match = where((ingroup eq ugroup[igroup]),nmatch)
                if (nmatch gt 1) then begin
                   cat[match].flag_split = 2^0 ; set the bit
                   struct_print, /no_head, struct_trimtags(cat[match],$
                     select=['XWIN_WORLD','YWIN_WORLD','XWIN_IMAGE',$
                     'YWIN_IMAGE','MAG_AUTO','MAGERR_AUTO','FLUX_RADIUS'])
                   print
                   if (nmatch gt 1L) then begin
;                     cat[match].nmulti = nmatch ; keep the brightest object
                      minmag = min(cat[match].mag_auto,minmagindx)
                      remove, minmagindx, match & finalindx[match] = -1L
                   endif
;                  case nmatch of
;                     2L: begin
;                        cat[match].nmulti = nmatch ; case (1) [potentially dumb]
;                        minmag = min(cat[match].mag_auto,minmagindx)
;                        remove, minmagindx, match & finalindx[match] = -1L
;                     end
;                     3L: begin
;                        cat[match].nmulti = nmatch ; case (1) [potentially dumb]
;                        minmag = min(cat[match].mag_auto,minmagindx) 
;                        remove, minmagindx, match & finalindx[match] = -1L
;                     end
;                  endcase
                endif 
             endfor             ; close UGROUP loop

             mult = where((finalindx eq -1L),nmult,comp=keep,ncomp=nkeep)
;            finalcat = cat[finalindx[keep]]
;            finalcat.number = lindgen(nkeep)+1L
             finalcat = cat
             
             splog, 'Identified '+string(nmult,format='(I0)')+' shredded sources'

; SE flags:

;  *   1 - The object has neighbours, bright and close enough to
;          significantly bias the MAG AUTO photometry, or bad pixels
;          (more than 10% of the integrated area affected)
;  *   2 - The object was originally blended with another one
;  *   4 - At least one pixel of the object is saturated
;  *   8 - The object is truncated (too close to an image boundary) 
;  *  16 - Object's aperture data are incomplete or corrupted
;  *  32 - Object's isophotal data are incomplete or corrupted
;  *  64 - A memory overfow occurred during deblending
;  * 128 - A memory overfow occurred during extraction

;         crap = where($
;           (finalcat.flags and 2L^0L) or $ ; close neighbors/many bad pixels
;           (finalcat.flags and 2L^2L) or $ ; saturated pixel
;           (finalcat.flags and 2L^3L),$    ; truncated
;           ncrap,comp=good,ncomp=ngood)

; my flags: (we will likely want to add flags for scattered light and
; the halos around bright stars)

;  * 1 - at least one pixel within the isophotal aperture has a weight
;        map value of zero
             
;            crap = where($
;              (finalcat.imaflags_iso gt 0) or $ ; masked pixel
;              (finalcat.flags and 2L^2L),$      ; saturated pixel
;              ncrap,comp=good,ncomp=ngood)
;
;            finalcat = finalcat[good]
;            splog, 'Identified and removed '+string(ncrap,$
;              format='(I0)')+' objects with specified SE flags'
             
             splog, '   Writing catalog '+catnames[icat]
             wsex, finalcat, outfile=catnames[icat]
             splog, '   Witing DS9 region file'+radecregnames[icat]
             write_ds9_regionfile, finalcat.xwin_world, finalcat.ywin_world, $
               filename=radecregnames[icat], color='red', symbol='box'

          endfor                ; close CATNAMES loop
       endfor                   ; close SUFFIX loop
    endif

; ---------------------------------------------------------------------------    
; derive the PSF in each bandpasses

    if keyword_set(psf) then begin

; NOT SUPPORTED!!  (code is left over from VIMOS_CATALOGS)
stop       

       catnames = catpath+'sg1120_'+['B','V','R']+'_'+version+'.chi2.cat'
       splog, 'Reading '+catnames[0] & bcat = rsex(catnames[0])
       splog, 'Reading '+catnames[1] & vcat = rsex(catnames[1])
       splog, 'Reading '+catnames[2] & rcat = rsex(catnames[2],use_row=1000)
       help, bcat, vcat, rcat
       
       stars = where((rcat.class_star gt 0.9))
;      im_plothist, bcat[stars].fwhm_image*0.205
       
;      plot, bcat.mag_best, bcat.class_star, xthick=2.0, ythick=2.0, $
;        charsize=2.0, charthick=2.0, ps=4, sym=0.5, xr=[14,25], $
;        yrange=[0,1]

stop
    endif
       
; ---------------------------------------------------------------------------    
; generate some QA plots

    if keyword_set(qaplots) then begin

; NOT SUPPORTED!!  (code is left over from VIMOS_CATALOGS)
stop       
       
; get the reddening
       hdr = headfits(mosaicpath+'sg1120_r.fits')
       extast, hdr, astr
       glactc, astr.crval[0], astr.crval[1], 2000.0, $
         gl, gb, 1, /deg
       ebv = dust_getval(gl,gb,/interp)
       splog, 'SG1120 E(B-V) = ', ebv

       allbands = ['g0','r0']
       filtinfo = im_filterspecs(filterlist='sdss_'+$
         allbands+'.par',/ver)
       kl = k_lambda(filtinfo.weff,/odon)
       niceprint, allbands, kl*ebv
       
; now make some plots       
       pagemaker, nx=1, ny=1, xpage=8.5, ypage=8.2, $
         xspace=0.0, yspace=0.0, xmargin=[1.5,0.3], $
         ymargin=[0.4,1.1], width=6.7, height=6.7, $
         /normal, position=pos

       psfile = catpath+'ldss3_gr_stellar_locus.ps'
       dfpsplot, psfile, /color, xsize=8.5, ysize=8.2
       
; identify stars

       minmag = 19.0 & maxmag = 22.5
       minr = 0.3 & maxr = 0.60
       
       mag = rcat_chi2.mag_auto
       rhalf = rcat_chi2.flux_radius1*rcat_chi2.a_image*pixscale ; [arcsec]
       stars = where((mag gt minmag) and (mag lt maxmag) and $
         (rhalf gt minr) and (rhalf lt maxr),nstars)
       
;      stars = where((bcat_chi2.flags eq 0) and (vcat_chi2.flags eq 0) and $
;        (rcat_chi2.flags eq 0) and (rcat_chi2.class_star gt 0.9) and $
;        (rcat_chi2.mag_auto lt 21.0),nstars)

       djs_plot, [0], [0], /nodata, /xsty, /ysty, xthick=4.0, ythick=4.0, $
         charthick=4.0, charsize=2.0, xrange=[16.5,28], yrange=[0,2.5], $
         xtitle='R_{Vega}', ytitle='Half-light Radius (arcsec)', $
         position=pos
       djs_oplot, mag, rhalf, ps=6, sym=0.4, thick=2.0
       djs_oplot, mag[stars], rhalf[stars], ps=6, sym=0.4, $
         color='red', thick=4.0

; synthesize the stellar locus

       p98 = read_98pickles()
       filterlist = 'sdss_'+['g0','r0']+'.par'
       vega2ab = k_vega2ab(filterlist=filterlist,/kurucz,/silent)
       lambda = k_lambda_to_edges(p98[0].wave)
       maggies = fltarr(n_elements(filterlist),n_elements(p98))
       for ii = 0L, n_elements(p98)-1L do maggies[*,ii] = $
         k_project_filters(lambda,p98[ii].flux,filterlist=filterlist,/silent)
       bv_p98 = -2.5*alog10(maggies[0,*]/maggies[1,*]) - (vega2ab[0]-vega2ab[1]) ; Vega
       vr_p98 = -2.5*alog10(maggies[1,*]/maggies[2,*]) - (vega2ab[1]-vega2ab[2]) ; Vega
;      plot, p98.teff, p98.bv-bv, ps=4, xsty=3, ysty=3, yr=[-0.1,0.1]
;      plot, p98.teff, p98.vr-vr, ps=4, xsty=3, ysty=3, yr=[-0.1,0.1]

; and finally make a color-color diagram

       bv = bcat_chi2[stars].mag_auto-vcat_chi2[stars].mag_auto - (kl[0]-kl[1])*ebv
       vr = vcat_chi2[stars].mag_auto-rcat_chi2[stars].mag_auto - (kl[1]-kl[2])*ebv
       
       djs_plot, [0], [0], /nodata, /xsty, /ysty, xthick=4.0, ythick=4.0, $
         charthick=4.0, charsize=2.0, xrange=[-0.3,1.6], yrange=[-0.6,2.3], $
         xtitle='V - R', ytitle='B - V', position=pos
       djs_oplot, vr, bv, ps=6, thick=4.0
       djs_oplot, vr_p98, bv_p98, ps=6, color='red', sym=2.0, thick=4.0

       legend, ['LDSS3','Pickles'], /left, /top, box=0, charsize=2.0, $
         charthick=4.0, color=djs_icolor(['default','red']), psym=[6,6], $
         thick=4.0

; in this one add 0.03 mag to the V-band magnitudes       

       bv = bcat_chi2[stars].mag_auto-(vcat_chi2[stars].mag_auto-0.00) - (kl[0]-kl[1])*ebv
       vr = (vcat_chi2[stars].mag_auto-0.05)-rcat_chi2[stars].mag_auto - (kl[1]-kl[2])*ebv
       
       djs_plot, [0], [0], /nodata, /xsty, /ysty, xthick=4.0, ythick=4.0, $
         charthick=4.0, charsize=2.0, xrange=[-0.3,1.6], yrange=[-0.6,2.3], $
         xtitle='(V-0.05) - R', ytitle='(B-0.05) - (V-0.05)', position=pos
       djs_oplot, vr, bv, ps=6, thick=4.0
       djs_oplot, vr_p98, bv_p98, ps=6, color='red', sym=2.0, thick=4.0

       legend, ['LDSS3','Pickles'], /left, /top, box=0, charsize=2.0, $
         charthick=4.0, color=djs_icolor(['default','red']), psym=[6,6], $
         thick=4.0

       dfpsclose
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf')
stop       
    endif
       
; ---------------------------------------------------------------------------    
; match against the redshift catalog
; ---------------------------------------------------------------------------    

    if keyword_set(zmatch) then begin
    
       zcat = sg1120_read_zcat(silent=silent)

       detect_suffix = ['','.chi2']
       for idetect = 0L, 1L do begin

          catnames = catpath+'sg1120_'+band+'_'+$ ; note clean!
            version+detect_suffix[idetect]+'.cat'
          outcatnames = repstr(catnames,'.cat','.withz.cat')
          radecregname = repstr(outcatnames,'.cat','.reg')
          ncat = n_elements(catnames)

          for icat = 0L, ncat-1L do begin

             splog, 'Reading '+catnames[icat]
             cat = rsex(catnames[icat])
             ngalaxy = n_elements(cat)
             
             splog, 'Matching '+file_basename(catnames[icat])+' to the redshift catalog.'
             spherematch, cat.xwin_world, cat.ywin_world, zcat.zcat_ra, zcat.zcat_dec, $
               match_radius_zcat/3600.0, match, zcatmatch, dist12, maxmatch=1
             splog, '  Found '+strtrim(n_elements(zcatmatch),2)+' matches'

; --------------------------------------------------
; code to show which sources don't match
             nomatch = lindgen(n_elements(zcat))
             remove, zcatmatch, nomatch

             these = where((zcat[nomatch].zcat_ra gt -900.0))
;            these = where((zcat[nomatch].zcat_ra gt -900.0) and $
;              (strmatch(zcat[nomatch].zcat_catalog,'*ldss3*',/fold)))
             nomatchfile = catpath+'ldss3_nomatch_'+band[icat]+'_'+$
               version+detect_suffix[idetect]
             write_ds9_regionfile, zcat[nomatch[these]].zcat_ra, $
               zcat[nomatch[these]].zcat_dec, $
               filename=nomatchfile+'.reg', color='green'
             
             struct_print, zcat[nomatch[these]], file=nomatchfile+'.dat'
; --------------------------------------------------

             zcatmatched = im_empty_structure(zcat,ncopies=ngalaxy,$
               empty_value=-999.0,empty_string='...')
             zcatmatched[match] = zcat[zcatmatch]
             catmatched = struct_addtags(cat,zcatmatched)

; write out the catalog and region files
             
             splog, 'Writing '+outcatnames[icat]
             wsex, catmatched, outfile=outcatnames[icat]

             withz = where(catmatched.zcat_z ge 0.0,nwithz)
             if (nwithz ne 0L) then begin
                splog, 'Writing '+radecregname[icat]
                write_ds9_regionfile, catmatched[withz].xwin_world, $
                  catmatched[withz].ywin_world, filename=radecregname[icat], $
                  color='red', symbol='box'
             endif
             splog, '###########################################################################'

          endfor
       endfor
       
    endif
       
return
end
