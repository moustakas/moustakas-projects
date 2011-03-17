pro ldss3_mosaics, band=band, update_headers=update_headers, sextractor=sextractor, $
  scamp=scamp, swarp=swarp, jpeg=jpeg, test=test
; jm07aug20 - generate LDSS3 g- and r-band mosaics
; jm09mar23nyu - major update to match the latest procedures developed
;   for VIMOS

    sexpath = sg1120_path(/sex)
    mosaicpath = ldss3_path(/mosaics)
    datapath = ldss3_path(/feb06)+'sg1120/'

    if (not keyword_set(band)) then band = '[g,r]'
    inbandpass = strsplit(repstr(repstr(band,'[',''),']',''),',',/extract)
    outbandpass = ['gprime','rprime']
    nband = n_elements(inbandpass)

    stiffconfig = sexpath+'default.stiff'
    scampconfig = sexpath+'default.scamp'
    swarpconfig = sexpath+'default.swarp'
    sexconfig = sexpath+'default.sex'
    sexparam = sexpath+'sg1120.sex.param'
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'

    http = 'http://sdss.physics.nyu.edu/ioannis/research/sg1120/sex/'
    xslsex = http+'sex.xsl'
    xslscamp = http+'scamp.xsl'
    xslswarp = http+'swarp.xsl'
    
; ---------------------------------------------------------------------------    
    if keyword_set(update_headers) then begin
; add the mean atmospheric extinction coefficients, the magnitude
; zero-points, and some other photometry keywords needed by scamp to
; the header; note that all the observations were effectively
; photometric, but the *most* photometric observations (determined by
; looking at the scamp output plots) are assigned PHOTOMFLAG=T; the
; magnitude zeropoints were iteratively determined using the SDSS (see
; LDSS3_SDSS_CALIBRATE)
;
;   readcol, getenv('ISPEC_DIR')+'/etc/ctioextinct.dat', ww, aa, /silent
;   niceprint, interpol(a,w,k_lambda_eff())
       gextcoeff = 0.221 & rextcoeff = 0.142
;      gmagzero = 28.346 & rmagzero = 28.236
       gmagzero = 28.000 & rmagzero = 28.000
       imagelist = file_search(datapath+'ra.????_sg1120_?_'+$
         band+'.fits',count=nimage)
       for ii = 0L, nimage-1 do begin
          for jj = 0, 1 do begin
             hdr = headfits(imagelist[ii],ext=jj+1)
             filter = sxpar(hdr,'FILTER')
             if strlowcase(strtrim(filter,2)) eq 'g' then begin
                magzero = gmagzero
                phot_k = gextcoeff 
             endif else begin
                magzero = rmagzero
                phot_k = rextcoeff
             endelse
             sxaddpar, hdr, 'MAGZERO', float(magzero), $
               ' magnitude zero-point (counts/s)', before='HISTORY'
             sxaddpar, hdr, 'PHOT_K', phot_k, ' atmospheric '+$
               'extinction coeff (mag/airmass)', after='MAGZERO'
             if strmatch(imagelist[ii],'*ra.2011_sg1120_1_g.fits*') or $
               strmatch(imagelist[ii],'*ra.2010_sg1120_1_r.fits*') then $
                 sxaddpar, hdr, 'PHOTFLAG', 'T', ' photometric? (T/F)', $
               after='PHOT_K' else $
                   sxaddpar, hdr, 'PHOTFLAG', 'F', ' photometric? (T/F)', $
               after='PHOT_K'
             modfits, imagelist[ii], 0, hdr, ext=jj+1
             print, file_basename(imagelist[ii]), sxpar(hdr,'PHOTFLAG')
          endfor
       endfor
    endif

; ---------------------------------------------------------------------------    
; generate SE catalogs
    
    if keyword_set(sextractor) then begin
       for ib = 0, nband-1 do begin
          imagelist = file_search(datapath+'ra.????_sg1120_?_'+$
            inbandpass[ib]+'.fits',count=nimage)
          weightlist = repstr(imagelist,'.fits','.weight.fits')
          rmslist = repstr(imagelist,'.fits','.rms.fits')
          flaglist = repstr(imagelist,'.fits','.flag.fits')

          catlist = repstr(imagelist,'.fits','.cat')
          radecreglist = repstr(catlist,'.cat','.reg')
          seglist = repstr(imagelist,'.fits','.seg.fits')

          info = im_headerforage(imagelist,ext=1) ; grab header info

; initialize the SE configuration parameters
          config = init_sex_config(nimage)
          configfile = datapath+'sex.config'
          
          config.catalog_name = catlist
          config.weight_image = rmslist
          config.flag_image = flaglist
          config.parameters_name = sexparam
          config.filter_name = sexconv
          config.starnnw_name = sexnnw

          config.catalog_type = 'FITS_LDAC'
          config.detect_thresh = 1.5
          config.analysis_thresh = 1.5
          config.weight_type = 'MAP_RMS'
          config.weight_gain = 'N'
          config.interp_type = 'NONE'
          config.nthreads = 4

          config.seeing_fwhm = '1.0'
          config.mag_zeropoint = 28.0
;         config.mag_zeropoint = info.mag0 ; use the magnitude zeropoint in scamp, not here
          
          config.checkimage_type = 'NONE' ; SEGMENTATION
          config.checkimage_name = seglist

          mwrfits, config, configfile+'.fits', /create

; do it!
          t0 = systime(1)
          im_sex, imagelist, config, silent=silent ; do not pass CONFIGFILE
          splog, 'Total time to generate SE catalogs = ', $
            (systime(1)-t0)/60.0, ' minutes'

; make a QAplot and write out region files       
          sexqaplot = mosaicpath+'qaplot_'+outbandpass[ib]+'_sex.ps'
          dfpsplot, sexqaplot
          for ii = 0L, nimage-1L do begin
             for jj = 0, 1 do begin
                cat1 = mrdfits(catlist[ii],2*jj+2,/silent)
                if (jj eq 0L) then cat = cat1 else cat = [cat,cat1]
             endfor
             xr = [10,24] & yrange1 = [-0.02,1.0] & yrange2 = [0.0,4.9]
             plot, [0], [0], /nodata, position=[0.14,0.55,0.96,0.93], xsty=1, ysty=1, $
               xrange=xr, yrange=yrange1, xtitle='', xtickname=replicate(' ',10), $
               ytitle='Class Star', charsize=1.8, xthick=5.0, $
               ythick=5.0, charthick=5.0, title=file_basename(repstr(imagelist[ii],'.fits',''))
             djs_oplot, cat1.mag_auto, cat1.class_star, psym=symcat(16)
             djs_plot, [0], [0], /nodata, /noerase, position=[0.14,0.12,0.96,0.55], $
               xsty=1, ysty=1, xrange=xr, yrange=yrange2, xtitle='Instrumental Magnitude', $
               ytitle='r_{eff} (arcsec)', charsize=1.8, xthick=5.0, $
               ythick=5.0, charthick=5.0
             djs_oplot, cat1.mag_auto, cat1.flux_radius*cat1.awin_image*0.188, $
               psym=symcat(16), symsize=0.8

             splog, 'Witing DS9 region file '+file_basename(radecreglist[ii])
             write_ds9_regionfile, cat.xwin_world, cat.ywin_world, $
               filename=radecreglist[ii], color='red'
          endfor
          dfpsclose
       endfor
    endif

; ---------------------------------------------------------------------------
; SCAMP
    if keyword_set(scamp) then begin
       catlist = file_search(datapath+'ra.????_sg1120_?_'+$
         band+'.cat')

; initialize the scamp configuration parameters
       config = init_scamp_config()
       configfile = mosaicpath+'scamp.config'

; this uses the USNO catalog built for me by Dustin       
       config.astref_catalog = 'FILE'
       config.astrefcat_name = sexpath+'sg1120_sdss_dr7_refcat.cat'
;      config.astrefcat_name = sexpath+'sg1120_usnob_refcat.cat'
       config.astrefcent_keys = 'RA,DEC'
       config.astreferr_keys = 'ERR_A,ERR_B'
       config.astrefmag_key = 'MAG'

; boost the weights
       config.astref_weight = '100.0'

;      config.astref_catalog = 'NOMAD-1' & config.astref_band = 'B'
;      config.astref_catalog = '2MASS'   & config.astref_band = 'Ks' 
;      config.astref_catalog = 'USNO-B1' & config.astref_band = 'Bj' 

;      config.crossid_radius = 1.5 ; this is important

       config.astrinstru_key = 'FILTER'
       config.photinstru_key = 'FILTER'
       config.magzero_key = 'MAGZERO'
       config.extinct_key = 'PHOT_K'
       config.photomflag_key = 'PHOTFLAG'
       config.save_refcatalog = 'N'
       config.refout_catpath = mosaicpath
       config.mergedoutcat_type = 'NONE'

       config.checkplot_type = strjoin(['ASTR_CHI2','ASTR_INTERROR1D','ASTR_INTERROR2D',$
         'ASTR_REFERROR1D','ASTR_REFERROR2D','DISTORTION','FGROUPS','PHOT_ERROR','PHOT_ZPCORR',$
         'PHOT_ZPCORR3D'],',')
       config.checkplot_name = strjoin('ldss3_'+['astr_chi2','astr_interror1d',$
         'astr_interror2d','astr_referror1d','referror2d','distort','fgroups',$
         'psphot_error','phot_zpcorr','phot_zpcorr3d'],',')

       config.xml_name = mosaicpath+'ldss3.scamp.xml'
       config.xsl_url = xslscamp

       t0 = systime(1)
       maxiter = 2

; think carefully before increasing DEGREE      

       for iter = 0, maxiter do begin 
          case iter of
             0: begin
                config.distort_degrees = '1'
                config.mosaic_type = 'LOOSE'
                config.pixscale_maxerr = '1.1'
                config.position_maxerr = '5.0'
                config.posangle_maxerr = '3.0'
                config.aheader_suffix = '.ahead'
             end
             1: begin
                config.distort_degrees = '3,3'
                config.mosaic_type = 'FIX_FOCALPLANE'
                config.pixscale_maxerr = '1.1'
                config.position_maxerr = '3.0'
                config.posangle_maxerr = '1.0'
                config.aheader_suffix = '.head'
             end
             else: begin
                config.distort_degrees = '5,5' ; '4,4'
                config.mosaic_type = 'FIX_FOCALPLANE'
                config.pixscale_maxerr = '1.1'
                config.position_maxerr = '1.0'
                config.posangle_maxerr = '1.0'
                config.aheader_suffix = '.head'
             end
          endcase

          mwrfits, config, configfile+'.fits', /create
          im_scamp, catlist, config, configfile=configfile, silent=silent

       endfor 
       splog, 'Total time to run scamp = ', (systime(1)-t0)/60.0, ' minutes.'
    endif

; ---------------------------------------------------------------------------
; SWARP

    if keyword_set(swarp) then begin
; initialize the swarp configuration parameters
       config = init_swarp_config()

       if keyword_set(test) then begin
          config.center = '00:00:00.0,+00:00:00.0'
          config.center_type = 'ALL' 
          config.pixelscale_type = 'MANUAL' 
          config.pixel_scale = '0.188'
          config.image_size = '0'
          config.header_only = 'Y'
          suffix = '_test'
       endif else begin
          config.center = '11:20:14.0,-12:03:12.0'
          config.center_type = 'MANUAL' 
          config.pixelscale_type = 'MANUAL' 
          config.pixel_scale = '0.188'
          config.image_size = '5520,7220'
          config.header_only = 'N'
          suffix = ''
       endelse

       config.blank_badpixels = 'Y'
       config.interpolate = 'N'
       config.write_fileinfo = 'Y'
       config.celestial_type = 'EQUATORIAL'
       config.copy_keywords = 'OBJECT,FILTER'
       config.xml_name = mosaicpath+'ldss3.swarp.xml'
       config.xsl_url = xslswarp

; build the individual mosaics; note that using the RMS images as
; weight maps produces some funky results
       config.write_xml = 'Y'
       config.subtract_back = 'Y'
       config.resampling_type = 'LANCZOS3'
       config.weight_type = 'MAP_WEIGHT'

       config.combine_type = 'WEIGHTED'

       for ib = 0L, nband-1L do begin

          configfile = mosaicpath+'swarp_'+outbandpass[ib]+'.config'

          config.imageout_name = mosaicpath+'sg1120_'+outbandpass[ib]+'.fits'
          config.weightout_name = repstr(config.imageout_name,'.fits','.weight.fits')

          imagelist = file_search(datapath+'ra.????_sg1120_?_'+$
            inbandpass[ib]+'.fits')
          rmslist = repstr(imagelist,'.fits','.rms.fits')
          weightlist = repstr(imagelist,'.fits','.weight.fits')

          config.weight_image = strjoin(weightlist,',') ; strjoin(rmslist,',') ; 

          mwrfits, config, configfile+'.fits', /create

          splog, 'Building '+strtrim(config.imageout_name,2)
          t0 = systime(1)
          im_swarp, imagelist, config, configfile=configfile, silent=silent
          splog, 'Total time to build '+config.imageout_name+' = ', $
            (systime(1)-t0)/60.0, ' minutes.'

          if (not keyword_set(test)) then begin
             fixme = [config.imageout_name,config.weightout_name]
             for ff = 0L, n_elements(fixme)-1L do begin
                hdr = headfits(fixme[ff])
                sxaddpar, hdr, 'OBJECT', 'SG1120 '+inbandpass[ib]
                modfits, fixme[ff], 0, hdr
             endfor
          endif 
       endfor 
          
; now build the chi2 image using all the filters 
       config.combine_type = 'CHI2'
       config.copy_keywords = 'OBJECT'
          
       config.imageout_name = mosaicpath+'sg1120_gr_chi2.fits'
       config.weightout_name = repstr(config.imageout_name,'.fits','.weight.fits')

       imagelist = file_search(datapath+'ra.????_sg1120_?_'+band+'.fits')
       weightlist = repstr(imagelist,'.fits','.weight.fits')
       config.weight_image = strjoin(weightlist,',')
       
       configfile = mosaicpath+'swarp_gr_chi2.config'
       mwrfits, config, configfile+'.fits', /create
       
       splog, 'Building '+strtrim(config.imageout_name,2)
       t0 = systime(1)
       im_swarp, imagelist, config, configfile=configfile, silent=silent
       splog, 'Total time to build '+config.imageout_name+' = ', $
         (systime(1)-t0)/60.0, ' minutes.'
       
       if (not keyword_set(test)) then begin
          fixme = [config.imageout_name,config.weightout_name]
          for ff = 0L, n_elements(fixme)-1L do begin
             hdr = headfits(fixme[ff])
             sxaddpar, hdr, 'OBJECT', 'SG1120 gr chi2'
             modfits, fixme[ff], 0, hdr
          endfor
       endif
          
    endif 
       
; ---------------------------------------------------------------------------
; color mosaic
    if keyword_set(jpeg) then begin
       tifffile = mosaicpath+'sg1120_gr.tiff'
       jpegfile = mosaicpath+'sg1120_gr.jpeg'
; JPEG       
       hdr = headfits(mosaicpath+'sg1120_gprime.fits')
       imsize = [sxpar(hdr,'NAXIS1'),sxpar(hdr,'NAXIS2')]
       nx = 512*(min(imsize)/512)

       splog, 'Reading '+mosaicpath+'sg1120_gprime.fits'
       gim = (mrdfits(mosaicpath+'sg1120_gprime.fits',0,/silent))$
         [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L,$
         imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
       splog, 'Reading '+mosaicpath+'sg1120_rprime.fits'
       rim = (mrdfits(mosaicpath+'sg1120_rprime.fits',0,/silent))$
         [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L,$
         imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]

       splog, 'Writing '+jpegfile
       scales = [1.5,1.3,3.0]*1D12
       nw_rgb_make, rim, rim, gim, name=jpegfile, scales=scales, $
         nonlinearity=3.0, rebinfactor=2, quality=75
; TIFF
;      splog, 'Writing '+tifffile
;      mosaic_file = strjoin(mosaicpath+'sg1120_'+['g','r']+'.fits',' ')
;      spawn, 'stiff -c '+stiffconfig+' '+mosaic_file+' -OUTFILE_NAME '+tifffile+$
;        ' -BINNING 2 -GAMMA_FAC 0.9 -COLOUR_SAT 1.8'
    endif

return
end
