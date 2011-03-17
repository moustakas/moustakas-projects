pro sg1120_tweak_hst_astrom, sextractor=sextractor, $
  scamp=scamp, swarp=swarp
; jm10jan29ucsd - tweak the astrometry on some HST mosaics
    
    datapath = sg1120_path()+'projects/mips/'
    allclust = file_search(datapath+'*-f???.fits',count=nclust)
;   allclust = allclust[0] & nclust = 1
    allroot = repstr(file_basename(allclust),'.fits','')
    outfile = repstr(allclust,'.fits','_tweak.fits')

    sexpath = sg1120_path(/sex)
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
; generate SE catalogs
    if keyword_set(sextractor) then begin
       for ii = 0, nclust-1 do begin
          imagelist = allclust[ii]
          catlist = repstr(imagelist,'.fits','.cat')
          radecreglist = repstr(catlist,'.cat','.reg')

;         info = im_headerforage(imagelist,ext=0) ; grab header info
       
; initialize the SE configuration parameters
          config = init_sex_config(1)
          configfile = datapath+'sex.config'
       
          config.catalog_name = catlist
          config.parameters_name = sexparam
          config.filter_name = sexconv
          config.starnnw_name = sexnnw

          config.catalog_type = 'FITS_LDAC'
          config.detect_thresh = 5.0
          config.analysis_thresh = 5.0
;         config.weight_type = 'MAP_RMS'
;         config.weight_gain = 'N'
          config.interp_type = 'NONE'
          config.nthreads = 4

          config.seeing_fwhm = 0.3
;         config.mag_zeropoint = 0.0
       
          config.checkimage_type = 'NONE' ; SEGMENTATION

          mwrfits, config, configfile+'.fits', /create

; do it!       
          t0 = systime(1)
          im_sex, imagelist, config, silent=silent ; do not pass CONFIGFILE
          splog, 'Total time to generate SE catalog = ', $
            (systime(1)-t0)/60.0, ' minutes'

; write a regions file          
          splog, 'Witing DS9 region file '+file_basename(radecreglist)
          cat = mrdfits(catlist,2)
          write_ds9_regionfile, cat.xwin_world, cat.ywin_world, $
            filename=radecreglist, color='red'
       endfor 
    endif

; ---------------------------------------------------------------------------
; SCAMP
    if keyword_set(scamp) then begin
       for ii = 3, 3 do begin
;      for ii = 0, nclust-1 do begin
          imagelist = allclust[ii]
          catlist = repstr(imagelist,'.fits','.cat')

; initialize the scamp configuration parameters
          config = init_scamp_config()
          configfile = datapath+'scamp'+allroot[ii]+'.config'

          config.astref_catalog = 'NOMAD-1'
          config.astref_band = 'B'
;         config.astref_catalog = '2MASS' 
;         config.astref_band = 'Ks' 
;         config.astref_catalog = 'USNO-B1'
;         config.astref_band = 'Bj' 

; boost the weights
;         config.astref_weight = '100.0'
;         config.crossid_radius = 1.5 ; this is important

; the OBJECT names are different between 2003 & 2006, so use this
; header tag to force scamp to define a different *astrometric*
; instrument between the two epochs; use the FILTER keyword to define
; the three (BVR) photometric instruments
;         config.astrinstru_key = 'FILTER,OBJECT'
;         config.photinstru_key = 'FILTER'
;         config.magzero_key = 'MAG0' ; 'MAGZERO'
;         config.save_refcatalog = 'N'
;         config.refout_catpath = datapath
          config.mergedoutcat_type = 'NONE'

          config.checkplot_type = strjoin(['ASTR_CHI2','ASTR_INTERROR1D','ASTR_INTERROR2D',$
            'ASTR_REFERROR1D','ASTR_REFERROR2D','DISTORTION','FGROUPS','PHOT_ERROR','PHOT_ZPCORR',$
            'PHOT_ZPCORR3D'],',')
          config.checkplot_name = strjoin(allroot[ii]+'_'+['astr_chi2','astr_interror1d',$
            'astr_interror2d','astr_referror1d','referror2d','distort','fgroups',$
            'psphot_error','phot_zpcorr','phot_zpcorr3d'],',')

          config.xml_name = datapath+allroot[ii]+'.scamp.xml'
          config.xsl_url = xslscamp

          t0 = systime(1)
          maxiter = 2

; think carefully before increasing DEGREE; for example DEGREE>3 when
; making a single mosaic results in DOOM!       
          for iter = 0, 1 do begin
             config.distort_degrees = '1'
             config.mosaic_type = 'UNCHANGED'
             config.pixscale_maxerr = '1.05'
             config.position_maxerr = '2.0'
             config.posangle_maxerr = '0.5'
             config.aheader_suffix = '.ahead'
             
             mwrfits, config, configfile+'.fits', /create
             im_scamp, catlist, config, configfile=configfile, silent=silent
          endfor  
          splog, 'Total time to run scamp = ', (systime(1)-t0)/60.0, ' minutes.'
       endfor
    endif

stop    

return
end
    
