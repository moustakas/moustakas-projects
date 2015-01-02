function get_hizlrg, rz=rz, rW1=rW1, zmag=zmag, magcut=magcut
    if n_elements(magcut) eq 0 then magcut = 99.0

    rzmin = 1.6                                          
    slope1 = 1.8
    int1 = -1.0
    hiz = where(zmag lt magcut and rz gt rzmin and $
      rW1 gt poly(rz,[int1,slope1]))

;; from Jeff Newman:
;    whsel_new = where( (z lt 20.56) $
;      AND (i gt 19.9) AND $
;      ((r-z) gt 1.6) AND $
;      ((r-W1) gt (1.8*(r-z)-1)) )
return, hiz
end

function get_hizelg, gr=gr, rz=rz, rmag=rmag, fluxoii=fluxoii, $
  magcut=magcut, oiifluxcut=oiifluxcut

    if n_elements(magcut) eq 0 then magcut = 99.0
    if n_elements(oiifluxcut) eq 0 then oiifluxcut = 5D-17 ; [erg/s/cm^2]

    rzmin = 0.3                                                              
    slope1 = 1.0                                                             
    slope2 = -1.0                                                            
    int1 = -0.2                                                              
    int2 = 1.2                                                               
    hiz = where(rmag lt magcut and fluxoii gt oiifluxcut and $
      rz gt rzmin and gr lt (poly(rz,[int1,slope1]) < poly(rz,[int2,slope2])))

return, hiz
end

pro mock_desi, elgs=elgs, lrgs=lrgs, stdstars=stdstars, qaplots=qaplots
; jm14dec01siena - build a mock set of objects for the DESI
; Winter/2014 data challenge

; ToDo:
; 1) varying doublet ratio; varying line-ratios
; 2) more realistic n(z,m) distribution
; 3) including contaminating galaxies

    common com_elg, elg_info, elg_restflux, elg_restwave, elg_zvals, elg_gr, elg_rz
    common com_star, star_info, star_restflux, star_restwave, star_ugriz
    common com_lrg, lrg_info, lrg_restflux, lrg_restwave, lrg_zvals, lrg_rz, lrg_rW1

    outdir = getenv('DESI_ROOT')+'/simulations/'
    light = im_light(/km)
    
    nelg = 5000
    nlrg = 5000
    nstd = 500

    elg_zmin = 0.6
    elg_zmax = 1.6
    elg_nz = 30
    elg_rbright = 21.0
    elg_rfaint = 23.5

    lrg_zmin = 0.5
    lrg_zmax = 1.1
    lrg_nz = 30
    lrg_zbright = 19.5
    lrg_zfaint = 20.6

    std_rbright = 15.0
    std_rfaint = 19.0

    elg_filterlist = decals_filterlist()
    lrg_filterlist = ['decam_r','decam_z','wise_w1']+'.par'

; --------------------------------------------------
    if keyword_set(elgs) then begin

       version = desi_elg_templates_version()
       templatepath = getenv('DESI_ROOT')+'/spectro/templates/'+$
         'elg_templates/'+version+'/'
       
       if n_elements(elg_restflux) eq 0L then begin
          obsfile = templatepath+'elg_templates_obs_'+version+'.fits.gz'
          restfile = templatepath+'elg_templates_'+version+'.fits.gz'
          
          elg_restflux = mrdfits(restfile,0,resthdr);,range=[0,999])
          elg_restwave = 10D^make_wave(resthdr)
;         obstemplates = mrdfits(obsfile,0,obshdr)
;         obswave = make_wave(obshdr)
          elg_info = mrdfits(obsfile,1);,range=[0,999])

; K-corrections are expensive, so compute them once on a regular
; redshift grid
          k_projection_table, grzmaggies, elg_restflux, k_lambda_to_edges(elg_restwave), $
            elg_zvals, elg_filterlist, zmin=elg_zmin, zmax=elg_zmax, nz=elg_nz
          elg_gr = -2.5*alog10(grzmaggies[*,*,0]/grzmaggies[*,*,1])
          elg_rz = -2.5*alog10(grzmaggies[*,*,1]/grzmaggies[*,*,2])
       endif
       ntemplate = n_elements(elg_info)

; build the output truth table
       out_template = {$
         objtype:  'ELG',$
         z:          0.0,$
         decam_r:    0.0,$
         decam_gr:   0.0,$
         decam_rz:   0.0,$
         oii_3727:   0.0,$
         templateid: 0L}

; randomly pick NELG templates at a time until we build up the
; requisite number of templates which satisfy our grz color-color and
; [OII] flux cuts
       delvarx, outinfo
       while n_elements(outinfo) lt nelg do begin
          these = long(randomu(seed,nelg)*ntemplate)

; choose a uniform redshift grid and r-magnitude distribution 
          zobj = randomu(seed,nelg)*(elg_zmax-elg_zmin)+elg_zmin
          rmag = randomu(seed,nelg)*(elg_rbright-elg_rfaint)+elg_rfaint

          gr_these = interpolate(elg_gr,findex(elg_zvals,zobj),these)
          rz_these = interpolate(elg_rz,findex(elg_zvals,zobj),these)
          
; scale the measured [OII] flux by the r-band (re)normalization factor 
          oii_these = elg_info[these].oii_3727*$ ; [erg/s/cm^2]
            10^(-0.4*(rmag-elg_info[these].decam_r)) 

; select all objects in the grz color box with sufficient [OII] flux
          keep = get_hizelg(gr=gr_these,rz=rz_these,$
            rmag=rmag,fluxoii=oii_these)
;         djs_plot, rz_these, gr_these, psym=8
;         djs_oplot, rz_these[keep], gr_these[keep], psym=8, color='orange'

          outinfo1 = replicate(out_template,n_elements(keep))
          outinfo1.z = zobj[keep]
          outinfo1.decam_r = rmag[keep]
          outinfo1.decam_gr = gr_these[keep]
          outinfo1.decam_rz = rz_these[keep]
          outinfo1.oii_3727 = oii_these[keep]
          outinfo1.templateid = these[keep]

          if n_elements(outinfo) eq 0L then outinfo = outinfo1 else $
            outinfo = [outinfo,outinfo1]
       endwhile

; just keep the first NELG templates
       outinfo = outinfo[0L:nelg-1]

       djs_plot, elg_rz, elg_gr, psym=3
       djs_oplot, outinfo.decam_rz, outinfo.decam_gr, psym=8, color='orange'

       djs_plot, outinfo.z, 1D17*outinfo.oii_3727, psym=8, /ylog, $
         xsty=3, ysty=3, yrange=[1,500]
       djs_oplot, elg_info.z, 1D17*elg_info.oii_3727, psym=8, color='red'
       
       djs_plot, outinfo.z, outinfo.decam_r, $
         psym=8, yr=[20,24.0], xsty=3, ysty=3
       djs_oplot, elg_info.z, elg_info.decam_r, $
         psym=8, color='red'       

; now loop through each template: redshift, normalize, and write out 
       npix = n_elements(elg_restwave)
       outflux = fltarr(npix,nelg)
;      for ii = 0L, 50 do begin
        for ii = 0L, nelg-1 do begin
          outwave = elg_restwave*(1+outinfo[ii].z)
          grzobj = reform(k_project_filters(k_lambda_to_edges(outwave),$
            reform(elg_restflux[*,outinfo[ii].templateid]),filterlist=elg_filterlist))
          grobj = -2.5*alog10(grzobj[0]/grzobj[1])
          rzobj = -2.5*alog10(grzobj[1]/grzobj[2])

          outflux1 = reform(elg_restflux[*,outinfo[ii].templateid])*$
            10.0^(-0.4*outinfo[ii].decam_r)/grzobj[1]
          outflux[*,ii] = outflux1

; check that we're doing this right - yes!
          print, ii, (-2.5*alog10(k_project_filters(k_lambda_to_edges(outwave),$
            outflux1,filterlist='decam_r.par')))[0], outinfo[ii].decam_r, grobj, $
            outinfo[ii].decam_gr, rzobj, outinfo[ii].decam_rz
;         plot, outwave, outflux1, xrange=3727*(1+outinfo[ii].z)+100*[-1,1], xsty=1

; update the output structure with the actual grz colors although
; these values should match to within ~one percent
          outinfo[ii].decam_gr = grobj
          outinfo[ii].decam_rz = rzobj
       endfor

; write out
        outfile = outdir+'elg_templates.fits'

        velpixsize_hires = 20D                      ; [km/s]
        pixsize_hires = velpixsize_hires/light/alog(10) ; [pixel size in log-10 A]
 
        mkhdr, outhdr, outflux, /extend
        sxdelpar, outhdr, 'DATE'
        sxdelpar, outhdr, 'COMMENT'
;       sxaddpar, outhdr, 'VERSION', version, ' template version number'
        sxaddpar, outhdr, 'OBJTYPE', 'ELG-DEEP2', ' object type'
        sxaddpar, outhdr, 'DISPAXIS', 1, ' dispersion axis'
        sxaddpar, outhdr, 'CTYPE1', 'WAVE-WAV', ' meaning of wavelength array'
        sxaddpar, outhdr, 'CUNIT1', 'Angstrom', ' units of wavelength array'
        sxaddpar, outhdr, 'CRPIX1', 1, ' reference pixel number'
        sxaddpar, outhdr, 'CRVAL1', min(alog10(elg_restwave)), ' reference log10(Angstrom)'
        sxaddpar, outhdr, 'CDELT1', pixsize_hires, ' delta log10(Angstrom)'
        sxaddpar, outhdr, 'LOGLAM', 1, ' log10 spaced wavelengths?'
        sxaddpar, outhdr, 'WAVEMIN', min(elg_restwave), ' minimum wavelength (Angstrom)'
        sxaddpar, outhdr, 'WAVEMAX', max(elg_restwave), ' maximum wavelength (Angstrom)'
        sxaddpar, outhdr, 'WAVEUNIT', 'Angstrom', ' wavelength units'
        sxaddpar, outhdr, 'AIRORVAC', 'vac', ' wavelengths in vacuum (vac) or air'
        sxaddpar, outhdr, 'VELSCALE', velpixsize_hires, ' pixel size in km/s'
        sxaddpar, outhdr, 'BUNIT', 'erg/s/cm2/A', ' spectrum flux units'
        sxaddpar, outhdr, 'FLUXUNIT', 'erg/s/cm2/A', ' spectrum flux units'
        
; we need to do some MWRFITS jujitsu to get the metadata header right 
        units = [$
          'OBJTYPE,,object type',$
          'Z,,object redshift',$
          'DECAM_R,,r-band magnitude',$
          'DECAM_GR,,g-r color',$
          'DECAM_RZ,,r-z color',$
          'OII_3727,,',$
          'TEMPLATEID,,unique template ID number (0-indexed)']
        newheader = [$
          'EXTNAME'+" = '"+string('METADATA','(A-8)')+"'           /"]

        mwrfits, 0, outfile, /create
        mwrfits, outinfo, outfile, /silent
        metahdr =  im_update_header(outfile,units,newheader)
    
        im_mwrfits, outflux, outfile, outhdr, /clobber, /nogzip
        im_mwrfits, outinfo, outfile, metahdr, /append, /nogzip
    endif

; --------------------------------------------------
; build LRG template models
    if keyword_set(lrgs) then begin

       version = desi_lrg_templates_version()
       templatepath = getenv('DESI_ROOT')+'/spectro/templates/'+$
         'lrg_templates/'+version+'/'
       
       if n_elements(lrg_restflux) eq 0L then begin
          restfile = templatepath+'lrg_templates_'+version+'.fits.gz'
          
          lrg_restflux1 = mrdfits(restfile,0,resthdr)
          lrg_restwave = 10D^make_wave(resthdr)
          lrg_info1 = mrdfits(restfile,1)

; choose just the super-solar, >1 Gyr models
          super = where(lrg_info1.zmetal ge 0.019 and lrg_info1.age gt 1.0)
          lrg_restflux = lrg_restflux1[*,super]
          lrg_info = lrg_info1[super]

; K-corrections are expensive, so compute them once on a regular
; redshift grid
          k_projection_table, rzW1maggies, lrg_restflux, k_lambda_to_edges(lrg_restwave), $
            lrg_zvals, lrg_filterlist, zmin=lrg_zmin, zmax=lrg_zmax, nz=lrg_nz

          lrg_rz = -2.5*alog10(rzW1maggies[*,*,0]/rzW1maggies[*,*,1])
          lrg_rW1 = -2.5*alog10(rzW1maggies[*,*,0]/rzW1maggies[*,*,2])
       endif
       ntemplate = n_elements(lrg_info)

; build the output truth table
       out_template = {$
         objtype:  'LRG',$
         z:          0.0,$
         decam_z:    0.0,$
         decam_rz:   0.0,$
         decam_rW1:  0.0,$
         templateid: 0L}

; randomly pick NELG templates at a time until we build up the
; requisite number of templates which satisfy our grz color-color and
; [OII] flux cuts
       delvarx, outinfo
       while n_elements(outinfo) lt nlrg do begin
          these = long(randomu(seed,nlrg)*ntemplate)

; choose a uniform redshift grid and r-magnitude distribution 
          zobj = randomu(seed,nlrg)*(lrg_zmax-lrg_zmin)+lrg_zmin
          zmag = randomu(seed,nlrg)*(lrg_zbright-lrg_zfaint)+lrg_zfaint

          rz_these = interpolate(lrg_rz,findex(lrg_zvals,zobj),these)
          rW1_these = interpolate(lrg_rW1,findex(lrg_zvals,zobj),these)

; select all objects in the rzW1 color box
          keep = get_hizlrg(rz=rz_these,rW1=rW1_these,zmag=zmag)
;          test = get_hizlrg(rz=lrg_rz,rW1=lrg_rW1,zmag=lrg_rW1*0+20)

;          djs_plot, lrg_rz, lrg_rW1, psym=8, xrange=[0,2.5], yrange=[-2,6]
;;         djs_plot, rz_these, rW1_these, psym=8
;;         djs_oplot, lrg_rz[test], lrg_rW1[test], psym=8, color='green'
;          djs_oplot, rz_these[keep], rW1_these[keep], psym=8, color='orange'
;          rzaxis = range(1.6,2.5,50)
;          djs_oplot, rzaxis, poly(rzaxis,[-1.5,2.0]), color='blue'
;          djs_oplot, [1.6,1.6], [poly(1.6,[-1.5,2.0]),!y.crange[1]], color='blue'

          outinfo1 = replicate(out_template,n_elements(keep))
          outinfo1.z = zobj[keep]
          outinfo1.decam_z = zmag[keep]
          outinfo1.decam_rz = rz_these[keep]
          outinfo1.decam_rW1 = rW1_these[keep]
          outinfo1.templateid = these[keep]

          if n_elements(outinfo) eq 0L then outinfo = outinfo1 else $
            outinfo = [outinfo,outinfo1]
       endwhile

; just keep the first NLRG templates
       outinfo = outinfo[0L:nlrg-1]

; now loop through each template: redshift, normalize, and write out 
       npix = n_elements(lrg_restwave)
       outflux = fltarr(npix,nlrg)
        for ii = 0L, nlrg-1 do begin
          outwave = lrg_restwave*(1+outinfo[ii].z)
          rzW1obj = reform(k_project_filters(k_lambda_to_edges(outwave),$
            reform(lrg_restflux[*,outinfo[ii].templateid]),filterlist=lrg_filterlist))
          rzobj = -2.5*alog10(rzW1obj[0]/rzW1obj[1])
          rW1obj = -2.5*alog10(rzW1obj[0]/rzW1obj[2])

          outflux1 = reform(lrg_restflux[*,outinfo[ii].templateid])*$
            10.0^(-0.4*outinfo[ii].decam_z)/rzW1obj[1]
          outflux[*,ii] = outflux1

; check that we're doing this right - yes!
          print, ii, (-2.5*alog10(k_project_filters(k_lambda_to_edges(outwave),$
            outflux1,filterlist='decam_z.par')))[0], outinfo[ii].decam_z, rzobj, $
            outinfo[ii].decam_rz, rW1obj, outinfo[ii].decam_rW1
;         plot, outwave, outflux1, xsty=3, ysty=3, /xlog
;         cc = get_kbrd(1)

; update the output structure with the actual grz colors although
; these values should match to within ~one percent
          outinfo[ii].decam_rz = rzobj
          outinfo[ii].decam_rW1 = rW1obj
       endfor

; write out
        outfile = outdir+'lrg_templates.fits'

        velpixsize = 10D                      ; [km/s]
        pixsize= velpixsize/light/alog(10) ; [pixel size in log-10 A]
 
        mkhdr, outhdr, outflux, /extend
        sxdelpar, outhdr, 'DATE'
        sxdelpar, outhdr, 'COMMENT'
;       sxaddpar, outhdr, 'VERSION', version, ' template version number'
        sxaddpar, outhdr, 'OBJTYPE', 'LRG', ' object type'
        sxaddpar, outhdr, 'DISPAXIS', 1, ' dispersion axis'
        sxaddpar, outhdr, 'CTYPE1', 'WAVE-WAV', ' meaning of wavelength array'
        sxaddpar, outhdr, 'CUNIT1', 'Angstrom', ' units of wavelength array'
        sxaddpar, outhdr, 'CRPIX1', 1, ' reference pixel number'
        sxaddpar, outhdr, 'CRVAL1', min(alog10(lrg_restwave)), ' reference log10(Angstrom)'
        sxaddpar, outhdr, 'CDELT1', pixsize, ' delta log10(Angstrom)'
        sxaddpar, outhdr, 'LOGLAM', 1, ' log10 spaced wavelengths?'
        sxaddpar, outhdr, 'WAVEMIN', min(lrg_restwave), ' minimum wavelength (Angstrom)'
        sxaddpar, outhdr, 'WAVEMAX', max(lrg_restwave), ' maximum wavelength (Angstrom)'
        sxaddpar, outhdr, 'WAVEUNIT', 'Angstrom', ' wavelength units'
        sxaddpar, outhdr, 'AIRORVAC', 'vac', ' wavelengths in vacuum (vac) or air'
        sxaddpar, outhdr, 'VELSCALE', velpixsize, ' pixel size in km/s'
        sxaddpar, outhdr, 'BUNIT', 'erg/s/cm2/A', ' spectrum flux units'
        sxaddpar, outhdr, 'FLUXUNIT', 'erg/s/cm2/A', ' spectrum flux units'
        
; we need to do some MWRFITS jujitsu to get the metadata header right 
        units = [$
          'OBJTYPE,,object type',$
          'Z,,object redshift',$
          'DECAM_z,,z-band magnitude',$
          'DECAM_RZ,,r-z color',$
          'DECAM_RW1,,r-W1 color',$
          'TEMPLATEID,,unique template ID number (0-indexed)']
        newheader = [$
          'EXTNAME'+" = '"+string('METADATA','(A-8)')+"'           /"]

        mwrfits, 0, outfile, /create
        mwrfits, outinfo, outfile, /silent
        metahdr =  im_update_header(outfile,units,newheader)
    
        im_mwrfits, outflux, outfile, outhdr, /clobber, /nogzip
        im_mwrfits, outinfo, outfile, metahdr, /append, /nogzip
     endif 

; --------------------------------------------------
; build standard-star templates
    if keyword_set(stdstars) then begin

       version = desi_stellar_templates_version()
       templatepath = getenv('DESI_ROOT')+'/spectro/templates/'+$
         'stellar_templates/'+version+'/'
       
       if n_elements(star_restflux) eq 0L then begin
          restfile = templatepath+'stellar_templates_'+version+'.fits.gz'
          
          star_restflux = mrdfits(restfile,0,resthdr)
          star_restwave = 10D^make_wave(resthdr)
          star_info = mrdfits(restfile,1)

; synthesize ugriz photometry; this isn't quite right because
; we should be just considering DECam/grz photometry
          star_ugriz = fltarr(5,n_elements(star_info))
          for ii = 0, n_elements(star_info)-1 do star_ugriz[*,ii] = $
            reform(k_project_filters(k_lambda_to_edges(star_restwave),$
            star_restflux[*,ii],filterlist=sdss_filterlist))
          star_ugriz = -2.5*alog10(star_ugriz)
       endif
       ntemplate = n_elements(star_info)
       
; define a distance in color space; see equation 4 in Dawson+13:
       ug = reform(star_ugriz[0,*]-star_ugriz[1,*])
       gr = reform(star_ugriz[1,*]-star_ugriz[2,*])
       ri = reform(star_ugriz[2,*]-star_ugriz[3,*])
       iz = reform(star_ugriz[3,*]-star_ugriz[4,*])

       mdist = sqrt((ug-0.82)^2 + (gr-0.30)^2 + (ri-0.09)^2 + (iz-0.02)^2)

; just 12 stars satisfy this criterion!  but they roughly match the
; distribution of physical properties of the BOSS standards
       keep = where(mdist lt 0.08,nkeep)

       boss = mrdfits(outdir+'boss-stds.fits',1)
       splog, median(boss.teff), median(boss.logg), median(boss.feh)
       struct_print, star_info[keep]

; assign r-band magnitudes uniformly
       rmag = randomu(seed,nstd)*(std_rbright-std_rfaint)+std_rfaint
       these = long(randomu(seed,nstd)*nkeep)

; build the output truth table
       outinfo = {$
         objtype:  'STD',$
         sdss_r:     0.0,$
         sdss_ug:    0.0,$
         sdss_gr:    0.0,$
         sdss_ri:    0.0,$
         sdss_iz:    0.0,$
         templateid: 0L}
       outinfo = replicate(outinfo,nstd)
       outinfo.sdss_r = rmag
       outinfo.sdss_ug = ug[keep[these]]
       outinfo.sdss_gr = gr[keep[these]]
       outinfo.sdss_ri = ri[keep[these]]
       outinfo.sdss_iz = iz[keep[these]]
       outinfo.templateid = keep[these]

; normalize each spectrum to the desired r-magnitude and write out
       restwave = k_lambda_to_edges(star_restwave)
       npix = n_elements(star_restwave)
       outflux = fltarr(npix,nstd)
       for ii = 0L, nstd-1 do begin
          robj = (reform(k_project_filters(restwave,$
            star_restflux[*,outinfo[ii].templateid],$
            filterlist='decam_r.par')))[0]

          outflux1 = reform(star_restflux[*,outinfo[ii].templateid])*$
            10.0^(-0.4*outinfo[ii].sdss_r)/robj
          outflux[*,ii] = outflux1

;         djs_plot, star_restwave, outflux1, xsty=1, ysty=1 ;, /xlog
;         struct_print, star_info[outinfo[ii].templateid], /no_head
;         cc = get_kbrd(1)
       endfor

; write out
        outfile = outdir+'std_templates.fits'

        cdelt1 = 1.4341626783D-05          ; [pixel size in log-10 A]
        velpixsize = cdelt1*light*alog(10) ; [km/s]
 
        mkhdr, outhdr, outflux, /extend
        sxdelpar, outhdr, 'DATE'
        sxdelpar, outhdr, 'COMMENT'
;       sxaddpar, outhdr, 'VERSION', version, ' template version number'
        sxaddpar, outhdr, 'OBJTYPE', 'STD', ' object type'
        sxaddpar, outhdr, 'DISPAXIS', 1, ' dispersion axis'
        sxaddpar, outhdr, 'CTYPE1', 'WAVE-WAV', ' meaning of wavelength array'
        sxaddpar, outhdr, 'CUNIT1', 'Angstrom', ' units of wavelength array'
        sxaddpar, outhdr, 'CRPIX1', 1, ' reference pixel number'
        sxaddpar, outhdr, 'CRVAL1', min(alog10(star_restwave)), ' reference log10(Angstrom)'
        sxaddpar, outhdr, 'CDELT1', cdelt1, ' delta log10(Angstrom)'
        sxaddpar, outhdr, 'LOGLAM', 1, ' log10 spaced wavelengths?'
        sxaddpar, outhdr, 'WAVEMIN', min(star_restwave), ' minimum wavelength (Angstrom)'
        sxaddpar, outhdr, 'WAVEMAX', max(star_restwave), ' maximum wavelength (Angstrom)'
        sxaddpar, outhdr, 'WAVEUNIT', 'Angstrom', ' wavelength units'
        sxaddpar, outhdr, 'AIRORVAC', 'vac', ' wavelengths in vacuum (vac) or air'
        sxaddpar, outhdr, 'VELSCALE', velpixsize, ' pixel size in km/s'
        sxaddpar, outhdr, 'BUNIT', 'erg/s/cm2/A', ' spectrum flux units'
        sxaddpar, outhdr, 'FLUXUNIT', 'erg/s/cm2/A', ' spectrum flux units'
        
; we need to do some MWRFITS jujitsu to get the metadata header right 
        units = [$
          'OBJTYPE,,object type',$
          'SDSS_R,,r-band magnitude',$
          'SDSS_UG,,u-g color',$
          'SDSS_GR,,g-r color',$
          'SDSS_RI,,r-i color',$
          'SDSS_IZ,,i-z color',$
          'TEMPLATEID,,unique template ID number (0-indexed)']
        newheader = [$
          'EXTNAME'+" = '"+string('METADATA','(A-8)')+"'           /"]

        mwrfits, 0, outfile, /create
        mwrfits, outinfo, outfile, /silent
        metahdr =  im_update_header(outfile,units,newheader)
    
        im_mwrfits, outflux, outfile, outhdr, /clobber, /nogzip
        im_mwrfits, outinfo, outfile, metahdr, /append, /nogzip
    endif 

; --------------------------------------------------
; make some QAplots
    if keyword_set(qaplots) then begin

; ####################
; LRGs
       if n_elements(lrg_info) eq 0L then begin
          version = desi_lrg_templates_version()
          templatepath = getenv('DESI_ROOT')+'/spectro/templates/'+$
            'lrg_templates/'+version+'/'

          restfile = templatepath+'lrg_templates_'+version+'.fits.gz'
          
          lrg_restflux1 = mrdfits(restfile,0,resthdr)
          lrg_restwave = 10D^make_wave(resthdr)
          lrg_info1 = mrdfits(restfile,1)

;; for the plots don't cut on age & metallicity!
;          super = where(lrg_info1.zmetal ge 0.019 and lrg_info1.age gt 1.0)
;          lrg_restflux = lrg_restflux1[*,super]
;          lrg_info = lrg_info1[super]

          lrg_restflux = lrg_restflux1
          lrg_info = lrg_info1

; K-corrections are expensive, so compute them once on a regular
; redshift grid
          k_projection_table, rzW1maggies, lrg_restflux, k_lambda_to_edges(lrg_restwave), $
            lrg_zvals, lrg_filterlist, zmin=lrg_zmin, zmax=lrg_zmax, nz=lrg_nz

          lrg_rz = -2.5*alog10(rzW1maggies[*,*,0]/rzW1maggies[*,*,1])
          lrg_rW1 = -2.5*alog10(rzW1maggies[*,*,0]/rzW1maggies[*,*,2])
       endif

       lrg_tempinfo = mrdfits(outdir+'lrg_templates.fits',1)

       psfile = outdir+'qa_lrg.ps'
       im_plotconfig, 0, pos1, psfile=psfile, height=4.0, ymargin=[0.4,6.6], xmargin=[1.2,0.3]
       djs_plot, [0], [0], /nodata, position=pos1, $
         xsty=1, ysty=1, xrange=[0.0,2.5], yrange=[-2,6], $
         xtitle='r - z', ytitle='r - W1'
       djs_oplot, lrg_rz, lrg_rW1, psym=symcat(16), $
         symsize=0.4, color=cgcolor('grey')
       djs_oplot, lrg_tempinfo.decam_rz, lrg_tempinfo.decam_rW1, $
         psym=symcat(9), color=cgcolor('dodger blue'), symsize=0.3

       rzaxis = range(1.6,2.5,50)
       djs_oplot, rzaxis, poly(rzaxis,[-1.5,2.0])
       djs_oplot, [1.6,1.6], [poly(1.6,[-1.5,2.0]),!y.crange[1]]

;      im_legend, ['Z/Z_{\odot}\ge1','SSP Age < 1 Gyr'], /left, /top, box=0, $
;        margin=0, charsize=1.8
;      im_legend, ['Z/Z_{\odot}\ge1','SSP Age < 1 Gyr'], /left, /top, box=0, $
;        margin=0, charsize=1.8

; magnitude vs redshift
       im_plotconfig, 0, pos2, ymargin=[6.0,1.1], height=4.5, xmargin=[1.2,0.3]
;      im_plotconfig, 12, pos2, ymargin=[6.0,1.1], height=4.5, xspace=1.0
       djs_plot, [0], [0], /nodata, /noerase, position=pos2[*,0], $
         xsty=1, ysty=1, xrange=[lrg_zmin*0.95,lrg_zmax*1.05], $
         yrange=[lrg_zbright-0.2,lrg_zfaint+0.2], $
         xtitle='Redshift', ytitle='z (AB mag)', ytickinterval=0.5
       djs_oplot, lrg_tempinfo.z, lrg_tempinfo.decam_z, psym=symcat(16), $
         symsize=0.3, color=cgcolor('dodger blue')

       im_plotconfig, psfile=psfile, /psclose, /pdf

stop

; ####################
; ELGs
       if n_elements(elg_info) eq 0L then begin
          version = desi_elg_templates_version()
          templatepath = getenv('DESI_ROOT')+'/spectro/templates/'+$
         'elg_templates/'+version+'/'
          obsfile = templatepath+'elg_templates_obs_'+version+'.fits.gz'
          elg_info = mrdfits(obsfile,1);,range=[0,999])
       endif

       elg_tempinfo = mrdfits(outdir+'elg_templates.fits',1)

       psfile = outdir+'qa_elg.ps'
       im_plotconfig, 0, pos1, psfile=psfile, height=4.0, $
         ymargin=[0.4,6.6], xmargin=[1.2,0.3]
       djs_plot, [0], [0], /nodata, position=pos1, $
         xsty=1, ysty=1, xrange=[-0.5,2], yrange=[-0.3,2], $
         xtitle='r - z', ytitle='g - r'
       djs_oplot, elg_info.decam_r-elg_info.decam_z, $
         elg_info.decam_g-elg_info.decam_r, psym=symcat(16), $
         symsize=0.4, color=cgcolor('grey')
       djs_oplot, elg_tempinfo.decam_rz, elg_tempinfo.decam_gr, $
         psym=symcat(9), color=cgcolor('dodger blue'), symsize=0.3

       rzmin = 0.3
       rzaxis = range(rzmin,1.5,100)
       slope1 = 1.0
       slope2 = -1.0
       int1 = -0.2
       int2 = 1.2
       
       djs_oplot, rzmin*[1,1], [!y.crange[0],poly(rzmin,[int1,slope1])], thick=8
       djs_oplot, rzaxis, poly(rzaxis,[int1,slope1])<poly(rzaxis,[int2,slope2]), thick=8

; magnitude vs redshift & [OII] flux vs redshift
       im_plotconfig, 12, pos2, ymargin=[6.0,1.1], height=4.5, $
         xspace=1.0, xmargin=[1.2,0.3]

       djs_plot, [0], [0], /nodata, /noerase, position=pos2[*,0], $
         xsty=1, ysty=1, xrange=[0.5,1.7], yrange=[0.3,20], $
         xtickinterval=0.4, /ylog, xtitle='Redshift', $
         ytitle='[O II] \lambda\lambda3726,29 (10^{-16} erg s^{-1} cm^{-2})'
       djs_oplot, elg_tempinfo.z, 1D16*elg_tempinfo.oii_3727, psym=symcat(16), $
         symsize=0.3, color=cgcolor('dodger blue')

       djs_plot, [0], [0], /nodata, /noerase, position=pos2[*,1], $
         xsty=1, ysty=1, xrange=[0.5,1.7], yrange=[elg_rbright-0.4,elg_rfaint+0.2], $
         xtitle='Redshift', ytitle='r (AB mag)', ytickinterval=1, $
         xtickinterval=0.4
       djs_oplot, elg_tempinfo.z, elg_tempinfo.decam_r, psym=symcat(16), $
         symsize=0.3, color=cgcolor('dodger blue')
       
       im_plotconfig, psfile=psfile, /psclose, /pdf

stop

; ####################
; standard stars
       if n_elements(star_restflux) eq 0L then begin
          version = desi_stellar_templates_version()
          templatepath = getenv('DESI_ROOT')+'/spectro/templates/'+$
            'stellar_templates/'+version+'/'
       
          restfile = templatepath+'stellar_templates_'+version+'.fits.gz'
          
          star_restflux = mrdfits(restfile,0,resthdr)
          star_restwave = 10D^make_wave(resthdr)
          star_info = mrdfits(restfile,1)

; synthesize ugriz photometry; this isn't quite right because
; we should be just considering DECam/grz photometry
          star_ugriz = fltarr(5,n_elements(star_info))
          for ii = 0, n_elements(star_info)-1 do star_ugriz[*,ii] = $
            reform(k_project_filters(k_lambda_to_edges(star_restwave),$
            star_restflux[*,ii],filterlist=sdss_filterlist))
          star_ugriz = -2.5*alog10(star_ugriz)
       endif

       std_tempinfo = mrdfits(outdir+'std_templates.fits',1)
       boss = mrdfits(outdir+'boss-stds.fits',1)

       psfile = outdir+'qa_std.ps'
       im_plotconfig, 0, pos1, psfile=psfile, height=4.0, $
         ymargin=[0.4,6.6], xmargin=[1.2,0.3]
       djs_plot, [0], [0], /nodata, position=pos1, $
         xsty=1, ysty=1, xrange=[-1,3], yrange=[-0.7,2], $
         xtitle='r - z', ytitle='g - r'
       djs_oplot, star_ugriz[2,*]-star_ugriz[4,*], $
         star_ugriz[1,*]-star_ugriz[2,*], psym=symcat(16), $
         symsize=0.7, color=cgcolor('grey')
       djs_oplot, std_tempinfo.sdss_ri+std_tempinfo.sdss_iz, $
         std_tempinfo.sdss_gr, psym=symcat(9), $
         color=cgcolor('dodger blue'), symsize=0.7

       im_legend, '[(u-g)-0.82]^2 + [(g-r)-0.30]^2 + !c!c '+$
         '[(r-i)-0.09]^2 + [(i-z)-0.02]^2 < 0.08', $
         /right, /bottom, box=0, charsize=1.5, margin=2

; magnitude histogram
       im_plotconfig, 12, pos2, ymargin=[6.0,1.1], height=4.5, $
         xspace=1.0, xmargin=[1.2,0.3]
       djs_plot, [0], [0], /nodata, /noerase, position=pos2[*,0], $
         xsty=1, ysty=1, yrange=[0,1.05], $
         xrange=[std_rbright-0.4,std_rfaint+0.2], $
         ytitle='Relative Number', xtitle='r (AB mag)'
       im_plothist, std_tempinfo.sdss_r, bin=0.2, /overplot, /peak

       djs_plot, [0], [0], /nodata, /noerase, position=pos2[*,1], $
         xsty=1, ysty=1, yrange=[4.5,8], xrange=[-3.0,0.5], $
         ytitle='T_{eff} (1000 K)', xtitle='[Fe/H]', $
         xtickinterval=1.0, ytickinterval=1.0
       djs_oplot, boss.feh, boss.teff/1000, psym=6, symsize=1.0
       djs_oplot, star_info[std_tempinfo.templateid].feh, $
         star_info[std_tempinfo.templateid].teff/1000.0, $
         psym=symcat(16), color=cgcolor('dodger blue'), symsize=1.0
       im_legend, ['BOSS Standards','STD Templates'], /left, /top, box=0, $
         psym=[6,16], color=['black','dodger blue'], charsize=1.4, $
         margin=0

       im_plotconfig, psfile=psfile, /psclose, /pdf
       
stop

    endif


stop

return
end
