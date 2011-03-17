pro unpack_ages, fluxed=fluxed, cleanfits=cleanfits, debug=debug, clobber=clobber
; jm06jan19uofa - prepare the SDSS-reduced spectra for sky subtraction
;                 using Wild's PCA algorithm; write FITS structures
;                 for each plate to SDSS_before_skysubpca
; jm06sep18nyu - updated
    
    rawpath = ages_path(/rawdata)
    outpath = ages_path(/spec1d)
    if keyword_set(fluxed) then begin
       outpath = outpath+'fluxed/before_skysubpca/'
       sproot = 'spHect'
       scale = 1d-17
    endif else begin
       outpath = outpath+'unfluxed/before_skysubpca/'
       sproot = 'spObs'
       scale = 1.0d
    endelse
    analysis_path = ages_path(/analysis)
    catalogs_path = ages_path(/catalogs)

    vredux_11 = '0301'
    vredux_20 = '0051'
    vredux = [vredux_11,vredux_20]

    datapath_11 = rawpath+vredux_11+'/'
    datapath_20 = rawpath+vredux_20+'/'
    datapath = [datapath_11,datapath_20]

    if keyword_set(cleanfits) then begin
       if (file_test(outpath,/directory) ne 0l) then begin
          splog, 'Remove all FITS files from '+outpath+' [y/n]?'
          cc = get_kbrd(1)
          if (strupcase(cc) eq 'y') then spawn, ['/bin/rm '+outpath+'*.fits*'], /sh
       endif
    endif

    badplates = [$
      '106',$ ; not fluxed
      '110',$ ; not fluxed
      '209',$ ; not fluxed
      '310',$ ; not fluxed
      '311' $ ; not fluxed
      ]

    badfluxplates = [$
      '104',$ ; "average" flux vector; not correct in detail
      '106',$ ; not fluxed
      '110',$ ; not fluxed
      '209',$ ; not fluxed
      '301',$ ; shape of standards "looks wrong"
      '302',$ ; adc not turned on and am=1.300
      '304',$ ; very poor fluxing ("king of bad plates")
      '310',$ ; not fluxed
      '311',$ ; not fluxed
      '313',$ ; adc not turned on (am=?)
      '314',$ ; adc not turned on (am=?)
      '315' $ ; adc not turned on (am=1.153)
      ]
    
; read the zmerge catalog to get PASS and APER, and the SPECTWEIGHT
; catalog to get the official set of redshifts, which includes
; redshifts from the SDSS
    zmerge = mrdfits(analysis_path+'catalog.zmerge.fits.gz',1,/silent)
;   weight = mrdfits(analysis_path+'catalog.spectweight.fits.gz',1,/silent)
;   zmerge.z = weight.redshift

; read the airmass file
;   readcol, catalogs_path+'airmass.dat', apass, afield, am, /silent
;   am_pass = string(apass,format='(i0)')+string(afield,format='(i2.2)')

    target_info_template = {$ ; output data structure
      ages_id:          0l,$
      galaxy:           '',$
      ra:             0.0d,$
      dec:            0.0d,$
      pass:              0,$
      aper:              0,$
;     airmass:        -1.0,$
      class:            '',$
      z:               0.0,$
      spzbest:         0.0,$
;     fluxing_problem:   0,$
      zmerge_multiple:   0,$
      minwave:         0.0,$
      maxwave:         0.0}

    if keyword_set(debug) then im_window, 0, xratio=0.5, /square
    if (keyword_set(debug) eq 0) then splog, filename=outpath+'unpack_ages.log'

    t0 = systime(1)
;   for ipath = 0, 0 do begin
;   for ipath = 1, 1 do begin
    for ipath = 0, 1 do begin
       
       sphectfits = file_basename(file_search(datapath[ipath]+sproot+'-???-'+$
         vredux[ipath]+'.fits*',count=npassfield))
       spzbestfits = file_basename(file_search(datapath[ipath]+'spZbest-???-'+$
         vredux[ipath]+'.fits*'))
       passfield = strmid(sphectfits,strlen(sproot)+1,3) ; NOT GENERAL

;      for jpass = 1, 1 do begin
       for jpass = 0, npassfield-1 do begin
          t1 = systime(1)

;         thisoutpath = outpath+passfield[jpass]
;         if (file_test(thisoutpath,/dir) eq 0) then $
;           spawn, 'mkdir '+thisoutpath, /sh
          
          splog, 'Unpacking PASSFIELD '+passfield[jpass]
          pass = long(passfield[jpass])

; read the data          
          bigwave   = mrdfits(datapath[ipath]+sphectfits[jpass],0,/silent)
          bigflux   = mrdfits(datapath[ipath]+sphectfits[jpass],1,/silent)
          biginvvar = mrdfits(datapath[ipath]+sphectfits[jpass],2,/silent)
          plugmap   = mrdfits(datapath[ipath]+sphectfits[jpass],5,/silent)
          spzbest   = mrdfits(datapath[ipath]+spzbestfits[jpass],1,/silent)

; identify TARGET and SKY fibers
          targetfibers = where((strtrim(plugmap.objtype,2) eq 'TARGET'),ntargetfibers)
          splog, '   Identified '+string(ntargetfibers,format='(I0)')+' target fibers.'
          skyfibers = where(strtrim(plugmap.objtype,2) eq 'SKY',nskyfibers)
          splog, '   Identified '+string(nskyfibers,format='(I0)')+' sky fibers.'

          spzbest_target = spzbest[targetfibers]

          target_info = replicate(target_info_template,ntargetfibers)
          target_info.pass            = pass
          target_info.aper            = targetfibers+1L
          target_info.galaxy          = 'ages_'+string(target_info.pass,for='(I3.3)')+'/'+$
            string(target_info.aper,for='(I3.3)')
          target_info.class           = spzbest_target.class
          target_info.spzbest         = spzbest_target.z
;         target_info.fluxing_problem = total(badfluxplates eq passfield[jpass])

; select the minimum and maximum wavelengths and the dispersion here;
; always shave off around 10 Angstroms (~8 pixels) from the red end
          maxwave = 9190.0D     ; 9150.0D ; 9199.6D ; max(wave)
          dwave = 1.2D
          case passfield[jpass] of
             '203': minwave = 4050.0D
             '205': minwave = 3850.0D
             '208': minwave = 3950.0D
             '213': minwave = 3850.0D
             '304': minwave = 3850.0D
             '309': minwave = 4150.0D
             else: minwave = 3700.0D ; 3740.0D
          endcase
          newwave = dindgen((maxwave-minwave)/dwave+1)*dwave+minwave ; output wavelength array
          
; loop on each TARGET fiber
;         for ifiber = 133, 134 do begin
          for ifiber = 0L, ntargetfibers-1L do begin

             wave = bigwave[*,targetfibers[ifiber]]
             flux = bigflux[*,targetfibers[ifiber]]
             invvar = biginvvar[*,targetfibers[ifiber]]
             npix = n_elements(wave)

             zero = where(invvar le 0.0,nzero)
             if (nzero eq npix) then begin ; special case of no data!
                splog, '      No data in TARGET fiber '+$
                  string(targetfibers[ifiber],format='(I0)')+'!'
                invvar = invvar*0.0
             endif

; grab the official AGES redshift from the ZMERGE catalog
             zmatch = where((pass eq zmerge.pass) and (targetfibers[ifiber]+1L eq zmerge.aper),nzmatch)
             if (nzmatch ne 0L) then begin
                if (nzmatch gt 1L) then begin
                   splog, '      Multiple matches for fiber '+string(targetfibers[ifiber]+1L,format='(I0)')+'.'
                   target_info[ifiber].zmerge_multiple = nzmatch
;                  print, zmatch & struct_print, struct_trimtags(zmerge[zmatch],except='*_N*'), /no_head
                   zdmin = min(zmerge[zmatch].dmin,zdminindx)
                   zmatch = zmatch[zdminindx]
                endif
                target_info[ifiber].ages_id = zmatch
                target_info[ifiber].z = zmerge[zmatch].z
                target_info[ifiber].ra = zmerge[zmatch].ra
                target_info[ifiber].dec = zmerge[zmatch].dec
             endif

; linearly resample; some pixels end up with zero error if the input
; INVVAR=0; set the error in those pixels (including extrapolated
; ones) to a large number
             var = 1.0/(invvar+(invvar le 0.0))*(invvar gt 0.0)
             newflux = rebin_spectrum(flux,wave,newwave)
             newvar = rebin_spectrum(var,wave,newwave)
;            x_specrebin, wave, flux, newwave, newflux, $
;              var=var, nwvar=newvar, /silent
             newflux = scale*newflux
             newvar = newvar*(newvar gt 0.0) ; enforce positivity
             newferr = scale*sqrt(newvar)

             zero = where(newferr eq 0.0,nzero,comp=good,ncomp=ngood)
             if (nzero ne 0) then newferr[zero] = 1E16
             if (ngood ne 0) then begin
                target_info[ifiber].minwave = min(newwave[good])
                target_info[ifiber].maxwave = max(newwave[good])
             endif

             inf = where((newvar lt 0.0) or (finite(invvar) eq 0) or $
               (finite(flux) eq 0),ninf)
             if (ninf ne 0L) then message, 'Help me!'

             if (ifiber eq 0L) then begin
                targetbigwave = newwave
                targetbigflux = newflux
                targetbigferr = newferr
             endif else begin
;               targetbigwave = [ [targetbigwave], [newwave] ]
                targetbigflux = [ [targetbigflux], [newflux] ]
                targetbigferr = [ [targetbigferr], [newferr] ]
             endelse

             if keyword_set(debug) then begin
                noextrap = where(newvar ne 0.0,nnoextrap)
                if (nnoextrap ne 0L) then begin
                   nsmooth = 7L
                   djs_plot, targetbigwave[noextrap], (1.0/scale)*smooth(targetbigflux[noextrap,ifiber],nsmooth), $
                     ps=10, xsty=3, ysty=3, xthick=2.0, ythick=2.0, charsize=1.5, $
                     charthick=2.0, position=[0.1,0.6,0.95,0.95], xtickname=replicate(' ',10), $
                     title='Pass = '+string(pass,format='(I0)')+', Fiber = '+string(ifiber+1L,format='(I3.3)')
                   djs_plot, targetbigwave[noextrap], (1.0/scale)*smooth(targetbigferr[noextrap,ifiber],nsmooth), $
                     ps=10, xsty=3, ysty=3, xthick=2.0, ythick=2.0, charsize=1.5, $
                     charthick=2.0, position=[0.1,0.1,0.95,0.6], /noerase
                   print, minmax(targetbigwave[noextrap])
                   cc = get_kbrd(1)
                endif
             endif
          endfor ; close FIBER loop

; loop on each SKY fiber
          for ifiber = 0L, nskyfibers-1L do begin
             wave = bigwave[*,skyfibers[ifiber]]
             flux = bigflux[*,skyfibers[ifiber]]
             invvar = biginvvar[*,skyfibers[ifiber]]
             npix = n_elements(wave)

             zero = where(invvar le 0.0,nzero)
             if (nzero eq npix) then begin ; special case of no data!
                splog, '      No data in SKY fiber '+$
                  string(skyfibers[ifiber],format='(I0)')+'!'
                invvar = invvar*0.0
             endif

; linearly resample; set the error of extrapolated pixels to zero
             var = 1.0/(invvar+(invvar le 0.0))*(invvar gt 0.0)
             newflux = rebin_spectrum(flux,wave,newwave)
             newvar = rebin_spectrum(var,wave,newwave)
;            x_specrebin, wave, flux, newwave, newflux, $
;              var=var, nwvar=newvar, /silent
             newflux = scale*newflux
             newvar = newvar*(newvar gt 0.0) ; enforce positivity
             newferr = scale*sqrt(newvar)
       
             zero = where(newferr eq 0.0,nzero)
             if (nzero ne 0) then newferr[zero] = 1E16

             inf = where((finite(invvar) eq 0) or (finite(flux) eq 0),ninf)
             if (ninf ne 0L) then message, 'Help me!'
             
             if (ifiber eq 0L) then begin
                skybigwave = newwave
                skybigflux = newflux
                skybigferr = newferr
             endif else begin
;               skybigwave = [ [skybigwave], [newwave] ]
                skybigflux = [ [skybigflux], [newflux] ]
                skybigferr = [ [skybigferr], [newferr] ]
             endelse

             if keyword_set(debug) then begin
                djs_plot, skybigwave, smooth(skybigflux[*,ifiber],1), ps=10, xsty=3, ysty=3, $
                  xthick=2.0, ythick=2.0, charsize=1.5, charthick=2.0
                cc = get_kbrd(1)
             endif
          endfor 
          
; only retain TARGETS with well-measured redshifts; exclude fibers
; with no data (TOTAL(FERR)=0.0)
          keep = where((target_info.z gt 0.0) and (total(targetbigferr,1) gt 0.0),nkeep)
          splog, '   Retaining '+string(nkeep,format='(I0)')+'/'+string(ntargetfibers,format='(I0)')+$
            ' target spectra with good redshifts.'
          spzbest_target = spzbest_target[keep]
          target_info = target_info[keep]

; note - keep WAVE in double precision!          
          final_data = {wave: double(targetbigwave), flux: float(targetbigflux[*,keep]), $
            ferr: float(targetbigferr[*,keep]), skyflux: float(skybigflux), $
            skyferr: float(skybigferr)}

; write out
          outfile = 'spectra_'+passfield[jpass]+'.fits'
          im_mwrfits, final_data, outpath+outfile, clobber=clobber
             
          outfile = 'target_info_'+passfield[jpass]+'.fits'
          im_mwrfits, target_info, outpath+outfile, clobber=clobber

          splog, 'Time for this pass = '+string((systime(1)-t1)/60.0,$
            format='(G0.0)')+' minutes'
          print

       endfor ; close PASS loop

    endfor ; close PATH loop
    print
    splog, 'Time to unpack spectra = '+string((systime(1)-t0)/60.0,$
      format='(G0.0)')+' minutes'

    infofiles = file_search(outpath+'target_info_???.fits.gz',count=npass)
    for ipass=0L, npass-1L do if (n_elements(allinfo) eq 0L) then $
      allinfo = mrdfits(infofiles[ipass],1,/silent) else $
      allinfo = struct_append(allinfo,mrdfits(infofiles[ipass],1,/silent))
    help, allinfo

    splog, /close
    
return
end    

