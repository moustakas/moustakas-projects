pro unpack_skysubpca, debug=debug, save=save
; jm05dec20uofa - write out some spectra to test Wild's PCA
;                 sky-subtraction code
    
    datapath = ages_path()+'skysubpca/'
    vredux = '0051'

    for ipath = 0L, n_elements(datapath)-1L do begin
       
       sphectfits = file_basename(file_search(datapath[ipath]+'spHect-???-'+vredux[ipath]+'.fits',count=npassfield))
       spzbestfits = file_basename(file_search(datapath[ipath]+'spZbest-???-'+vredux[ipath]+'.fits'))
       passfield = strmid(sphectfits,7,3) ; NOT GENERAL

       for jpass = 0L, npassfield-1L do begin

          splog, 'Unpacking PASSFIELD '+passfield[jpass]+'.'

; read the data          
          
          bigwave   = mrdfits(datapath[ipath]+sphectfits[jpass],0,/silent)
          bigflux   = mrdfits(datapath[ipath]+sphectfits[jpass],1,/silent)
          biginvvar = mrdfits(datapath[ipath]+sphectfits[jpass],2,/silent)
          plugmap   = mrdfits(datapath[ipath]+sphectfits[jpass],5,/silent)
          spzbest   = mrdfits(datapath[ipath]+spzbestfits[jpass],1,/silent)

; identify TARGET and SKY fibers

          skyfibers = where(strtrim(plugmap.objtype,2) eq 'SKY',nskyfibers)
          targetfibers = where((strtrim(plugmap.objtype,2) eq 'TARGET') and $
            (spzbest.zwarning eq 0L),ntargetfibers)

          zobj = spzbest[targetfibers].z
          class = spzbest[targetfibers].class
          subclass = spzbest[targetfibers].subclass
          
; loop on each TARGET fiber

          for ifiber = 0L, ntargetfibers-1L do begin

             wave = bigwave[*,targetfibers[ifiber]]
             flux = bigflux[*,targetfibers[ifiber]]
             invvar = biginvvar[*,targetfibers[ifiber]]
             npix = n_elements(wave)

             flux = djs_maskinterp(flux,invvar le 0,wave,/const)
             invvar = djs_maskinterp(invvar,invvar le 0,wave,/const)
             ferr = 1.0/sqrt(invvar)

; resample the wavelength array

             minwave = 3700.0 ; min(wave)
             maxwave = 9200.0 ; max(wave)
;            dwave = (max(wave)-min(wave)) / (npix-1.0) ; mean dispersion
             dwave = 1.2

             newwave = findgen((maxwave-minwave)/dwave+1)*dwave+minwave
             newflux = 1E-17*interpol(flux,wave,newwave)
             newferr = 1E-17*sqrt(interpol(ferr^2.0,wave,newwave))

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
                djs_plot, targetbigwave, smooth(targetbigflux[*,ifiber],1), ps=10, xsty=3, ysty=3, $
                  xthick=2.0, ythick=2.0, charsize=1.5, charthick=2.0
                cc = get_kbrd(1)
             endif
             
          endfor

; loop on each SKY fiber

          for ifiber = 0L, nskyfibers-1L do begin

             wave = bigwave[*,skyfibers[ifiber]]
             flux = bigflux[*,skyfibers[ifiber]]
             invvar = biginvvar[*,skyfibers[ifiber]]
             npix = n_elements(wave)

             flux = djs_maskinterp(flux,invvar le 0,wave,/const)
             invvar = djs_maskinterp(invvar,invvar le 0,wave,/const)
             ferr = 1.0/sqrt(invvar)

; resample the wavelength array

             minwave = 3700.0 ; min(wave)
             maxwave = 9200.0 ; max(wave)
;            dwave = (max(wave)-min(wave)) / (npix-1.0) ; mean dispersion
             dwave = 1.2

             newwave = findgen((maxwave-minwave)/dwave+1)*dwave+minwave
             newflux = 1E-17*interpol(flux,wave,newwave)
             newferr = 1E-17*sqrt(interpol(ferr^2.0,wave,newwave))

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

; write out

          if keyword_set(save) then begin
             
             cmsave, filename=datapath+'ages_skysubpca.idlsave', zobj, class, subclass, $
               skybigwave, skybigflux, skybigferr, targetbigwave, targetbigflux, $
               targetbigferr, /verbose

          endif

       endfor 

    endfor 
          
    stop
    
return
end    
