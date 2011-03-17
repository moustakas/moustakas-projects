pro fit_mzlzlocal, clobber=clobber
; jm09mar25nyu - fit and write out the local LZ and MZ relations
;   (no plotting)
; jm10aug16ucsd - major rewrite

; read the data    
    sdssancillary = read_mz_sample(/mzhii_ancillary,/sdss)
    sdssmass = read_mz_sample(/mzhii_mass,/sdss)
    sdssohdust = read_mz_sample(/sdss,/mzhii_log12oh)
    sdssohnodust = read_mz_sample(/sdss,/mzhii_log12oh,/nodust)

    agesancillary = read_mz_sample(/mzhii_ancillary)
    agesmass = read_mz_sample(/mzhii_mass)
    agesohdust = read_mz_sample(/mzhii_log12oh)
    agesohnodust = read_mz_sample(/mzhii_log12oh,/nodust)
;   mockages_noevol = read_mockages_sample(/zbin1)
;   mockages_levol = read_mockages_sample(/evolve,/zbin1)

    mzpath = ages_path(/projects)+'mz/'

; --------------------------------------------------    
; fit the SDSS MZ relation using three different models, all three
; abundance calibrations, and using both EWs and reddening-corrected
; fluxes
    fit_minmass = 9.0
    minmass = 8.7
    maxmass = 11.0
    binsize = 0.1
;   splog, 'testing mingal!'
    mingal = 100 ; 1000 ; 
    
    for jj = 0, 1 do begin
       case jj of
          0: begin
             suffix = '' ; EWs, the default
             flux = 0
             sdssoh = sdssohdust
          end
          1: begin
             suffix = '_fluxcor' ; reddening-corrected fluxes
             flux = 1
             sdssoh = sdssohnodust
          end
       endcase 
       
       closed_file = mzpath+'mzlocal_sdss'+suffix+'_closedbox.fits'
       poly_file = mzpath+'mzlocal_sdss'+suffix+'_poly.fits'
       doublepl_file = mzpath+'mzlocal_sdss'+suffix+'_doublepl.fits'
       brokenpl_file = mzpath+'mzlocal_sdss'+suffix+'_brokenpl.fits'

       for ii = 0, 2 do begin ; loop on each calibration
          t04 = 0 & m91 = 0 & kk04 = 0
          case ii of
             0: t04 = 1
             1: m91 = 1
             2: kk04 = 1
          endcase
          if keyword_set(t04) then calib = 't04'
          if keyword_set(m91) then calib = 'm91'
          if keyword_set(kk04) then calib = 'kk04'
          info = mzlz_grab_info(sdssoh,sdssancillary,sdssmass,$
            flux=flux,t04=t04,m91=m91,kk04=kk04)
          
          closed1 = fit_mz_closedbox(info.mass,info.oh,info.weight,$
            oh_err=info.oh_err,binsize=binsize,minmass=minmass,$
            maxmass=maxmass,mingal=mingal,fit_minmass=fit_minmass,$
            verbose=verbose)

;; test fit_minmass          
;          mm = range(8.0,12.0,30)
;          fit_minmass = [8.8,9.2,9.6,10.0]
;          djs_plot, info.mass, info.oh, psym=3, ysty=3, yr=[8.4,9.3]
;;         for bb = 0, 0 do begin
;          for bb = 0, n_elements(fit_minmass)-1 do begin
;             closed1 = fit_mz_closedbox(info.mass,info.oh,info.weight,$
;               oh_err=info.oh_err,binsize=binsize,minmass=minmass,$
;               maxmass=maxmass,mingal=mingal,fit_minmass=fit_minmass[bb],$
;               verbose=0)
;             if (bb eq 0) then djs_oplot, closed1.bin_mass, $
;               closed1.bin_oh, psym=6, color='red', sym=2
;             djs_oplot, mm, mz_closedbox(mm,closed1.coeff), color='blue'
;             splog, closed1.coeff, closed1.chi2, closed1.scatter
;          endfor
              
; fit a simple polynomial
          polyfit1 = fit_mz_poly(info.mass,info.oh,info.weight,$
            oh_err=info.oh_err,binsize=binsize,minmass=minmass,$
            maxmass=maxmass,mingal=mingal,fit_minmass=fit_minmass,$
            verbose=verbose)
; fit a smooth double-powerlaw
          doublepl1 = fit_mz_doublepl(info.mass,info.oh,info.weight,$
            oh_err=info.oh_err,binsize=binsize,minmass=minmass,$
            maxmass=maxmass,mingal=mingal,fit_minmass=fit_minmass,$
            verbose=0)
; fit a broken powerlaw
          brokenpl1 = fit_mz_brokenpl(info.mass,info.oh,info.weight,$
            oh_err=info.oh_err,binsize=binsize,minmass=minmass,$
            maxmass=maxmass,mingal=mingal,fit_minmass=fit_minmass,$
            verbose=0)

          closed1 = struct_addtags({calib: calib},closed1)
          polyfit1 = struct_addtags({calib: calib},polyfit1)
          doublepl1 = struct_addtags({calib: calib},doublepl1)
          brokenpl1 = struct_addtags({calib: calib},brokenpl1)
          if (ii eq 0) then begin
             closed = closed1
             polyfit = polyfit1
             doublepl = doublepl1
             brokenpl = brokenpl1
          endif else begin
             closed = [closed,closed1]
             polyfit = [polyfit,polyfit1]
             doublepl = [doublepl,doublepl1]
             brokenpl = [brokenpl,brokenpl1]
          endelse
       endfor          
       im_mwrfits, closed, closed_file, /clobber
       im_mwrfits, polyfit, poly_file, /clobber
       im_mwrfits, doublepl, doublepl_file, /clobber
       im_mwrfits, brokenpl, brokenpl_file, /clobber
    endfor

stop    
    
; --------------------------------------------------    
; fit the SDSS and low-redshift AGES LZ relations using all three
; abundance calibrations  
    faintmag = -18.0
    brightmag = -23.0
    binsize = 0.1
    mingal = 100

    for kk = 0, 1 do begin ; loop on the sample
       case kk of
          0: ssuffix = 'sdss'
          1: ssuffix = 'ages'
       endcase 

       for jj = 0, 1 do begin   ; loop on EWs and fluxes
          case jj of
             0: begin
                suffix = ''     ; EWs, the default
                flux = 0
                if (kk eq 0) then stroh = sdssohdust else stroh = agesohdust
             end
             1: begin
                suffix = '_fluxcor' ; reddening-corrected fluxes
                flux = 1
                if (kk eq 0) then stroh = sdssohnodust else stroh = agesohnodust
             end 
          endcase 
          lzfile = mzpath+'lzlocal_'+ssuffix+suffix+'.fits'

          for ii = 0, 2 do begin ; loop on each calibration
             t04 = 0 & m91 = 0 & kk04 = 0
             case ii of
                0: t04 = 1
                1: m91 = 1
                2: kk04 = 1
             endcase
             if keyword_set(t04) then calib = 't04'
             if keyword_set(m91) then calib = 'm91'
             if keyword_set(kk04) then calib = 'kk04'
             case kk of
                0: info = mzlz_grab_info(stroh,sdssancillary,sdssmass,$
                  flux=flux,t04=t04,m91=m91,kk04=kk04)
                1: info = mzlz_grab_info(stroh,agesancillary,agesmass,$
                  flux=flux,t04=t04,m91=m91,kk04=kk04,zmin=0.05,zmax=0.15)
             endcase
             lzfit1 = fit_lz(info.mb_ab,info.oh,info.weight,binsize=binsize, $
               brightmag=brightmag,faintmag=faintmag,mingal=mingal,verbose=verbose)
             lzfit1 = struct_addtags({calib: calib},lzfit1)
             if (ii eq 0) then lzfit = lzfit1 else lzfit = [lzfit,lzfit1]
          endfor          
          im_mwrfits, lzfit, lzfile, /clobber
       endfor
    endfor

return
end
    
    
