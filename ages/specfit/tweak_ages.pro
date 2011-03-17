function ages_smooth, fluxin
return, smooth(medsmooth(fluxin,51),151,/nan,/edge_truncate)
end

pro tweak_ages, derive_tweak=derive_tweak, apply_tweak=apply_tweak
; jm08jan10nyu - the flux-calibration of some of the AGES spectra at
;                the blue wavelengths is poor; however, we can use our
;                spectral synthesis fits to estimate a multiplicative
;                correction vector; this routine constructs a QA plot
;                comparing the median ratio of the data (observed
;                spectrum) to the best-fitting BC03 model, so that we
;                can determine the wavelength below which the
;                correction function should be applied; it also writes
;                out the correction function for each plate; note that
;                the discrepancy at red wavelengths is due to the red
;                leak, which is an *additive* effect, so it is
;                accounted for by subtracting the median-smooth
;                function in IFITSPEC

; jm08aug - major update!
    
; note that ideally we would derive this function as a function of
; galaxy color    

    version = ages_version(/ispec_specfit)
    spec1dpath = ages_path(/spec1d)+'fluxed/after_skysubpca/'
    spec1dpath_out = ages_path(/spec1d)+'fluxed/tweak/'
    specfitpath = ages_path(/specfit)
    notweakpath = specfitpath+'fluxed/notweak/'+version+'/'
    tweakpath = specfitpath+'fluxed/tweak/'+version+'/'

    psfile = specfitpath+'tweak_ages_'+version+'.ps'
    fitsfile = specfitpath+'tweak_ages_'+version+'.fits'

    badplates = [106,110,209,310,311] ; unfluxed plates -- ignore

    ages = read_ages(/ancillary)
    allpass = ages.pass
    pass = allpass[uniq(allpass,sort(allpass))]
;   pass = cmset_op(pass,'and',/not2,badplates)
    npass = n_elements(pass)

; ##################################################    
; derive the tweak
; ##################################################    

    if keyword_set(derive_tweak) then begin
       
       dfpsplot, psfile, /color, /landscape
;      im_window, 0, xratio=0.9, yratio=0.7

; initialize the output tweak structure       
       
       result = {pass: 0, unfluxed: 0, ngalaxy: 0, waveratio: fltarr(4650), meanratio: fltarr(4650), $
         sigratio: fltarr(4650), medratio: fltarr(4650), smoothmeanratio: fltarr(4650), $
         smoothsigratio: fltarr(4650), smoothmedratio: fltarr(4650), smoothsmoothmeanratio: fltarr(4650), $
         smoothsmoothsigratio: fltarr(4650), smoothsmoothmedratio: fltarr(4650)}
       result = replicate(result,npass)

       for ipass = 0L, npass-1L do begin

          bad = where(pass[ipass] eq badplates,nbad)
          if (nbad ne 0L) then result[ipass].unfluxed = 1 ; flag
          
          specfitfile = file_search(notweakpath+$
            '?????_'+string(pass[ipass],format='(I0)')+'_specfit.fits.gz')
          specdatafile = repstr(specfitfile,'_specfit','_specdata')

          splog, 'Reading '+specdatafile
          specdata = mrdfits(strtrim(specdatafile,2),1,/silent)
          splog, 'Reading '+specfitfile
          specfit = mrdfits(specfitfile,0,/silent)

; only use high S/N galaxies at 0.01<z<1
          
          keep = where((specdata.z_obj gt 0.01) and (specdata.z_obj lt 1.0) and $
            (specdata.continuum_snr gt 5.0) and (specdata.continuum_chi2 lt 5.0),ngalaxy)
          
          result[ipass].pass = pass[ipass]
          result[ipass].ngalaxy = ngalaxy

          npix = (size(specfit,/dim))[0]
          fluxratio = fltarr(npix,ngalaxy)+1.0
          smoothfluxratio = fltarr(npix,ngalaxy)
          waveratio = reform(specfit[*,0,0])*(1.0+specdata[0].z_abs) ; they are all the same

          keep = lindgen(ngalaxy)
          for igal = 0L, ngalaxy-1L do begin ; observed spectrum divided by BC03 model
             good = where((specfit[*,5,igal] gt 0.0) and (specfit[*,2,igal] gt 0.0),ngood)
             if (ngood eq 0L) then keep[igal] = -1L else begin ; shouldn't happen, but catch anyway
                fluxratio[good,igal] = specfit[good,1,igal]/specfit[good,2,igal]
                smoothfluxratio[good,igal] = ages_smooth(specfit[good,1,igal])/$
                  ages_smooth(specfit[good,2,igal])
             endelse
          endfor
          keep = keep[where(keep ne -1L)]
          fluxratio = fluxratio[*,keep]
          smoothfluxratio = smoothfluxratio[*,keep]
          
          meanratio = fltarr(npix)
          medratio = fltarr(npix)
          sigratio = fltarr(npix)

          smoothmeanratio = fltarr(npix)
          smoothmedratio = fltarr(npix)
          smoothsigratio = fltarr(npix)

          for ipix = 0L, npix-1L do begin
             djs_iterstat, fluxratio[ipix,*], median=med, $
               mean=mn, sigma=sig, sigrej=3.0
             meanratio[ipix] = mn
             medratio[ipix] = med
             sigratio[ipix] = sig

             djs_iterstat, smoothfluxratio[ipix,*], median=med, $
               mean=mn, sigma=sig, sigrej=3.0
             smoothmeanratio[ipix] = mn
             smoothmedratio[ipix] = med
             smoothsigratio[ipix] = sig
          endfor

;         meanratio = total(fluxratio,2)/float(ngalaxy)
;         medratio = djs_median(fluxratio,2)
;         smoothmeanratio = total(smoothfluxratio,2)/float(ngalaxy)
;         smoothmedratio = djs_median(smoothfluxratio,2)

          smoothsmoothmeanratio = ages_smooth(smoothmeanratio)
          smoothsmoothmedratio = ages_smooth(smoothmedratio)
          smoothsmoothsigratio = ages_smooth(smoothsigratio)

; store the results

          result[ipass].waveratio[0L:npix-1L] = waveratio
          result[ipass].meanratio[0L:npix-1L] = meanratio
          result[ipass].sigratio[0L:npix-1L] = sigratio
          result[ipass].medratio[0L:npix-1L] = medratio
          result[ipass].smoothmeanratio[0L:npix-1L] = smoothmeanratio
          result[ipass].smoothsigratio[0L:npix-1L] = smoothsigratio
          result[ipass].smoothmedratio[0L:npix-1L] = smoothmedratio
          result[ipass].smoothsmoothmeanratio[0L:npix-1L] = smoothsmoothmeanratio
          result[ipass].smoothsmoothmedratio[0L:npix-1L] = smoothsmoothmedratio
          result[ipass].smoothsmoothsigratio[0L:npix-1L] = smoothsmoothsigratio
          
; generate the QA plot       
          
          djs_plot, waveratio, meanratio, ps=10, xsty=3, ysty=1, charsize=2.0, $
            charthick=2.0, yrange=[0.6,1.6], title='Pass '+string(pass[ipass],format='(I0)'), $
            xtitle='Observed Wavelength (\AA)', ytitle='Data / Model', color='gray', $
            xthick=2.0, ythick=2.0
          djs_oplot, !x.crange, [1,1]
;         djs_oplot, waveratio, smoothmedratio, color='dark green'
          djs_oplot, waveratio, smoothsmoothmedratio, color='red', thick=2
          djs_oplot, waveratio, smoothsmoothmeanratio, color='blue', thick=3
          djs_oplot, waveratio, smoothsmoothmeanratio+smoothsmoothsigratio, color='blue', thick=3, line=2
          djs_oplot, waveratio, smoothsmoothmeanratio-smoothsmoothsigratio, color='blue', thick=3, line=2
          legend, 'Ngal = '+string(ngalaxy,format='(I0)'), /right, /top, box=0, $
            charsize=2.0, charthick=2.0, clear=keyword_set(ps)
          
;         if (not keyword_set(ps)) then cc = get_kbrd(1)

       endfor

       splog, 'Writing '+fitsfile
       mwrfits, result, fitsfile, /create
       spawn, 'gzip -f '+fitsfile, /sh

       splog, 'Writing '+psfile
       dfpsclose
       spawn, 'gzip -f '+psfile, /sh

    endif

; ##################################################    
; apply the tweak
; ##################################################    

; the following spectra have elements of crappiness and should be
; truncated and/or removed:
;   312/081
;   
;   
;   
;       
    
    if keyword_set(apply_tweak) then begin

       splog, 'Reading '+fitsfile+'.gz'
       tweak = mrdfits(fitsfile+'.gz',1,/silent)

;      for ipass = 12, n_elements(tweak)-1L do begin
       for ipass = 0L, n_elements(tweak)-1L do begin

          thispass = spec1dpath+'ages_'+string(pass[ipass],$
            format='(I0)')+'.fits.gz'
          splog, 'Reading '+thispass
          specpass = mrdfits(thispass,1,/silent)

          specfitfile = file_search(notweakpath+$
            '?????_'+string(pass[ipass],format='(I0)')+'_specfit.fits.gz')
          specdatafile = repstr(specfitfile,'_specfit','_specdata')
          splog, 'Reading '+specdatafile
          specdata = mrdfits(strtrim(specdatafile,2),1,/silent)
          splog, 'Reading '+specfitfile
          specfit = mrdfits(specfitfile,0,/silent)

          npix = (size(specfit,/dim))[0]

; flux tweak

          fluxcor = specpass.wave*0.0+1.0
          fixflux = 1
          case pass[ipass] of
             105: maxwave = 4600.0
             111: maxwave = 4600.0
             112: maxwave = 4500.0
             113: maxwave = 4900.0
             115: maxwave = 4700.0
             201: maxwave = 4800.0
             203: maxwave = 5300.0
             205: maxwave = 4800.0
             206: maxwave = 4300.0
             207: maxwave = 4300.0
             208: maxwave = 5400.0
             210: maxwave = 4800.0
             213: maxwave = 5100.0
             301: maxwave = 5200.0
             303: maxwave = 4600.0
             304: maxwave = 5200.0
             307: maxwave = 5200.0
             308: maxwave = 5200.0
             309: maxwave = 5200.0
             313: maxwave = 4200.0
             314: maxwave = 4500.0
             422: maxwave = 4400.0
             else: fixflux = 0
          endcase             

          if fixflux then begin
             these = where((specpass.wave le maxwave),comp=notthese)
             fluxcor[these] = (tweak[ipass].smoothsmoothmeanratio[0:npix-1])[these]
;            djs_plot, specpass.wave, fluxcor, xsty=3, ysty=3
          endif

; now loop through each spectrum, subtract the red leak, and apply the
; flux tweak; don't do anything to the unfluxed plates

          if (tweak[ipass].unfluxed eq 0) then begin
             ngalaxy = n_elements(specpass.z)
             for igal = 0L, ngalaxy-1L do begin
;               djs_plot, specpass.wave, specpass.flux[*,igal], xsty=3, ysty=3
                these = where(specpass.wave ge 8500.0,nthese)
                specpass.flux[these,igal] = specpass.flux[these,igal] - $
                  reform(specfit[these,4,igal])/(1.0+specdata[igal].z_abs)
                specpass.flux[*,igal] = specpass.flux[*,igal]/fluxcor
                specpass.ferr[*,igal] = specpass.ferr[*,igal]/fluxcor
;               djs_oplot, specpass.wave, specpass.flux[*,igal], color='red'
;               cc = get_kbrd(1)
             endfor
          endif
             
; write out

          outfile = spec1dpath_out+'ages_'+string(pass[ipass],$
            format='(I0)')+'.fits'
          splog, 'Writing '+outfile
          mwrfits, specpass, outfile, /create
          spawn, 'gzip -f '+outfile

       endfor  

    endif
    
stop
    
return
end
    
