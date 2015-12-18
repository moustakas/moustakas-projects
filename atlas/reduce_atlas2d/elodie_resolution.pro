;+
; NAME:
;       ELODIE_RESOLUTION
;
; PURPOSE:
;       Derive the spectral resolution of the integrated
;       spectrophotometric atlas by comparing drift-scanned
;       observations of stars that appear in the ELODIE spectral
;       database.  
;
; CALLING SEQUENCE:
;       elodie_resolution, /doplot
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       doplot - generate a plot
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS: 
;       The ELODIE data are of sufficiently high resolution that we
;       can assume that the instrumental resolution is zero.
;
; PROCEDURES USED:
;       IFORAGE(), RD1DSPEC(), GET_ELEMENT, COMBINE1FIBER,
;       GCONVOLVE(), ELODIE_TABLE(), IM_READ_ELODIE()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 May 13 U of A
;-

pro elodie_resolution, doplot=doplot, postscript=postscript

    datapath = atlas_path(/templates1d)
    pushd, datapath

; read the 1D stellar spectra

;   A5 G8III K0III G8III G0 G2III F0V A0 B9V F3V F9V B6V B6V K2 A2 A3 A0 A0 A0 A0 A0 K0 G5III K3III K1III B5

    speclist = findfile('*.ms.fits')
    forage = iforage(speclist)

; read the ELODIE database and grab the appropriate stars

    elodie = elodie_table()
    
    doit = match_string(strcompress(forage.galaxy,/remove),elodie.star,findex=match,/exact)
    good = where(match ne -1L,nstar)
    elodie = elodie[match[good]]

    keep = [1,5,9,10,13,15,21,23,24]
;   keep = [1,2,5,9,10,13,21,22,23]
    nstar = n_elements(keep)
    elodie = elodie[keep]
    good = good[keep]
    
;   print_struct, elodie

    forage = forage[good]
    speclist = speclist[good]
    star = strcompress(forage.galaxy,/remove)
    niceprint, elodie.star, forage.galaxy, speclist, star

    splog, 'Matching '+strn(nstar)+' ELODIE stars.'

; initialize the Gaussian widths

    fwhm2sig = 2.0*sqrt(2.0*alog(2.0)) ; convert FWHM to sigma [= 2.35]

    fwhmmin = 5.0  ; [Angstrom]
    fwhmmax = 15.0 ; [Angstrom]
    dfwhm = 1.0    ; [Angstrom]
    
    sigma = (findgen((fwhmmax-fwhmmin)/dfwhm+1)*dfwhm+fwhmmin)/fwhm2sig ; sigma [Angstrom]
    nsigma = n_elements(sigma)

    elodie_cd1_1 = 0.2                 ; ELODIE dispersion [Angstrom/pixel]
    sigpix = float(sigma/elodie_cd1_1) ; Gaussian sigma width [pixels]
       
; initialize the wavelength bins

    wave1 = 4300.0
    wave2 = 6500.0
    dwave = 200.0

    nbins = long((wave2-wave1)/dwave)
    wavebins = findgen(nbins)*dwave+wave1+dwave/2.0

    chi2arr = fltarr(nbins,nsigma)
    resolution = fltarr(nbins,nstar)

;   if keyword_set(doplot) then window, 2, xs=300, ys=300
    for i = 0L, nstar-1L do begin

; read the ELODIE data and the ATLAS drift-scanned star

       id = string(elodie[i].elodie_id,format='(I5.5)')
       splog, 'Reading star '+elodie[i].star+', ID='+id+', Type='+strn(elodie[i].sp_type)
       elodieflux = im_read_elodie(id+'.fits.gz',wave=elodiewave)
       
       scube = rd1dspec(speclist[i],/silent)
       atlasflux = scube.spec
       atlasferr = scube.sigspec
       atlaswave = scube.wave
       icleanup, scube

; solve for the radial velocity shift and apply it using the full
; wavelength range.  maybe first rebin the data to have constant S/N
; per pixel?

       temp = gconvolve(elodieflux,9.0/fwhm2sig/elodie_cd1_1,/edge_truncate)
       zans = im_ztweak(temp,elodiewave,atlasflux,atlaswave,vmin=-100,$
         vmax=+100,/normalize,minwave=wave1,maxwave=wave2,doplot=0)

       if zans.pshift_err ne -1.0 then begin
          
          elodieflux = interpol(elodieflux,elodiewave,elodiewave*(1+zans.zshift))
          elodiewave = elodiewave*(1+zans.zshift)

       endif

; loop on each wavelength bin       
       
       for k = 0L, nbins-1L do begin

          mnwave = wavebins[k]-dwave/2.0
          mxwave = wavebins[k]+dwave/2.0
          splog, 'Wavelength range=['+strn(mnwave)+','+strn(mxwave)+']'
          
          get_element, elodiewave, [mnwave,mxwave]*[-1,1], xx
          efl = im_normalize(elodieflux[xx[0]:xx[1]],/max)
          eflwave = elodiewave[xx[0]:xx[1]]

          get_element, atlaswave, [mnwave,mxwave]*[-1,1], xx
          flux = im_normalize(atlasflux[xx[0]:xx[1]],/max,const=const)
          ferr = atlasferr[xx[0]:xx[1]]/const
          wave = atlaswave[xx[0]:xx[1]]

; loop on each Gaussian sigma width.  the first vector may have a
; small radial velocity offset, which we remove after the first
; iteration
          
          efl_broad = efl*0.0
;         efl_in = efl
          for j = 0L, nsigma-1L do begin

             print, format='("Sigma width (FWHM)=",F5.2," Angstrom.",A1,$)', sigma[j]*fwhm2sig, string(13b)

; broaden the ELODIE spectrum
             
             efl_broad = gconvolve(efl,sigpix[j],/edge_truncate)
;            efl_broad = gconvolve(efl_in,sigpix[j],/edge_truncate)
             
; resample the ELODIE spectrum onto the wavelength spacing of our data 

             combine1fiber, alog10(eflwave), efl_broad, efl_broad*0.0+1.0, $
               newloglam=alog10(wave), newflux=efl_new
             good = where(efl_new gt 0.0,npix)

             efl_new = efl_new[good]
             wave = wave[good]
             flux = flux[good]             
             ferr = ferr[good]             
             
; solve for the smoothly varying change in continuum shape and divide
; this curve into our data

             binflux = im_binspec(flux,wave,binsize=(mxwave-mnwave)*0.05,binwave=binwave)
             binefl = im_binspec(efl_new,wave,binsize=(mxwave-mnwave)*0.05,binwave=binwave)

             smooth = binflux / binefl
             smooth_big = interpol(smooth,binwave,wave)
             binflux = binflux / smooth
             flux = flux / smooth_big

;             if keyword_set(doplot) and (j eq 0L) then begin
;                wset, 2
;                plot, wave, flux, ps=10, xsty=3, ysty=3, xthick=2.0, ythick=2.0, $
;                  charsize=1.2, charthick=2.0, thick=2.0
;;               oplot, binwave, binflux, ps=4, syms=2.0
;                djs_oplot, wave, efl_new, ps=10, color='red'
;;               oplot, binwave, binefl, ps=2, syms=2.0
;             endif

; compute the chi2 statistic
             
;            chi2arr[k,j] = total((flux-efl_new)^2)
             chi2arr[k,j] = total((flux-efl_new)^2/ferr^2.0)/npix

          endfor
          print
;         chi2arr[k,*] = chi2arr[k,*]/max(chi2arr[k,*])
          
          findchi2min, sigma*fwhm2sig, reform(chi2arr[k,*]), minchi2, best_sigma
          resolution[k,i] = best_sigma
          splog, 'Best-fitting width (FWHM)='+string(best_sigma,format='(F5.2)')+' Angstrom.'

          if keyword_set(doplot) then begin
             window, 0, xs=300, ys=300
             djs_plot, sigma*fwhm2sig, chi2arr[k,*], xsty=3, ysty=3, thick=2.0, line=0, $
               xthick=2.0, ythick=2.0, charsize=1.2, charthick=2.0, $
               xtitle='\sigma(FWHM) ['+angstrom()+']', ytitle='\chi^{2}'
;            cc = get_kbrd(1)
          endif
          
       endfor 
          
    endfor

; write out

    cmsave, filename='elodie_resolution.idlsave'

; remove any measurements that floored against the minimum or maximum
; allowable resolution    

stop
       

return
end    
