pro deep2_specfit, zcat, thismask=thismask, firstmask=firstmask, $
  lastmask=lastmask, test=test, doplot=doplot, index=index1
; jm07sep26nyu - based on DEEP2_SPECFIT
; jm08sep04nyu - major rewrite based on the latest AGES_SPECFIT; now
;   uses ISPECLINEFIT_UNFLUXED() to do the continuum + line-fitting 

; ssh -X prism
; echo "deep2_specfit" | idl > & deep2_specfit.log.01 &

    fwhm2sig = 2D0*sqrt(2.0*alog(2.0))
    if keyword_set(doplot) or keyword_set(test) then $
      silent = 0 else silent = 1
    
    if (n_elements(zcat) eq 0L) then zcat = read_deep2_zcat(/good)

    version = deep2_version(/ispec)
    spec1dpath = deep2_path(/dr3)
    base_specfitpath = deep2_path(/specfit,/dr3)
    specfitpath = deep2_path(/specfit)+version+'/'

    linefile = base_specfitpath+'elinelist_'+version+'.dat'

    if (n_elements(thismask) eq 0L) then thismask = fix(zcat[uniq(zcat.maskname,$
      sort(zcat.maskname))].maskname)
    nmask = n_elements(thismask)

    if (n_elements(firstmask) eq 0L) then firstmask = 0L
    if (n_elements(lastmask) eq 0L) then lastmask = nmask-1L
    
; fit each mask separately

    t0 = systime(1)
    for imask = firstmask, lastmask do begin

; read the spectra

       these = where(thismask[imask] eq fix(zcat.maskname),nthese)
       if (nthese eq 0L) then begin
          splog, 'No spectra for mask '+strtrim(thismask[imask],2)
          continue
       endif
       zcat_mask = zcat[these]
       
       if (n_elements(index1) ne 0L) then zcat_mask = zcat[these[index1]]
       nobj = n_elements(zcat_mask)

       splog, 'Reading spectra from mask '+strtrim(thismask[imask],2)

       npix = 2*4096L
       wave = dblarr(npix,nobj)
       flux = wave*0.0D
       invvar = wave*0.0D

; possibly interpolate where the spectra do not overlap and set the
; inverse variance there equal to zero; also possibly normalize the
; spectra to be the same where they (nearly) overlap 
       for iobj = 0L, nobj-1L do begin

          spec1 = mrdfits(strtrim(spec1dpath+zcat_mask[iobj].file,2),1,/silent)
          spec2 = mrdfits(strtrim(spec1dpath+zcat_mask[iobj].file,2),2,/silent)
          
          wave[*,iobj] = [spec1.lambda,spec2.lambda]
          flux[*,iobj] = [spec1.spec,spec2.spec]
          invvar[*,iobj] = [spec1.ivar,spec2.ivar]

          srt = sort(wave[*,iobj])
          wave[*,iobj] = wave[srt,iobj]
          flux[*,iobj] = flux[srt,iobj]
          invvar[*,iobj] = invvar[srt,iobj]
          
;         djs_plot, [0], [0], xr=minmax(wave[*,iobj]), yr=minmax(flux[*,iobj])
;         djs_oplot, spec2.lambda, spec2.spec, color='red', ps=10
;         djs_oplot, spec1.lambda, spec1.spec, color='blue', ps=10
;         cc = get_kbrd(1)

       endfor
       
       suffix = string(thismask[imask],format='(I0)')       
       if keyword_set(test) then suffix = suffix+'_test'
       
       specres = 0.56*fwhm2sig  ; [FWHM, from Willmer et al.]
       vlinemaxtweak = 500.0    ; [km/s]
       sigmax = replicate(5000.0,nobj)

       specdata = ispeclinefit_unfluxed(wave,flux,invvar,zobj=zcat_mask.z,$
         specres=specres,sigmax=sigmax,galaxy=zcat_mask.galaxy,outpath=specfitpath,$
         suffix=suffix,linefile=linefile,specfit=specfit,/nologfile,$
         vlinemaxtweak=vlinemaxtweak,/clobber,silent=silent,doplot=doplot,$
         specdatafile=specdatafile)

; prepend the ZCAT structure to and overwrite SPECDATAFILE

       mwrfits, struct_addtags($
         struct_trimtags(zcat_mask,except=['MINWAVE','MAXWAVE']),$
         struct_trimtags(specdata,except=['GALAXY'])), specdatafile, /create
       spawn, 'gzip -f '+specdatafile, /sh

    endfor 
    splog, 'Total '+string((systime(1)-t0)/3600.0,format='(G0.0)')+' hours.'

; merge all the masks and write out

stop

    deep2_merge_specfit, alldata, test=test, /write
    
stop
    
return
end    
