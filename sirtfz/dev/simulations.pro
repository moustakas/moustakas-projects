pro simulations, snr=snr, mindetect=mindetect, nzsub=nzsub, skyfrac=skyfrac, image=image
;+
; NAME:
;	SIMULATIONS	
;
; PURPOSE:
;	Generate realistic photometric catalogs for observations in
;	user-selected bandpasses to user-selected limiting depths.
;
; CALLING SEQUENCE:
;	simulations, [snr=, mindetect=, nzsub=, image=, skyfrac=]
;
; OPTIONAL INPUTS:
;	snr       - signal-to-noise ratio for a galaxy to be detected
;	mindetect - the galaxy must be detected in at least mindetect
;                   filters to be written to the catalog
;	skyfrac   - fraction of the sky to observe in arcminutes
;                   (e.g. skyfrac=32 means 32 sq arcmin)
;	image     - simulate an image in a MIPS band (e.g. "MIPS 70") 
;	
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	IDL save sets.  Also number counts and redshift distribution
;	plots. 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; PROCEDURES USED:
;	FILTER_MATCH(), 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 May 12, U of A, written
;-

; useful notes:  1 deg^2 = 3.046175D-4 sr; 1 arcmin^2 = 2.77777D-4 deg^2

    common sirtf_simulations
    common cosmology

    if strupcase(sirtf.templates) ne 'DEVRIENDT' then message, 'These models are not usable!'
    if not keyword_set(skyfrac) then skyfrac = 60D    ; sky coverage in arcmin
    if not keyword_set(snr) then snr = 5.0    ; signal-to-noise threshold for detection
    if not keyword_set(nzsub) then nzsub = 5L ; number of redshift subsections

    alpha = sirtf.evolution.alpha & beta = sirtf.evolution.beta

; match the currently-selected filters with the available filters

    obands = filter_match(*sirtf.filters)
    nbands = n_elements(obands)
    filters = sirtf.bandcube[obands].bandnames ; filter names

    path = (sirtf_datapath())[0]+'catalogs/SIMULATIONS/'
    idlname = path+'simulations_'+strn(skyfrac,format='(G0.0)')+'_'+$
      strn(alpha,format='(G0.0)')+strn(beta,format='(G0.0)')+'.idlsave'

    if keyword_set(image) then begin
       imindx = where(obands eq (filter_match(image))[0],bmatch) & imindx = imindx[0]
       if bmatch eq 0L then message, 'Filter mismatch!'
       imfname = path+'imsim_'+sirtf.bandcube[obands[imindx]].mininames+'_'+$
         strn(skyfrac,format='(G0.0)')+'_'+strn(alpha,format='(G0.0)')+strn(beta,format='(G0.0)')+'.idlsave'
    endif
       
    if not keyword_set(mindetect) then mindetect = (5B<byte(nbands))

    dlum = alog(10D)*sirtf.lf.lum*sirtf.lf.logbinsz/2.5D ; luminosity function width

    zarray = *sirtf.redshift.zarray ; redshift array and resolution
    dz = sirtf.redshift.dz
    nztot = n_elements(zarray)
    nlf = sirtf.lf.nlf
    
    flim = 10D^(-findgen(401)/50.)*1000.0 ; limiting flux array (mJy)
    nflim = n_elements(flim)
    dnds = fltarr(nflim,nbands) ; number counts
    dndz = fltarr(nztot,nbands) ; redshift distribution
    
    nsources = 0L
    for i = 0L, nztot/nzsub-1L do begin

       if i gt 0L then append = 1L

       zsub = zarray[i*nzsub:(i+1L)*nzsub-1L]  ; redshift subarray
       nz = n_elements(zsub)

       dvdz = dvcomoving(zsub)*3.04617D-4*2.77777D-4*skyfrac^2/1D18 ; comoving volume (Mpc^3/skyfrac^2)

       savenames = ['gflux_'+strn(i+1),'ferror_'+strn(i+1),'catinfo_'+strn(i+1)]
       
       gflux = ptr_new()              ; catalog fluxes
       ferror = fltarr(nbands,nlf*nz) ; flux error
       catinfo = fltarr(3,nlf*nz)     ; catalog info (ngalaxies,type,redshift)
    
       l_evol = (1.0+zsub)^alpha ; luminosity evolution
       d_evol = (1.0+zsub)^beta  ; density evolution

       if keyword_set(image) then begin

          imflux = ptr_new()
          imfluxname = ['imflux_'+strn(i+1)]

       endif
    
       for j = 0L, nz-1L do begin ; loop on redshift
       
; gaussian perturb the number of galaxies in each LF bin
       
          ngals_true = reform(d_evol[j] * sirtf.lf.phi * dlum * dvdz[j] * dz)
          ngals = round(ngals_true + sqrt(ngals_true)*randomn(seed,nlf)) > 0L

; observe a random SED in each LF bin
       
          observe_each_sed, zsub[j], l_evol[j], sirtf.sedcube, sirtf.lf, sirtf.bandcube[obands], flux, type

          for k = 0L, nlf-1L do begin ; loop on the luminosity function
       
             if ngals[k] gt 0L then begin
       
                dflux = fltarr(ngals[k],nbands)
                detect = bytarr(ngals[k],nbands)

                for x = 0L, nbands-1L do begin

                   ferror[x,j*nlf+k] = sqrt(sirtf.bandcube[obands[x]].flimit^2 + $ ; photon/confusion noise
                                        (0.05*flux[k,x])^2)                        ; flux error floor (5%)
          
                   dflux[*,x] = flux[k,x] + randomn(seed,ngals[k])*ferror[x,j*nlf+k] ; Gaussian-perturbed flux
                   detect[*,x] = dflux[*,x] gt snr*sirtf.bandcube[obands[x]].flimit  ; 1B means snr-sigma detection

                   bad = where(detect[*,x] eq 0B,nbad) ; set undetected fluxes to zero
                   if nbad ne 0L then dflux[bad,x] = 0.0

                   for y = 0L, nflim-1L do if (flux[k,x] gt flim[y]) then $ ; number counts (no noise)
                     dnds[y,x] = dnds[y,x] + ngals_true[k]

                   if (flux[k,x] gt snr*sirtf.bandcube[obands[x]].flimit) then $ ; z-distribution (no noise)
                     dndz[j+i*nz,x] = dndz[j+i*nz,x] + ngals_true[k]
                   
                endfor 

                good = where(byte(total(detect,2)) ge mindetect,ngood) ; detected in at least mindetect bands
                if ngood ne 0L then begin

                   goodflux = transpose(dflux[good,*]) ; catalog galaxy fluxes
                   if ptr_valid(gflux) then $
                     gflux = ptr_new([[*gflux],[goodflux]]) else $ 
                     gflux = ptr_new(goodflux)

                endif

                nsources = nsources + ngood
                catinfo[*,j*nlf+k] = [ngood,type[k],zsub[j]] ; catalog information

                if keyword_set(image) then begin

                   tempflux = flux[k,imindx]+(randomu(seed,ngals[k])-0.5)*flux[k,imindx]
                   if ptr_valid(imflux) then $
                     imflux = ptr_new([*imflux,tempflux]) else imflux = ptr_new(tempflux)
                   delvarx, tempflux
                   
                endif
                
             endif else catinfo[*,j*nlf+k] = [0,type[k],zsub[j]]

          endfor ; LF loop

       endfor    ; redshift loop

       if ptr_valid(gflux) then begin
          print, 'Writing ', strupcase(savenames)
          cmsave, filename=idlname, *gflux, ferror, catinfo, names=savenames, append=append
          ptr_free, gflux
       endif
          
       if keyword_set(image) then cmsave, filename=imfname, *imflux, names=imfluxname, append=append
       if ptr_valid(imflux) then ptr_free, imflux
       
    endfor       ; redshift subarray loop

    cmsave, filename=idlname, nsources, filters, /append

    filtids = sirtf.bandcube[obands[0]].mininames
    mininames = sirtf.bandcube[obands].mininames
    for i = 1L, nbands-1L do filtids = filtids+'_'+mininames[i]
    flimits = sirtf.bandcube[obands].flimit
    
; write the number counts data to disk

    ncountsfile = path+'NCOUNTS/counts_'+strn(skyfrac,format='(G0.0)')+'_'+$
      strn(alpha,format='(G0.0)')+strn(beta,format='(G0.0)')+'_'+filtids+'.idlsave'
    cmsave, filename=ncountsfile, dnds, flim, filters, mininames
    
; write the redshift distribution data to disk

    zdistfile = path+'ZDIST/zdist_'+strn(skyfrac,format='(G0.0)')+'_'+$
      strn(alpha,format='(G0.0)')+strn(beta,format='(G0.0)')+'_'+filtids+'.idlsave'
    cmsave, filename=zdistfile, dndz, zarray, filters, mininames, flimits

return
end
