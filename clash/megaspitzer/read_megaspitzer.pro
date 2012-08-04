function read_megaspitzer
; jm12jul25siena - read the photometry of the dropouts

    catpath = clash_path(/megaspitzer)
    filt = megaspitzer_filterlist(short=short,weff=weff,zpt=zpt,nice=nice)
    nfilt = n_elements(filt)

; fluxes are in nJy and have *not* been corrected for neither Galactic
; extinction nor aperture losses
    cat = rsex(catpath+'alldropouts.cat')
    ngal = n_elements(cat)
    cat = struct_addtags(cat,replicate({ch1_flux: 0.0, ch1_fluxerr: 1000.0, $
      ch2_flux: 0.0, ch2_fluxerr: 1000.0},ngal))
    
; fluxes are in nJy and have been corrected for Galactic extinction
; and aperture losses    
;   cat = rsex(catpath+'dropouts.cat')
;   ngal = n_elements(cat)

; assign rough redshifts
    cat = struct_addtags(cat,replicate({z: 0.0, hmag: 0.0},ngal))
    cat.hmag = -2.5*alog10(cat.f160w_flux)+31.4
    
    drop = ['B','V','i','z','Y']
    zz = [3.8, 4.9, 5.9, 6.8, 7.8]
    for ii = 0, n_elements(drop)-1 do begin
       match = where(strmatch(cat.galaxy,'*'+drop[ii]+'-*',/fold),nmatch)
       if nmatch ne 0 then cat[match].z = zz[ii]
    endfor

; add WeiOne
    cat1 = cat[0]
    wei = read_santorini()
    santorini_to_maggies, wei, maggies, ivar, filterlist=filt

    fact = 10D^(+0.4D*31.4) ; conversion to nJy
    for ii = 0, n_elements(filt)-1 do begin
       wei.(2*ii+4) = maggies[ii]*fact
       wei.(2*ii+1+4) = fact/sqrt(ivar[ii])
    endfor
    
    cat1 = im_struct_assign(wei,cat1)
    cat1.galaxy = 'WeiOne'
    cat = [cat,cat1]
    
return, cat
end
