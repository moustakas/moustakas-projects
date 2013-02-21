pro nsa_ir_wise
; jm13feb11siena - compute the IR properties of the NSA sample using
; the WISE photometry

    common nsa_wise, nsa, wise
    if n_elements(nsa) eq 0L then nsa = read_nsa(wise=wise)
    ngal = n_elements(nsa)

    filters = [galex_filterlist(),sdss_filterlist(),wise_filterlist()]
    nfilt = n_elements(filters)
    
    path = getenv('IM_DATA_DIR')+'/nsa/'
    outfile = path+'nsa_ir_wise.fits'

; nuLnu is the rest-frame luminosity 
    out = struct_addtags(struct_trimtags(nsa,select=['iauname','ra','dec','zdist']),$
      replicate({maggies: fltarr(nfilt), ivarmaggies: fltarr(nfilt), $
      lir_chary: -999.0, lir_dale: -999.0, $
;     lir_rieke: -999.0, $
      nuLnu_chary: dblarr(3)-999.0, nuLnu_dale: dblarr(3)-999.0, $
;     nuLnu_rieke: dblarr(3)-999.0, $
      ivar_nulnu_chary: dblarr(3)-999.0, ivar_nulnu_dale: dblarr(3)-999.0},ngal));, $
;     ivar_nulnu_rieke: dblarr(3)-999.0},ngal))
;     sfr_chary: -999.0, sfr_dale: -999.0, sfr_rieke: -999.0, $
;     uvsfr_chary: -999.0, uvsfr_dale: -999.0, uvsfr_rieke: -999.0, $
;     indx_chary: -999.0, indx_dale: -999.0, indx_rieke: -999.0, $
;     irx_chary: -999.0, irx_dale: -999.0, irx_rieke: -999.0, $
;     a1500_chary: -999.0, a1500_dale: -999.0, a1500_rieke: -999.0},ngal))

; add the photometry
    wise_to_maggies, wise, wmm, wiv
    out.maggies = [nsa.nmgy*1D-9,wmm]
    out.ivarmaggies = [nsa.nmgy_ivar/1D-9^2,wiv]
    
; throw out upper limits to make sure the K-corrections in
; im_simple_kcorrect get computed from the band where the object was
; detected
    these = where(filters eq 'wise_w3.par' or filters eq 'wise_w4.par')
    good = where(total(out.maggies[these] gt 0,1) gt 0)
    wmaggies = out[good].maggies[these]
    wivarmaggies = out[good].ivarmaggies[these]*(wmaggies gt 0)

    out[good].lir_chary = alog10(im_wise2lir(out[good].zdist,wmaggies,wivarmaggies,$
      /chary,model_indx=indx_chary,nulnu=nulnu_chary,$
      ivar_nulnu=ivar_nulnu_chary,modelwave=modelwave_chary,modelmab=modelmab_chary))
    out[good].lir_dale = alog10(im_wise2lir(out[good].zdist,wmaggies,wivarmaggies,$
      /dale,model_indx=indx_dale,nulnu=nulnu_dale,$
      ivar_nulnu=ivar_nulnu_dale,modelwave=modelwave_dale,modelmab=modelmab_dale))
;   out[good].lir_rieke = alog10(im_wise2lir(out[good].zdist,wmaggies,wivarmaggies,$
;     /rieke,model_indx=indx_rieke,nulnu=nulnu_rieke,$
;     ivar_nulnu=ivar_nulnu_rieke,modelwave=modelwave_rieke,modelmab=modelmab_rieke))

    out[good].nulnu_chary = nulnu_chary
    out[good].nulnu_dale = nulnu_dale
;   out[good].nulnu_rieke = nulnu_rieke

    out[good].ivar_nulnu_chary = ivar_nulnu_chary
    out[good].ivar_nulnu_dale = ivar_nulnu_dale
;   out[good].ivar_nulnu_rieke = ivar_nulnu_rieke

; compute upper limits using just the 22-micron flux (ignore the
; 12-micron flux, if any)
    crap = where((out.maggies[these[1]] eq 0) and (out.ivarmaggies[these[1]] gt 0))
;   crap = where((total(out.maggies[these] gt 0,1) eq 0) and (out.ivarmaggies[these[1]] gt 0))
    wmaggies = 1.0/sqrt(out[crap].ivarmaggies[these])
    wivarmaggies = 1.0/(0.05*wmaggies)^2
    wivarmaggies[0,*] = 0

    lir_chary = alog10(im_wise2lir(out[crap].z,wmaggies,wivarmaggies,$
      /chary,model_indx=indx_chary,nulnu=nulnu_chary,$
      ivar_nulnu=ivar_nulnu_chary,modelwave=modelwave_chary,modelmab=modelmab_chary))
    lir_dale = alog10(im_wise2lir(out[crap].z,wmaggies,wivarmaggies,$
      /dale,model_indx=indx_dale,nulnu=nulnu_dale,$
      ivar_nulnu=ivar_nulnu_dale,modelwave=modelwave_dale,modelmab=modelmab_dale))
;    lir_rieke = alog10(im_wise2lir(out[crap].z,wmaggies,wivarmaggies,$
;      /rieke,model_indx=indx_rieke,nulnu=nulnu_rieke,$
;      ivar_nulnu=ivar_nulnu_rieke,modelwave=modelwave_rieke,modelmab=modelmab_rieke))

    out[crap].lir_chary = -lir_chary
    out[crap].lir_dale = -lir_dale
;   out[crap].lir_rieke = -lir_rieke

    out[crap].nulnu_chary = -nulnu_chary
    out[crap].nulnu_dale = -nulnu_dale
;   out[crap].nulnu_rieke = -nulnu_rieke

;   niceprint, cat.galaxy, out.lir_chary, out.nulnu_chary, cat.maggies[these[1]], cat.ivarmaggies[these[1]]
    
; write out    
    im_mwrfits, out, outfile, clobber=clobber

    
    
stop
return
end
