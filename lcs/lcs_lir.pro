function wise2lir, phot, chary=chary, dale=dale, rieke=rieke, $
  filters=filters, debug=debug
    lsun = im_lsun()
    ngal = n_elements(phot)
;   ir = replicate({lir: -999.0, nuLnu: dblarr(3)-999.0, $
;     ivar_nuLnu: dblarr(3)-999.0},ngal)
    ir = replicate({lir: -999.0, nuLnu24: -999.0, nuLnu24_err: -999.0},ngal)

    these = where(filters eq 'wise_w3.par' or filters eq 'wise_w4.par')
    good = where(total(phot.maggies[these] gt 0,1) gt 0)
    wmaggies = phot[good].maggies[these]
    wivarmaggies = phot[good].ivarmaggies[these]*(wmaggies gt 0)
    
    ir[good].lir = alog10(im_wise2lir(phot[good].zdist,wmaggies,wivarmaggies,$
      chary=chary,dale=dale,rieke=rieke,model_indx=indx,nulnu=nulnu,$
      ivar_nulnu=ivar_nulnu,debug=debug))
    ir[good].nuLnu24 = reform(alog10(nuLnu[2,*]/lsun))
    ir[good].nuLnu24_err = reform(1.0/sqrt(ivar_nuLnu[2,*])/nuLnu[2,*]/alog(10))

; compute upper limits using just the 22-micron flux (ignore the
; 12-micron flux, if any)
    crap = where((phot.maggies[these[1]] eq 0) and (phot.ivarmaggies[these[1]] gt 0))
    wmaggies = 1.0/sqrt(phot[crap].ivarmaggies[these])
    wivarmaggies = 1.0/(0.05*wmaggies)^2
    wivarmaggies[0,*] = 0

    ir[crap].lir = -alog10(im_wise2lir(phot[crap].zdist,wmaggies,wivarmaggies,$
      chary=chary,dale=dale,rieke=rieke,model_indx=indx,nulnu=nulnu))
    ir[crap].nuLnu24 = reform(-alog10(nuLnu[2,*]/lsun))

; rearrange the structure tag names
    if keyword_set(chary) then suffix = 'chary'
    if keyword_set(dale) then suffix = 'dale'
    if keyword_set(rieke) then suffix = 'rieke'
    tags = tag_names(ir)
    ir = im_struct_trimtags(ir,select=tags,newtags=tags+'_'+suffix)
return, ir
end

pro lcs_lir, clobber=clobber
; jm13feb26siena - read the output of LCS_MEMBERSHIP and compute IR 
; properties for the sample

    common lcs_wise, allwise
    if n_elements(allwise) eq 0L then allwise = mrdfits(getenv('IM_DATA_DIR')+$
      '/nsa/nsa_v0_1_2_wise.fits.gz',1)

    lcspath = getenv('LCS_DATA')
    filters = [galex_filterlist(),sdss_filterlist(),wise_filterlist()]
    nfilters = n_elements(filters)

    cl = rsex(getenv('LCS_DIR')+'/lcs_sample.cat')
    ncl = n_elements(cl)
    struct_print, cl

    for ii = 0, ncl-1 do begin
       nsafile = lcspath+strlowcase(cl[ii].cluster)+'_nsa.fits.gz'
       nsa = mrdfits(nsafile,1)
       wise = allwise[nsa.nsaid] ; pick the ones we care about
       ngal = n_elements(nsa)

; build the output structure, including the complete photometry
       wise_to_maggies, wise, wmm, wiv
       maggies = [nsa.nmgy*1D-9,wmm]
       ivarmaggies = [nsa.nmgy_ivar/1D-9^2,wiv]
       
       out = struct_addtags(struct_trimtags(nsa,select=['nsaid','ra','dec','zdist']),$
         replicate({maggies: fltarr(nfilters), ivarmaggies: fltarr(nfilters)},ngal))
       out.maggies = maggies
       out.ivarmaggies = ivarmaggies

stop       
       
       ir_chary = wise2lir(out,filters=filters,/chary,debug=0)
       ir_dale = wise2lir(out,filters=filters,/dale,debug=0)
       out = struct_addtags(out,ir_chary)
       out = struct_addtags(out,ir_dale)

       gd = where(out.lir_chary gt 0.0)
       lim = where(out.lir_chary lt 0.0 and out.lir_chary gt -900)
       djs_plot, out[gd].zdist, out[gd].lir_chary, psym=8, xsty=3, ysty=3, yr=[8.5,11.5]
       djs_oplot, out[lim].zdist, -out[lim].lir_chary, psym=8, color='orange'
;      djs_oplot, out.zdist, out.lir_dale, psym=6, color='cyan'

;; correlate SFR/M with EW(Ha)
;       sfrm = (alog10(4.5D-44)+alog10(3.826D33)+out.lir_chary)-$
;         alog10(nsa.mass)+9
;       plot, nsa.haew, sfrm, psym=6, xr=[0.1,100], yr=[-3,1], /xlog
       
; write out
       lirfile = lcspath+strlowcase(cl[ii].cluster)+'_lir.fits'
       im_mwrfits, out, lirfile, clobber=clobber
    endfor

return
end
    

    
