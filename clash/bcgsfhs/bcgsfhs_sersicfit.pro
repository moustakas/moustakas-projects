function get_radius, radius_sb, nrad=nrad, inrad=inrad, outrad=outrad
; everything in kpc

    radius = radius_sb
    nrad = n_elements(radius)
    
    deltar = (shift(radius,-1)-radius)/2.0
    deltar[nrad-1] = (radius[nrad-1]-radius[nrad-2])/2.0
    outrad = radius+deltar
    inrad = radius-shift(deltar,1)
    inrad[0] = 0.0

;   niceprint, inrad, radius, outrad    
    
; --------------------
;; old but functioning code:    
;    if n_elements(nrad) eq 0 then nrad = 20
;    nrad1 = nrad-1
;    inrad = range(min(radius_sb),rmax,nrad1,/asinh) ; inner radius [kpc]
;    outrad = inrad+(shift(inrad,-1)-inrad)           ; outer radius [kpc]
;
;; adjust the edges
;    outrad[nrad1-1] = (inrad[nrad1-1]-inrad[nrad1-2])+inrad[nrad1-1]
;    outrad = [inrad[0],outrad]
;    inrad = [0D,inrad]
;    radius = (outrad-inrad)/2.0+inrad
; --------------------

;   niceprint, inrad, radius, outrad

return, radius
end

function get_sersic_variance, radius_kpc, sersic=sersic, debug=debug
; variance in the SB profile of a Sersic fit at fixed radius,
; accounting for the covariance in the paramaeters

    nran = 100
    
    rand = mrandomn(seed,sersic.sersic_covar,nran)
    rand[*,0] += sersic.sersic_sb0
    rand[*,1] += sersic.sersic_k
    rand[*,2] += sersic.sersic_n

    sb = fltarr(n_elements(radius_kpc),nran)

    for ii = 0, nran-1 do sb[*,ii] = bcgsfhs_sersic_func(radius_kpc,$
      [rand[ii,0],rand[ii,1],rand[ii,2]])

    sb_var = radius_kpc*0.0
    for ii = 0, n_elements(radius_kpc)-1 do sb_var[ii] = stddev(sb[ii,*])^2
    
; debugging plot    
    if keyword_set(debug) then begin
       dfpsclose & im_plotfaves
       djs_plot, radius_kpc, sb[*,0], xsty=3, ysty=3, /xlog, $
         /ylog, xrange=[min(radius_kpc)>1E-4,max(radius_kpc)], yrange=minmax(sb)
       for jj = 1, nran-1 do djs_oplot, radius_kpc, sb[*,jj]
       oploterror, radius_kpc, bcgsfhs_sersic_func(radius_kpc,params=sersic), $
         sqrt(sb_var), color='blue', errcolor='blue', psym=8
       cc = get_kbrd(1)

; S/N    
       djs_plot, radius_kpc, bcgsfhs_sersic_func(radius_kpc,$
         params=sersic)/sqrt(sb_var), xsty=3, ysty=3
       cc = get_kbrd(1)
    endif
    
return, sb_var
end

pro bcgsfhs_sersic2, rr, sb, scoeff, sb_ivar=sb_ivar, $
  init_params=init_params, fixed=fixed, sersicfit=sersicfit, $
  fixdevac=fixdevac, verbose=verbose
; fit a double Sersic function

; parse the data; take the log
    good = where(sb_ivar gt 0.0 and finite(sb),nn)
    xx = rr[good]
    xxsb = sb[good]
    xxsb_ivar = sb_ivar[good] 

    nparam = 6
    parinfo = replicate({value: 0D, fixed: 0, limited: [0,0], $
      limits: [0D,0D]},nparam)

; sbe1 > 0 
    parinfo[0].limited = [1,0]
    parinfo[0].limits = [1D-12,0D]
; re1>0
    parinfo[1].limited = [1,1]
    parinfo[1].limits = [0.01D,500D]
; 1<n1<10
    parinfo[2].limited = [1,1]
    parinfo[2].limits = [0.1D,10D]

; sbe2 > 0
    parinfo[3].limited = [1,0]
    parinfo[3].limits = [1D-12,0D]
; re2>0
    parinfo[4].limited = [1,1]
    parinfo[4].limits = [0.01D,500D]
; 1<n2<10
    parinfo[5].limited = [1,1]
    parinfo[5].limits = [0.1D,10D]        
    
    if n_elements(init_params) eq nparam then $
      parinfo.value = init_params else $
      parinfo.value = [median(xxsb),1D,1D,0.1*median(xxsb),1D,4D]
    if n_elements(fixed) eq nparam then parinfo.fixed = fixed

    if n_elements(fixed) ne 0 then parinfo.fixed = fixed

    if keyword_set(fixdevac) then begin
       parinfo[5].value = 4D
       parinfo[5].fixed = 1
    endif
    
    params = mpfitfun('bcgsfhs_sersic2_func',xx,xxsb,parinfo=parinfo,$
      yfit=sersicfit,perror=perror,covar=covar,weights=xxsb_ivar,$
      dof=dof,bestnorm=chi2,status=status,quiet=keyword_set(verbose) eq 0,$
      functargs={parinfo: parinfo})

;   djs_plot, rr, -2.5*alog10(xxsb), /xlog, psym=8, xsty=3, ysty=3, $
;     yr=-2.5*alog10(minmax(xxsb))
;   djs_oplot, rr, -2.5*alog10(sersicfit), color='red'
;   djs_oplot, rr, bcgsfhs_sersic_func(rr,[-2.5*alog10(params[0]),params[1],params[2]]), $
;     color='cyan', line=5
;   djs_oplot, rr, bcgsfhs_sersic_func(rr,[-2.5*alog10(params[3]),params[4],params[5]]), $
;     color='orange', line=5
;   cc = get_kbrd(1)
    
    scoeff = {$
      sersic2_status: status,$
      sersic2_chi2:     chi2,$
      sersic2_dof:       dof,$
      sersic2_covar:   covar,$

      sersic2_sbe1: params[0], $
      sersic2_re1:  params[1], $
      sersic2_n1:   params[2], $
      sersic2_sbe1_err: perror[0], $
      sersic2_re1_err:  perror[1], $
      sersic2_n1_err:   perror[2], $

      sersic2_sbe2: params[3], $
      sersic2_re2:  params[4], $
      sersic2_n2:   params[5], $
      sersic2_sbe2_err: perror[3], $
      sersic2_re2_err:  perror[4], $
      sersic2_n2_err:   perror[5]}

;     sersic2_total1:         0.0, $
;     sersic2_total2: 0.0}
;   scoeff.sersic2_total1 = cumsersic_total(params[0:2])
;   scoeff.sersic2_total2 = cumsersic_total(params[3:5])
    
return
end

pro bcgsfhs_sersic, rr, sb, scoeff, sb_ivar=sb_ivar, $
  init_params=init_params, fixed=fixed, sersicfit=sersicfit, $
  verbose=verbose
; fit a single Sersic function; RR should be in kpc; the model
; returned is in magnitudes 

; parse the data; take the log
    good = where(sb_ivar gt 0.0 and finite(sb),nn)
    xx = rr[good]
    xxsb = -2.5*alog10(sb[good])
    xxsb_ivar = sb_ivar[good]*(alog(10)*sb[good]/2.5)^2.0
;   xxsb_ivar = sb_ivar[good]*sb[good]^2.0
;   xxsb = sb[good]
;   xxsb_ivar = sb_ivar[good] 

    nparam = 3
    parinfo = replicate({value: 0D, fixed: 0, limited: [0,0], $
      limits: [0D,0D]},nparam)

; sbe - surface brightness at re
    parinfo[0].limited = [1,1]
    parinfo[0].limits = [10D,35D]
; re>0
    parinfo[1].limited = [1,1]
    parinfo[1].limits = [0.01D,500D]
; 0<n<10
    parinfo[2].limited = [1,1]
    parinfo[2].limits = [0.1D,10D]

    if n_elements(init_params) eq nparam then $
      parinfo.value = init_params else $
      parinfo.value = [median(xxsb), 10D, 4D]
    if n_elements(fixed) eq nparam then parinfo.fixed = fixed
    
    params = mpfitfun('bcgsfhs_sersic_func',xx,xxsb,parinfo=parinfo,$
      perror=perror,covar=covar,weights=xxsb_ivar,bestnorm=chi2,dof=dof,$
      status=status,quiet=keyword_set(verbose) eq 0,yfit=sersicfit,functargs={parinfo: parinfo})

    factor = sqrt(chi2/dof)
    scoeff = {$
      sersic_status: status,$
      sersic_chi2:     chi2,$
      sersic_dof:       dof,$
      sersic_covar:   covar,$

      sersic_sbe: params[0],$
      sersic_re:  params[1],$
      sersic_n:   params[2],$
      sersic_sbe_err: perror[0],$ ; *factor,$
      sersic_re_err:  perror[1],$ ; *factor,$
      sersic_n_err:   perror[2]}  ; *factor,$

;     sersic_total: cumsersic_total(params)}

return
end

pro bcgsfhs_sersic2_multiband, rr, sb, wave, scoeff, sb_ivar=sb_ivar, $
  init_params=init_params, fixed=fixed, sersicfit=sersicfit, $
  fixdevac=fixdevac, verbose=verbose
; fit a double-Sersic function to multiple bands simultaneously; solve
; for the best-fit half-light radius and Sersic n parameter, but allow
; the surface brightness at re to vary 

; parse the data
    good = where(sb_ivar gt 0.0 and finite(sb),nn)
    xx = rr[good]
    xxsb = sb[good]
    xxsb_ivar = sb_ivar[good] 

; figure out how many bands we need to fit and then set up the
; parameters 
    nband = n_elements(uniq(wave,sort(wave)))
    nparam = 2*nband+4

    parinfo = replicate({value: 0D, fixed: 0, $
      limited: [0,0], limits: [0D,0D]},nparam)

; re1>0
    parinfo[0].limited = [1,1]
    parinfo[0].limits = [0.01D,500D]
; 1<n1<10
    parinfo[1].limited = [1,1]
    parinfo[1].limits = [0.1D,10D]
; re2>0
    parinfo[2].limited = [1,1]
    parinfo[2].limits = [0.01D,500D]
; 1<n2<10
    parinfo[3].limited = [1,1]
    parinfo[3].limits = [0.1D,10D]

; sbe1 - surface brightness at re
    for ib = 0, nband-1 do begin
       parinfo[4+ib].limited = [1,0]
;      parinfo[4+ib].limits = [0D,0D]
       parinfo[4+ib].limits = [1D-12,0D]
    endfor

; sbe2 - surface brightness at re
    for ib = 0, nband-1 do begin
       parinfo[4+nband+ib].limited = [1,0]
;      parinfo[4+nband+ib].limits = [0D,0D]
       parinfo[4+nband+ib].limits = [1D-12,0D]
    endfor

    if n_elements(fixed) eq nparam then parinfo.fixed = fixed
    if n_elements(init_params) eq nparam then $
      parinfo.value = init_params else $
        parinfo.value = [1D,1D,30D,4D,replicate(median(xxsb),nband),$
      0.1*replicate(median(xxsb),nband)]

    if keyword_set(fixdevac) then begin
       parinfo[3].value = 4D
       parinfo[3].fixed = 1
    endif

;   struct_print, parinfo
    params = mpfitfun('bcgsfhs_sersic2_multiband_func',xx,xxsb,parinfo=parinfo,$
      yfit=sersicfit,perror=perror,covar=covar,weights=xxsb_ivar,dof=dof,$
      bestnorm=chi2,status=status,quiet=keyword_set(verbose) eq 0,$
      functargs={parinfo: parinfo, wave: wave})

    factor = sqrt(chi2/dof)
    for ib = 0, nband-1 do begin
       params1 = [params[4+ib],params[0],params[1]]
       params2 = [params[4+nband+ib],params[2],params[3]]

       perror1 = [perror[4+ib],perror[0],perror[1]]
       perror2 = [perror[4+nband+ib],perror[2],perror[3]]
       
       scoeff1 = {$
         sersic2_status: status,$
         sersic2_chi2:     chi2,$
         sersic2_dof:       dof,$
         sersic2_covar:   covar,$

         sersic2_sbe1: params1[0],$
         sersic2_re1:  params1[1],$
         sersic2_n1:   params1[2],$
         sersic2_sbe1_err: perror1[0],$ ; *factor,$
         sersic2_re1_err:  perror1[1],$ ; *factor,$
         sersic2_n1_err:   perror1[2],$ ; *factor,$

         sersic2_sbe2: params2[0],$
         sersic2_re2:  params2[1],$
         sersic2_n2:   params2[2],$
         sersic2_sbe2_err: perror2[0],$ ; *factor,$
         sersic2_re2_err:  perror2[1],$ ; *factor,$
         sersic2_n2_err:   perror2[2]}  ; *factor,$
       if ib eq 0 then scoeff = scoeff1 else scoeff = [scoeff,scoeff1]
    endfor

return
end

pro bcgsfhs_sersic_multiband, rr, sb, wave, scoeff, sb_ivar=sb_ivar, $
  init_params=init_params, fixed=fixed, sersicfit=sersicfit, verbose=verbose
; fit multiple bands simultaneously; solve for the best-fit half-light
; radius and Sersic n parameter, but allow the surface brightness at
; re to vary 

; parse the data; take the log
    good = where(sb_ivar gt 0.0 and finite(sb),nn)
    xx = rr[good]
    xxsb = -2.5*alog10(sb[good])
    xxsb_ivar = sb_ivar[good]*(alog(10)*sb[good]/2.5)^2.0
    
; figure out how many bands we need to fit and then set up the
; parameters 
    nband = n_elements(uniq(wave,sort(wave)))
    nparam = nband+2

    parinfo = replicate({value: 0D, fixed: 0, $
      limited: [0,0], limits: [0D,0D]},nparam)

; re>0
    parinfo[0].limited = [1,1]
    parinfo[0].limits = [0.1D,500D]
; 0<n<10
    parinfo[1].limited = [1,1]
    parinfo[1].limits = [0.1D,10D]
; sbe - surface brightness at re
    for ib = 0, nband-1 do begin
       parinfo[2+ib].limited = [1,1]
       parinfo[2+ib].limits = [10D,35D]
    endfor
    
    if n_elements(init_params) eq nparam then $
      parinfo.value = init_params else $
        parinfo.value = [10D, 4D, replicate(median(xxsb),nband)]
    if n_elements(fixed) eq nparam then parinfo.fixed = fixed
    
    params = mpfitfun('bcgsfhs_sersic_multiband_func',xx,xxsb,parinfo=parinfo,$
      perror=perror,covar=covar,weights=xxsb_ivar,bestnorm=chi2,dof=dof,$
      status=status,quiet=keyword_set(verbose) eq 0,yfit=sersicfit,$
      functargs={parinfo: parinfo, wave: wave})

    factor = sqrt(chi2/dof)
    for ib = 0, nband-1 do begin
       params1 = [params[2+ib],params[0],params[1]]
;      covar1 = fltarr(3,3)
;      for ii = 0, 2 do covar1[*,ii] = covar[[0,1,2+ib]]
       
       scoeff1 = {$
         sersic_status: status,$
         sersic_chi2:     chi2,$
         sersic_dof:       dof,$
         sersic_covar:   covar,$

         sersic_sbe: params[2+ib],$
         sersic_re:  params[0],$
         sersic_n:   params[1],$
         sersic_sbe_err: perror[2+ib],$ ; *factor,$
         sersic_re_err:  perror[0],$ ; *factor,$
         sersic_n_err:   perror[1]}  ; *factor,$

;        sersic_total: cumsersic_total(params1)}
       if ib eq 0 then scoeff = scoeff1 else scoeff = [scoeff,scoeff1]
    endfor

return
end

pro bcgsfhs_sersicfit, dofit=dofit, dophot=dophot, clobber=clobber, $
  qaplot_seds=qaplot_seds, qaplot_sbprofiles=qaplot_sbprofiles, $
  qaplot_colorradius=qaplot_colorradius, verbose=verbose
; jm13oct22siena - fit various Sersic models to the output of
; BCGSFHS_ELLIPSE 

    ellpath = bcgsfhs_path(/ellipse)
    sersicpath = bcgsfhs_path(/sersic)

; read the sample
    sample = read_bcgsfhs_sample()

    splog, 'Come back to A2261!'    
;   sample = sample[13]
    ncl = n_elements(sample)

    pixscale = 0.065D           ; [arcsec/pixel]
    errfloor = 0.0D ; 0.02      ; magnitude error floor on my SB measurements 
    nphotradius = 15            ; number of radial bins
    photradius_multiplier = range(0.03,3.0,nphotradius,/log) ; search for "snippet", below

; do 2-component Sersic fitting?    
    dosersic2 = 1
    fixdevac = 1 ; fix the 2nd Sersic model to n=4?

; ##################################################
; fit single and double Sersic models to every band
    if keyword_set(dofit) then begin
; wrap on each cluster    
;      for ic = 2, 2 do begin
       for ic = 0, ncl-1 do begin
          arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
          
          cluster = strtrim(sample[ic].shortname,2)
          splog, 'Sersic fitting cluster '+cluster
          
; read the data and Marc's SB profiles to determine the last radius at
; which the models are reliable
          modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
          allband = strtrim(modphot.band,2)
          
          pp = read_bcg_profiles(cluster,these_filters=allband)
          nfilt = n_elements(modphot)

; first fit the (rest-frame) near-IR bands simultaneously, then the
; optical bands together, and finally the bluest bands; we do this
; because to help constrain the low S/N bands
          delvarx, nir, opt, blue, nir_multiband
          nir_multiband = 0 ; by default fit the near-IR bands individually
          case cluster of
             'a209': begin
                blue = where(allband eq 'f475w' or allband eq 'f435w' or allband eq 'f390w')
                nir = lindgen(nfilt) & remove, [blue], nir ; everything else
             end
             'a383': begin
                opt = where(allband eq 'f775w' or allband eq 'f814w')
                blue = where(allband eq 'f475w' or allband eq 'f435w')
                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
             end
             'macs0744': begin
                blue = where(allband eq 'f475w' or allband eq 'f555w')
                nir = lindgen(nfilt) & remove, [blue], nir ; everything else
             end
             'a611': begin
                opt = where(allband eq 'f850lp' or allband eq 'f606w' or $
                  allband eq 'f814w' or allband eq 'f775w')
                blue = where(allband eq 'f475w' or allband eq 'f435w')
                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
                nir_multiband = 1 ; constrain!
             end
             'macs1149': begin
                opt = where(allband eq 'f850lp' or allband eq 'f606w' or $
                  allband eq 'f814w' or allband eq 'f775w')
                blue = where(allband eq 'f475w' or allband eq 'f435w')
                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
                nir_multiband = 1 ; constrain!
             end
             'a1423': begin
                opt = where(allband eq 'f850lp' or allband eq 'f814w' or $
                  allband eq 'f775w' or allband eq 'f625w')
                blue = where(allband eq 'f606w' or allband eq 'f475w' or $
                  allband eq 'f435w' or allband eq 'f390w')
                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
             end
             'macs1206': begin
                opt = where(allband eq 'f850lp' or allband eq 'f814w' or $
                  allband eq 'f775w' or allband eq 'f625w')
                blue = where(allband eq 'f606w' or allband eq 'f475w' or $
                  allband eq 'f435w' or allband eq 'f390w')
                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
                nir_multiband = 1 ; constrain!
             end
             'clj1226': begin
                blue = -1
                opt = where(allband eq 'f850lp' or allband eq 'f814w' or $
                  allband eq 'f775w' or allband eq 'f625w' or allband eq 'f606w')
                nir = lindgen(nfilt) & remove, [opt], nir ; everything else
                nir_multiband = 1 ; constrain!
             end
             'macs1311': begin
                opt = where(allband eq 'f775w' or allband eq 'f625w' or allband eq 'f606w')
                blue = where(allband eq 'f475w' or allband eq 'f435w')
                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
             end
             'macs1720': begin
                opt = where(allband eq 'f850lp' or allband eq 'f814w' or $
                  allband eq 'f775w' or allband eq 'f625w')
                blue = where(allband eq 'f606w' or allband eq 'f475w')
                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
             end
             'macs2129': begin
                opt = where(allband eq 'f850lp' or allband eq 'f814w')
                blue = where(allband eq 'f475w' or allband eq 'f435w')
                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
             end
             'rxj2129': begin
                opt = where(allband eq 'f850lp' or allband eq 'f814w' or $
                  allband eq 'f775w' or allband eq 'f625w' or allband eq 'f606w')
                blue = where(allband eq 'f475w' or allband eq 'f435w')
                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
                nir_multiband = 1 ; constrain!
             end
; the blue bands for ms2137 don't fit well!
             'ms2137': begin
                opt = where(allband eq 'f850lp' or allband eq 'f814w' or $
                  allband eq 'f775w' or allband eq 'f625w')
                blue = where(allband eq 'f475w' or allband eq 'f435w' or allband eq 'f390w')
                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
                nir_multiband = 1 ; constrain!
             end
             'rxj2248': begin
                opt = where(allband eq 'f850lp' or allband eq 'f814w' or $
                  allband eq 'f775w' or allband eq 'f625w' or allband eq 'f606w')
                blue = where(allband eq 'f475w' or allband eq 'f435w')
                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
                nir_multiband = 1 ; constrain!
             end
             else: begin
                opt = -1
                blue = where(modphot.weff/1D4 lt 0.5,comp=nir)
             end
          endcase
          nnir = n_elements(nir)
          if opt[0] eq -1 then nopt = 0 else nopt = n_elements(opt)
          if blue[0] eq -1 then nblue = 0 else nblue = n_elements(blue)

;         nopt = 0
;         nir = where(modphot.weff/(1+sample[ic].z)/1D4 ge 0.4,nnir,comp=blue,ncomp=nblue)
;         splog, 'Tying blue bands: '+strjoin(modphot[blue].band,', ')
          
;         nir = where(modphot.weff/(1+sample[ic].z)/1D4 ge 0.6,nnir)
;         opt = where(modphot.weff/(1+sample[ic].z)/1D4 lt 0.6 and $
;           modphot.weff/(1+sample[ic].z)/1D4 gt 0.4,nopt)
;         blue = where(modphot.weff/(1+sample[ic].z)/1D4 le 0.4,nblue)

; output structure
          out = struct_addtags(struct_trimtags(modphot,$
            select=['file','band','weff','sblimit','ra','dec','mge_*']),$
            replicate({amin_kpc: 0.0, amax_kpc: 0.0, rmin_kpc: 0.0, $
            rmax_kpc: 0.0, $
            sersic_covar_nir:  fltarr(2+nnir,2+nnir),$
            sersic_covar_opt:  fltarr(2+nopt,2+nopt),$
            sersic_covar_blue: fltarr(2+nblue,2+nblue)},nfilt))
          if dosersic2 then begin
             out = struct_addtags(out,replicate({$
               sersic2_covar_nir:  fltarr(2*nnir+4,2*nnir+4),$
               sersic2_covar_opt:  fltarr(2*nopt+4,2*nopt+4),$
               sersic2_covar_blue: fltarr(2*nblue+4,2*nblue+4)},nfilt))
          endif

; near-IR, fitted individually
          if nir_multiband eq 0 then begin
             for ib = 0, nnir-1 do begin

; if AMAX_KPC occurs when the SB profile is below the formal 1-sigma
; surface brightness limit of the data, then cut it off
                amax_kpc = max(pp[nir[ib]].sma,mxindx) ; [kpc]
                if pp[nir[ib]].mu[mxindx] gt modphot[nir[ib]].sblimit then $
                  amax_kpc = pp[nir[ib]].sma[mxindx-1]
                
                modgood = where(modphot[nir[ib]].majora*pixscale*arcsec2kpc le amax_kpc and $
                  modphot[nir[ib]].sb0fit gt 0 and modphot[nir[ib]].sb0fit_ivar gt 0,nmodgood)
                amin_kpc = min(modphot[nir[ib]].majora[modgood])*pixscale*arcsec2kpc

; the equivalent radius needs to be sorted (note: the ellipse
; parameters in BCGSFHS_ELLIPSE are monotonic in *semi-major axis*,
; not equivalent radius!)  also add a minimum error floor to the
; surface brightnesses
                sb = modphot[nir[ib]].sb0fit[modgood]*1D
                sb_var_floor = (sb*errfloor)^2.0
                sb_ivar = 1D/(1D/modphot[nir[ib]].sb0fit_ivar[modgood]+sb_var_floor)
                radius_kpc = modphot[nir[ib]].radius_kpc[modgood] ; [kpc]
                
                rmin_kpc = min(radius_kpc)
                rmax_kpc = max(radius_kpc)
                
                out[nir[ib]].amin_kpc = amin_kpc
                out[nir[ib]].amax_kpc = amax_kpc
                out[nir[ib]].rmin_kpc = rmin_kpc
                out[nir[ib]].rmax_kpc = rmax_kpc
                
                srt = sort(radius_kpc)
                radius_kpc = radius_kpc[srt]
                sb = sb[srt]
                sb_ivar = sb_ivar[srt]
             
; fit with a single-Sersic and then a double-Sersic
                bcgsfhs_sersic, radius_kpc, sb, sersic, sb_ivar=sb_ivar, verbose=verbose
                if dosersic2 then bcgsfhs_sersic2, radius_kpc, sb, sersic2, $
                  sb_ivar=sb_ivar, fixdevac=fixdevac, verbose=verbose

; expand the output structure
                if ib eq 0 then begin
                   out = struct_addtags(out,im_empty_structure(sersic[0],ncopies=nfilt))
                   if dosersic2 then out = struct_addtags(out,im_empty_structure(sersic2[0],ncopies=nfilt))
                endif 
             
                out[nir[ib]] = im_struct_assign(sersic,out[nir[ib]],/nozero)
                if dosersic2 then out[nir[ib]] = im_struct_assign(sersic2,out[nir[ib]],/nozero)
             endfor
          endif

; near-IR, fitted together
          if nir_multiband then begin
             for ib = 0, nnir-1 do begin
                amax_kpc = max(pp[nir[ib]].sma,mxindx) ; [kpc]
                if pp[nir[ib]].mu[mxindx] gt modphot[nir[ib]].sblimit then $
                  amax_kpc = pp[nir[ib]].sma[mxindx-1]

                modgood = where(modphot[nir[ib]].majora*pixscale*arcsec2kpc le amax_kpc and $
                  modphot[nir[ib]].sb0fit gt 0 and modphot[nir[ib]].sb0fit_ivar gt 0,nmodgood)
                amin_kpc = min(modphot[nir[ib]].majora[modgood])*pixscale*arcsec2kpc
                
                sb = modphot[nir[ib]].sb0fit[modgood]*1D
                sb_var_floor = (sb*errfloor)^2.0
                sb_ivar = 1D/(1D/modphot[nir[ib]].sb0fit_ivar[modgood]+sb_var_floor)
                radius_kpc = modphot[nir[ib]].radius_kpc[modgood] ; [kpc]

                rmin_kpc = min(radius_kpc)
                rmax_kpc = max(radius_kpc)
                
                out[nir[ib]].amin_kpc = amin_kpc
                out[nir[ib]].amax_kpc = amax_kpc
                out[nir[ib]].rmin_kpc = rmin_kpc
                out[nir[ib]].rmax_kpc = rmax_kpc
             
                srt = sort(radius_kpc)
                radius_kpc = radius_kpc[srt]
                sb = sb[srt]
                sb_ivar = sb_ivar[srt]
             
                if ib eq 0 then begin
                   fit_sb = sb
                   fit_sb_ivar = sb_ivar
                   fit_radius_kpc = radius_kpc
                   fit_wave = replicate(modphot[nir[ib]].weff,nmodgood)
                endif else begin
                   fit_sb = [fit_sb,sb]
                   fit_sb_ivar = [fit_sb_ivar,sb_ivar]
                   fit_radius_kpc = [fit_radius_kpc,radius_kpc]
                   fit_wave = [fit_wave,replicate(modphot[nir[ib]].weff,nmodgood)]
                endelse
             endfor
             
             bcgsfhs_sersic_multiband, fit_radius_kpc, fit_sb, fit_wave, $
               multisersic, sb_ivar=fit_sb_ivar, verbose=verbose
             out = struct_addtags(out,im_empty_structure(multisersic[0],ncopies=nfilt))
             
             for ib = 0, nnir-1 do begin
                out[nir[ib]] = im_struct_assign(multisersic[ib],out[nir[ib]],/nozero)
                out[nir[ib]].sersic_covar = 0
                out[nir[ib]].sersic_covar_nir = multisersic[ib].sersic_covar
             endfor
             
             if dosersic2 then begin
                bcgsfhs_sersic2_multiband, fit_radius_kpc, fit_sb, fit_wave, $
                  multisersic2, sb_ivar=fit_sb_ivar, sersicfit=sersicfit, $
                  fixdevac=fixdevac, verbose=verbose
                out = struct_addtags(out,im_empty_structure(multisersic2[0],ncopies=nfilt))
                for ib = 0, nnir-1 do begin
                   out[nir[ib]] = im_struct_assign(multisersic2[ib],out[nir[ib]],/nozero)
                   out[nir[ib]].sersic2_covar = 0
                   out[nir[ib]].sersic2_covar_nir = multisersic2[ib].sersic2_covar
                endfor
             endif
          endif

; optical
          if nopt gt 0 then begin
             for ib = 0, nopt-1 do begin
                amax_kpc = max(pp[opt[ib]].sma,mxindx) ; [kpc]
                if pp[opt[ib]].mu[mxindx] gt modphot[opt[ib]].sblimit then $
                  amax_kpc = pp[opt[ib]].sma[mxindx-1]
                
                modgood = where(modphot[opt[ib]].majora*pixscale*arcsec2kpc le amax_kpc and $
                  modphot[opt[ib]].sb0fit gt 0 and modphot[opt[ib]].sb0fit_ivar gt 0,nmodgood)
                amin_kpc = min(modphot[opt[ib]].majora[modgood])*pixscale*arcsec2kpc
                
                sb = modphot[opt[ib]].sb0fit[modgood]*1D
                sb_var_floor = (sb*errfloor)^2.0
                sb_ivar = 1D/(1D/modphot[opt[ib]].sb0fit_ivar[modgood]+sb_var_floor)
                radius_kpc = modphot[opt[ib]].radius_kpc[modgood] ; [kpc]

                rmin_kpc = min(radius_kpc)
                rmax_kpc = max(radius_kpc)

                out[opt[ib]].amin_kpc = amin_kpc
                out[opt[ib]].amax_kpc = amax_kpc
                out[opt[ib]].rmin_kpc = rmin_kpc
                out[opt[ib]].rmax_kpc = rmax_kpc
             
                srt = sort(radius_kpc)
                radius_kpc = radius_kpc[srt]
                sb = sb[srt]
                sb_ivar = sb_ivar[srt]
             
                if ib eq 0 then begin
                   fit_sb = sb
                   fit_sb_ivar = sb_ivar
                   fit_radius_kpc = radius_kpc
                   fit_wave = replicate(modphot[opt[ib]].weff,nmodgood)
                endif else begin
                   fit_sb = [fit_sb,sb]
                   fit_sb_ivar = [fit_sb_ivar,sb_ivar]
                   fit_radius_kpc = [fit_radius_kpc,radius_kpc]
                   fit_wave = [fit_wave,replicate(modphot[opt[ib]].weff,nmodgood)]
                endelse
             endfor
             
             bcgsfhs_sersic_multiband, fit_radius_kpc, fit_sb, fit_wave, $
               multisersic, sb_ivar=fit_sb_ivar, verbose=verbose
             for ib = 0, nopt-1 do begin
                out[opt[ib]] = im_struct_assign(multisersic[ib],out[opt[ib]],/nozero)
                out[opt[ib]].sersic_covar = 0
                out[opt[ib]].sersic_covar_opt = multisersic[ib].sersic_covar
             endfor
          
             if dosersic2 then begin
                bcgsfhs_sersic2_multiband, fit_radius_kpc, fit_sb, fit_wave, $
                  multisersic2, sb_ivar=fit_sb_ivar, sersicfit=sersicfit, $
                  fixdevac=fixdevac, verbose=verbose
                for ib = 0, nopt-1 do begin
                   out[opt[ib]] = im_struct_assign(multisersic2[ib],out[opt[ib]],/nozero)
                   out[opt[ib]].sersic2_covar = 0
                   out[opt[ib]].sersic2_covar_opt = multisersic2[ib].sersic2_covar
                endfor
             endif
          endif

; blue
          if nblue ne 0 then begin
             for ib = 0, nblue-1 do begin
                amax_kpc = max(pp[blue[ib]].sma,mxindx) ; [kpc]
                if pp[blue[ib]].mu[mxindx] gt modphot[blue[ib]].sblimit then $
                  amax_kpc = pp[blue[ib]].sma[mxindx-1]
                
                modgood = where(modphot[blue[ib]].majora*pixscale*arcsec2kpc le amax_kpc and $
                  modphot[blue[ib]].sb0fit gt 0 and modphot[blue[ib]].sb0fit_ivar gt 0,nmodgood)
                amin_kpc = min(modphot[blue[ib]].majora[modgood])*pixscale*arcsec2kpc
                
                sb = modphot[blue[ib]].sb0fit[modgood]*1D
                sb_var_floor = (sb*errfloor)^2.0
                sb_ivar = 1D/(1D/modphot[blue[ib]].sb0fit_ivar[modgood]+sb_var_floor)
                radius_kpc = modphot[blue[ib]].radius_kpc[modgood] ; [kpc]

                rmin_kpc = min(radius_kpc)
                rmax_kpc = max(radius_kpc)

                out[blue[ib]].amin_kpc = amin_kpc
                out[blue[ib]].amax_kpc = amax_kpc
                out[blue[ib]].rmin_kpc = rmin_kpc
                out[blue[ib]].rmax_kpc = rmax_kpc
             
                srt = sort(radius_kpc)
                radius_kpc = radius_kpc[srt]
                sb = sb[srt]
                sb_ivar = sb_ivar[srt]
             
                if ib eq 0 then begin
                   fit_sb = sb
                   fit_sb_ivar = sb_ivar
                   fit_radius_kpc = radius_kpc
                   fit_wave = replicate(modphot[blue[ib]].weff,nmodgood)
                endif else begin
                   fit_sb = [fit_sb,sb]
                   fit_sb_ivar = [fit_sb_ivar,sb_ivar]
                   fit_radius_kpc = [fit_radius_kpc,radius_kpc]
                   fit_wave = [fit_wave,replicate(modphot[blue[ib]].weff,nmodgood)]
                endelse
             endfor
             
             bcgsfhs_sersic_multiband, fit_radius_kpc, fit_sb, fit_wave, $
               multisersic, sb_ivar=fit_sb_ivar, verbose=verbose
             for ib = 0, nblue-1 do begin
                out[blue[ib]] = im_struct_assign(multisersic[ib],out[blue[ib]],/nozero)
                out[blue[ib]].sersic_covar = 0
                out[blue[ib]].sersic_covar_blue = multisersic[ib].sersic_covar
             endfor
             
             if dosersic2 then begin
                bcgsfhs_sersic2_multiband, fit_radius_kpc, fit_sb, fit_wave, $
                  multisersic2, sb_ivar=fit_sb_ivar, sersicfit=sersicfit, $
                  fixdevac=fixdevac, verbose=verbose
                for ib = 0, nblue-1 do begin
                   out[blue[ib]] = im_struct_assign(multisersic2[ib],out[blue[ib]],/nozero)
                   out[blue[ib]].sersic2_covar = 0
                   out[blue[ib]].sersic2_covar_blue = multisersic2[ib].sersic2_covar
                endfor
             endif
          endif

; write out
          im_mwrfits, out, sersicpath+cluster+'-sersic.fits', clobber=clobber
       endfor
    endif
    
; ##################################################
; build QAplots
    ncol = 3 ; number of columns
    rr = [0,range(0.01,200,500,/log)] ; equivalent radius [kpc]

; -------------------------    
; QAplot: SEDs
    if keyword_set(qaplot_seds) then begin
       psfile = sersicpath+'qa_seds.ps'
       im_plotconfig, 0, pos, psfile=psfile, charsize=1.8, $
         height=5.0
    
       xrange = [0.3,2.0]
       yrange = [26,10]
    
       for ic = 0, ncl-1 do begin
          cluster = strtrim(sample[ic].shortname,2)
          phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
          nrad = n_elements(phot[0].photradius_kpc)
          nfilt = n_elements(phot)

          djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
            xrange=xrange, yrange=yrange, xtitle='Wavelength \lambda (\mu'+'m)', $
            ytitle='Magnitude (AB)', /xlog, title=strupcase(cluster)

; integrated light       
          mab = maggies2mag(phot.maggies_int,ivarmaggies=phot.ivarmaggies_int,$
            magerr=maberr,lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
          used = where(mab gt -90.0,nused)
          upper = where(mab lt -90.0 and mabupper gt -90,nupper)
          
          oploterror, phot[used].weff/1D4, mab[used], mabhierr[used], $
            psym=-symcat(16), symsize=1.5, color=cgcolor('firebrick'), $
            /hibar, errcolor=cgcolor('firebrick')
          oploterror, phot[used].weff/1D4, mab[used], mabloerr[used], psym=3, $
            color=cgcolor('firebrick'), /lobar, errcolor=cgcolor('firebrick')

; radial bins       
          for ir = 0, nrad-1 do begin
             mab = maggies2mag(phot.maggies[ir],ivarmaggies=phot.ivarmaggies[ir],$
               magerr=maberr,lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
             used = where(mab gt -90.0,nused)
             upper = where(mab lt -90.0 and mabupper gt -90,nupper)
             if (nused ne 0) then begin
                oploterror, phot[used].weff/1D4, mab[used], mabhierr[used], $
                  psym=-symcat(16), symsize=1.5, color=cgcolor('dodger blue'), $
                  /hibar, errcolor=cgcolor('dodger blue')
                oploterror, phot[used].weff/1D4, mab[used], mabloerr[used], psym=3, $
                  color=cgcolor('dodger blue'), /lobar, errcolor=cgcolor('dodger blue')
             endif
             if (nupper ne 0) then begin
                djs_oplot, [phot[upper].weff/1D4], [mabupper[upper]], $
                  psym=symcat(11,thick=6), symsize=2.0, color=cgcolor('forest green')
             endif
          endfor 
       endfor 
       im_plotconfig, psfile=psfile, /psclose, /pdf
    endif

; -------------------------    
; QAplot: color-radius plots, relative to F160W
    if keyword_set(qaplot_colorradius) then begin
       psfile = sersicpath+'qa_color_sersic.ps'
       im_plotconfig, 0, pos, psfile=psfile, charsize=1.1

       for ic = 0, ncl-1 do begin
          cluster = strtrim(sample[ic].shortname,2)
          splog, cluster

          arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
          
          sersic = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent)
          modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
          nfilt = n_elements(modphot)-1 ; relative to F160W
          
          nrow = ceil(nfilt/float(ncol))
          pos = im_getposition(nx=ncol,ny=nrow,yspace=0.1,xspace=0.8*[1,1],$
            xmargin=[1.0,0.2],width=1.9)
          for ib = 1, nfilt-1 do begin
             band = strtrim(strupcase(modphot[ib].band),2)

             modgood = where($
               modphot[ib].majora*pixscale*arcsec2kpc le sersic[ib].amax_kpc and $
               modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0 and $
               modphot[0].majora*pixscale*arcsec2kpc le sersic[0].amax_kpc and $
               modphot[0].sb0fit gt 0 and modphot[0].sb0fit_ivar gt 0)
             modbad = where($
               (modphot[ib].majora*pixscale*arcsec2kpc gt sersic[ib].amax_kpc or $
               modphot[0].majora*pixscale*arcsec2kpc gt sersic[0].amax_kpc) and $
               modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0 and $
               modphot[0].sb0fit gt 0 and modphot[0].sb0fit_ivar gt 0)
             
             radius_kpc = modphot[ib].radius_kpc[modgood]  ; [kpc]
             c1 = -2.5*alog10(modphot[ib].sb0fit[modgood])
             c2 = -2.5*alog10(modphot[0].sb0fit[modgood])
             c1err = 2.5/(alog(10)*modphot[ib].sb0fit[modgood]*sqrt(modphot[ib].sb0fit_ivar[modgood]))
             c2err = 2.5/(alog(10)*modphot[0].sb0fit[modgood]*sqrt(modphot[0].sb0fit_ivar[modgood]))
;            c1err = modphot[ib].sb0fit_ivar[modgood]*(modphot[ib].sb0fit[modgood]*alog(10)/2.5)^2
             
             color = c1 - c2
             colorerr = sqrt(c1err^2+c2err^2)

             colorbad = -2.5*alog10(modphot[ib].sb0fit[modbad]/modphot[0].sb0fit[modbad])

             yrange = [min(color-colorerr),max(color+colorerr)]
;            splog, band, yrange
             
             sersiccolor = bcgsfhs_sersic_func(rr,params=sersic[ib])-$
               bcgsfhs_sersic_func(rr,params=sersic[0])
             if dosersic2 then sersic2color = -2.5*alog10(bcgsfhs_sersic2_func(rr,params=sersic[ib])/$
               bcgsfhs_sersic2_func(rr,params=sersic[0]))
             
             if ib-1 eq 1 then title = strupcase(cluster) else delvarx, title
             if ib-1 ge nfilt-4 then begin
                xtitle = 'Equivalent Radius (kpc)'
                delvarx, xtickname
             endif else begin
                xtitle = ''
                xtickname = replicate(' ',10)
             endelse
;            if (ib-1 mod 3) eq 0 then begin
;               delvarx, ytickname
;            endif else begin
;               ytickname = replicate(' ',10)
;            endelse
             
             djs_plot, [0], [0], /nodata, /xlog, noerase=ib-1 gt 0, $
               xrange=[0.3,200], xsty=1, yrange=yrange, position=pos[*,ib-1], $
               xtickname=xtickname, ytickname=ytickname, title=title, $
               symsize=0.5, ysty=1, ytitle=band+'-F160W', xtitle=xtitle
             oploterror, radius_kpc, color, colorerr, psym=3, /nohat, $
               errcolor=cgcolor('light grey')
             djs_oplot, radius_kpc, color, psym=symcat(16), symsize=0.5

;            djs_oplot, modphot[ib].radius_kpc[modbad], colorbad, $
;              psym=symcat(6), symsize=0.3; color=cgcolor('medium grey')
             
;            im_legend, band+'-F160W', /right, /top, box=0, margin=0, charsize=1.0
             djs_oplot, rr, sersiccolor, color=cgcolor('firebrick')
             if dosersic2 then djs_oplot, rr, sersic2color, color=cgcolor('dodger blue')
          endfor
          
;         xyouts, min(pos[0,*])-0.06, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
;           textoidl(''), orientation=90, align=0.5, charsize=1.4, /norm
;         xyouts, (max(pos[2,*])-min(pos[0,*]))/2.0+min(pos[0,*]), min(pos[1,*])-0.06, $
;           textoidl('Equivalent Radius (kpc)'), align=0.5, charsize=1.1, /norm
       endfor
       im_plotconfig, psfile=psfile, /psclose, /pdf
    endif

; -------------------------    
; QAplot: SB profiles and the Sersic fits
    if keyword_set(qaplot_sbprofiles) then begin
       psfile = sersicpath+'qa_sersic.ps'
       im_plotconfig, 0, pos, psfile=psfile, charsize=1.3

       for ic = 0, ncl-1 do begin
          cluster = strtrim(sample[ic].shortname,2)
          splog, cluster

          arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]

          sersic = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent)
          modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
          nfilt = n_elements(modphot)

          nrow = ceil(nfilt/float(ncol))
          pos = im_getposition(nx=ncol,ny=nrow,yspace=0.0,xspace=0.0,$
            xmargin=[0.9,0.4],width=2.4)
          for ib = 0, nfilt-1 do begin
             band = strtrim(strupcase(modphot[ib].band),2)

             modgood = where(modphot[ib].majora*pixscale*arcsec2kpc le sersic[ib].amax_kpc and $
               modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0)
             modbad = where(modphot[ib].majora*pixscale*arcsec2kpc gt sersic[ib].amax_kpc and $
               modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0)
             
             sb = modphot[ib].sb0fit[modgood]*1D
             sb_var_floor = (sb*errfloor)^2.0
             sb_ivar = 1D/(1D/modphot[ib].sb0fit_ivar[modgood]+sb_var_floor)

             radius_kpc = modphot[ib].radius_kpc[modgood]  ; [kpc]

             if ib eq 1 then title = strupcase(cluster) else delvarx, title
             if ib ge nfilt-3 then begin
                delvarx, xtickname
             endif else begin
                xtickname = replicate(' ',10)
             endelse
             if (ib mod 3) eq 0 then begin
                delvarx, ytickname
             endif else begin
                ytickname = replicate(' ',10)
             endelse
             
             djs_plot, radius_kpc, -2.5*alog10(sb), psym=symcat(16), /xlog, noerase=ib gt 0, $
               xrange=[0.3,200], xsty=1, yrange=[29,16], position=pos[*,ib], $
               xtickname=xtickname, ytickname=ytickname, title=title, $
               symsize=0.5, ytickinterval=3, ysty=1

             label = [$
               '\chi^{2}_{\nu, single}='+$
               strtrim(string(sersic[ib].sersic_chi2/sersic[ib].sersic_dof,format='(F12.2)'),2),$
               '\mu_{e}='+strtrim(string(sersic[ib].sersic_sbe,format='(F12.1)'),2)+','+$
               'n='+strtrim(string(sersic[ib].sersic_n,format='(F12.2)'),2)+','+$
               'r_{e}='+strtrim(string(sersic[ib].sersic_re,format='(F12.1)'),2)+' kpc']
             
;            label = [$
;              '\mu_{e}='+strtrim(string(sersic[ib].sersic_sbe,format='(F12.1)'),2),$
;              'r_{e}='+strtrim(string(sersic[ib].sersic_re,format='(F12.1)'),2)+' kpc',$
;              'n='+strtrim(string(sersic[ib].sersic_n,format='(F12.2)'),2),$
;              '\chi^{2}_{\nu, single}='+strtrim(string(sersic[ib].sersic_chi2/$
;              sersic[ib].sersic_dof,format='(F12.2)'),2)]
             if dosersic2 then begin
                if sersic[ib].sersic2_sbe1 eq 0.0 or sersic[ib].sersic2_sbe2 eq 0.0 then begin
                   label = [label,'Sersic-2 dropped']
                endif else begin
                   label = [label,$
                     '\chi^{2}_{\nu, double}='+strtrim(string(sersic[ib].sersic2_chi2/$
                     sersic[ib].sersic2_dof,format='(F12.2)'),2),$
                     '\mu_{e1}='+strtrim(string(-2.5*alog10(sersic[ib].sersic2_sbe1),format='(F12.1)'),2)+','+$
                     'n_{1}='+strtrim(string(sersic[ib].sersic2_n1,format='(F12.2)'),2)+','+$
                     'r_{e1}='+strtrim(string(sersic[ib].sersic2_re1,format='(F12.2)'),2)+' kpc',$
                     '\mu_{e2}='+strtrim(string(-2.5*alog10(sersic[ib].sersic2_sbe2),format='(F12.1)'),2)+','+$
                     'n_{2}='+strtrim(string(sersic[ib].sersic2_n2,format='(F12.2)'),2)+','+$
                     'r_{e2}='+strtrim(string(sersic[ib].sersic2_re2,format='(F12.1)'),2)+' kpc']
                endelse
             endif
             im_legend, label, /left, /bottom, box=0, margin=0, charsize=0.7, charthick=1.8

             djs_oplot, modphot[ib].radius_kpc[modbad], -2.5*alog10(modphot[ib].sb0fit[modbad]), $
               psym=symcat(9), color=cgcolor('medium grey'), symsize=0.5
             
             im_legend, band, /right, /top, box=0, margin=0, charsize=1.0
             
             djs_oplot, rr, bcgsfhs_sersic_func(rr,params=sersic[ib]), $
               color=cgcolor('firebrick')
             if ib gt 0 then djs_oplot, rr, bcgsfhs_sersic_func(rr,$
               params=sersic[0]), color=cgcolor('forest green')
             
             if dosersic2 then begin
                djs_oplot, rr, -2.5*alog10(bcgsfhs_sersic2_func(rr,params=sersic[ib])), $
                  color=cgcolor('dodger blue')
                if sersic[ib].sersic2_sbe1 eq 0.0 or sersic[ib].sersic2_sbe2 eq 0.0 then begin
                   splog, '  '+band+': second Sersic dropped!'
                endif else begin
                   djs_oplot, rr, bcgsfhs_sersic_func(rr,[-2.5*alog10(sersic[ib].sersic2_sbe1),$
                     sersic[ib].sersic2_re1,sersic[ib].sersic2_n1]), color=cgcolor('orange'), line=2
                   djs_oplot, rr, bcgsfhs_sersic_func(rr,[-2.5*alog10(sersic[ib].sersic2_sbe2),$
                     sersic[ib].sersic2_re2,sersic[ib].sersic2_n2]), color=cgcolor('orange'), line=2
;                  if band eq 'F435W' then stop
                endelse
             endif

             djs_oplot, [70.0,10^!x.crange[1]], modphot[ib].sblimit*[1,1], line=0
;            djs_oplot, [6.0,10^!x.crange[1]], modphot[ib].sblimit*[1,1], line=0
;            djs_oplot, sersic[ib].rmax_kpc*[1,1], [!y.crange[0]-0.2,!y.crange[0]-6], line=0
          endfor
          
          xyouts, min(pos[0,*])-0.06, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
            textoidl('\mu (mag arcsec^{-2})'), orientation=90, align=0.5, charsize=1.4, /norm
          xyouts, (max(pos[2,*])-min(pos[0,*]))/2.0+min(pos[0,*]), min(pos[1,*])-0.06, $
            textoidl('Equivalent Radius (kpc)'), align=0.5, charsize=1.4, /norm
       endfor
       im_plotconfig, psfile=psfile, /psclose, /pdf

    endif
       
; ##################################################
; do radial aperture photometry in each band, using the Sersic models
; to extrapolate inward and outward
    if keyword_set(dophot) then begin

;; -------------------------
;; from this snippet of code, I conclude that we should measure photometry
;; in NRAD radial apertures from approximately 0.05-3 times Re in the
;; F160W band 
;       for ic = 0, ncl-1 do begin
;          cluster = strtrim(sample[ic].shortname,2)
;          jj = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent,row=0)
;          print, jj.sersic_re, jj.rmin_kpc, jj.rmax_kpc, jj.rmin_kpc/jj.sersic_re, $
;            jj.rmax_kpc/jj.sersic_re, '   '+cluster
;       endfor
;; -------------------------

; wrap on each cluster    
;      for ic = 3, 3 do begin
       for ic = 0, ncl-1 do begin
          arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]

          cluster = strtrim(sample[ic].shortname,2)
          splog, 'Measuring aperture photometry for cluster '+cluster

          sersic = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent)
          modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
          nfilt = n_elements(modphot)

; initialize the radial photometry structure using F160W as the
; reference band
          photradius_kpc = get_radius(photradius_multiplier*30.0,$ ; fixed radii for all clusters
            nrad=nphotradius,inrad=photradius_kpc_in,outrad=photradius_kpc_out)
;         photradius_kpc = get_radius(photradius_multiplier*sersic[0].sersic2_re2,$ ; use 2nd component!
;           nrad=nphotradius,inrad=photradius_kpc_in,outrad=photradius_kpc_out)
;         photradius_kpc = get_radius(photradius_multiplier*sersic[0].sersic_re,$
;           nrad=nphotradius,inrad=photradius_kpc_in,outrad=photradius_kpc_out)
;         niceprint, photradius_kpc_in, photradius_kpc, photradius_kpc_out
             
          phot = replicate({$
            photradius_kpc:          photradius_kpc,$
            photradius_kpc_in:    photradius_kpc_in,$
            photradius_kpc_out:  photradius_kpc_out,$
            maggies:            fltarr(nphotradius),$
            ivarmaggies:        fltarr(nphotradius),$
            maggies_int_obs:                    0.0,$ ; observed integrated flux
            ivarmaggies_int_obs:                0.0,$
            maggies_int:                        0.0,$ ; integrated flux (Sersic extrapolated)
            ivarmaggies_int:                    0.0,$
            dabmag_in:                          0.0,$ ; magnitude correction from extrapolating inward
            dabmag_out:                         0.0,$ ; magnitude correction from extrapolating outward
            dabmag_int: 0.0},nfilt)                   ; AB magnitude difference between _int and _int_obs,
                                                      ;  or the sum of DABMAG_IN and DABMAG_OUT
          phot = struct_addtags(sersic,phot)
          
; now loop on each band          
          for ib = 0, nfilt-1 do begin
             band = strtrim(strupcase(modphot[ib].band),2)

             modgood = where(modphot[ib].majora*pixscale*arcsec2kpc le sersic[ib].amax_kpc and $
               modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0,nmodgood)

             sb = modphot[ib].sb0fit[modgood]*1D
             sb_var_floor = (sb*errfloor)^2.0
             sb_ivar = 1D/(1D/modphot[ib].sb0fit_ivar[modgood]+sb_var_floor)
             radius_kpc = modphot[ib].radius_kpc[modgood] ; [kpc]
             
             srt = sort(radius_kpc)
             radius_kpc = radius_kpc[srt]
             sb = sb[srt]
             sb_ivar = sb_ivar[srt]
             sb_var = 1.0/sb_ivar

; use the Sersic model to extrapolate the SB profile inward and
; outward and then integrate to get the total magnitude; our radii
; need to be in arcseconds to make the units work out when we
; integrate
             radius_arcsec = radius_kpc/arcsec2kpc
             radius_arcsec_extrap_in = [0,range(min(radius_arcsec)*1E-3,min(radius_arcsec)*0.999,20,/log)]
             radius_arcsec_extrap_out = range(max(radius_arcsec)*1.001,500.0/arcsec2kpc,30,/log)
             
             int_radius_arcsec = [radius_arcsec_extrap_in,radius_arcsec,radius_arcsec_extrap_out]

; single and 2-component Sersic fluxes
;            int_sb = [10.0^(-0.4*bcgsfhs_sersic_func(radius_arcsec_extrap_in*arcsec2kpc,params=sersic[ib])),$
;              sb,10.0^(-0.4*bcgsfhs_sersic_func(radius_arcsec_extrap_out*arcsec2kpc,params=sersic[ib]))]
             int_sb = [bcgsfhs_sersic2_func(radius_arcsec_extrap_in*arcsec2kpc,params=sersic[ib]),$
               sb,bcgsfhs_sersic2_func(radius_arcsec_extrap_out*arcsec2kpc,params=sersic[ib])]

; the code below is fancier because it takes into account the covariance in
; the Sersic parameters, but for now just extrapolate the last
; measured variance
             int_sb_var = [interpolate(sb_var,findex(radius_arcsec,radius_arcsec_extrap_in)),$
               sb_var,interpolate(sb_var,findex(radius_arcsec,radius_arcsec_extrap_out))]
             
;; get the variance of the SB profile accounting for the covariance in
;; the Sersic parameters
;             if sersic.sersic_sb0_err eq 0.0 or sersic.sersic_k_err eq 0.0 or $
;               sersic.sersic_n_err eq 0.0 then begin
;                splog, band+' - ERROR IS ZERO!'
;                int_sb_var = [$
;                  radius_arcsec_extrap_in*0.0+median(sb_var),$
;                  sb_var,$
;                  radius_arcsec_extrap_out*0.0+median(sb_var)]
;             endif else begin
;                int_sb_var = [$
;                  get_sersic_variance(radius_arcsec_extrap_in,sersic=sersic)>sb_var[0],$
;                  sb_var,$
;                  get_sersic_variance(radius_arcsec_extrap_out,sersic=sersic)>sb_var[n_elements(sb_var)-1]]
;             endelse

; get the total galaxy magnitude as well as the size of the magnitude
; correction due to the Sersic extrapolation
             phot[ib].maggies_int = 2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*int_sb)
             phot[ib].ivarmaggies_int = 1.0/(2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*int_sb_var))
             
             phot[ib].maggies_int_obs = 2.0*!pi*im_integral(radius_arcsec,radius_arcsec*sb)
             phot[ib].ivarmaggies_int_obs = 1.0/(2.0*!pi*im_integral(radius_arcsec,radius_arcsec*sb_var))

;; single Sersic             
;             phot[ib].dabmag_in = -2.5*alog10(1+2.0*!pi*im_integral(radius_arcsec_extrap_in,$
;               radius_arcsec_extrap_in*10.0^(-0.4*bcgsfhs_sersic_func($
;               radius_arcsec_extrap_in*arcsec2kpc,params=sersic[ib])))/phot[ib].maggies_int_obs)
;             phot[ib].dabmag_out = -2.5*alog10(1+2.0*!pi*im_integral(radius_arcsec_extrap_out,$
;               radius_arcsec_extrap_out*10.0^(-0.4*bcgsfhs_sersic_func($
;               radius_arcsec_extrap_out*arcsec2kpc,params=sersic[ib])))/phot[ib].maggies_int_obs)
; double Sersic             
             phot[ib].dabmag_in = -2.5*alog10(1+2.0*!pi*im_integral(radius_arcsec_extrap_in,$
               radius_arcsec_extrap_in*bcgsfhs_sersic2_func($
               radius_arcsec_extrap_in*arcsec2kpc,params=sersic[ib]))/phot[ib].maggies_int_obs)
             phot[ib].dabmag_out = -2.5*alog10(1+2.0*!pi*im_integral(radius_arcsec_extrap_out,$
               radius_arcsec_extrap_out*bcgsfhs_sersic2_func($
               radius_arcsec_extrap_out*arcsec2kpc,params=sersic[ib]))/phot[ib].maggies_int_obs)

             phot[ib].dabmag_int = -2.5*alog10(phot[ib].maggies_int/phot[ib].maggies_int_obs)

;            print, phot[ib].band, 2.0*!pi*im_integral(radius_arcsec_extrap_out,$
;              radius_arcsec_extrap_out*10.0^(-0.4*bcgsfhs_sersic_func($
;              radius_arcsec_extrap_out*arcsec2kpc,params=sersic[ib])))

; now do photometry in each aperture
             for ir = 0, nphotradius-1 do begin
                phot[ib].maggies[ir] = 2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*int_sb,$
                  phot[ib].photradius_kpc_in[ir]/arcsec2kpc,phot[ib].photradius_kpc_out[ir]/arcsec2kpc)
                phot[ib].ivarmaggies[ir] = 1.0/(2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*$
                  int_sb_var,phot[ib].photradius_kpc_in[ir]/arcsec2kpc,phot[ib].photradius_kpc_out[ir]/arcsec2kpc))
             endfor
             
;            for ir = 0, nphotradius-1 do begin
;               if max(radius_kpc) gt phot[ib].photradius_kpc[ir] then begin
;                  phot[ib].maggies[ir] = 2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*int_sb,$
;                    phot[ib].photradius_in_kpc[ir],phot[ib].photradius_out_kpc[ir])
;                  phot[ib].ivarmaggies[ir] = 1.0/(2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*$
;                    int_sb_var,phot[ib].photradius_in_kpc[ir],phot[ib].photradius_out_kpc[ir]))
;               endif else begin ; extrapolation
;                  phot[ib].maggies[ir] = 0.0
;                  phot[ib].ivarmaggies[ir] = 1.0/(10.0^(-0.4*modphot[ib].sblimit))^2
;               endelse
;            endfor
;            if total(int_sb_var eq 0) ne 0.0 then stop
             
          endfor 

; write out
          im_mwrfits, phot, sersicpath+cluster+'-phot.fits', clobber=clobber
       endfor
    endif

return
end



;; this (working) code is to use the F160W to constrain the free
;; parameters of the Sersic model
;          if ib gt 0 then begin
;             bcgsfhs_sersic, radius_kpc, sb, sersic, sb_ivar=sb_ivar, $
;               init_params=[out[0].sersic_sb0,out[0].sersic_k,$
;               out[0].sersic_n], fixed=[0,1,1]
;          endif else begin
;             bcgsfhs_sersic, radius_kpc, sb, sersic, sb_ivar=sb_ivar
;          endelse

