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

function get_radius, radius_sb, nrad=nrad, rmax=rmax, $
  inrad=inrad, outrad=outrad
; everything in kpc
    if n_elements(nrad) eq 0 then nrad = 20
    nrad1 = nrad-1
    inrad = range(min(radius_sb),rmax,nrad1,/asinh) ; inner radius [kpc]
    outrad = inrad+(shift(inrad,-1)-inrad)           ; outer radius [kpc]

; adjust the edges
    outrad[nrad1-1] = (inrad[nrad1-1]-inrad[nrad1-2])+inrad[nrad1-1]
    outrad = [inrad[0],outrad]
    inrad = [0D,inrad]
    radius = (outrad-inrad)/2.0+inrad
    
;   niceprint, inrad, radius, outrad
return, radius
end

pro bcgsfhs_sersic2, rr, sb, scoeff, sb_ivar=sb_ivar, $
  init_params=init_params, fixed=fixed, sersicfit=sersicfit
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
    parinfo[0].limits = [1D-15,0D]
; re1>0
    parinfo[1].limited = [1,1]
    parinfo[1].limits = [1D-3,500D]
; 1<n1<10
    parinfo[2].limited = [1,1]
    parinfo[2].limits = [1D,10D]

; sbe2 > 0
    parinfo[3].limited = [1,0]
    parinfo[3].limits = [1D-15,0D]
; re2>0
    parinfo[4].limited = [1,1]
    parinfo[4].limits = [1D-3,50D]
; 1<n2<10
    parinfo[5].limited = [1,1]
    parinfo[5].limits = [1D,10D]

    if n_elements(init_params) eq nparam then $
      parinfo.value = init_params else $
      parinfo.value = [median(xxsb),1D,1D,0.1*median(xxsb),4D,10D]
    if n_elements(fixed) eq nparam then parinfo.fixed = fixed

    if n_elements(fixed) ne 0 then parinfo.fixed = fixed

    quiet = 1
    params = mpfitfun('bcgsfhs_sersic2_func',xx,xxsb,parinfo=parinfo,$
      yfit=sersicfit,perror=perror,covar=covar,weights=xxsb_ivar,$
      dof=dof,bestnorm=chi2,status=status,quiet=quiet,$
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
  init_params=init_params, fixed=fixed, sersicfit=sersicfit
; fit a single Sersic function

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
    parinfo[1].limits = [1D-3,500D]
; 0<n<10
    parinfo[2].limited = [1,1]
    parinfo[2].limits = [1D,10D]

    if n_elements(init_params) eq nparam then $
      parinfo.value = init_params else $
      parinfo.value = [median(xxsb), 10D, 4D]
    if n_elements(fixed) eq nparam then parinfo.fixed = fixed
    
    quiet = 1
    params = mpfitfun('bcgsfhs_sersic_func',xx,xxsb,parinfo=parinfo,$
      perror=perror,covar=covar,weights=xxsb_ivar,bestnorm=chi2,dof=dof,$
      status=status,quiet=quiet,yfit=sersicfit,functargs={parinfo: parinfo})

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
  init_params=init_params, fixed=fixed, sersicfit=sersicfit
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
    parinfo[0].limits = [1D-3,500D]
; 1<n1<10
    parinfo[1].limited = [1,1]
    parinfo[1].limits = [0.1D,10D]
; re2>0
    parinfo[2].limited = [1,1]
    parinfo[2].limits = [1D-3,500D]
; 1<n2<10
    parinfo[3].limited = [1,1]
    parinfo[3].limits = [0.1D,10D]

; sbe1 - surface brightness at re
    for ib = 0, nband-1 do begin
       parinfo[4+ib].limited = [1,0]
       parinfo[4+ib].limits = [0D,0D]
    endfor

; sbe2 - surface brightness at re
    for ib = 0, nband-1 do begin
       parinfo[4+nband+ib].limited = [1,0]
       parinfo[4+nband+ib].limits = [0D,0D]
    endfor

    if n_elements(fixed) eq nparam then parinfo.fixed = fixed
    if n_elements(init_params) eq nparam then $
      parinfo.value = init_params else $
        parinfo.value = [0.1D,1.5D,30D,4D,0.1*replicate(median(xxsb),nband),$
      replicate(median(xxsb),nband)]

    quiet = 0
    params = mpfitfun('bcgsfhs_sersic2_multiband_func',xx,xxsb,parinfo=parinfo,$
      yfit=sersicfit,perror=perror,covar=covar,weights=xxsb_ivar,dof=dof,$
      bestnorm=chi2,status=status,quiet=quiet,functargs={parinfo: parinfo, wave: wave})

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
  init_params=init_params, fixed=fixed, sersicfit=sersicfit
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
    parinfo[0].limits = [1D-3,300D]
; 0<n<10
    parinfo[1].limited = [1,1]
    parinfo[1].limits = [1D,10D]
; sbe - surface brightness at re
    for ib = 0, nband-1 do begin
       parinfo[2+ib].limited = [1,1]
       parinfo[2+ib].limits = [10D,35D]
    endfor
    
    if n_elements(init_params) eq nparam then $
      parinfo.value = init_params else $
        parinfo.value = [10D, 4D, replicate(median(xxsb),nband)]
    if n_elements(fixed) eq nparam then parinfo.fixed = fixed
    
    quiet = 1
    params = mpfitfun('bcgsfhs_sersic_multiband_func',xx,xxsb,parinfo=parinfo,$
      perror=perror,covar=covar,weights=xxsb_ivar,bestnorm=chi2,dof=dof,$
      status=status,quiet=quiet,yfit=sersicfit,functargs={parinfo: parinfo, wave: wave})

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

pro bcgsfhs_sersicfit, dofit=dofit, dophot=dophot, qaplot=qaplot, clobber=clobber
; jm13oct22siena - fit various Sersic models to the output of
; BCGSFHS_ELLIPSE 

    ellpath = bcgsfhs_path()+'ellipse/'
    sersicpath = bcgsfhs_path()+'sersic/'

; read the sample
    sample = read_bcgsfhs_sample()
    ncl = n_elements(sample)

    nphotradius = 10            ; number of radial bins
;   nphotradius = 15            ; number of radial bins
    rmax_kpc = 100.0            ; [kpc]
    pixscale = 0.065D           ; [arcsec/pixel]
    errfloor = 0.0D             ; error floor on my SB measurements 
;   errfloor = 0.02D            ; error floor on my SB measurements 

; build the radius vector
    
    

; do 2-component Sersic fitting?    
    dosersic2 = 1

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
          pp = read_bcg_profiles(cluster,these_filters=strtrim(modphot.band,2))
          nfilt = n_elements(modphot)
       
; fit the red bands independent, but the bluest bands simultaneously
; because of the low S/N
          if cluster eq 'macs2129' or cluster eq 'macs0744' then wcut = 6000.0 else wcut = 5000.0
          blue = where(modphot.weff lt wcut,nblue,comp=red,ncomp=nred)
          for ib = 0, nred-1 do begin
             amax = max(pp[red[ib]].sma,mxindx) ; [kpc]
             if pp[red[ib]].mu[mxindx] gt modphot[red[ib]].sblimit then $
               amax = pp[red[ib]].sma[mxindx-1]

             modgood = where(modphot[red[ib]].majora*pixscale*arcsec2kpc le amax and $
               modphot[red[ib]].sb0fit gt 0 and modphot[red[ib]].sb0fit_ivar gt 0,nmodgood)
          
; the equivalent radius needs to be sorted (note: the ellipse
; parameters in BCGSFHS_ELLIPSE are monotonic in semi-major axis, not
; equivalent radius!)  also add a minimum error floor to the surface
; brightnesses
             sb = modphot[red[ib]].sb0fit[modgood]*1D
             sb_var_floor = (sb*errfloor)^2.0
             sb_ivar = 1D/(1D/modphot[red[ib]].sb0fit_ivar[modgood]+sb_var_floor)
             radius_kpc = modphot[red[ib]].radius_kpc[modgood] ; [kpc]
             
             srt = sort(radius_kpc)
             radius_kpc = radius_kpc[srt]
             sb = sb[srt]
             sb_ivar = sb_ivar[srt]
             
; fit with a single-Sersic and then a double-Sersic
             bcgsfhs_sersic, radius_kpc, sb, sersic, sb_ivar=sb_ivar
             if dosersic2 then bcgsfhs_sersic2, radius_kpc, sb, sersic2, sb_ivar=sb_ivar

; pack into a structure
             if ib eq 0 then begin
                out = struct_addtags(struct_trimtags(modphot,select=[$
                  'file','band','weff','sblimit','ra','dec','mge_*']),$
                  im_empty_structure(sersic,ncopies=nfilt))
                out = struct_addtags(out,replicate({sersic_covar_multiband: $
                  fltarr(2+nblue,2+nblue)},nfilt))
                if dosersic2 then begin
                   out = struct_addtags(out,im_empty_structure(sersic2,ncopies=nfilt))
                   out = struct_addtags(out,replicate({sersic2_covar_multiband: $
                     fltarr(2*nblue+4,2*nblue+4)},nfilt))
                endif
             endif 
             
             out[red[ib]] = im_struct_assign(sersic,out[red[ib]],/nozero)
             if dosersic2 then out[red[ib]] = im_struct_assign(sersic2,out[red[ib]],/nozero)
          endfor 
          
; now fit the blue bands simultaneously 
          for ib = 0, nblue-1 do begin
             amax = max(pp[blue[ib]].sma,mxindx) ; [kpc]
             if pp[blue[ib]].mu[mxindx] gt modphot[blue[ib]].sblimit then amax = pp[blue[ib]].sma[mxindx-1]

             modgood = where(modphot[blue[ib]].majora*pixscale*arcsec2kpc le amax and $
               modphot[blue[ib]].sb0fit gt 0 and modphot[blue[ib]].sb0fit_ivar gt 0,nmodgood)
             
             sb = modphot[blue[ib]].sb0fit[modgood]*1D
             sb_var_floor = (sb*errfloor)^2.0
             sb_ivar = 1D/(1D/modphot[blue[ib]].sb0fit_ivar[modgood]+sb_var_floor)
             radius_kpc = modphot[blue[ib]].radius_kpc[modgood] ; [kpc]
             
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
            multisersic, sb_ivar=fit_sb_ivar
          for ib = 0, nblue-1 do begin
             out[blue[ib]] = im_struct_assign(multisersic[ib],out[blue[ib]],/nozero)
             out[blue[ib]].sersic_covar = 0
             out[blue[ib]].sersic_covar_multiband = multisersic[ib].sersic_covar
          endfor
          niceprint, modphot.band, out.sersic_n, out.sersic_n_err, $
            out.sersic_re, out.sersic_re_err, out.sersic_sbe, out.sersic_sbe_err
          
          if dosersic2 then begin
             bcgsfhs_sersic2_multiband, fit_radius_kpc, fit_sb, fit_wave, $
               multisersic2, sb_ivar=fit_sb_ivar, sersicfit=sersicfit
             for ib = 0, nblue-1 do begin
                out[blue[ib]] = im_struct_assign(multisersic2[ib],out[blue[ib]],/nozero)
                out[blue[ib]].sersic2_covar = 0
                out[blue[ib]].sersic2_covar_multiband = multisersic2[ib].sersic2_covar
             endfor
          endif
; write out
          im_mwrfits, out, sersicpath+cluster+'-sersic.fits', clobber=clobber
       endfor
    endif

    
; ##################################################
; build a QAplot
    if keyword_set(qaplot) then begin
       ncol = 3
       rr = [0,range(0.01,200,500,/log)]
       
; -------------------------
; second QAplot: color-radius plots, relative to F160W
       psfile = sersicpath+'qa_color_sersic.ps'
       im_plotconfig, 0, pos, psfile=psfile, charsize=1.1
;      for ic = 0, 0 do begin
       for ic = 0, ncl-1 do begin
          cluster = strtrim(sample[ic].shortname,2)
          splog, cluster

          arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
          
          sersic = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent)
          modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
          pp = read_bcg_profiles(cluster,these_filters=strtrim(modphot.band,2))
          nfilt = n_elements(modphot)-1 ; relative to F160W
          
          nrow = ceil(nfilt/float(ncol))
          pos = im_getposition(nx=ncol,ny=nrow,yspace=0.1,xspace=0.8*[1,1],$
            xmargin=[1.0,0.2],width=1.9)
          for ib = 1, nfilt-1 do begin
             band = strtrim(strupcase(modphot[ib].band),2)

             amax = max(pp[ib].sma,mxindx)<max(pp[0].sma,mxindx) ; [kpc]
             if pp[ib].mu[mxindx] gt modphot[ib].sblimit then amax = pp[ib].sma[mxindx-1]

             modgood = where(modphot[ib].majora*pixscale*arcsec2kpc le amax and $
               modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0 and $
               modphot[0].majora*pixscale*arcsec2kpc le amax and $
               modphot[0].sb0fit gt 0 and modphot[0].sb0fit_ivar gt 0)
             modbad = where((modphot[ib].majora*pixscale*arcsec2kpc gt amax or $
               modphot[0].majora*pixscale*arcsec2kpc gt amax) and $
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

stop       

; -------------------------
; first QAplot: SB profiles and the Sersic fits
       psfile = sersicpath+'qa_sersic.ps'
       im_plotconfig, 0, pos, psfile=psfile, charsize=1.3
;      for ic = 2, 2 do begin
       for ic = 0, ncl-1 do begin
          cluster = strtrim(sample[ic].shortname,2)
          splog, cluster

          arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]

          sersic = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent)
          modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
          pp = read_bcg_profiles(cluster,these_filters=strtrim(modphot.band,2))
          nfilt = n_elements(modphot)

          nrow = ceil(nfilt/float(ncol))
          pos = im_getposition(nx=ncol,ny=nrow,yspace=0.0,xspace=0.0,$
            xmargin=[0.9,0.4],width=2.4)
          for ib = 0, nfilt-1 do begin
;         for ib = 0, nfilt-1 do begin
             band = strtrim(strupcase(modphot[ib].band),2)

             amax = max(pp[ib].sma,mxindx) ; [kpc]
             if pp[ib].mu[mxindx] gt modphot[ib].sblimit then amax = pp[ib].sma[mxindx-1]

             modgood = where(modphot[ib].majora*pixscale*arcsec2kpc le amax and $
               modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0)
             modbad = where(modphot[ib].majora*pixscale*arcsec2kpc gt amax and $
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
               xrange=[0.3,200], xsty=1, yrange=[28,16], position=pos[*,ib], $
               xtickname=xtickname, ytickname=ytickname, title=title, $
               symsize=0.5, ytickinterval=3, ysty=1

             label = [$
               '\mu_{e}='+strtrim(string(sersic[ib].sersic_sbe,format='(F12.1)'),2),$
               'r_{e}='+strtrim(string(sersic[ib].sersic_re,format='(F12.1)'),2)+' kpc',$
;              'r_{e}='+strtrim(string(sersic[ib].sersic_re,format='(F12.1)'),2)+' kpc',$
               'n='+strtrim(string(sersic[ib].sersic_n,format='(F12.2)'),2),$
               '\chi^{2}_{\nu, single}='+strtrim(string(sersic[ib].sersic_chi2/$
               sersic[ib].sersic_dof,format='(F12.2)'),2)]
             if dosersic2 then begin
                if sersic[ib].sersic2_sbe1 eq 0.0 or sersic[ib].sersic2_sbe2 eq 0.0 then begin
                   label = [label,'Sersic-2 dropped']
                endif else begin
                   label = [label,'\chi^{2}_{\nu, double}='+strtrim(string(sersic[ib].sersic2_chi2/$
                     sersic[ib].sersic2_dof,format='(F12.2)'),2)]
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
                endelse
             endif

             djs_oplot, [10.0,10^!x.crange[1]], modphot[ib].sblimit*[1,1], $
               line=0                              ;, color=cgcolor('grey')
;            djs_oplot, amax*[1,1], !y.crange, line=0 ;, color=cgcolor('grey')
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
; wrap on each cluster    
       for ic = 0, 0 do begin
;      for ic = 0, ncl-1 do begin
          cluster = strtrim(sample[ic].shortname,2)
          splog, cluster

          sersic = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent)
          modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
          reffilt = where(strtrim(modphot.band,2) eq 'f160w')
          nfilt = n_elements(modphot)

          for ib = 0, nfilt-1 do begin
             band = strtrim(strupcase(modphot[ib].band),2)

; define the radial photometry structure
             if ib eq reffilt then begin ; not general!
                photradius_kpc = get_radius(radius_kpc,nrad=nphotradius,rmax=rmax_kpc,$
                  inrad=photradius_in_kpc,outrad=photradius_out_kpc)
                niceprint, photradius_kpc_in, photradius_kpc, photradius_kpc_out
                
                phot_template = {$
                  amax:                               0.0,$ ; [kpc]
                  photradius_kpc:          photradius_kpc,$
                  photradius_in_kpc:    photradius_in_kpc,$
                  photradius_out_kpc:  photradius_out_kpc,$
;                 photisdata:         intarr(nphotradius),$ ; 1=observed photometry; 0=extrapolated photometry 
                  maggies:            fltarr(nphotradius),$
                  ivarmaggies:        fltarr(nphotradius),$
                  maggies_int_obs:                    0.0,$ ; observed integrated flux
                  ivarmaggies_int_obs:                0.0,$
                  maggies_int:                        0.0,$ ; integrated flux (Sersic extrapolated)
                  ivarmaggies_int:                    0.0,$
                  dabmag_int:         0.0} ; AB magnitude difference between _int and _int_obs
             endif 
             phot1 = phot_template
             phot1.amax = amax
stop             
; use the Sersic model to extrapolate the SB profile inward and outward
             radius_arcsec = radius_kpc/arcsec2kpc ; the data [arcsec]
             radius_arcsec_extrap_in = [0,range(min(radius_arcsec)*1E-3,min(radius_arcsec)*0.95,20,/log)]
             radius_arcsec_extrap_out = range(max(radius_arcsec)*1.05,500.0,50,/log)
             
             int_radius_arcsec = [radius_arcsec_extrap_in,radius_arcsec,radius_arcsec_extrap_out]
             int_sb = [bcgsfhs_sersic_func(radius_arcsec_extrap_in*arcsec2kpc,params=sersic),sb,$
               bcgsfhs_sersic_func(radius_arcsec_extrap_out*arcsec2kpc,params=sersic)]
             
; get the variance of the SB profile accounting for the covariance in
; the Sersic parameters
             if sersic.sersic_sb0_err eq 0.0 or sersic.sersic_k_err eq 0.0 or $
               sersic.sersic_n_err eq 0.0 then begin
                splog, band+' - ERROR IS ZERO!'
                int_sb_var = [$
                  radius_arcsec_extrap_in*0.0+median(sb_var),$
                  sb_var,$
                  radius_arcsec_extrap_out*0.0+median(sb_var)]
             endif else begin
                int_sb_var = [$
                  get_sersic_variance(radius_arcsec_extrap_in,sersic=sersic)>sb_var[0],$
                  sb_var,$
                  get_sersic_variance(radius_arcsec_extrap_out,sersic=sersic)>sb_var[n_elements(sb_var)-1]]
             endelse
             
             if total(int_sb_var eq 0) ne 0.0 then stop
          
; do the photometry           
             phot1.maggies_int = 2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*int_sb)
             phot1.ivarmaggies_int = 1.0/(2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*int_sb_var))
             
             phot1.maggies_int_obs = 2.0*!pi*im_integral(radius_arcsec,radius_arcsec*sb)
             phot1.ivarmaggies_int_obs = 1.0/(2.0*!pi*im_integral(radius_arcsec,radius_arcsec*sb_var))
             
             phot1.dabmag_int = -2.5*alog10(phot1.maggies_int_obs/phot1.maggies_int)
             
             for ir = 0, nphotradius-1 do begin
                if max(radius_kpc) gt phot1.photradius_kpc[ir] then begin
                   phot1.maggies[ir] = 2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*int_sb,$
                     phot1.photradius_in_kpc[ir],phot1.photradius_out_kpc[ir])
                   phot1.ivarmaggies[ir] = 1.0/(2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*$
                     int_sb_var,phot1.photradius_in_kpc[ir],phot1.photradius_out_kpc[ir]))
                endif else begin ; extrapolation
                   phot1.maggies[ir] = 0.0
                   phot1.ivarmaggies[ir] = 1.0/(10.0^(-0.4*modphot[ib].sblimit))^2
                endelse
             endfor
             
; pack everything into a structure
             out1 = struct_addtags(struct_trimtags(modphot[ib],select=[$
               'file','band','weff','sblimit','ra','dec','mge_*']),phot1)
             if dosersic2 eq 0 then out1 = struct_addtags(out1,sersic) else $
               out1 = struct_addtags(struct_addtags(out1,sersic),sersic2)
             if ib eq 0 then out = out1 else out = [out,out1]
          endfor
          
; write out
          im_mwrfits, out, sersicpath+cluster+'-phot.fits', clobber=clobber
       endfor 
    endif

stop    

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

