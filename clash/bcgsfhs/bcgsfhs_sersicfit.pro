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

pro bcgsfhs_sersic2, rr, sb, scoeff, sb_err=sb_err, $
  sb_ivar=sb_ivar, init_params=init_params, fixed=fixed, $
  covar=covar, sersicfit=sersicfit
; fit a double Sersic function

    parinfo = replicate({value: 0D, fixed: 0, limited: [0,0], $
      limits: [0D,0D]},6)

; sb01 > 0 
    parinfo[0].limited = [1,0]
    parinfo[0].limits = [0D,0D]
; k1>0
    parinfo[1].limited = [1,0]
    parinfo[1].limits = [0D,0D]
; 1<n1<10
    parinfo[2].limited = [1,1]
    parinfo[2].limits = [1D,10D]
; sb02 > 0
    parinfo[3].limited = [1,0]
    parinfo[3].limits = [0D,0D]
; k2>0
    parinfo[4].limited = [1,0]
    parinfo[4].limits = [0D,0D]
; 1<n2<10
    parinfo[5].limited = [1,1]
    parinfo[5].limits = [1D,10D]

    if n_elements(init_params) eq 0 then begin
       parinfo.value = [max(sb), 4.0, 4.0, 0.1*max(sb), 4.0, 4.0]
    endif else parinfo.value = init_params

    if n_elements(fixed) ne 0 then parinfo.fixed = fixed

    good = where(sb_ivar gt 0.0 and finite(sb) eq 1, nn)
    xx = rr[good]
    xxsb = sb[good]
    xxsb_ivar = sb_ivar[good]
    
    params = mpfitfun('bcgsfhs_sersic2_func',xx,xxsb,parinfo=parinfo,$
      yfit=sersicfit,perror=perror,covar=covar,weights=xxsb_ivar,$
      status=status,quiet=1)
    
    scoeff = {$
      sersic2_sb01: max(sb), sersic2_k1: 0D, sersic2_n1: 0D, $
      sersic2_sb01_err: 0.0, sersic2_k1_err: 0.0, sersic2_n1_err: 0.0, $
      sersic2_sb02: max(sb), sersic2_k2: 0D, sersic2_n2: 0D, $
      sersic2_sb02_err: 0.0, sersic2_k2_err: 0.0, sersic2_n2_err: 0.0, $
      sersic2_covar: fltarr(6,6), sersic2_status: status, $
      sersic2_total1: 0.0, sersic2_total2: 0.0}
    
    scoeff.sersic2_sb01 = params[0]
    scoeff.sersic2_k1 = params[1]
    scoeff.sersic2_n1 = params[2]
    scoeff.sersic2_sb02 = params[3]
    scoeff.sersic2_k2 = params[4]
    scoeff.sersic2_n2 = params[5]
    
    scoeff.sersic2_sb01_err = perror[0]
    scoeff.sersic2_k1_err = perror[1]
    scoeff.sersic2_n1_err = perror[2]
    scoeff.sersic2_sb02_err = perror[3]
    scoeff.sersic2_k2_err = perror[4]
    scoeff.sersic2_n2_err = perror[5]
    
    scoeff.sersic2_covar = covar
    scoeff.sersic2_total1 = cumsersic_total(params[0:2])
    scoeff.sersic2_total2 = cumsersic_total(params[3:5])
    
return
end

pro bcgsfhs_sersic, rr, sb, scoeff, sb_err=sb_err, $
  sb_ivar=sb_ivar, init_params=init_params, fixed=fixed, $
  covar=covar, sersicfit=sersicfit
; fit a single Sersic function
    
    parinfo = replicate({value: 0D, fixed: 0, limited: [0,0], $
      limits: [0D,0D]},3)

; sb0 > max(sb0)
    parinfo[0].limited = [1,0]
    parinfo[0].limits = [0D,0D]
;   parinfo[0].limited = [1,0]
;   parinfo[0].limits = [max(sb),0D]
; k>0
    parinfo[1].limited = [1,1]
    parinfo[1].limits = [0D,1000D]
; n>0, <10
    parinfo[2].limited = [1,1]
    parinfo[2].limited = [1D,20D]

    if n_elements(init_params) eq 0 then begin
       parinfo.value = [max(sb), 10D, 4D]
;      parinfo.value = [max(sb), 4.0, 4.0]
    endif else parinfo.value = init_params 

    if n_elements(fixed) eq 3 then parinfo.fixed = fixed
    
    good = where(sb_ivar gt 0. and finite(sb) eq 1, nn)
    xx = rr[good]
    xxsb = sb[good]
    xxsb_ivar = sb_ivar[good] 

    quiet = 1
    params = mpfitfun('bcgsfhs_sersic_func',xx,xxsb,parinfo=parinfo,$
      perror=perror,covar=covar,weights=xxsb_ivar,bestnorm=chi2,dof=dof,$
      status=status,quiet=quiet,yfit=sersicfit)
;   if total(perror eq 0) ne 0.0 then stop

    factor = sqrt(chi2/dof)
    scoeff = {$
      sersic_status: status,$
      sersic_chi2:     chi2,$
      sersic_dof:       dof,$
      sersic_sb0: params[0],$
      sersic_k:   params[1],$
      sersic_n:   params[2],$
      sersic_sb0_err: perror[0],$ ; *factor,$
      sersic_k_err:   perror[1],$ ; *factor,$
      sersic_n_err:   perror[2],$ ; *factor,$
      sersic_covar:       covar,$
      sersic_total: cumsersic_total(params)}

return
end

pro bcgsfhs_sersic_wave, rr, sb, wave, scoeff, sb_err=sb_err, $
  sb_ivar=sb_ivar, init_params=init_params, fixed=fixed, $
  sersicfit=sersicfit
; fit for the wavelength dependence of the central surface brightness,
; Sersic n, and half-light radius

;     
    
    parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
      limits:[0.D,0.D]},3)

    parinfo[0].limited = [0,0]
    parinfo[0].limits = [0D,0D]
;   parinfo[0].limited = [1,0]
;   parinfo[0].limits = [max(sb),0D]
; k>0
    parinfo[1].limited = [1,1]
    parinfo[1].limits = [0D,100D]
; n>0, <10
    parinfo[2].limited = [1,1]
    parinfo[2].limited = [1D,10D]
   
    if n_elements(init_params) eq 0 then begin
       parinfo.value = [max(sb), 4.0, 4.0]
    endif else parinfo.value = init_params 

    if n_elements(fixed) ne 0 then parinfo.fixed = fixed

    good = where(sb_ivar gt 0. and finite(sb) eq 1, nn)
    xx = rr[good]
    xxsb = sb[good]
    xxsb_ivar = sb_ivar[good] 

    params = mpfitfun('bcgsfhs_sersic_wave_func',xx,xxsb,parinfo=parinfo,$
      perror=perror,covar=covar,weights=xxsb_ivar,$
      status=status,quiet=1,yfit=sersicfit)

    scoeff = {sersic_sb0: max(sb), sersic_k:0.d, sersic_n: 0D, $
      sersic_sb0_err:0., sersic_k_err:0., sersic_n_err: 0.0, $
      sersic_covar:fltarr(3,3), sersic_status: status, $
      sersic_total: 0.}

    scoeff.sersic_sb0 = params[0]
    scoeff.sersic_k = params[1]
    scoeff.sersic_n = params[2]
    scoeff.sersic_sb0_err = perror[0]
    scoeff.sersic_k_err = perror[1]
    scoeff.sersic_n_err = perror[2]
    scoeff.sersic_covar = covar
    scoeff.sersic_total = cumsersic_total(params)

return
end

pro bcgsfhs_sersicfit, debug=debug, clobber=clobber
; jm13oct22siena - fit various Sersic models to the output of
; BCGSFHS_ELLIPSE 

; read the sample
    sample = read_bcgsfhs_sample()
    ncl = n_elements(sample)

    nphotradius = 10            ; number of radial bins
;   nphotradius = 15            ; number of radial bins
    rmax_kpc = 100.0            ; [kpc]
    pixscale = 0.065D           ; [arcsec/pixel]
    errfloor = 0.05D            ; error floor on my SB measurements 
    
    ellpath = bcgsfhs_path()+'ellipse/'
    sersicpath = bcgsfhs_path()+'sersic/'

    dosersic2 = 0
    
    psfile = sersicpath+'qa_sersic.ps'

    ncol = 3
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.3
    
; wrap on each cluster    
    for ic = 4, 6 do begin
;   for ic = 4, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Working on cluster '+cluster

       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       
       modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
       reffilt = where(strtrim(modphot.band,2) eq 'f160w') ; reference filter
       nfilt = n_elements(modphot)

; read Marc's SB profiles to determine the last radius at which
; the models are reliable 
       pp = read_bcg_profiles(cluster,these_filters=strtrim(modphot.band,2))
       
; test a wavelength-dependent Sersic functional fit
;       for ib = 0, nfilt-1 do begin
;
;          modphot[ib].sb0fit_ivar gt 0 and modphot[ib].sb0fit gt $
;            10^(-0.4*modphot[ib].sblimit))

; build a QAplot
       nrow = ceil(nfilt/float(ncol))
       pos = im_getposition(nx=ncol,ny=nrow,yspace=0.0,xspace=0.0,$
         xmargin=[0.9,0.4],width=2.4)
       
       for ib = 0, nfilt-1 do begin
;      for ib = 11, 11 do begin
          band = strtrim(strupcase(modphot[ib].band),2)

          amax = max(pp[ib].sma) ; [kpc]
          modgood = where(modphot[ib].majora*pixscale*arcsec2kpc le amax and $
            modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0)
          modbad = where(modphot[ib].majora*pixscale*arcsec2kpc gt amax and $
            modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0)
          
; the equivalent radius needs to be sorted (note: the ellipse
; parameters in BCGSFHS_ELLIPSE are monotonic in semi-major axis, not
; equivalent radius!)  also add a 2% error floor to the surface
; brightnesses           
          sb = modphot[ib].sb0fit[modgood]*1D
          sb_var_floor = (sb*errfloor)^2.0
          sb_ivar = 1D/(1D/modphot[ib].sb0fit_ivar[modgood]+sb_var_floor)
          
          radius = modphot[ib].radius[modgood]*pixscale ; [arcsec]
          radius_kpc = modphot[ib].radius_kpc[modgood] ; [kpc]

          srt = sort(radius)
          radius = radius[srt]
          radius_kpc = radius_kpc[srt]
          sb = sb[srt]
          sb_ivar = sb_ivar[srt]
          sb_var = 1.0/sb_ivar ; should never be zero
          
; fit with a single-Sersic and then a double-Sersic
          bcgsfhs_sersic, radius_kpc, sb, sersic, sb_ivar=sb_ivar, $
            covar=sersic_covar
          
;; this (working) code is to use the F160W to constrain the free
;; parameters of the Sersic model
;          if ib gt 0 then begin
;             bcgsfhs_sersic, radius_kpc, sb, sersic, sb_ivar=sb_ivar, $
;               init_params=[out[0].sersic_sb0,out[0].sersic_k,$
;               out[0].sersic_n], fixed=[0,1,1]
;          endif else begin
;             bcgsfhs_sersic, radius_kpc, sb, sersic, sb_ivar=sb_ivar
;          endelse

          if dosersic2 then begin
             bcgsfhs_sersic2, radius_kpc, sb, sersic2, sb_ivar=sb_ivar, $
               covar=sersic2_covar
          endif
          
; define the radial photometry structure
          if ib eq reffilt then begin ; not general!
             photradius_kpc = get_radius(radius_kpc,nrad=nphotradius,rmax=rmax_kpc,$
               inrad=photradius_in_kpc,outrad=photradius_out_kpc)
             niceprint, photradius_kpc_in, photradius_kpc, photradius_kpc_out

;            isdata = photradius_in_kpc
             phot_template = {$
               amax:                               0.0,$ ; [kpc]
               photradius_kpc:          photradius_kpc,$
               photradius_in_kpc:    photradius_in_kpc,$
               photradius_out_kpc:  photradius_out_kpc,$
;              photisdata:         intarr(nphotradius),$ ; 1=observed photometry; 0=extrapolated photometry 
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

;         help, sersic, sersic2, /str
;         splog, band, sersic.sersic_n, sersic2.sersic2_n1, sersic2.sersic2_n2
          
; make the QAplot          
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

          print, band, sersic.sersic_k, 1.0/sersic.sersic_k
          
          rr = [0,range(0.01,200,500,/log)]
          djs_plot, radius_kpc, -2.5*alog10(sb), psym=symcat(16), /xlog, noerase=ib gt 0, $
            xrange=[0.3,200], xsty=1, yrange=[28,16], position=pos[*,ib], $
            xtickname=xtickname, ytickname=ytickname, title=title, $
            symsize=0.5, ytickinterval=3, ysty=1
          label = [$
            '\mu_{0}='+strtrim(string(-2.5*alog10(sersic.sersic_sb0),format='(F12.2)'),2),$
;           '\mu_{0}='+strtrim(string(-2.5*alog10(sersic.sersic_sb0),format='(F12.2)'),2)+' mag arcsec^{-2}',$
            'r_{e}='+strtrim(string(sersic.sersic_k,format='(F12.1)'),2)+' kpc',$
            'n='+strtrim(string(sersic.sersic_n,format='(F12.2)'),2)]
          im_legend, label, /left, /bottom, box=0, margin=0, charsize=0.8
          
          djs_oplot, modphot[ib].radius_kpc[modbad], -2.5*alog10(modphot[ib].sb0fit[modbad]), $
            psym=symcat(9), color=cgcolor('medium grey'), symsize=0.5
          
          im_legend, band, /right, /top, box=0, margin=0, charsize=1.0

          djs_oplot, rr, -2.5*alog10(bcgsfhs_sersic_func(rr,params=out[ib])), $
            color=cgcolor('firebrick')
          if ib gt 0 then djs_oplot, rr, -2.5*alog10(bcgsfhs_sersic_func(rr,$
            params=out[0])), color=cgcolor('forest green')

          if dosersic2 then begin
             djs_oplot, rr, -2.5*alog10(bcgsfhs_sersic2_func(rr,params=out[ib])), $
               color=cgcolor('red')
             djs_oplot, rr, -2.5*alog10(bcgsfhs_sersic_func(rr,[out[ib].sersic2_sb01,$
               out[ib].sersic2_k1,out[ib].sersic2_n1])), color=cgcolor('orange'), line=2
             djs_oplot, rr, -2.5*alog10(bcgsfhs_sersic_func(rr,[out[ib].sersic2_sb02,$
               out[ib].sersic2_k2,out[ib].sersic2_n2])), color=cgcolor('orange'), line=2
          endif

          djs_oplot, [10.0,10^!x.crange[1]], modphot[ib].sblimit*[1,1], $
            line=0;, color=cgcolor('grey')
          djs_oplot, amax*[1,1], !y.crange, line=0;, color=cgcolor('grey')

       endfor

       xyouts, min(pos[0,*])-0.06, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
         textoidl('\mu (mag arcsec^{-2})'), orientation=90, align=0.5, charsize=1.4, /norm
       xyouts, (max(pos[2,*])-min(pos[0,*]))/2.0+min(pos[0,*]), min(pos[1,*])-0.06, $
         textoidl('Equivalent Radius (kpc)'), align=0.5, charsize=1.4, /norm
       
; write out
       im_mwrfits, out, sersicpath+cluster+'-sersic.fits', clobber=clobber
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf

;;         bgt_ellipse_sersicradius, ellipse, outradius=outrad
;;         bgt_ellipse_radius() ; get the half-light radius

return
end
