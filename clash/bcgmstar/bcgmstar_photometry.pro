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

function get_sersic_variance, radius_kpc, sersic=sersic, covar=covar, debug=debug
; variance in the SB profile of a Sersic fit at fixed radius,
; accounting for the covariance in the paramaeters

; this code assumes that alpha,beta are both zero!    
    
    nran = 200
    rand = mrandomn(seed,covar,nran)
    
    rand[*,0] += sersic.sersic_sb0
    rand[*,1] += sersic.sersic_k
    rand[*,2] += sersic.sersic_n

    sb = fltarr(n_elements(radius_kpc),nran)

    for ii = 0, nran-1 do sb[*,ii] = bcgmstar_sersic_func(radius_kpc,$
      [rand[ii,0],rand[ii,1],rand[ii,2]])

    sb_var = radius_kpc*0.0
    for ii = 0, n_elements(radius_kpc)-1 do sb_var[ii] = stddev(sb[ii,*])^2
    
; debugging plot    
    if keyword_set(debug) then begin
       dfpsclose & im_plotfaves
       djs_plot, radius_kpc, sb[*,0], xsty=3, ysty=3, /xlog, $
         /ylog, xrange=[min(radius_kpc)>1E-4,max(radius_kpc)], yrange=minmax(sb)
       for jj = 1, nran-1 do djs_oplot, radius_kpc, sb[*,jj]
       oploterror, radius_kpc, bcgmstar_sersic_func(radius_kpc,params=sersic), $
         sqrt(sb_var), color='blue', errcolor='blue', psym=8
       cc = get_kbrd(1)

; S/N    
       djs_plot, radius_kpc, bcgmstar_sersic_func(radius_kpc,$
         params=sersic)/sqrt(sb_var), xsty=3, ysty=3
       cc = get_kbrd(1)
    endif
    
return, sb_var
end

pro bcgmstar_photometry, debug=debug, clobber=clobber
; jm14aug12siena - do radial aperture photometry in each band, using
; the Sersic models to extrapolate inward and outward 

    ellpath = bcgmstar_path(/ellipse)
    sersicpath = bcgmstar_path(/sersic)

; read the sample
    sample = read_bcgmstar_sample()
;   sample = sample[0]
    ncl = n_elements(sample)
    
    pixscale = 0.065D           ; [arcsec/pixel]
    errfloor = 0.0D ; 0.02      ; magnitude error floor on my SB measurements 
    nphotradius = 20            ; number of radial bins
    photradius_multiplier = range(0.03,3.0,nphotradius,/log) ; search for "snippet", below

    nrr = 1000
    modelrr = [0D,range(1E-4,300,nrr,/log)] ; equivalent radius [kpc]
    
    nran = 200

;; -------------------------
;; from this snippet of code, I conclude that we should measure photometry
;; in NRAD radial apertures from approximately 0.05-3 times Re in the
;; F160W band 
;    for ic = 0, ncl-1 do begin
;       cluster = strtrim(sample[ic].shortname,2)
;       jj = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent,row=0)
;       print, jj.sersic_re, jj.rmin_kpc, jj.rmax_kpc, jj.rmin_kpc/jj.sersic_re, $
;         jj.rmax_kpc/jj.sersic_re, '   '+cluster
;    endfor
;; -------------------------

; wrap on each cluster    
;   for ic = 3, 3 do begin
    for ic = 0, ncl-1 do begin
       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Measuring aperture photometry for cluster '+cluster
       
       sersic = read_bcgmstar_sersic(cluster,radius=modelrr,$
         model=modelsb,results=results,covar=covar)
       modelsb = 10D^(-0.4*modelsb) ; note!
       
       modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
       nfilt = n_elements(modphot)
       
; initialize the radial photometry structure using F160W as the
; reference band
       photradius_kpc = get_radius(photradius_multiplier*30.0,$ ; fixed radii for all clusters
         nrad=nphotradius,inrad=photradius_kpc_in,outrad=photradius_kpc_out)
;      photradius_kpc = get_radius(photradius_multiplier*sersic[0].sersic2_re2,$ ; use 2nd component!
;        nrad=nphotradius,inrad=photradius_kpc_in,outrad=photradius_kpc_out)
;      photradius_kpc = get_radius(photradius_multiplier*sersic[0].sersic_re,$
;        nrad=nphotradius,inrad=photradius_kpc_in,outrad=photradius_kpc_out)
;      niceprint, photradius_kpc_in, photradius_kpc, photradius_kpc_out
             
       phot = replicate({$
         sbe_mean:                           0.0,$ ; mean surface brightness (within re)
         sbe_mean_err:                       0.0,$
         re_mean:                            0.0,$ ; mean half-light radius (should match sersic_all_re) kpc
;        re_mean_err:                        0.0,$ ; in the clusters that were fitted with single Sersics

         photradius_kpc:          photradius_kpc,$
         photradius_kpc_in:    photradius_kpc_in,$
         photradius_kpc_out:  photradius_kpc_out,$
         maggies:            fltarr(nphotradius),$
         ivarmaggies:        fltarr(nphotradius),$
         maggies_int_obs:                    0.0,$ ; observed integrated flux
         ivarmaggies_int_obs:                0.0,$
         maggies_int:                        0.0,$ ; integrated flux (Sersic extrapolated)
         ivarmaggies_int:                    0.0,$
         maggies_30:                         0.0,$ ; flux within 30 kpc (Sersic extrapolated)
         ivarmaggies_30:                     0.0,$
         maggies_50:                         0.0,$ ; flux within 50 kpc (Sersic extrapolated)
         ivarmaggies_50:                     0.0,$
         maggies_70:                         0.0,$ ; flux within 70 kpc (Sersic extrapolated)
         ivarmaggies_70:                     0.0,$
         dabmag_in:                          0.0,$ ; magnitude correction from extrapolating inward
         dabmag_out:                         0.0,$ ; magnitude correction from extrapolating outward
         dabmag_int: 0.0},nfilt) ; AB magnitude difference between _int and _int_obs,
                                ;  or the sum of DABMAG_IN and DABMAG_OUT
       phot = struct_addtags(sersic,phot)
       
; now loop on each band          
       for ib = 0, nfilt-1 do begin
          band = strtrim(strupcase(modphot[ib].band),2)
          
          modgood = where(modphot[ib].majora*pixscale*arcsec2kpc le sersic[ib].amax_kpc and $
            modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0,nmodgood)
          
          radius_kpc = modphot[ib].radius_kpc[modgood] ; [kpc]
          sb = modphot[ib].sb0fit[modgood]
          sbvar = 1.0/modphot[ib].sb0fit_ivar[modgood]

;         sb = -2.5*alog10(modphot[ib].sb0fit[modgood])
;         sberr = 2.5/(modphot[ib].sb0fit[modgood]*sqrt(modphot[ib].sb0fit_ivar[modgood])*alog(10))
;         sbvar = sberr^2.0
          
          srt = sort(radius_kpc)
          radius_kpc = radius_kpc[srt]
          sb = sb[srt]
          sbvar = sbvar[srt]

; use the Sersic model to extrapolate the SB profile inward and
; outward and then integrate to get the total magnitude; our radii
; need to be in arcseconds to make the units work out when we
; integrate
          radius_arcsec = radius_kpc/arcsec2kpc
          radius_arcsec_extrap_in = [0,range(min(radius_arcsec)*1E-3,min(radius_arcsec)*0.999,50,/log)]
          radius_arcsec_extrap_out = range(max(radius_arcsec)*1.001,1000.0/arcsec2kpc,200,/log)
          int_radius_arcsec = [radius_arcsec_extrap_in,radius_arcsec,radius_arcsec_extrap_out]

          if cluster eq 'a209' or cluster eq 'a2261' then begin
             int_sb_in = 10D^(-0.4*bcgmstar_sersic2_func(radius_arcsec_extrap_in*$
               arcsec2kpc,params=sersic[ib],/allbands))
             int_sb_out = 10D^(-0.4*bcgmstar_sersic2_func(radius_arcsec_extrap_out*$
               arcsec2kpc,params=sersic[ib],/allbands))
          endif else begin
             int_sb_in = 10D^(-0.4*bcgmstar_sersic_func(radius_arcsec_extrap_in*$
               arcsec2kpc,params=sersic[ib],/allbands))
             int_sb_out = 10D^(-0.4*bcgmstar_sersic_func(radius_arcsec_extrap_out*$
               arcsec2kpc,params=sersic[ib],/allbands))
          endelse 
          int_sb = [int_sb_in,sb,int_sb_out]

;; -------------------------
;; test plots          
;          plot, int_radius_arcsec, int_sb, psym=4, /xlog, /ylog, xr=[0.01,100], symsize=0.4, ysty=3
;          djs_oplot, int_radius_arcsec, interpol(modelsb[*,ib],modelrr,int_radius_arcsec*arcsec2kpc), color='orange'
;
;          resid = interpol(modelsb[*,ib],modelrr,radius_kpc)/sb-1
;          djs_plot, radius_kpc, -2.5*alog10(sb), psym=8, yr=[26,16], /xlog, xrange=[0.3,120], xsty=1, ysty=1, yrange=[-1,1]
;          djs_oplot, modelrr, modelsb[*,ib], color='orange'
;          djs_plot, radius_kpc, resid, psym=8, yr=[-0.7,0.7]
;; -------------------------

;; single and 2-component Sersic fluxes
;;         int_sb = [10.0^(-0.4*bcgmstar_sersic_func(radius_arcsec_extrap_in*arcsec2kpc,params=sersic[ib])),$
;;           sb,10.0^(-0.4*bcgmstar_sersic_func(radius_arcsec_extrap_out*arcsec2kpc,params=sersic[ib]))]
;          int_sb = [bcgmstar_sersic2_func(radius_arcsec_extrap_in*arcsec2kpc,params=sersic[ib]),$
;            sb,bcgmstar_sersic2_func(radius_arcsec_extrap_out*arcsec2kpc,params=sersic[ib])]

; get the variance of the SB profile by just extrapolating the
; last/first measured variance
;         int_sbvar = [interpolate(sbvar,findex(radius_arcsec,radius_arcsec_extrap_in)),$
;           sbvar,interpolate(sbvar,findex(radius_arcsec,radius_arcsec_extrap_out))]

; get the variance of the SB profile accounting for the covariance in
; the Sersic parameters; note that the variance is generally very
; small in the outer regions!
          rand = mrandomn(seed,covar,nran)

          n = reform(rand[*,0]+sersic[ib].sersic_all_n)
          re = reform(rand[*,1]+sersic[ib].sersic_all_re)
          sbe = reform(rand[*,ib+2]+sersic[ib].sersic_all_sbe)
          
          sbmodel = fltarr(n_elements(int_radius_arcsec),nran)
          for ii = 0, nran-1 do sbmodel[*,ii] = 10D^(-0.4*bcgmstar_sersic_func($
            int_radius_arcsec*arcsec2kpc,[sbe[ii],re[ii],n[ii]]))
          sbmodel_var = stddev(sbmodel,dim=2)^2

          int_sbvar = [$
            interpolate(sbmodel_var,findex(int_radius_arcsec,radius_arcsec_extrap_in)),$
            sbvar,$
            interpolate(sbmodel_var,findex(int_radius_arcsec,radius_arcsec_extrap_out))]
          
; debugging plot    
          if keyword_set(debug) then begin
             djs_plot, int_radius_arcsec*arcsec2kpc, int_sb, xr=[0.1,100], $
               xsty=1, ysty=1, /xlog, /ylog
             djs_oplot, int_radius_arcsec*arcsec2kpc, int_sb+sqrt(int_sbvar), line=5
             djs_oplot, int_radius_arcsec*arcsec2kpc, int_sb-sqrt(int_sbvar), line=5
             djs_oplot, radius_kpc, sb, psym=symcat(16), color='blue'
             djs_oplot, modelrr, modelsb[*,ib], color='orange'
             for ii = 0, nran-1 do djs_oplot, int_radius_arcsec*arcsec2kpc, sbmodel[*,ii], color='red'
             cc = get_kbrd(1)
          endif
    
; get the total galaxy magnitude in multiple physical apertures, as
; well as the size of the magnitude correction due to the Sersic
; extrapolation
          phot[ib].maggies_int = 2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*int_sb)
          phot[ib].ivarmaggies_int = 1.0/(2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*int_sbvar))
             
          phot[ib].maggies_int_obs = 2.0*!pi*im_integral(radius_arcsec,radius_arcsec*sb)
          phot[ib].ivarmaggies_int_obs = 1.0/(2.0*!pi*im_integral(radius_arcsec,radius_arcsec*sbvar))

; within 30, 50, 70 kpc
          phot[ib].maggies_30 = 2.0*!pi*im_integral(int_radius_arcsec,$
            int_radius_arcsec*int_sb,0D,30D/arcsec2kpc)
          phot[ib].maggies_50 = 2.0*!pi*im_integral(int_radius_arcsec,$
            int_radius_arcsec*int_sb,0D,50D/arcsec2kpc)
          phot[ib].maggies_70 = 2.0*!pi*im_integral(int_radius_arcsec,$
            int_radius_arcsec*int_sb,0D,70D/arcsec2kpc)

          phot[ib].ivarmaggies_30 = 1.0/(2.0*!pi*im_integral(int_radius_arcsec,$
            int_radius_arcsec*int_sbvar,0D,30D/arcsec2kpc))
          phot[ib].ivarmaggies_50 = 1.0/(2.0*!pi*im_integral(int_radius_arcsec,$
            int_radius_arcsec*int_sbvar,0D,50D/arcsec2kpc))
          phot[ib].ivarmaggies_70 = 1.0/(2.0*!pi*im_integral(int_radius_arcsec,$
            int_radius_arcsec*int_sbvar,0D,70D/arcsec2kpc))

; calculate some magnitude differences of convenience
          phot[ib].dabmag_in = -2.5*alog10(1+2.0*!pi*im_integral(radius_arcsec_extrap_in,$
            radius_arcsec_extrap_in*int_sb_in)/phot[ib].maggies_int_obs)

          phot[ib].dabmag_out = -2.5*alog10(1+2.0*!pi*im_integral(radius_arcsec_extrap_out,$
            radius_arcsec_extrap_out*int_sb_out)/phot[ib].maggies_int_obs)
          phot[ib].dabmag_int = -2.5*alog10(phot[ib].maggies_int/phot[ib].maggies_int_obs)

; now do photometry in each aperture
          for ir = 0, nphotradius-1 do begin
             phot[ib].maggies[ir] = 2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*int_sb,$
               phot[ib].photradius_kpc_in[ir]/arcsec2kpc,phot[ib].photradius_kpc_out[ir]/arcsec2kpc)
             phot[ib].ivarmaggies[ir] = 1.0/(2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*$
               int_sbvar,phot[ib].photradius_kpc_in[ir]/arcsec2kpc,phot[ib].photradius_kpc_out[ir]/arcsec2kpc))
          endfor
             
;         for ir = 0, nphotradius-1 do begin
;            if max(radius_kpc) gt phot[ib].photradius_kpc[ir] then begin
;               phot[ib].maggies[ir] = 2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*int_sb,$
;                 phot[ib].photradius_in_kpc[ir],phot[ib].photradius_out_kpc[ir])
;               phot[ib].ivarmaggies[ir] = 1.0/(2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*$
;                 int_sbvar,phot[ib].photradius_in_kpc[ir],phot[ib].photradius_out_kpc[ir]))
;            endif else begin ; extrapolation
;               phot[ib].maggies[ir] = 0.0
;               phot[ib].ivarmaggies[ir] = 1.0/(10.0^(-0.4*modphot[ib].sblimit))^2
;            endelse
;         endfor
;         if total(int_sbvar eq 0) ne 0.0 then stop
          
; get the mean surface brightness within Re (see Binney & Merrifield,
; pg 186 for some equations which I used to check that I'm
; going the calculation correctly); get Re numerically because some
; clusters were fitted with two Sersic functions
          light = [0,fltarr(nrr)]
          for ir = 1, nrr do light[ir] = 2.0*!pi*im_integral(int_radius_arcsec,$
            int_radius_arcsec*int_sb,0D,modelrr[ir]/arcsec2kpc)
          phot[ib].re_mean = interpol(modelrr,light/phot[ib].maggies_int,0.5)

          phot[ib].sbe_mean = 2.0*!pi*im_integral(int_radius_arcsec,int_radius_arcsec*int_sb,$
            0D,phot[ib].re_mean/arcsec2kpc)/(!pi*(phot[ib].re_mean/arcsec2kpc)^2)
       endfor 

; write out
       im_mwrfits, phot, sersicpath+cluster+'-phot.fits', clobber=clobber
    endfor

return
end
    
