pro mzplot_apbias, ps=ps
; jm10jul27ucsd - model the systematic effects of aperture bias

    mzpath = mz_path()
    qapath = mzpath+'qaplots/'
    if keyword_set(ps) then begin
       pspath = qapath
       suffix = '.ps'
    endif else begin
       pspath = mz_path(/paper)
       suffix = '.eps'
    endelse

    ages = read_mz_sample(/mzhii_ancillary)
    mass = read_mz_sample(/mzhii_mass)

; ---------------------------------------------------------------------------    
; mean gradient properties from SINGS
;   singshii = mrdfits(sings_path(/projects)+'log12oh/'+$
;     'sings_log12oh_hiiregions_'+sings_log12oh_version()+'.fits.gz',1)
;   sings = mrdfits(sings_path(/projects)+'log12oh/'+$
;     'sings_log12oh_'+sings_log12oh_version()+'.fits.gz',1)
;   keep = where(singshii.gradient_flag)
;   singshii = singshii[keep]
;   sings = sings[keep]
;   factor = sings.r25*!dtor/60.0*sings.distance*1E3 ; [kpc]
;   slope = singshii.hii_kk04_slope[0]/factor ; [dex/kpc]
;   slopeerr = singshii.hii_kk04_slope[1]/factor
;   slopemn = djs_mean(slope)
;   slopesig = djsig(slope)
;   slopemn = sings_weighted_mean(slope,slopeerr,wsigma=slopesig)
;   im_plothist, slope, bin=0.01
    slopemn = -0.03 
    slopesig = 0.02
    
; mean galaxy sizes (r50) and sersic indices from Blanton &
; Moustakas'09 
    r50mn = 0.35            ; =2.24 kpc for h=0.7
    r50sig = 0.25           ; =0.81 kpc for h=0.7
    sersicmn = 1.8
    sersicsig = 0.6
;   r50mn = 0.193288-alog10(0.7)            ; =2.23 kpc for h=0.7
;   r50sig = 0.251851                       ; =0.81 kpc for h=0.7
;   sersicmn = 1.86
;   sersicsig = 0.52
    
; simulation parameters
    rmin = 0.1D
    rmax = 30D
    nr = 50
    nbigr = 150
    raxis = range(rmin,rmax,nr,/log)
    bigraxis = range(1D-3,100D,nbigr,/log) ; for the integrals

    nmodel = 1000
    sim = replicate({r50: 0.0, sersic: 0.0, $
      slope: 0.0, scatter: 0.0, infiber: fltarr(nr), $
      intoh: 0.0, radoh: fltarr(nr), $
      deltaoh: fltarr(nr)},nmodel)
    sim.r50 = 10^(randomn(seed,nmodel)*r50sig+r50mn)
    sim.sersic = (randomn(seed,nmodel)*sersicsig+sersicmn)>0.5
    sim.slope = (randomn(seed,nmodel)*slopesig+slopemn)<0
;   sim.slope = randomu(seed,nmodel)*slopesig+slopemn
;   sim.scatter = 0.0
    sim.scatter = (randomn(seed,nmodel)*0.02+0.06)>0.01

;   djs_plot, [0], [0], /nodata, xrange=[0.001,rmax], $
;     yrange=[0.1,1], xsty=1, ysty=1, /ylog
    for ii = 0L, nmodel-1 do begin
; construct the Sersic model; see Graham & Driver (2005), equations
; 1&4
       bigsb = exp(-get_sersicb(sim[ii].sersic)*((bigraxis/sim[ii].r50)^(1.0/sim[ii].sersic)-1.0))
       bigsb = bigsb/max(bigsb)
       biggradient = poly(bigraxis,[9.0,sim[ii].slope])+$
         randomn(seed,nbigr)*sim[ii].scatter
;      plot, bigraxis, biggradient, psym=6, ysty=3, /xlog
;      djs_oplot, bigraxis, poly(bigraxis,[9.0,sim[ii].slope]), line=0

;      djs_oplot, bigraxis, bigsb
;      sb = exp(-get_sersicb(sim[ii].sersic)*((raxis/sim[ii].r50)^(1.0/sim[ii].sersic)-1.0))
;      gradient = poly(raxis,[9.0,sim[ii].slope])

       totlight = im_integral(bigraxis,bigraxis*bigsb)
       sim[ii].intoh = im_integral(bigraxis,bigraxis*bigsb*biggradient)/totlight ; integrated abundance
;      print, sim[ii].intoh, im_integral(raxis,sb*gradient)/im_integral(raxis,sb)
       for jj = 0, nr-1 do begin
          sim[ii].infiber[jj] = im_integral(bigraxis,bigraxis*bigsb,$
            0.0,raxis[jj])/totlight ; light fraction
          sim[ii].radoh[jj] = im_integral(bigraxis,bigraxis*bigsb*biggradient,0.0,raxis[jj])/$
            im_integral(bigraxis,bigraxis*bigsb,0D,raxis[jj])
       endfor
       sim[ii].deltaoh = sim[ii].radoh-sim[ii].intoh
;      djs_plot, raxis, sim[ii].dlogoh, xsty=3, ysty=3
    endfor

; compute the running median and make the plot
    fracbin = im_medxbin(sim.infiber,sim.deltaoh,0.1,$
      minx=0.0,maxx=1.0,minpts=10,/verbose)
    simfile = mzpath+'mzerror_apbias_simulation.fits'
    splog, 'Writing '+simfile
    mwrfits, sim, simfile, /create
    mwrfits, fracbin, simfile
    spawn, 'gzip -f '+simfile, /sh
    
;   ytitle=textoidl('\Delta<Z(r)> (dex)'), $
;   ytitle=textoidl('\Delta<Z(r)>=<Z(r)>-<Z_{int}> (dex)'), $

    im_plotconfig, 0, pos, psfile=pspath+'mzerror_apbias'+suffix, height=5.0
    hogg_scatterplot, sim.infiber, sim.deltaoh, position=pos, $
      xsty=1, ysty=1, yrange=[-0.05,0.5], xrange=[0,1], /nogrey, $
      /outliers, levels=[0.5,0.75,0.95], xtitle='F(r) (Light-Fraction)', $
      ytitle='Metallicity Bias (dex)', /internal, xnpix=51, ynpix=51, $
      outcolor=im_color('gray50',10), ccolor=im_color('gray50',10)
    oploterror, fracbin.xbin, fracbin.medy, fracbin.quant75-fracbin.medy, $
      psym=-symcat(6,thick=6), symsize=3, errthick=6, /hibar, $
      color=im_color('midnight blue',11), errcolor=im_color('midnight blue',11), thick=6
    oploterror, fracbin.xbin, fracbin.medy, fracbin.medy-fracbin.quant25, $
      psym=-symcat(6,thick=6), symsize=3, errthick=6, /lobar, $
      color=im_color('midnight blue',103), errcolor=im_color('midnight blue',103), thick=6
    im_plotconfig, /psclose

; ---------------------------------------------------------------------------    
; light fraction vs redshift in bins of stellar mass 

    xtitle = 'Redshift'
    ytitle = 'I-band Light Fraction'
    zrange = [0.03,0.75]
    fracrange = [-0.001,0.53]
    levels = [0.1,0.25,0.5,0.75,0.90]

    fracbin = im_medxbin(ages.z,ages.infiber_i,0.05,$
      minx=0.05,minpts=20,/ver)

; split by mass
    minmass = [9.0, 9.5,10.0,10.5]
    maxmass = [9.5,10.0,10.5,11.0]
    masslabel = string(minmass,format='(F4.1)')+$
      '-'+strtrim(string(maxmass,format='(F4.1)'),2)
    masslabel[0] = 'log (M/M_{\odot})='+masslabel[0]

    color1 = ['dodger blue','tan','purple','forest green']
    psym1 = [16,6,15,5]
    symsize1 = [0.6,0.3,0.6,0.3]
    nbin = n_elements(minmass)
    
    psfile = pspath+'z_vs_infiber'+suffix
    im_plotconfig, 0, pos, psfile=psfile, thick=6, $
      xmargin=[1.3,0.4], width=6.8, height=5.0

    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=zrange, yrange=fracrange, xtitle=xtitle, $
      ytitle=ytitle
    im_legend, masslabel, color=color1, psym=psym1, charsize=1.5, $
      /right, /top, box=0, spacing=1.8, margin=0

    ngal = n_elements(ages)
    bigpsym = intarr(ngal)
    bigsymsize = fltarr(ngal)
    bigcolor = strarr(ngal)
    for ii = 0, nbin-1 do begin
       indx = where((mass.mass_avg gt minmass[ii]) and $
         (mass.mass_avg lt maxmass[ii]))
       bigpsym[indx] = psym1[ii]
       bigsymsize[indx] = symsize1[ii]
       bigcolor[indx] = color1[ii]
;      djs_oplot, ages[indx].z, ages[indx].infiber_i, psym=symcat(psym1[ii],thick=8), $
;        symsize=symsize1[ii], color=im_color(color1[ii],101)
    endfor
    indx = shuffle_indx(ngal)
    for ii = 0, ngal-1 do plots, ages[indx[ii]].z, ages[indx[ii]].infiber_i, $
      psym=symcat(bigpsym[indx[ii]],thick=8), symsize=bigsymsize[indx[ii]], $
      color=im_color(bigcolor[indx[ii]]), noclip=0
    
;   djs_oplot, ages.z, ages.infiber_i, psym=symcat(16), $
;     symsize=0.4, color='grey'
;   mzages_hogg_scatterplot, ages.z, ages.infiber_i, position=pos, $
;     xsty=1, ysty=1, xrange=zrange, yrange=fracrange, xtitle=xtitle, $
;     ytitle=ytitle, levels=levels
    oploterror, fracbin.xbin, fracbin.medy, fracbin.quant75-fracbin.medy, $
      psym=-symcat(15), symsize=2.0, errthick=6, /hibar;, $
;     color=im_color('firebrick',101), errcolor=im_color('firebrick',101)
    oploterror, fracbin.xbin, fracbin.medy, fracbin.medy-fracbin.quant25, $
      psym=-symcat(15), symsize=2.0, errthick=6, /lobar;, $
;     color=im_color('firebrick',101), errcolor=im_color('firebrick',101)

    im_plotconfig, /psclose

stop
    
return
end
    
