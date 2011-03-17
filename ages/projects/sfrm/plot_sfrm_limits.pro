pro plot_sfrm_limits, ps=ps
; jm10feb05ucsd - plot the output from COMPUTE_SFRM_LIMITS

    sfrmpath = ages_path(/projects)+'sfrm/'
    paperpath = ages_path(/papers)+'sfrm/'

    parent = read_sfrm_sample()
    zbins = sfrm_zbins(nzbins)
    
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

; --------------------------------------------------
; limiting stellar mass versus redshift: all, quiescent, and
; star-forming 
    zrange = [0.0,0.8]
    massrange = [7.7,12.3]
    ztitle = 'Redshift'
    masstitle = sfrm_masstitle()
    jv2ab = k_vega2ab(filterlist='twomass_J.par',/kurucz,/silent)

    nuvmr = parent.galex_absmag[1]-parent.ugriz_absmag[2]
    rmj = parent.ugriz_absmag[2]-(parent.ubvrijhk_absmag[5]+jv2ab)
    qq = select_quiescent(nuvmr,rmj,active=aa)
    
    psfile = paperpath+'mass_limits'+suffix
    im_plotconfig, 4, pos, psfile=psfile

; all galaxies    
    djs_plot, [0], [0], /nodata, position=pos[*,0], $
      xrange=zrange, yrange=massrange, xsty=1, ysty=1, $
      xtitle='', ytitle='', xtickname=replicate(' ',10)
    mips = where(parent.phot_mips24 gt 0.0,comp=notmips)
    djs_oplot, parent[notmips].z, parent[notmips].mass, $
      psym=symcat(16), symsize=0.3, color=fsc_color('dodger blue',100)
    djs_oplot, parent[mips].z, parent[mips].mass, $
      psym=symcat(6,thick=2), symsize=0.3, color=fsc_color('firebrick',101)
    legend, 'All', /right, /bottom, box=0, margin=0, charsize=1.7

; quiescent galaxies    
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], $
      xrange=zrange, yrange=massrange, xsty=1, ysty=1, $
      xtitle='', ytitle=masstitle, xtickname=replicate(' ',10)
    mips = where(parent[qq].phot_mips24 gt 0.0,comp=notmips)
    djs_oplot, parent[qq[notmips]].z, parent[qq[notmips]].mass, $
      psym=symcat(16), symsize=0.3, color=fsc_color('dodger blue',100)
    djs_oplot, parent[qq[mips]].z, parent[qq[mips]].mass, $
      psym=symcat(6,thick=2), symsize=0.3, color=fsc_color('firebrick',101)
    legend, 'Quiescent', /right, /bottom, box=0, margin=0, charsize=1.7

; quiescent galaxies    
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,2], $
      xrange=zrange, yrange=massrange, xsty=1, ysty=1, $
      xtitle=ztitle, ytitle=''
    mips = where(parent[aa].phot_mips24 gt 0.0,comp=notmips)
    djs_oplot, parent[aa[notmips]].z, parent[aa[notmips]].mass, $
      psym=symcat(16), symsize=0.3, color=fsc_color('dodger blue',100)
    djs_oplot, parent[aa[mips]].z, parent[aa[mips]].mass, $
      psym=symcat(6,thick=2), symsize=0.3, color=fsc_color('firebrick',101)
    legend, 'Star-Forming', /right, /bottom, box=0, margin=0, charsize=1.7

    im_plotconfig, /psclose

stop    
    
; --------------------------------------------------
; limiting stellar mass based on SSP/CSP models

; desired redshift grid
    sspzform = 5.0
    tauzform = 2.0
    csfzform = 2.0
    dz = 0.01 & minz = 0.04 & maxz = 0.78
    zz = (findgen((maxz-minz)/dz+1)*dz+minz)>0.01
    nz = n_elements(zz)

; specify the models    
    modelspath = getenv('ISEDFIT_SFHGRID_DIR')+'/basemodels/bc03/'
    sspfits = 'chab_Z0.02_tau_00.0Gyr.fits.gz'
    taufits = 'chab_Z0.02_tau_03.0Gyr.fits.gz'
    csffits = 'chab_Z0.02_tau_100Gyr.fits.gz'

    ssp = mrdfits(modelspath+sspfits,1,/silent)
    tau = mrdfits(modelspath+taufits,1,/silent)
    csf = mrdfits(modelspath+csffits,1,/silent)
    wave = csf.wave

; interpolate the models at the appropriate ages, given ZFORM and the
; desired redshift grid
    sspflux = interpolate(ssp.flux,findex(ssp.age/1D9,getage(zz)-getage(sspzform)))
    tauflux = interpolate(tau.flux,findex(tau.age/1D9,getage(zz)-getage(tauzform)))
    csfflux = interpolate(csf.flux,findex(csf.age/1D9,getage(zz)-getage(csfzform)))

; now given the apparent magnitude (I=20), compute the rest-frame
; luminosity and stellar mass at each redshift
    in_filterlist = 'ndwfs_I.par'
    out_filterlist = 'sdss_r0.par'
    sdss_band_shift = 0.1
    solarmag = k_solar_magnitudes(filterlist=out_filterlist,$
      band_shift=sdss_band_shift,/silent)

    I20 = 19.95
    Ivega2ab = k_vega2ab(filterlist=in_filterlist,/kurucz,/silent)
    maggies = reform(10^(-0.4*(I20+Ivega2ab)),1,1)

    massvz = replicate({z: 0.0, ssp_Mr: 0.0, tau_Mr: 0.0, csf_Mr: 0.0, $
      ssp_mass: 0.0, tau_mass: 0.0, csf_mass: 0.0},nz)
    massvz.z = zz
    for ii = 0L, nz-1L do begin
; SSP
       kk = im_simple_kcorrect(zz[ii],maggies,maggies*0.0+1.0,$
         in_filterlist,out_filterlist,wave,sspflux[*,ii],$
         band_shift=sdss_band_shift,absmag=abs)
       massvz[ii].ssp_Mr = abs
       massvz[ii].ssp_mass = 10^(-0.4*(abs-solarmag))
; CSF
       kk = im_simple_kcorrect(zz[ii],maggies,maggies*0.0+1.0,$
         in_filterlist,out_filterlist,wave,csfflux[*,ii],$
         band_shift=sdss_band_shift,absmag=abs)
       massvz[ii].csf_Mr = abs
       massvz[ii].csf_mass = 10^(-0.4*(abs-solarmag))
; TAU
       kk = im_simple_kcorrect(zz[ii],maggies,maggies*0.0+1.0,$
         in_filterlist,out_filterlist,wave,tauflux[*,ii],$
         band_shift=sdss_band_shift,absmag=abs)
       massvz[ii].tau_Mr = abs
       massvz[ii].tau_mass = 10^(-0.4*(abs-solarmag))
    endfor

; --------------------------------------------------
; make a second plot showing M_{0.1g} *and* Mass versus redshift

    limits = mrdfits(sfrmpath+'sfrm_limits_all.fits.gz',1)

    colors = ['red','dark green','blue']
    zrange = [0.0,0.8]
    mgrange = [-15.4,-23.6]
    massrange = [7.7,12.3]
    ztitle = 'Redshift'
    mgtitle = textoidl('M_{0.1g} - 5 log (h_{70})')
    masstitle = sfrm_masstitle()

    psfile = paperpath+'mg_mass_limits'+suffix
    im_plotconfig, 6, pos, psfile=psfile, xmargin=[1.2,0.3]

    djs_plot, [0], [0], /nodata, position=pos[*,0], $
      xrange=zrange, yrange=mgrange, xsty=1, ysty=1, $
      xtitle='', ytitle=mgtitle, xtickname=replicate(' ',10)
    djs_oplot, parent.z, parent.ugriz_absmag[1], psym=symcat(16), $
      symsize=0.3, color='grey'

    notzero = where(limits.mglim_75 ne 0.0)
    djs_oplot, limits.zaxis[notzero], limits.mglim_75[notzero], $
      color='red', psym=symcat(6,thick=8)

;   djs_oplot, limits.zaxis, limits.mglim_50, color=colors[0], line=0, psym=-8
;   djs_oplot, limits.zaxis, limits.mglim_75, color=colors[1], line=0, psym=-8
;   djs_oplot, limits.zaxis, limits.mglim_95, color=colors[2], line=0, psym=-8
    djs_oplot, limits.zaxis, limits.mglim_poly, line=0, thick=6

    for ii = 0, nzbins-1 do begin
       if (ii eq 0) then djs_oplot, zbins[ii].zlo*[1,1], $
         [limits.mglim[ii],!y.crange[1]], line=1, thick=10
       djs_oplot, [zbins[ii].zlo,zbins[ii].zup], limits.mglim[ii]*[1,1], line=1, thick=10
       djs_oplot, zbins[ii].zup*[1,1], [limits.mglim[ii],!y.crange[1]], line=1, thick=10
    endfor
    
;   djs_oplot, limits.zbin, limits.mglim, psym=symcat(6,thick=8), $
;     color='orange', symsize=4.0
;   im_legend, ['50%','75%','95%'], /right, /bottom, $
;     box=0, color=colors, line=[0,0,0], pspacing=1.2

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], $
      xrange=zrange, yrange=massrange, xsty=1, ysty=1, $
      xtitle=ztitle, ytitle=masstitle
    djs_oplot, parent.z, parent.mass, psym=symcat(16), $
      symsize=0.3, color='grey'

    notzero = where(limits.minmass_75 ne 0.0)
    djs_oplot, limits.zaxis[notzero], limits.minmass_75[notzero], $
      color='red', psym=symcat(6,thick=8)

;   djs_oplot, limits.zaxis, limits.minmass_50, color=colors[0], line=0, psym=-8
;   djs_oplot, limits.zaxis, limits.minmass_75, color=colors[1], line=0, psym=-8
;   djs_oplot, limits.zaxis, limits.minmass_95, color=colors[2], line=0, psym=-8
    djs_oplot, limits.zaxis, limits.minmass_poly, line=0, thick=6

    for ii = 0, nzbins-1 do begin
       if (ii eq 0) then djs_oplot, zbins[ii].zlo*[1,1], $
         [limits.minmass[ii],!y.crange[1]], line=5, thick=5
       djs_oplot, [zbins[ii].zlo,zbins[ii].zup], limits.minmass[ii]*[1,1], line=5, thick=5
       djs_oplot, zbins[ii].zup*[1,1], [limits.minmass[ii],!y.crange[1]], line=5, thick=5
    endfor
    
;   djs_oplot, limits.zbin, limits.minmass, psym=symcat(6,thick=8), $
;     color='orange', symsize=3.0
;   im_legend, ['50%','75%','95%'], /right, /bottom, $
;     box=0, color=colors, line=[0,0,0], pspacing=1.2

    im_plotconfig, /psclose
    
; --------------------------------------------------
; limiting stellar mass versus redshift just for all galaxies
    zrange = [0.0,0.8]
    massrange = [7.7,12.3]
    ztitle = 'Redshift'
    masstitle = sfrm_masstitle()

    psfile = paperpath+'mass_limits'+suffix
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.2,0.3], height=5.0

    djs_plot, [0], [0], /nodata, position=pos, $
      xrange=zrange, yrange=massrange, xsty=1, ysty=1, $
      xtitle=ztitle, ytitle=masstitle
    mips = where(parent.phot_mips24 gt 0.0,comp=notmips)
;   djs_oplot, parent.z, parent.mass, psym=symcat(16), $
;     symsize=0.3, color='grey'
    djs_oplot, parent[notmips].z, parent[notmips].mass, $
      psym=symcat(16), symsize=0.3, color=fsc_color('dodger blue',100)
    djs_oplot, parent[mips].z, parent[mips].mass, $
      psym=symcat(6,thick=2), symsize=0.3, color=fsc_color('firebrick',101)

    notzero = where(limits.minmass_75 ne 0.0)
;   djs_oplot, limits.zaxis[notzero], limits.minmass_75[notzero], $
;     color='red', psym=symcat(6,thick=8)
;   djs_oplot, limits.zaxis, limits.minmass_poly, line=0, $
;     thick=4;, color=fsc_color('goldenrod',101)

; models    
;   djs_oplot, massvz.z, alog10(massvz.csf_mass), line=0, $
;     color=fsc_color('firebrick',101), thick=8
;   djs_oplot, massvz.z, alog10(massvz.ssp_mass), line=5, $
;     color=fsc_color('royal blue',101), thick=8

    ymax = 12.0 ; !y.crange[1]
    for ii = 0, nzbins-1 do begin
       if (ii eq 0) then djs_oplot, zbins[ii].zlo*[1,1], $
         [limits.minmass[ii],ymax], line=5, thick=4
       djs_oplot, [zbins[ii].zlo,zbins[ii].zup], limits.minmass[ii]*[1,1], line=5, thick=4
       djs_oplot, zbins[ii].zup*[1,1], [limits.minmass[ii],ymax], line=5, thick=4
    endfor
    
    im_plotconfig, /psclose

stop    
    
return
end
    
