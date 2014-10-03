pro clash_hizsearch

    path = clash_path()+'projects/hizsearch/'
    
;   filt = clash_filterlist()
    filt = ['clash_wfc3_f110w.par','clash_wfc3_f125w.par','clash_wfc3_f160w.par']
    ii = im_filterspecs(filterlist=filt)
    
    fsps = im_read_fsps(/flam)
    ssp = im_convolve_sfh(fsps,tau=0D)
    cont = im_convolve_sfh(fsps,tau=20D)

    zform = 15.0
    zz = findgen(8)+1
    nzz = n_elements(zz)

    av = 3.0
    klam = k_lambda(fsps.wave,/calz,r_v=rv)
    alam = klam*(av/rv)

    igmgrid = mrdfits(getenv('IMPRO_DIR')+'/etc/igmtau_grid.fits.gz',1)
    
; F125W-F160W vs F110W-F125W
    out = replicate({z: 0.0, age: 0.0, ssp125_160: 0.0, ssp110_125: 0.0, $
      cont125_160: 0.0, cont110_125: 0.0, dcont125_160: 0.0, dcont110_125: 0.0},nzz)
    out.z = zz

    for jj = 0, nzz-1 do begin
       out[jj].age = getage(zz[jj])-getage(zform) ; age at z=zobs
       ageindx = findex(fsps.age/1D9,out[jj].age)
       wave = fsps.wave*(1+zz[jj])
       edgewave = k_lambda_to_edges(wave)
; deal with the IGM
       windx = findex(igmgrid.wave,wave)
       zindx = findex(igmgrid.zgrid,zz[jj])
       igm = interpolate(igmgrid.igm,windx,zindx,/grid,missing=1.0)
; SSP
       sflux = interpolate(ssp,ageindx)*igm
       maggies = k_project_filters(edgewave,sflux,filterlist=filt)
       out[jj].ssp110_125 = -2.5*alog10(maggies[0]/maggies[1])
       out[jj].ssp125_160 = -2.5*alog10(maggies[1]/maggies[2])
; continous
       cflux = interpolate(cont,ageindx)*igm
       maggies = k_project_filters(edgewave,cflux,filterlist=filt)
       out[jj].cont110_125 = -2.5*alog10(maggies[0]/maggies[1])
       out[jj].cont125_160 = -2.5*alog10(maggies[1]/maggies[2])
; dusty continous
       dcflux = cflux*10^(-0.4*alam)*igm
       maggies = k_project_filters(edgewave,dcflux,filterlist=filt)
       out[jj].dcont110_125 = -2.5*alog10(maggies[0]/maggies[1])
       out[jj].dcont125_160 = -2.5*alog10(maggies[1]/maggies[2])

;      djs_plot, [0], [0], /nodata, xrange=[0.5E4,2.5E4], yrange=[1E-6,1E2], /xlog, /ylog
;      djs_oplot, wave, sflux/interpol(sflux,wave,ii[2].weff), color='red'
;      djs_oplot, wave, cflux/interpol(cflux,wave,ii[2].weff), color='blue'
;      djs_oplot, wave, dcflux/interpol(dcflux,wave,ii[2].weff), color='green'
;      for mm = 0, n_elements(filt)-1 do begin
;         ff = ii[mm].filtf[0:ii[mm].filtn-1]
;         djs_oplot, ii[mm].filtw[0:ii[mm].filtn-1], ff/max(ff)*!y.crange[1]*0.2, $
;           line=mm
;      endfor
;      cc = get_kbrd(1)
    endfor

    psfile = path+'z10_and_allgals_models.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], width=6.8, height=5.5
    
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=[-0.5,2.0], $
      yrange=[-0.5,1.5], xtitle='F125W - F160W', ytitle='F110W - F125W', $
      position=pos
    im_legend, ['Passively Evolving','Continuous Star Formation, A_{V}=0',$
      'Continuous Star Formation, A_{V}=3 mag'], /left, /top, box=0, $
      charsize=1.4, color=['red','dodger blue','forest green'], $
      line=[0,3,5], pspacing=1.9, psym=-[16,15,17], spacing=2.0
    
    djs_oplot, out.ssp125_160, out.ssp110_125, psym=-symcat(16), $
      symsize=1.3, color=im_color('red'), line=0
    djs_oplot, out.cont125_160, out.cont110_125, psym=-symcat(15), $
      symsize=1.3, color=im_color('dodger blue'), line=3
    djs_oplot, out.dcont125_160, out.dcont110_125, psym=-symcat(17), $
      symsize=1.3, color=im_color('forest green'), line=5

    plots, out[0].ssp125_160, out[0].ssp110_125, psym=symcat(9), $
      symsize=2.0, color=im_color('red')
    plots, out[0].cont125_160, out[0].cont110_125, psym=symcat(6,thick=6), $
      symsize=2.0, color=im_color('dodger blue')
    plots, out[0].dcont125_160, out[0].dcont110_125, psym=symcat(5,thick=6), $
      symsize=2.5, color=im_color('forest green')

    xyouts, out[0].ssp125_160+0.02, out[0].ssp110_125-0.1, 'z=1', align=0.0, $
      /data, charsize=1.2, color=im_color('red')
    xyouts, out[nzz-1].ssp125_160, out[nzz-1].ssp110_125+0.05, 'z=8', align=0.0, $
      /data, charsize=1.2, color=im_color('red')

    xyouts, out[0].cont125_160+0.05, out[0].cont110_125, 'z=1', align=0.0, $
      /data, charsize=1.2, color=im_color('dodger blue')
    xyouts, out[nzz-1].cont125_160+0.05, out[nzz-1].cont110_125, 'z=8', align=0.0, $
      /data, charsize=1.2, color=im_color('dodger blue')
    
    xyouts, out[0].dcont125_160, out[0].dcont110_125-0.05, 'z=1', align=0.0, $
      /data, charsize=1.2, color=im_color('forest green')
    xyouts, out[nzz-1].dcont125_160+0.05, out[nzz-1].dcont110_125, 'z=8', align=0.0, $
      /data, charsize=1.2, color=im_color('forest green')
    
    readcol, path+'hiz_candidates.txt', $
      cl, id, ra, dec, xpos, ypos, hz110, hz125, hz160, $
      format='A,L,D,D,F,F,F,F,F', /silent
    djs_oplot, hz125-hz160, hz110-hz125, psym=symcat(9,thick=6), $
      symsize=1.3

    cc = read_santorini()
    plots, cc.f125w_mag-cc.f160w_mag, cc.f110w_mag-cc.f125w_mag, $
      psym=symcat(15), symsize=2.0
    xyouts, 1.2, 0.7, 'z=9.6 candidate', align=0.0, /data, charsize=1.5

;   im_legend, 'Tracks range from z=1-8', /left, /bottom, box=0, charsize=1.3
    
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
    
stop    
    
return
end
    
