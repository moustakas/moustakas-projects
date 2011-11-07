pro mzplot_literature, ps=ps
; jm11oct25ucsd - compare with the literature

    mzpath = mz_path()
    qapath = mzpath+'qaplots/'
    if keyword_set(ps) then begin
       pspath = qapath
       suffix = '.ps'
    endif else begin
       pspath = mz_path(/paper)
       suffix = '.eps'
    endelse

; read the super-important output of FIT_MZLZEVOL()    
    mzavg = mrdfits(mzpath+'mzevol_avg.fits.gz',1,/silent)
    zbins = mz_zbins(nz)
    massbins = mz_massbins(nmassbins)

; --------------------------------------------------
; MZ relation in two redshift bins compared to data from the
; literature 
    calib = 't04'
;   calib = 'kk04'
    mzlocal = mrdfits(mzpath+'mzlocal_sdss_ews_'+calib+'.fits.gz',1)
    mzevol = mrdfits(mzpath+'mzevol_'+calib+'.fits.gz',1,/silent)

    xrange = [9.0,11.5]
    yrange = [8.35,9.25]
    xtitle1 = mzplot_masstitle()
    ytitle1 = mzplot_ohtitle()+'_{T04}'
;   ytitle1 = textoidl('\Delta'+'<12 + log (O/H)>')

    maxis1 = range(9.5,xrange[1]-0.1,100)
    mzref = mzevol_func(maxis1,[mzlocal.coeff,mzavg.r0[1],mzavg.p0[1]],z=0.1,qz0=mzavg.qz0)
    
    psfile = pspath+'mzevol_literature'+suffix
    im_plotconfig, 0, pos, psfile=psfile, height=5.2;, height=4.7, width=6.5, xmargin=[1.6,0.4]
;   im_plotconfig, 6, pos, psfile=psfile;, height=4.7, width=6.5, xmargin=[1.6,0.4]

;; -------------------------
;; z=0.3-0.5
;;   yrange = [-0.3,+0.05]
;    zval = 0.4
;
;    djs_plot, [0], [0], /nodata, ysty=1, xsty=1, position=pos[*,0], $
;      xrange=xrange, yrange=yrange, xtitle='', ytitle=ytitle1, $
;      xtickname=replicate(' ',10)
;    djs_oplot, !x.crange, [0,0], line=0, color=im_color('grey')
;
;    good = where(mzevol.ohmean_bymass[*,3] gt -900.0) ; z=0.4
;    oploterror, mzevol.medmass_bymass[good,3], mzevol.ohmean_bymass[good,3], $
;      mzevol.ohmean_bymass_err[good,3], psym=6, symsize=2 ;, $
;;     psym=symcat(calibpsym[ii],thick=5), symsize=calibsymsize[ii], $
;;     color=im_color(calibcolor[ii],101), errcolor=im_color(calibcolor[ii],101)
;    djs_oplot, maxis1, mzref, line=0
;    djs_oplot, maxis1, mzevol_func(maxis1,[mzlocal.coeff,mzavg.r0[1],mzavg.p0[1]],$ ; =T04
;      z=0.3,qz0=mzavg.qz0), line=5
;    djs_oplot, maxis1, mzevol_func(maxis1,[mzlocal.coeff,mzavg.r0[1],mzavg.p0[1]],$ ; =T04
;      z=0.5,qz0=mzavg.qz0), line=5

; -------------------------
; z=0.7-0.9
;   yrange = [-0.5,+0.05]
    djs_plot, [0], [0], /nodata, ysty=1, xsty=1, position=pos, $
      xrange=xrange, yrange=yrange, xtitle=xtitle1, ytitle=ytitle1
    polyfill, [maxis1,reverse(maxis1)],[mzref-mzlocal.scatter/2,reverse(mzref+mzlocal.scatter/2)], $
      /data, color=im_color('light gray'), noclip=0, /fill
;   djs_oplot, maxis1, mzref, line=0
    good = where(mzevol.ohmean_bymass[*,5] gt -900.0) ; z=0.65
    oploterror, mzevol.medmass_bymass[good,5], mzevol.ohmean_bymass[good,5], $
      mzevol.ohmean_bymass_err[good,5], psym=symcat(15), symsize=2.0

    for ii = nz-2, nz-1 do begin
       good = where(mzevol.ohmean_bymass[*,ii] gt -900.0) ; z=0.65
       oploterror, mzevol.medmass_bymass[good,ii], mzevol.ohmean_bymass[good,ii], $
         mzevol.ohmean_bymass_err[good,ii], psym=symcat(6), symsize=2.0
    endfor

    djs_oplot, maxis1, mzevol_func(maxis1,[mzlocal.coeff,mzavg.r0[1],mzavg.p0[1]],$ ; =T04
      z=0.7,qz0=mzavg.qz0), line=5
    djs_oplot, maxis1, mzevol_func(maxis1,[mzlocal.coeff,mzavg.r0[1],mzavg.p0[1]],$ ; =T04
      z=0.9,qz0=mzavg.qz0), line=5

; ##########    
; savaglio+05: h=
    logt = alog10(getage(0.1))
    logm = 10.0
;   print, -7.5903+2.5315*logm-0.09649*logm^2+5.1733*logt-0.3944*logt^2-0.4030*logt*logm

; ##########    
; cowie & barger 08; h=0.7; T04 metallicities; Salpeter IMF

; from their Fig 63 derived using dexter
    data = [$ ; z=0.475-0.9 (median z=0.75)
      [10.148,8.7487,0.04,0.03],$
      [10.451,8.8196,0.05,0.05],$
      [10.747,8.9431,0.04,0.08],$
      [11.043,8.9719,0.07,0.21]]
    data[0,*] = data[0,*]-alog10(1.54)

; from their Table 3, converting back to Chabrier IMF by subtracting
; alog10(1.54) from the stellar mass     
;    data = [$                   ; z=0.05-0.475 (median z=0.44)
;      [djs_mean([ 9.30, 9.75]),8.77-0.10,0.08],$
;      [djs_mean([ 9.75,10.25]),8.91-0.10,0.07],$
;      [djs_mean([10.25,10.75]),9.02-0.10,0.10]]
;    data[0,*] = data[0,*]-alog10(1.54)
    
;    data = [$                   ; z=0.475-0.9 (median z=0.75)
;      [djs_mean([ 9.75,10.25]),8.91-0.21,0.04],$
;      [djs_mean([10.25,10.75]),9.02-0.20,0.04],$
;      [djs_mean([10.75,11.25]),9.10-0.12,0.05]]
;    data[0,*] = data[0,*]-alog10(1.54)

    oploterror, data[0,*], data[1,*], data[2,*], psym=symcat(14), $
      color=im_color('forest green'), symsize=1.5, /hibar
    oploterror, data[0,*], data[1,*], data[3,*], psym=symcat(14), $
      color=im_color('forest green'), symsize=1.5, /lobar
    maxis = range(min(data[0,*])-0.05,max(data[0,*])+0.05,50)
;   djs_oplot, maxis, poly(maxis-alog10(1.54)-10.0,[8.70,0.17]), $
;     line=0, color=im_color('forest green')
        
; ##########    
; Lama+09: h=0.7, Chabrier IMF; metallicities recalibrated to T04 

; from their Fig 15 derived using dexter
;   mm = [9.1242,9.5083,9.7807,10.160] ; galaxies at z=0.5-0.6
;   oh = [8.6260,8.6725,8.7739,8.8542]
;   oherr = oh*0+0.05
;   oploterror, mm, oh, oherr, psym=symcat(16), color=im_color('navy')

    data = [$ ; z=0.6-0.8
      [9.2702,8.5305],$
      [9.5950,8.5906],$
      [9.8931,8.6397],$
      [10.304,8.7194]]
    oploterror, data[0,*], data[1,*], data[1,*]*0+0.05, psym=symcat(16), $
      color=im_color('navy'), symsize=1.5
    maxis = range(min(data[0,*])-0.05,max(data[0,*])+0.05,50)
    djs_oplot, maxis, poly(maxis-10.0,[8.65,0.150]), $
      line=0, color=im_color('navy')
    
;   oploterror, [10.19], [-0.23], [0.19], psym=symcat(16), color=im_color('navy')
;   oploterror, [9.4], [-0.12], [0.20], psym=symcat(9), color=im_color('navy')

; ##########    
; zahid+09: h=0.7, Chabrier IMF; metallicities derived using KK04/R23

; difference between their MZ relation and their recalibrated SDSS MZ relation
;   poly(maxis-10.0,[8.923,0.24,-0.06])-poly(maxis-10,[9.051,0.151,-0.104])    
    
; this is their local MZ relation    
;   djs_oplot, maxis1, poly(maxis1-10,[9.051,0.151,-0.104]), color='green'
    
; from their Fig 6 derived using dexter
    data = [$ ; z=0.75-0.82
      [9.2483,8.6907],$
      [9.3165,8.7698],$
      [9.3848,8.7778],$
      [9.4399,8.7752],$
      [9.4885,8.7457],$
      [9.5633,8.8020],$
      [9.6404,8.8314],$
      [9.7242,8.8395],$
      [9.7903,8.8609],$
      [9.8673,8.9266],$
      [9.9665,8.9360],$
      [10.063,8.9319],$
      [10.178,8.9654],$
      [10.327,9.0043],$
      [10.585,9.0458]]

; convert to T04 using Table 3 in Kewley & Ellison 08
    data[1,*] = poly(data[1,*],[-461.2352D,158.4484D,-17.94607D,0.682817D]) ; [KK04-->T04]
    
    oploterror, data[0,*], data[1,*], data[1,*]*0+0.03, psym=symcat(6,thick=5), $
      color=im_color('orange red'), symsize=1.3
    maxis = range(min(data[0,*])-0.05,max(data[0,*])+0.05,50)

    mzfit = poly(maxis-10.0,[8.923,0.24,-0.06])
    mzfit = poly(mzfit,[-461.2352D,158.4484D,-17.94607D,0.682817D]) ; [KK04-->T04]
    djs_oplot, maxis, mzfit, line=0, color=im_color('orange red')

; make a legend
    im_legend, ['CB08','L09','ZKB11'], /right, /bottom, box=0, $
      psym=[14,16,6], color=['forest green','navy','orange red'], $
      charsize=1.5, symsize=[2.2,1.7,1.7]
    
    im_plotconfig, /psclose, psfile=psfile

stop        
    

return
end
