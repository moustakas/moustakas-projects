pro talk_qa_one_mask, pos=pos, keycolor=keycolor
; this is Alison's routine to make Fig 5 in the survey paper,
; but modified for Keynote
    
    survey='cfhtls_xmm'
    run='oct07a'
    mask='0210'

    runs=yanny_readone(getenv('PRIMUS_DIR')+'/data/primus_masks.par')
    isurvey=where(runs.survey eq survey and runs.run eq run)
    runs=runs[isurvey[0]]

    if(file_test(getenv('PRIMUS_DIR')+'/data/primus_masks_refrac_'+ $
      survey+'_'+run+'.par') gt 0) then $
        refrac=yanny_readone(getenv('PRIMUS_DIR')+'/data/primus_masks_refrac_'+ $
      survey+'_'+run+'.par')
    irefrac=(where(refrac.mask eq mask))[0]
    refrac=refrac[irefrac]

    masks=mrdfits(primus_name('masks', survey=survey, run=run),1)
    masks_stats=mrdfits(primus_name('masks_stats', survey=survey, run=run),1)
    imask=where(masks.mask eq mask)
    slits=mrdfits(primus_name('slits', survey=survey, run=run, mask=mask)+'.gz', 1)
    target=mrdfits(primus_name('targets', survey=survey), 1)
    read_mangle_polygons, primus_name('masks_poly',  run=run,survey=survey), poly
    if(runs.bsmask eq 'half') then $
      read_mangle_polygons, primus_name('bsmaskhalf', survey=survey), bsmask $
    else $
      read_mangle_polygons, primus_name('bsmask', survey=survey), bsmask 
    poly=poly[imask]

;openw,unit, primus_name('masks_dir', survey=survey, run=run)+ $
;  '/'+masks[imask].obsbase+'.sub', /get_lun
    iref=where((slits.primus_target AND $
      primus_flagval('PRIMUS_TARGET', 'REFERENCE')) gt 0 AND $
      slits.priority eq -99., nref)
;printf, unit, '1 1'
    for i=0L, nref-1L do begin
       radec2ccdnum,slits[iref[i]].ra, slits[iref[i]].dec, masks[imask].racen, $
         masks[imask].deccen, masks[imask].pa, ccdnum=ccdnum
       if(ccdnum eq -1) then message, 'reference star near edge'
       pa_mask= masks[imask].pa-90.
       baade_radec2xy, slits[iref[i]].ra, slits[iref[i]].dec, $
         masks[imask].racen, masks[imask].deccen, pa_mask, x=xf, y=yf
       imacs_xy2xy, xf, yf, xccd, yccd
       imacs_xy2pix, xccd, yccd, col, row, ccd=ccdnum
;    printf, unit, strtrim(string(ccdnum),2)+' '+ $
;      strtrim(string(col),2)+' '+ $
;      strtrim(string(row),2)+' 100 100'
    endfor
;free_lun, unit

    if(tag_indx(masks, 'NODEAST') eq -1) then $
      nodeast=8 $ 
    else $
      nodeast=masks[imask].nodeast
    if(tag_indx(masks, 'NODNORTH') eq -1) then $
      nodnorth=16 $
    else $
      nodnorth=masks[imask].nodnorth
    if(tag_indx(masks, 'SHUFFLE') eq -1) then $
      shuffle=8 $
    else $
      shuffle=masks[imask].shuffle

    corners=fltarr(2,4,8)
    for i=0L, 7L do $
      corners[*,*,i]=ccdoutline(masks[imask].racen, masks[imask].deccen, $
      masks[imask].pa, i+1)

    inz=where(runs.priority ne 0. and runs.priority ne 99., nnz)

;fsc_color:
    colstr=['forest green','yellow','dodger blue','magenta']

    sizes=[0.8,0.6,0.7,0.7]
    priority=[10000.,308.,208.,108.]

    mmra=minmax(corners[0,*,*])
    mmdec=minmax(corners[1,*,*])
    cosdec=cos(!DPI/180.*masks[imask].deccen)
    mdiff=max([abs(mmra-masks[imask].racen)*cosdec, $
      abs(mmdec-masks[imask].deccen)])
    plot, [0], [0], /nodata, xra=masks[imask].racen+mdiff/cosdec*[1.05, -1.05], $
      yra=masks[imask].deccen+mdiff*[-1.05, 1.05], position=pos, $
      xtitle=textoidl('\alpha_{J2000} (deg)'), $
      ytitle=textoidl('\delta_{J2000} (deg)'), color=keycolor, xsty=1, ysty=1
; background of mask region - grey
    plot_poly, poly, minside=40, dangle=0.01, color=fsc_color('grey',101), /fill, /over

; bright star mask regions - white
    plot_poly, bsmask, minside=40, dangle=0.01, color=fsc_color('white',101), /over, /fill

    iin=where(is_in_polygon(poly, ra=target.ra, dec=target.dec) gt 0 AND $
      (target.primus_target AND $
      primus_flagval('PRIMUS_TARGET', 'TARGET')) gt 0)

; this is the *untargeted* targets - don't show!  (little black dots)
;   djs_oplot, target[iin].ra, target[iin].dec, psym=3

    hogg_usersym, 10, /fill
;isort=sort(runs.priority)
;iprior=isort[uniq(runs.priority[isort])]
;priority=runs.priority[iprior]

    for i=0,3 do begin
       ii=where(slits.priority eq priority[i], nii)
       if(nii gt 0) then $
         djs_oplot,[slits[ii].ra], [slits[ii].dec],psym=8, $
         color=fsc_color(colstr[i],10), $
         symsize=sizes[i]
    endfor
    hogg_usersym, 4, /fill
    ii=where((slits.primus_target AND $
      primus_flagval('PRIMUS_TARGET', 'REFERENCE')) gt 0 AND $
      slits.priority eq -99.)
    djs_oplot,slits[ii].ra, slits[ii].dec,psym=8, $
      color='black', symsize=2.
    hogg_usersym, 10, /fill
    
    plot,[0,0],[0,0], xra=masks[imask].racen+mdiff/cosdec*[1.05, -1.05], $
      yra=masks[imask].deccen+mdiff*[-1.05, 1.05], position=pos,/noerase,$
;     ,xmargin=[8,2],ymargin=[6,2],/noerase, $
      color=keycolor, xsty=1, ysty=1
    for i=0, 7 do begin
       oplot, corners[0,[0,1,2,3,0],i], corners[1,[0,1,2,3,0],i], th=2, color=keycolor
       mra=mean(corners[0,*,i])
       mdec=mean(corners[1,*,i])
    endfor 

return
end

pro talk_10sep_caltech, keynote=keynote
; jm10sep21ucsd - miscellaneous plots for my 2010/Sep Caltech talk 

    common com_talk, sdss_parent
    
    mfpath = mf_path()
    isedpath = mf_path(/isedfit)
    
    rootpath = getenv('RESEARCHPATH')+'/talks/2010/10sep_caltech/'
    datapath = rootpath+'data/'
    if keyword_set(keynote) then talkpath = rootpath+'keynote/' else $
      talkpath = rootpath+'figures/'

; --------------------------------------------------
; zvz - read the data
    zerod = primus_read_zerod(field='zspec',zmerge=zmerge)
    zcut = where((zmerge.z_spec gt 0.0) and (zmerge.z_spec lt 1.2) and $
      (zmerge.good_spec eq 1) and (zerod.zprimus_zconf ge 3) and $
      (strtrim(zerod.zprimus_class,2) eq 'GALAXY') and $
      ((zerod.fieldnum lt 3 and zerod.targ_mag le 22.5) or $
      (zerod.fieldnum eq 4 and zerod.targ_mag le 22.8) or $
      (zerod.fieldnum eq 5 and zerod.targ_mag le 22.8) or $
      (zerod.fieldnum eq 6 and zerod.targ_mag le 22.5)),nzconf3)

    allfield = strtrim(zmerge[zcut].field,2)
    allzconf = zerod[zcut].zprimus_zconf
    allzprimus = zerod[zcut].zprimus
    allzspec = zmerge[zcut].z_spec

; use just a fraction of the deep2 sources
    isdeep = where(strmatch(allfield,'*deep*'),nisdeep,comp=notdeep)
    nkeep = long(nisdeep*0.2)
    keep = shuffle_indx(nisdeep,num_sub=nkeep,seed=seed)

    zconf = [allzconf[notdeep],allzconf[isdeep[keep]]]
    zprimus = [allzprimus[notdeep],allzprimus[isdeep[keep]]]
    zspec = [allzspec[notdeep],allzspec[isdeep[keep]]]
    field = [allfield[notdeep],allfield[isdeep[keep]]]

; shuffle the indices so that they don't overlap on the plot
    npts = n_elements(zprimus)
    indx = shuffle_indx(npts)
    zconf = zconf[indx]
    zprimus = zprimus[indx]
    zspec = zspec[indx]
    field = field[indx]

; make the plot    
    psfile = talkpath+'zvz.ps'
    im_plotconfig, 1, pos, psfile=psfile, keynote=keynote, charsize=1.5

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')
    
; plotting preferences    
    nfield = 4
    prefs = replicate({field: '', nicefield: '', color: '', $
      psym: 0, symsize: 1.0},nfield)
    prefs.field = ['cosmos','calib','deep2_02hr','deep2_23hr']
    prefs.nicefield = ['zCOSMOS','VVDS','DEEP2 02^{hr}','DEEP2 23^{hr}']
    if keyword_set(keynote) then $
      prefs.color = ['white','medium orchid','dodger blue','tan'] else $
        prefs.color = ['black','medium orchid','dodger blue','tan']
    prefs.psym = [14,9,15,4]
    prefs.symsize = [1.1,0.9,0.9,1.1]*0.7

    zrange = [-0.06,1.27]
    xtitle1 = 'Redshift from High-Resolution Spectroscopy'
    ytitle1 = 'PRIMUS Redshift'

    zaxis = findgen((2.0-(-1.0))/0.01)*0.01+(-1.0)
    dzcut_03 = 0.03*(1.0+zaxis)
    dzcut_10 = 0.10*(1.0+zaxis)

; left panel: Q>3; render each point one at a time to make sure
; they don't overlap 
    plot, [0], [0], /nodata, /xsty, /ysty, position=pos[*,0], $
      xrange=zrange, yrange=zrange, xtitle='', ytitle=ytitle1, $
      color=keycolor
    for ii = 0, npts-1 do begin
       this = where(field[ii] eq strtrim(prefs.field,2))
       plots, zspec[ii], zprimus[ii], psym=symcat(prefs[this].psym,thick=6), $
         color=fsc_color(prefs[this].color,101), symsize=prefs[this].symsize
    endfor
    oplot, zaxis, zaxis, line=0, thick=5, color=keycolor
    oplot, zaxis, zaxis-dzcut_03, line=5, thick=5, color=keycolor
    oplot, zaxis, zaxis+dzcut_03, line=5, thick=5, color=keycolor
    oplot, zaxis, zaxis-dzcut_10, line=1, thick=5, color=keycolor
    oplot, zaxis, zaxis+dzcut_10, line=1, thick=5, color=keycolor

    stats = primus_zvz_stats(zspec,zprimus)
    label = [$
;     'Q\geq3',$
      'F(>0.03)='+strtrim(string(100.0*stats.frac_bad_03,format='(F12.1)'),2)+'%',$
      'F(>0.10)='+strtrim(string(100.0*stats.frac_bad_10,format='(F12.1)'),2)+'%',$
      '\sigma_{\Delta'+'z/(1+z)}='+strtrim(string(stats.mad,format='(F12.4)'),2)]
    legend, textoidl(label), /left, /top, box=0, charsize=1.3, $
      margin=0, textcolor=keycolor
    xyouts, 0.43, 0.82, textoidl('Q>3'), align=0.5, /data, color=keycolor
    
; right panel: Q=4
    plot, [0], [0], /nodata, /xsty, /ysty, /noerase, position=pos[*,1], $
      xrange=zrange, yrange=zrange, xtitle='', ytickname=replicate(' ',10), $
      color=keycolor
    best = where(zconf eq 4,nbest)
    for ii = 0, nbest-1 do begin
       this = where(field[best[ii]] eq strtrim(prefs.field,2))
       plots, zspec[best[ii]], zprimus[best[ii]], psym=symcat(prefs[this].psym,thick=6), $
         color=fsc_color(prefs[this].color,101), symsize=prefs[this].symsize
    endfor
    oplot, zaxis, zaxis, line=0, thick=5, color=keycolor
    oplot, zaxis, zaxis-dzcut_03, line=5, thick=5, color=keycolor
    oplot, zaxis, zaxis+dzcut_03, line=5, thick=5, color=keycolor
    oplot, zaxis, zaxis-dzcut_10, line=1, thick=5, color=keycolor
    oplot, zaxis, zaxis+dzcut_10, line=1, thick=5, color=keycolor
    
    stats = primus_zvz_stats(zspec[best],zprimus[best])
    label = [$
;     'Q=4',$
      'F(>0.03)='+strtrim(string(100.0*stats.frac_bad_03,format='(F12.1)'),2)+'%',$
      'F(>0.10)='+strtrim(string(100.0*stats.frac_bad_10,format='(F12.1)'),2)+'%',$
;     '\sigma_{\Delta'+'z/(1+z)}='+strtrim(string(100.0*stats.mad,format='(F12.3)')+'%',2)]
      '\sigma_{\Delta'+'z/(1+z)}='+strtrim(string(stats.mad,format='(F12.4)'),2)]
    legend, textoidl(label), /left, /top, box=0, charsize=1.3, $
      margin=0, textcolor=keycolor
    xyouts, 0.43, 0.82, 'Q=4', align=0.5, /data, color=keycolor
    
; final legend
    im_legend, prefs.nicefield, textcolor=keycolor, $
      psym=prefs.psym, color=strtrim(prefs.color,2), $
      /right, /bottom, box=0, margin=0, charsize=1.2, $
      symsize=prefs.symsize*2.3, symthick=6
    
; x-title
    xyouts, pos[2,0], pos[1,0]-0.12, xtitle1, align=0.5, $
      /norm, color=keycolor
    
;    magcut = [22.5,22.5,22.8,22.8]
;    for ii = 0, nfield-1 do begin
;       zerod = primus_read_zerod(field=field[ii],/alltarg,zmerge=zmerge)
;       zcut = where((zerod.zprimus_zconf ge 4) and (strtrim(zerod.zprimus_class,2) eq 'GALAXY') and $
;         (zmerge.z_spec gt 0.0) and (zmerge.z_spec lt 1.2) and $
;         (zmerge.good_spec eq 1) and (zerod.targ_mag le magcut[ii]))
;       zerod = zerod[zcut]
;       zmerge = zmerge[zcut]
;
;       if (ii eq 0) then begin
;          allzerod = zerod
;          allzmerge = zmerge
;       endif else begin
;          allzerod = [allzerod,zerod]
;          allzmerge = [allzmerge,zmerge]
;       endelse
;    endfor       

    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; redshift success rate in bins of magnitude and color
; choose your fields    
    field = ['cdfs','cosmos','xmm','es1'] ; sort by i,R mag
    nfield = n_elements(field)

; plot parameters    
    nmagpix = 50
    ncolorpix = 50
    levels = [0.10,0.25,0.5,0.75,0.9]

    mincol = 50 ; 100 
    ndiv = 5
    loadct, 16 ; 15
    ticknames = string(range(0.0,1.0,ndiv+1),format='(F5.2)')

    magrange = [18.0,23.5]
    colorrange = [-0.5,2.0]
    nmock = 2000L
    mag = randomu(seed,nmock,nmock)*(magrange[1]-magrange[0])+magrange[0]
    color = randomu(seed,nmock,nmock)*(colorrange[1]-colorrange[0])+colorrange[0]
    
    psfile = talkpath+'completeness.ps'
    ccharsize = 1.5
    im_plotconfig, 2, pos, psfile=psfile, charsize=ccharsize, $
      xmargin=[1.2D,2.25D], ymargin=[0.2D,0.9D], xspace=1.0D, $
      yspace=1D, width=3.5D*[1,1], height=2.8D*[1,1], keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    for ii = 0, nfield-1 do begin
       faint = mf_maglimit(field[ii],bright=bright)
       params = get_completeness_params(field[ii])
       zsuccess = get_mf_zsuccess(mag,color,field=field[ii],/noextrap)
       
       good = where(zsuccess gt 0.0)
       mf_hogg_scatterplot, mag[good], color[good], weight=zsuccess[good], position=pos[*,ii], $
         noerase=(ii gt 0), xsty=1, ysty=1, xrange=magrange, $
         yrange=colorrange, xnpix=14, ynpix=14, $ 
         xtitle=params.nice_filter+' (AB mag)', ytitle=params.nice_color_filter, $
         color=keycolor, exp=2.0, darkest=mincol, $
         levels=levels;, /nocon;, cannotation=cannotation, /nocontours
       djs_oplot, faint*[1,1], !y.crange, line=0, thick=10, color=keycolor
; overplot the axes in black to make the tick marks visible 
       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,ii], $
         xsty=1, ysty=1, xrange=magrange, yrange=colorrange, $
         xtickname=replicate(' ',10), ytickname=replicate(' ',10)
       
       legend, strupcase(strtrim(field[ii],2)), position=[pos[0,ii]-0.01,pos[3,ii]], $
         /norm, box=0, charsize=1.4, charthick=3.5;textcolor=djs_icolor('brown')
    endfor

    off1 = 0.03 & hei = 0.3 & wid = 0.05 & off2 = 0.13
    cpos = [pos[2,1]+off1,(pos[1,1]-pos[3,3])/2+pos[3,3]-hei,pos[2,1]+off1+wid,(pos[1,1]-pos[3,3])/2+pos[3,3]+hei]
    primus_fsc_colorbar, position=cpos, range=[0.0,1.0], /vertical, $
      charsize=ccharsize, minor=4, division=ndiv, ticknames=ticknames, $
      /invert, bottom=mincol, /right, color=keycolor
    xyouts, cpos[0]+off2, (cpos[3]-cpos[1])/2.0+cpos[1], 'Redshift Success', $
      orientation=270, align=0.5, charsize=1.8, /norm, color=keycolor

    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; mock stellar mass function
    psfile = talkpath+'fake_mf.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, $
      xmargin=[1.5,0.2], width=6.8, height=5.5, charsize=2

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    mf_fit_sdss = read_binned_mf(/bestfit,/sdss,$
      active=active,quiescent=quiescent,noevol=noevol)

    maxis1 = range(9.1,12.5,100)
    xrange = [9.0,12.3]
    yrange = [5E-6,6E-2]

    plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, /ylog, $
      yrange=yrange, xrange=xrange, xtickinterval=1, $
      color=keycolor, xtitle='Stellar Mass', $
      ytitle='Number Density', xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10)
    djs_oplot, maxis1, mf_schechter(maxis1,mf_fit_sdss), $
      line=0, thick=10, color=djs_icolor('cyan')
;   x1 = 11
;   x2 = 11.5
;   y1 = 2E-3
;   y2 = 1E-2
;   arrow, x1, y1, x1, y2, /data, hsize=-0.25, hthick=8, thick=8
;   arrow, x1, y1, x2, y1, /data, hsize=-0.25, hthick=8, thick=8
    
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; Omega* vs redshift + SDSS + PRIMUS + literature    
    H0 = (100.0*mf_h100())*3.1556926D7*1D6/3.086D19 ; [Myr^-1]
    bigG = 4.4983D-21 ; [Mpc^3/Myr^2/M_sun]
    rhocrit = 3.0*H0^2.0/(8.0*!pi*bigG) ; [M_sun/Mpc^3, h=0.7]
;   rhocrit = 1.35D11 ; [M_sun/Mpc^3, h=0.7, from Pettini 2006]

    zbins = mf_zbins(nzbins)
    field = get_mf_fields(/nodeep,/nodls)
    prefs = mfplot_prefs(field)
    
    psfile = talkpath+'z_vs_omega.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.5,0.2], $
      width=6.8, height=5.8, charsize=2, keynote=keynote
    
    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      yrange=[-3.3,-2], xrange=[0,1.08], xtitle='Redshift', $
      ytitle=mfplot_omegatitle(), ytickinterval=0.5, color=keycolor

; wilkins+; note: Chabrier = Wilkins - 0.03 (ish)
    wilk = rsex(getenv('PAPERSPATH')+'/literature/data/08wilkins.sex')
    zwilk = total([[wilk.zmin],[wilk.zmax]],2)/2.0
    loz = where(zwilk lt 0.15,nloz)
    zwilk[loz] = zwilk[loz]+randomn(seed,nloz)*0.01
    
    oploterror, zwilk, alog10(wilk.omega)-0.03, wilk.omega_err/wilk.omega/alog(10), $
      psym=symcat(16), color=fsc_color('tan',101), errcolor=fsc_color('tan',101), $
      symsize=1.5

; PRIMUS - field-by-field and then the average
    for ii = 0, n_elements(field)-1 do begin
       rho = get_mf_param(field=field[ii],/rho,param_err=rho_err,noevol=noevol)
       omega = alog10(rho/rhocrit)
       omega_err = rho_err/rho/alog(10)
       oploterror, zbins.zbin, omega, omega_err, psym=symcat(prefs[ii].psym,thick=8), $
         symsize=3, errthick=8, color=fsc_color(prefs[ii].color,101), $
         errcolor=fsc_color(prefs[ii].color,101)
    endfor
    rho = get_mf_param(/rho,param_err=rho_err,noevol=noevol)
    omega = alog10(rho/rhocrit)
    omega_err = rho_err/rho/alog(10)
    oploterror, zbins.zbin, omega, omega_err, psym=symcat(6,thick=8), $
      symsize=3, errthick=8

; sdss    
    rho_sdss = get_mf_param(/rho,/sdss,param_err=rho_err,noevol=noevol)
    omega_sdss = alog10(rho_sdss/rhocrit)
    plots, 0.1, omega_sdss, psym=symcat(34,thick=10), symsize=5, $
      color=fsc_color('firebrick',101)

; put the legend last
    label = [strupcase(strtrim(prefs.field,2)),'SDSS','Wilkins+08']
    psym = [prefs.psym,34,16]
    color = [strtrim(prefs.color,2),'firebrick','tan']
    symsize = [prefs.symsize*1.8,2.0,1.5]

    im_legend, label[0:2], psym=psym[0:2], color=color[0:2], $
      /left, /bottom, box=0, margin=0, charsize=1.6, $
      symsize=symsize[0:2], textcolor=keycolor, symthick=6, $
      position=[0.03,!y.crange[0]+0.08], /data
    im_legend, label[3:5], psym=psym[3:5], color=color[3:5], $
      /left, /bottom, box=0, margin=0, charsize=1.6, $
      symsize=symsize[3:5], textcolor=keycolor, symthick=6, $
      position=[0.26,!y.crange[0]+0.08], /data
;   im_legend, label, psym=psym, color=color, $
;     /left, /bottom, box=0, margin=0, charsize=1.6, $
;     symsize=symsize, textcolor=keycolor, symthick=6
      
    
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; PRIMUS stellar mass functions for all galaxies
    fields = get_mf_fields(/nodeep2,/nodls)
    psfile = talkpath+'binned_mf_all_evol.ps'
    plot_binned_mf, psfile=psfile, field=fields, keynote=keynote
    
; --------------------------------------------------
; limiting stellar mass versus redshift
    field = get_mf_fields(/nodeep2,/nodls)
    nfield = n_elements(field)
    zbins = mf_zbins(nzbins,zmin=zmin,zmax=zmax)

    xrange = [0.16,1.1]
    yrange = [8.05,12.0]
    zaxis = range(zmin,zmax,100)
    
    psfile = talkpath+'mass_vs_redshift.ps'
    im_plotconfig, 2, pos, psfile=psfile, charsize=1.7, $
      width=4.2D*[1,1], height=2.9*[1,1], xmargin=[1.1D,0.3D], $
      keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

; loop on each field
    for ii = 0, nfield-1 do begin
       parent = read_mf_parent(field=field[ii])
       faint = mf_maglimit(field[ii],nice_filter=nice_filter)
       
       if (ii eq 0) or (ii eq 2) then begin
          ytitle = mfplot_masstitle()
          delvarx, ytickname
       endif else begin
          delvarx, ytitle
          ytickname = replicate(' ',10)
       endelse
       if (ii lt 2) then begin
          xtitle = '' & xtickname = replicate(' ',10)
       endif else begin
          xtitle = 'Redshift' & delvarx, xtickname 
       endelse

       mfplot_scatterplot, parent.zprimus, parent.mass_avg[1], $
         weight=parent.final_weight, noerase=(ii gt 0), $
         position=pos[*,ii], xrange=xrange, yrange=yrange, $
         xsty=1, ysty=1, ytitle=ytitle, xtickname=xtickname, $
         ytickname=ytickname, npix=20, cthick=2.0, xtitle=xtitle, $
         color=keycolor
; overplot the axes in black to make the tick marks visible 
       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,ii], $
         xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
         xtickname=replicate(' ',10), ytickname=replicate(' ',10)
       legend, strupcase(strtrim(field[ii],2))+' ('+nice_filter+'<'+$
         string(faint,format='(F4.1)')+')', /left, /top, $
         box=0, margin=0, charsize=1.3
    endfor
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; scaling between exposure time and resolution
    psfile = talkpath+'exptime_scaling.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, $
      xmargin=[1.5,0.2], width=6.8, height=5.5

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    bigr = range(1,200,50)
    exptime_photoz = bigr^2
    exptime_primus = bigr

    plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='Spectral Resolution, R', xrange=[1,max(bigr)], $
      yrange=minmax(exptime_photoz), $
      /xlog, /ylog, ytitle='Relative Exposure Time', color=keycolor
    djs_oplot, bigr, exptime_primus, line=0, color='blue', thick=10
    djs_oplot, bigr, exptime_photoz, line=5, color='orange', thick=10

;   djs_oplot, 15*[1,1], 10^!y.crange, line=
    
    xyouts, 30, 8.0, textoidl('\propto'+'R'), align=0.5, $
      color=keycolor, charsize=2
    xyouts, 10, 300.0, textoidl('\propto'+'R^{2}'), align=0.5, $
      color=keycolor, charsize=2
    legend, ['PRIMUS',"Photo-z's"], /left, /top, box=0, $
      pspacing=1.9, line=[0,5], color=djs_icolor(['blue','orange']), $
      textcolor=keycolor
    
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; SDSS MF for all galaxies
    mf_data = read_binned_mf(/sdss)
    mf_fit = read_binned_mf(/bestfit,/sdss)

    psfile = talkpath+'binned_mf_sdss_lit.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.5,0.4], $
      width=6.6, height=6.4, charsize=2, keynote=keynote

    xrange = [9.4,12.0]
    yrange = [4E-6,1.5E-2]
;   yrange = [1E-7,2E-2]
    maxis1 = range(8+0.04,12.5,100)
    
    cline = 3 & ccolor = 'orange'
    bline = 5 & bcolor = 'blue'
    pline = 0 & pcolor = 'red'
    
    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      /ylog, yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytitle=mfplot_phititle(), color=keycolor
;     ['log(M^{*}/M'+sunsymbol()+')='+$
;     strtrim(string(mf_fit.mstar,format='(F12.2)'),2)+,$
;     ['M^{*}='+strmid(strtrim(10^mf_fit.mstar,2),0,4)+'\pm'+$
;     strmid(strtrim(mf_fit.mstar_err*mf_fit.mstar*alog(10),2),0,4)+$
;     ' M'+sunsymbol(),$
    label = textoidl($
      ['M^{*}='+strmid(strtrim(10^mf_fit.mstar,2),0,4)+'\times10^{10} h_{70}^{-2}M'+sunsymbol(),$
      '\Phi^{*}='+strmid(strtrim(mf_fit.phistar,2),0,6)+' h_{70}^3 Mpc^{-3}',$
      '\alpha='+strtrim(string(mf_fit.alpha,format='(F12.2)'),2)])
    legend, label, /left, /top, box=0, position=[9.5,3E-4], $
      /data, spacing=2.9, textcolor=keycolor
    
;   legend, ['Cole+01','Bell+03','Panter+07','This Paper'], $
;     /left, /bottom, box=0, line=[cline,bline,pline,0], $
;     pspacing=1.7, charsize=1.6, thick=8, $
;     color=djs_icolor([ccolor,bcolor,pcolor,'default']);, margin=0
    
    full = where((mf_data.fullbin eq 1)); and (mf_data.phi gt 10^!y.crange[0]))
    oploterror, mf_data.mass[full], mf_data.phi[full], mf_data.phierr[full], $
      color=keycolor, errcolor=keycolor, thick=8, psym=10
;   djs_oplot, mf_data.mass[full], mf_data.phi[full], $
;     color=fsc_color('black',100), thick=8, psym=10

;; overplot fits from the literature
;    mfoplot_lit, maxis1, color=ccolor, line=cline, thick=8, /cole
;    mfoplot_lit, maxis1, color=bcolor, line=bline, thick=8, /bell
;    mfoplot_lit, maxis1, color=pcolor, line=pline, thick=8, /panter
; overplot our best fit    
    djs_oplot, maxis1, mf_schechter(maxis1,mf_fit), line=5, thick=10, $
      color=fsc_color('dodger blue',101)
    
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; SED-fitting example    
    field = 'xmm'
    params = read_isedfit_paramfile(isedpath+field+'_isedfit.par',sfhgrid='02')
    index = 6144
    galaxy = 'J344800.1-040108.4'
    model = isedfit_restore(paramfile,isedfit,params=params,index=index,iopath=isedpath)
    filters = strtrim(get_field_filters(field,nice=nice_filters),2)
    filtinfo = im_filterspecs(filterlist=filters)
;   ised = read_mf_isedfit(field=field,sfhgrid='02',params=params)

    xtitle1 = textoidl('Observed Wavelength (\AA)')
    xtitle2 = textoidl('Rest Wavelength (\AA)')
    ytitle1 = textoidl('m_{AB}')

    psfile = talkpath+'isedfit_example.ps'
    im_plotconfig, 8, pos, psfile=psfile, ymargin=[1.0,1.1], keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

; make the plot       
    yrange = [28.5,18.5]
    xrange1 = [1000.0,6E4]
    xrange2 = xrange1/(1.0+isedfit.zobj)
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, /xlog, $
      xtickinterval=1000, position=pos, color=keycolor
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, /xlog, color=keycolor
    oplot, model.wave, model.flux, line=0, color=keycolor; color='grey'

; overplot the filter names
    nice_filters = repstr(repstr(nice_filters,'ch1','[3.6]'),'ch2','[4.5]')
    for ii = 0, n_elements(filters)-3 do xyouts, filtinfo[ii].weff, $
      28.0, nice_filters[ii], /data, align=0.5, color=keycolor, $
      charsize=1.4
    
; overplot the observed and model photometry
    used = where((isedfit.maggies gt 0.0) and $ ; used in the fitting
      (isedfit.ivarmaggies gt 0.0),nused)
    notused = where((isedfit.maggies gt 0.0) and $ ; not used in the fitting
      (isedfit.ivarmaggies eq 0.0),nnotused)
    nodata = where((isedfit.maggies eq 0.0) and $ ; no measurement
      (isedfit.ivarmaggies eq 0.0),nnodata)

    djs_oplot, filtinfo[used].weff, -2.5*alog10(isedfit[used].bestmaggies), $
      psym=symcat(6,thick=6), symsize=2.5

; overplot galex, cfhtls, and irac
    galex = where(strmatch(filters,'*galex*',/fold))
    mab = maggies2mag(isedfit.maggies[galex],$
      ivar=isedfit.ivarmaggies[galex],magerr=mab_err)
    oploterror, filtinfo[galex].weff, mab, mab_err, psym=symcat(16), $
      symsize=2.0, color=fsc_color('dodger blue',101), $
      errcolor=fsc_color('dodger blue',101), errthick=!p.thick
    
    optical = where(strmatch(filters,'*capak*',/fold))
    mab = maggies2mag(isedfit.maggies[optical],$
      ivar=isedfit.ivarmaggies[optical],magerr=mab_err)
    oploterror, filtinfo[optical].weff, mab, mab_err, psym=symcat(16), $
      symsize=2.0, color=fsc_color('tan',101), $
      errcolor=fsc_color('tan',101), errthick=!p.thick
    
    irac = where(strmatch(filters,'*irac*',/fold))
    mab = maggies2mag(isedfit.maggies[irac],$
      ivar=isedfit.ivarmaggies[irac],magerr=mab_err)
    oploterror, filtinfo[irac].weff, mab, mab_err, psym=symcat(16), $
      symsize=2.0, color=fsc_color('firebrick',101), $
      errcolor=fsc_color('firebrick',101), errthick=!p.thick
    
    label = [$
      galaxy,'z = '+string(isedfit.zobj,format='(F6.4)'),$
      'log (M/M'+sunsymbol()+') = '+strtrim(string(isedfit.mass_avg,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(isedfit.mass_err,format='(F12.2)'),2)]

;    label = [$
;      'log (M/M'+sunsymbol()+') = '+strtrim(string(isedfit.mass,format='(F12.2)'),2)+$
;      ' ('+strtrim(string(isedfit.mass_avg,format='(F12.2)'),2)+')',$
;      'E(B-V) = '+strtrim(string(isedfit.ebv,format='(F12.3)'),2)+$
;      ' ('+strtrim(string(isedfit.ebv_avg,format='(F12.3)'),2)+')',$
;      'Z/Z'+sunsymbol()+' = '+strtrim(string(isedfit.Z/0.02,format='(F12.2)'),2)+$
;      ' ('+strtrim(string(isedfit.Z_avg/0.02,format='(F12.2)'),2)+')',$
;      '\tau = '+strtrim(string(isedfit.tau,format='(F12.1)'),2)+$
;      ' ('+strtrim(string(isedfit.tau_avg,format='(F12.2)'),2)+') Gyr',$
;      'Age = '+strtrim(string(isedfit.age,format='(F12.2)'),2)+$
;      ' ('+strtrim(string(isedfit.age_avg,format='(F12.2)'),2)+') Gyr',$
;      '\psi_{100} = '+strtrim(string(isedfit.sfr100,format='(F12.3)'),2)+$
;      ' ('+strtrim(string(isedfit.sfr100_avg,format='(F12.3)'),2)+') M'+sunsymbol()+' yr^{-1}',$
;      'b_{100} = '+strtrim(string(isedfit.birthrate,format='(F12.3)'),2)+$
;      ' ('+strtrim(string(isedfit.birthrate_avg,format='(F12.3)'),2)+')']
     legend, textoidl(label), /left, /top, box=0, charsize=1.5, $
       margin=0, textcolor=keycolor, spacing=2.2
;    label = textoidl([strtrim(repstr(galaxy,'_',' '),2),$
;      'z = '+string(isedfit.zobj,format='(F6.4)'),'\chi^{2} = '+$
;      strtrim(string(isedfit.chi2,format='(F12.2)'),2)])
;    legend, label, /right, /bottom, box=0, spacing=1.5, charsize=1.4, margin=0

    im_plotconfig, psfile=psfile, /psclose, keynote=keynote

stop    
    
; --------------------------------------------------
; SDSS - 0.5(u-z) vs 0.5z
    zbins = mf_zbins(/sdss)
    if (n_elements(sdss_parent) eq 0L) then sdss_parent = read_mf_parent(field='sdss')
    uv = where(sdss_parent.in_galex_window and (sdss_parent.fuv_good or sdss_parent.nuv_good))
    mz = sdss_parent[uv].k_ugriz_absmag_05[4]
    umz = sdss_parent[uv].k_ugriz_absmag_05[0]-sdss_parent[uv].k_ugriz_absmag_05[4]
    zmk = sdss_parent[uv].k_ugriz_absmag_05[4]-sdss_parent[uv].k_ubvrijhk_absmag_05[7]
    weight = sdss_parent[uv].galex_weight
    
    psfile = talkpath+'mz_vs_umz_sdss.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    xrange = [-18.3,-23.7]
    yrange = [0.3,6.9]
    xtitle = mfplot_mztitle()
    ytitle = textoidl('^{0.5}(u - z)')

    mfplot_scatterplot, mz, umz, weight=weight, position=pos, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, /sdss, color=keycolor;, $
;     ccolor=keycolor ; fsc_color('tan',101)
; overplot the axes in black to make the tick marks visible 
    djs_plot, [0], [0], /nodata, /noerase, position=pos, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    xyouts, -19.2, 6.3, 'Quiescent', align=0.5, charsize=1.8
    xyouts, -19.5, 0.8, 'Star-Forming', align=0.5, charsize=1.8
;   xyouts, -22.8, 0.8, 'Star-Forming', align=0.5, charsize=1.6

    junk = mf_select_quiescent(z=zbins.zbin,umz_cut=umz_cut)
    djs_oplot, [!x.crange[0]-0.1,!x.crange[1]+0.1], $
      umz_cut*[1,1], line=0, color='red', thick=8

; show a reddening vector
    ebv = 0.3
    mz_true = -23.3
    umz_true = 0.9

    uweff = k_lambda_eff(filterlist='sdss_u0.par',band_shift=0.5)
    zweff = k_lambda_eff(filterlist='sdss_z0.par',band_shift=0.5)
    
    mz_red = mz_true + 0.4*ebv*k_lambda(zweff,/calz)
    umz_red = umz_true + 0.4*ebv*(k_lambda(uweff,/calz)-k_lambda(zweff,/calz))

    arrow, mz_true, umz_true, mz_red, umz_red, /data, $
      hsize=-0.25, hthick=8, thick=8
    xyouts, (mz_red-mz_true)/2.0+mz_true+0.3, (umz_red-umz_true)/2.0+umz_true-0.22, $
      orientation=-45, /data, 'E(B-V) = '+string(ebv,format='(F3.1)'), $
      charsize=1.3, align=0.5

    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; remake Fig 7 in Ken's paper
    psfile = talkpath+'11wong_fig7.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, $
      height=5.2, xmargin=[1.7,0.2], width=6.6, charsize=2

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='Redshift', xrange=[0.27,0.71], $
      yrange=[0.02,-0.28], ytitle=textoidl('\Delta<NUV-r>'), $
      color=keycolor, ytickinterval=0.1
    oplot, !x.crange, [0,0], line=1, thick=4, color=keycolor
    oploterror, [mean([0.25,0.5]),mean([0.5,0.75])], $  ; 30 kpc
      [-0.137,-0.149], [0.068,0.099], psym=-symcat(6,thick=8), line=0, $
      color=fsc_color('tan',101), errcolor=fsc_color('tan',101), $
      thick=8, errthick=8, symsize=3
    oploterror, [mean([0.25,0.5]),mean([0.5,0.75])]-0.03, $ ; 50 kpc
      [-0.095,-0.073], [0.047,0.060], psym=-symcat(15), line=5, $
      color=fsc_color('dodger blue',101), errcolor=fsc_color('dodger blue',101), $
      thick=8, errthick=8, symsize=3
    im_legend, ['r_{p}<30 h^{-1} kpc','r_{p}<50 h^{-1} kpc'], $
      /left, /top, box=0, charsize=1.7, $
      pspacing=1.9, line=[0,5], color=['tan','dodger blue'], $
      textcolor=keycolor, psym=[6,15], symsize=2, thick=8, symthick=8
    
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; GALEX Emphot QAplot
    pixscale = 1.5 ; [arcsec/pixel]
    stamp = (15.0*60)/pixscale   ; 15 arcmin=600 pixels
    xcen = 1850
    ycen = 2175

    im = mrdfits(datapath+'CDFS_04-nd-int.fits.gz',0)
    imsize = size(im,/dimension)
;   xsize = imsize[0] & xcen = xsize/2.0
;   ysize = imsize[1] & ycen = ysize/2.0
;   xaxis = (findgen(xsize)-xcen)          ;*pixscale ; [arcsec]
;   yaxis = (findgen(ysize)-ycen)          ;*pixscale ; [arcsec]
    im = im[xcen-stamp/2.0:xcen+stamp/2.0,$
      ycen-stamp/2.0:ycen+stamp/2.0]
    
    omin = 50
    omax = 250
;   beta1 = 5.0 & alpha1 = 10.0
;   beta2 = 5.0 & alpha2 = 10.0
    quant1 = [0.01,0.99]
    quant2 = [0.01,0.99]
    
    psfile = talkpath+'galex_emphot_diff.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, bits=24
    diff = mrdfits(datapath+'CDFS_04-nd_v1-diff.fits.gz',0)
    diff = diff[xcen-stamp/2.0:xcen+stamp/2.0,$
      ycen-stamp/2.0:ycen+stamp/2.0]
    qq = weighted_quantile(diff,quant=quant2)
    img = asinhscl(diff,/negative,alpha=alpha2,beta=beta2,$
      omin=omin,omax=omax,min=qq[0],max=qq[1])
    plotimage, img, /normal, position=pos, margin=0, xsty=5, ysty=5, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10)
;     imgxrange=minmax(xaxis), imgyrange=minmax(yaxis)
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

    psfile = talkpath+'galex_emphot_raw.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, bits=24
    qq = weighted_quantile(im,quant=quant1)
    img = asinhscl(im,/negative,alpha=alpha1,beta=beta1,$
      omin=omin,omax=omax,min=qq[0],max=qq[1])
    plotimage, img, /normal, position=pos, margin=0, xsty=5, ysty=5, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10)
;     imgxrange=minmax(xaxis), imgyrange=minmax(yaxis)
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; survey comparison plot (from Alison's code)
    psfile = talkpath+'survey_comp.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, $
      xmargin=[1.5,0.2], width=6.8, height=5.5

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    readcol,datapath+'survey_info.dat',survey,nobj,zmin,$
      zmax,zmed,area,mag,depth,idepth,v,logv,$
      format='a,l,f,f,f,f,a,f,f,f', /silent
    spec=[0,1,2,3,6,7,8,9,10]   ; spec-z surveys
    phot=[4,5]                  ; photo-z surveys

    logv = logv + alog10(0.7^3) ; h=1-->0.7
    xr = [5.9,8.7] + alog10(0.7^3) ; h=1-->0.7
    
; have symsize scale with depth
    symm=(idepth-10)/4.

; color code:  primus=red (filled square), photz=grey (filled triangle), specz=blue (filled diamond)
; don't worry about different redshift ranges for this plot

    plot,logv[spec],nobj[spec],/ylog,xtit=textoidl('log Volume (h_{70}^{-3} Mpc^3)'),$
;   plot,logv[spec],nobj[spec],/ylog,xtit=textoidl('log Volume (!8h_{70}^{-3}!3 Mpc^3)'),$
      ytit='Number of Redshifts',psym=4,xrange=xr,yr=[1e3,2e6],$
      symsize=1.0, color=keycolor, xsty=1, ysty=1

    if keyword_set(keynote) then col = fsc_color('goldenrod',101) else col = fsc_color('charcoal',10)
    for i=0,n_elements(spec)-1 do begin
       oplot,[logv[spec[i]],logv[spec[i]]],[nobj[spec[i]],nobj[spec[i]]],$
         psym=symcat(14),symsize=symm[spec[i]],color=col
    endfor
    for i=0,n_elements(phot)-1 do begin
       oplot,[logv[phot[i]],logv[phot[i]]],[nobj[phot[i]],nobj[phot[i]]],$
         psym=symcat(17),symsize=symm[phot[i]]-0.4,color=fsc_color('royal blue',10)
    endfor

; PRIMUS:
    oplot,[logv[11],logv[11]],[nobj[11],nobj[11]],psym=symcat(15),$
      symsize=symm[11]-0.8,color=fsc_color('firebrick',10)

    off=0.07   
    xyouts,logv[0]+0.02,nobj[0]+5e3,survey[0],charsize=1.6, color=keycolor
    xyouts,logv[1]+off,nobj[1]-1.4e3,survey[1],charsize=1.6, color=keycolor
    xyouts,logv[2]+off,nobj[2],survey[2],charsize=1.6, color=keycolor
    xyouts,logv[3]+off,nobj[3]-7e2,survey[3],charsize=1.6, color=keycolor

    xyouts,logv[4]-0.3,nobj[4]-7e4,'COSMOS-30',charsize=1.6, color=keycolor
;   xyouts,6.1+ alog10(0.7^3),nobj[4]+1e4,'COSMOS-30',charsize=1.6, color=keycolor
;   xyouts,6.1+ alog10(0.7^3),nobj[4]-3e4,'phot-z',charsize=1.6, color=keycolor
    xyouts,logv[5]+off,nobj[5]-1e3,survey[5],charsize=1.6, color=keycolor
    xyouts,logv[6]+off,nobj[6],survey[6],charsize=1.6, color=keycolor
    xyouts,logv[7]+off+0.01,nobj[7]-2e3,survey[7],charsize=1.6, color=keycolor
    xyouts,logv[8]+off,nobj[8],survey[8],charsize=1.6, color=keycolor
    xyouts,logv[9]+off,nobj[9],survey[9],charsize=1.6, color=keycolor
    xyouts,logv[10]+off,nobj[10],survey[10],charsize=1.6, color=keycolor

    xyouts,logv[11]-0.25,nobj[11]+4e4,'PRIMUS',charsize=2.2, color=keycolor

    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

; --------------------------------------------------
; plot three example PRIMUS spectra    
    psfile = talkpath+'example_spectra.ps'
    im_plotconfig,4,pos,psfile=psfile, xmargin=[1,0.3],$
      thick=5,charsize=1.9,height=2.6*[1,1,1], $
      keynote=keynote

    readcol, datapath+'blue.example.1.pair1.dat',wb2a,fb2a, /silent
    readcol, datapath+'blue.example.1.pair2.dat',wb2b,fb2b, /silent
    readcol, datapath+'red.example.1.pair1.dat',wr2a,fr2a, /silent
    readcol, datapath+'red.example.1.pair2.dat',wr2b,fr2b, /silent
    readcol, datapath+'qso.example.1.pair1.dat',wq2a,fq2a, /silent
    readcol, datapath+'qso.example.1.pair2.dat',wq2b,fq2b, /silent

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    xr = [4300,9400]
    
; Red galaxy
    plot,wr2a,fr2a,position=pos[*,0],xr=xr,psym=10,yr=[0.1,9],$
      xtickname=replicate(' ',5),yminor=1, color=keycolor, xsty=1, ysty=1
    oplot,wr2b,fr2b*1.2,color=fsc_color('firebrick',101),psym=10
    xyouts,7300,1.6,'red galaxy z=0.67', color=keycolor

; blue galaxy
    plot,wb2a,fb2a,position=pos[*,1],/noerase,xr=xr,psym=10,yr=[0.1,9],$
      xtickname=replicate(' ',5),yminor=1, color=keycolor, xsty=1, ysty=1, $
      ytitle=textoidl('Flux (10^{-17} erg s^{-1} cm^{-2} \AA^{-1})')
    oplot,wb2b,fb2b*1.13,color=fsc_color('firebrick',101),psym=10
    xyouts,7050,1.2,'blue galaxy z=0.62', color=keycolor
    xyouts,6110,7.5,'[OII]',charsize=1.4, color=keycolor
    xyouts,8200,6.5,'H'+textoidl('\beta')+'+[OIII]',charsize=1.4, color=keycolor

; QSO
    plot,wq2a,fq2a,position=pos[*,2],/noerase,xr=xr,psym=10,yr=[5,27],$
      yminor=1, color=keycolor, xsty=1, ysty=1, $
      xtitle=textoidl('Observed Wavelength (\AA)')
    oplot,wq2b,fq2b,color=fsc_color('firebrick',101),psym=10
    xyouts,7800,20,'AGN z=1.62', color=keycolor
    xyouts,5150,20,'[CIII]',charsize=1.4, color=keycolor
    xyouts,7600,13,'[MgII]',charsize=1.4, color=keycolor

    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; single-mask layout visualization
    psfile = talkpath+'one_mask_layout.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.2,0.3], keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')
    talk_qa_one_mask, pos=pos, keycolor=keycolor    

    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    

return
end

;
;    psfile = talkpath+'galex_emphot.ps'
;    im_plotconfig, 6, pos, psfile=psfile, keynote=keynote, bits=24, $
;      width=4.0, height=[4.0,4.0], yspace=0.05
;
;    omin = 50
;    omax = 250
;    beta = 5.0
;    alpha = 10.0
;    quant1 = [0.1,0.9]
;    quant2 = [0.1,0.9]
;    
;; raw image
;;   im = alog10(im)
;;   zz = zscale_range(im,0.1)
;;   img = bytscl(im,min=zz[0],max=zz[1],top=omax)
;;   img = bytscl(omax-img,min=omin,top=omax)
;
;    qq = weighted_quantile(im,quant=quant1)
;;   img = logscl(im,/negative,omin=omin,omax=omax,min=qq[0],max=qq[1])
;    img = asinhscl(im,/negative,alpha=alpha,beta=beta,$
;      omin=omin,omax=omax,min=qq[0],max=qq[1])
;
;    plotimage, img, /normal, position=pos[*,0], margin=0, imgxrange=minmax(xaxis), $
;      imgyrange=minmax(yaxis), xtickname=replicate(' ',10), $
;      ytickname=replicate(' ',10), xsty=5, ysty=5
;
;;; difference image
;;   im = alog10(diff>min(diff))
;;   zz = zscale_range(im,0.1)
;;   img = bytscl(im,min=zz[0],max=zz[1],top=omax)
;;   img = bytscl(omax-img,min=omin,top=omax)
;
;    qq = weighted_quantile(diff,quant=quant2)
;    img = asinhscl(diff,/negative,alpha=alpha,beta=beta,$
;      omin=omin,omax=omax,min=qq[0],max=qq[1])
;    
;    plotimage, img, /normal, position=pos[*,1], margin=0, imgxrange=minmax(xaxis), /noerase, $
;      imgyrange=minmax(yaxis), xtickname=replicate(' ',10), $
;      ytickname=replicate(' ',10), xsty=5, ysty=5
;    
;    im_plotconfig, psfile=psfile, keynote=keynote, /psclose
;
