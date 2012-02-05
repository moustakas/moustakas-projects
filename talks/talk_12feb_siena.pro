pro pie_allsurveys, pos, keynote=keynote, keycolor=keycolor, $
  zrange=zrange, rot=rot, thaxes=thaxes

; SDSS    
    zz = read_vagc_garching(/postlss)
    ww = where(zz.dec gt -2 and zz.dec lt 2)
    pie_plot, zz[ww].z, zz[ww].ra, psym=3, symsize=0.1, rlabel='Redshift', $
      rrange=zrange, thaxes=thaxes, axiscolor=keycolor, $
      position=pos, rotate=rot, color=keycolor

; zcosmos
    zz = mrdfits(getenv('PRIMUS_DATA')+'/photo/zcosmos/ZCOSMOS_VIMOS_BRIGHT_DR2_TABLE.fits.gz',1)
    ww = where(primus_zmerge_flag(zz.cc,/zcosmos))
    pie_plot, zz[ww].z, zz[ww].alpha_j2000, psym=3, symsize=0.1, $
      /over, rrange=zrange, color=im_color('forest green'), rotate=rot

; AGES
    zz = read_ages(/phot)
    ww = where(zz.main_flag)
    pie_plot, zz[ww].z, zz[ww].ra, psym=3, symsize=0.1, $
      /over, rrange=zrange, color=im_color('tomato'), rotate=rot

; DEEP2
    zz = mrdfits(deep2_path(/analysis)+'zcat.dr3.v1_0.uniq.fits.gz',1)
    ww = where(zz.zquality ge 3 and zz.ra/15 lt 15.0 and zz.ra gt 14.0) ; aegis field
    pie_plot, zz[ww].z, zz[ww].ra, psym=3, symsize=0.1, $
      /over, rrange=zrange, color=im_color('dodger blue'), rotate=rot

    ww = where(zz.zquality ge 3 and (zz.ra/15 gt 15.0 or zz.ra lt 14.0) and $ ; other fields field
      zz.z gt 0.7)
    pie_plot, zz[ww].z, zz[ww].ra, psym=3, symsize=0.1, $
      /over, rrange=zrange, color=im_color('dodger blue'), rotate=rot

return
end    

pro talk_12feb_siena, keynote=keynote, noevol=noevol
; jm12jan28ucsd - miscellaneous plots for my 2012 Feb talk at Siena

    common com_siena, model, mstar, isedfit
    
    mfpath = mf_path()
    mfspath = mf_path(/mfs)
    isedpath = mf_path(/isedfit)
    mzpath = mz_path()
    
    datapath = getenv('IM_RESEARCH_DIR')+'/talks/2012/12feb_siena/'
    if keyword_set(keynote) then talkpath = datapath+'keynote/' else $
      talkpath = datapath

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('black')

; ---------------------------------------------------------------------------
; evolution of the sf MF in one panel
    zbins = mf_zbins(nzbins)
    xrange = [8.7,12.2]
    yrange = [-5,-1.5]

    ss = read_mf_vmax('sdss',noevol=noevol,/final,/log,/active)
    pp = read_mf_vmax(noevol=noevol,/avgfield,/final,/log,/active)

    if keyword_set(keynote) then col1 = 'white' else col1 = 'black'
    col = [col1,'coral','steel blue','forest green','orchid','orange','tan']
    for ii = -1, nzbins-1 do begin
       psfile = talkpath+'sf_evol_zbin'+strtrim(ii+1,2)+'.ps'
       im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, height=5.5

       plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
         ytitle=mfplot_phititle(), $
         xtickinterval=1, xminor=4, color=keycolor, ytickinterval=1

       good = where(ss.limit eq 1 and ss.number ge 3)
       phimass = ss.mass[good]
       polyfill, [phimass,reverse(phimass)],[ss.phi_lower_stat[good]-0.03,$
         reverse(ss.phi_upper_stat[good]+0.03)], $
         /data, color=im_color(col[0],10), noclip=0, /fill

       if ii ge 0 then for iz = 0, ii do begin
          good = where(pp[iz].limit eq 1 and pp[iz].number ge 3)
          phimass = pp[iz].mass[good]
          polyfill, [phimass,reverse(phimass)],[pp[iz].phi_lower_stat[good]-0.03,$
            reverse(pp[iz].phi_upper_stat[good]+0.03)], $
            /data, color=im_color(col[iz+1]), noclip=0, /fill
       endfor
       im_plotconfig, /psclose, psfile=psfile, keynote=keynote
    endfor

stop    
    
; ---------------------------------------------------------------------------
; evolution of the qq MF in one panel
    zbins = mf_zbins(nzbins)
    xrange = [8.7,12.2]
    yrange = [-5,-2]

    ss = read_mf_vmax('sdss',noevol=noevol,/final,/log,/quiescent)
    pp = read_mf_vmax(noevol=noevol,/avgfield,/final,/log,/quiescent)

    if keyword_set(keynote) then col1 = 'white' else col1 = 'black'
    col = [col1,'coral','steel blue','forest green','orchid','orange','tan']
    for ii = -1, nzbins-1 do begin
       psfile = talkpath+'qq_evol_zbin'+strtrim(ii+1,2)+'.ps'
       im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, height=5.5

       plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
         ytitle=mfplot_phititle(), $
         xtickinterval=1, xminor=4, color=keycolor, ytickinterval=1

       good = where(ss.limit eq 1 and ss.number ge 3)
       phimass = ss.mass[good]
       polyfill, [phimass,reverse(phimass)],[ss.phi_lower_stat[good]-0.03,$
         reverse(ss.phi_upper_stat[good]+0.03)], $
         /data, color=im_color(col[0],10), noclip=0, /fill

       if ii ge 0 then for iz = 0, ii do begin
          good = where(pp[iz].limit eq 1 and pp[iz].number ge 3)
          phimass = pp[iz].mass[good]
          polyfill, [phimass,reverse(phimass)],[pp[iz].phi_lower_stat[good]-0.03,$
            reverse(pp[iz].phi_upper_stat[good]+0.03)], $
            /data, color=im_color(col[iz+1]), noclip=0, /fill
       endfor
       im_plotconfig, /psclose, psfile=psfile, keynote=keynote
    endfor


stop    

; ---------------------------------------------------------------------------
; sdss+primus - quiescent and SF stellar mass function in three 
; redshift bins 
    supergrid = 1
    xrange = [8.7,12.2]
    yrange = [-5,-1]

    zbins_sdss = mf_zbins(/sdss)
    zbins = mf_zbins(nzbins)
    
    maxis1 = range(xrange[0]+0.2,12.5,100)
    scolor = 'black' & sline = 5

    qqcolor = 'firebrick'  & qqpolycolor = 'coral'
    qqfitcolor = 'dark red' & qqline = 0

    sfcolor = 'navy' & sfpolycolor = 'powder blue'
    sffitcolor = 'navy'     & sfline = 5
    
    psfile = talkpath+'mf_qq_sf_evol.ps'
    im_plotconfig, 3, pos, psfile=psfile, keynote=keynote, $
      ymargin=[1.0,1.5]

; -------------------------
; SDSS
    plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle='', $
      ytitle=textoidl('log (Number Density) (Mpc^{-3} dex^{-1})'), $
      xtickinterval=1, xminor=4, color=keycolor, ytickinterval=1
    im_legend, '<Age> = '+strtrim(string(getage(zbins_sdss.zbin),format='(F12.1)'),2)+' Gyr', $
      /right, /top, box=0, margin=0, charsize=1.5, textcolor=keycolor

    qq = read_mf_vmax('sdss',noevol=noevol,/final,/log,/quiescent)
    sf = read_mf_vmax('sdss',noevol=noevol,/final,/log,/active)
    qqfit = read_mf_vmax('sdss',/bestfit,noevol=noevol,/final,/quiescent)
    sffit = read_mf_vmax('sdss',/bestfit,noevol=noevol,/final,/active)

;   djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=qqfit)), $
;     line=qqline, thick=5, color=im_color(qqfitcolor)
;   djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=sffit)), $
;     line=sfline, thick=5, color=im_color(sffitcolor)

; quiescent    
    good = where(qq.limit eq 1 and qq.number ge 3)
    phimass = qq.mass[good]
    polyfill, [phimass,reverse(phimass)],[qq.phi_lower_stat[good]-0.03,$
      reverse(qq.phi_upper_stat[good]+0.03)], $
      /data, color=im_color(qqpolycolor), noclip=0, /fill

; star-forming        
    good = where(sf.limit eq 1 and sf.number ge 3)
    phimass = sf.mass[good]
    polyfill, [phimass,reverse(phimass)],[sf.phi_lower_stat[good]-0.03,$
      reverse(sf.phi_upper_stat[good]+0.03)], $
      /data, color=im_color(sfpolycolor), noclip=0, /fill

; -------------------------
; PRIMUS
    qq = read_mf_vmax(noevol=noevol,/avgfield,/final,/log,/quiescent)
    sf = read_mf_vmax(noevol=noevol,/avgfield,/final,/log,/active)

; zbin2
    iz = 1
    plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, ytickname=replicate(' ',10), $
      xtitle=textoidl('log (Stellar Mass) (M'+sunsymbol()+')'), $
      ytitle='', xtickinterval=1, xminor=4, color=keycolor
    im_legend, '<Age> = '+strtrim(string(getage(zbins[iz].zbin),format='(F12.1)'),2)+' Gyr', $
      /right, /top, box=0, margin=0, charsize=1.5, textcolor=keycolor

    good = where(sf[iz].limit eq 1 and sf[iz].number ge 3)
    phimass = sf[iz].mass[good]
    polyfill, [phimass,reverse(phimass)],[sf[iz].phi_lower_stat[good],$
      reverse(sf[iz].phi_upper_stat[good])], $
      /data, color=im_color(sfpolycolor), noclip=0, /fill
    
    good = where(qq[iz].limit eq 1 and qq[iz].number ge 3)
    phimass = qq[iz].mass[good]
    polyfill, [phimass,reverse(phimass)],[qq[iz].phi_lower_stat[good],$
      reverse(qq[iz].phi_upper_stat[good])], $
      /data, color=im_color(qqpolycolor), noclip=0, /fill

; zbin5
    iz = 4
    plot, [0], [0], /nodata, /noerase, position=pos[*,2], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, ytickname=replicate(' ',10), $
      xtitle='', ytitle='', xtickinterval=1, xminor=4, color=keycolor
    im_legend, '<Age> = '+strtrim(string(getage(zbins[iz].zbin),format='(F12.1)'),2)+' Gyr', $
      /right, /top, box=0, margin=0, charsize=1.5, textcolor=keycolor

    good = where(sf[iz].limit eq 1 and sf[iz].number ge 3)
    phimass = sf[iz].mass[good]
    polyfill, [phimass,reverse(phimass)],[sf[iz].phi_lower_stat[good],$
      reverse(sf[iz].phi_upper_stat[good])], $
      /data, color=im_color(sfpolycolor), noclip=0, /fill
    
    good = where(qq[iz].limit eq 1 and qq[iz].number ge 3)
    phimass = qq[iz].mass[good]
    polyfill, [phimass,reverse(phimass)],[qq[iz].phi_lower_stat[good],$
      reverse(qq[iz].phi_upper_stat[good])], $
      /data, color=im_color(qqpolycolor), noclip=0, /fill

    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

stop    


; ---------------------------------------------------------------------------
; sdss - compare quiescent and active samples
    psfile = talkpath+'mf_sdss.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.8, height=6.0, charsize=2.0, keynote=keynote

    xrange = [8.8,12.2]
    yrange = [-6,-2]
    maxis1 = range(xrange[0]+0.15,12.5,100)

    subsample = ['active','quiescent']
    mfcolor = ['powder blue','tomato']

    psymsize = [0.9,1.1,0.9]*1.2
    mfpsym = [16,14,15]

    plot, [0], [0], /nodata, position=pos, xsty=9, ysty=9, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytitle=mfplot_phititle(), xtickinterval=1, color=keycolor
;   im_legend, ['All','Quiescent','Star-Forming'], /left, /bottom, box=0, $
;     color=reverse(mfcolor), psym=reverse(mfpsym), symsize=[1.7,2.1,1.9], $
;     symthick=8, spacing=2.4, charsize=1.8
    
    for jj = 0, n_elements(subsample)-1 do begin
       case jj of
          0: begin
             mfdata = read_mf_vmax('sdss',/final,/active,/log)
          end
          1: begin
             mfdata = read_mf_vmax('sdss',/final,/quiescent,/log)
          end
       endcase
       
       good = where(mfdata.limit eq 1)
       few = where(mfdata.limit eq 1 and mfdata.number le 3,nfew)

       phimass = mfdata.mass[good]
       polyfill, [phimass,reverse(phimass)],[mfdata.phi_lower_stat[good],$
         reverse(mfdata.phi_upper_stat[good])], $
         /data, color=im_color(mfcolor[jj]), noclip=0, /fill
       
;      oploterror, mfdata.mass[good], mfdata.phi[good], mfdata.phierr_upper[good], $
;        psym=symcat(mfpsym[jj],thick=5), symsize=psymsize[jj], errthick=5, $
;        color=im_color(mfcolor[jj],100+jj), errcolor=im_color(mfcolor[jj],100+jj), /hibar, /nohat
;      oploterror, mfdata.mass[good], mfdata.phi[good], mfdata.phierr_lower[good], psym=3, $
;        symsize=psymsize[jj], errthick=5, color=im_color(mfcolor[jj],100+jj), $
;        errcolor=im_color(mfcolor[jj],100+jj), /lobar, /nohat
    endfor

    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

stop    
    
; --------------------------------------------------
; pie diagram for sdss + primus
    zrange = [0,1.2]
    rot = 270.0
    thaxes = 90.0
    
    psfile = talkpath+'pie_sdss_primus.ps'
    im_plotconfig, 0, pos, psfile=psfile, width=7.0, $
      height=7.0, charsize=2, keynote=keynote
    zz = read_vagc_garching(/postlss)
    ww = where(zz.dec gt -2 and zz.dec lt 2)
    pie_plot, zz[ww].z, zz[ww].ra, psym=3, symsize=0.1, rlabel='Redshift', $
      rrange=zrange, thaxes=thaxes, axiscolor=keycolor, $
      color=keycolor, position=pos, rotate=rot

    field = ['cdfs','xmm','es1','deep2_02hr','deep2_23hr','cosmos','dls']
    for ii = 0, n_elements(field)-1 do begin
       zz = primus_read_zerod(field=field[ii])
       gal = where(strtrim(zz.zprimus_class,2) eq 'GALAXY')
       pie_plot, zz[gal].zprimus, zz[gal].ra, psym=3, symsize=0.2, $
         /over, rrange=zrange, color='orange', rotate=rot
    endfor
    im_plotconfig, psfile=psfile, /psclose, keynote=keynote

; --------------------------------------------------
; survey comparison plot (from Alison's code)
    psfile = talkpath+'survey_comp.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, $
      xmargin=[1.5,0.2], width=6.8, height=5.5

    readcol,datapath+'survey_info.dat',survey,nobj,zmin,$
      zmax,zmed,area,mag,depth,idepth,v,logv, comment=';', $
      format='a,l,f,f,f,f,a,f,f,f', /silent
;   spec=lindgen(9)
    spec=[0,3,6,9]   ; spec-z surveys
    if keyword_set(keynote) then col = 'white' else col = 'black'
    color = ['dodger blue','forest green','tomato',col]
;   spec=[0,1,2,3,6,7,8,9,10]   ; spec-z surveys
    phot=[4,5]                  ; photo-z surveys

    logv = logv + alog10(0.7^3) ; h=1-->0.7
    xr = [5.9,8.7] + alog10(0.7^3) ; h=1-->0.7
    
; have symsize scale with depth
    symm=(idepth-10)/4.

; color code:  primus=red (filled square), photz=grey (filled triangle), specz=blue (filled diamond)
; don't worry about different redshift ranges for this plot

    plot,logv[spec],nobj[spec],/ylog,xtit=textoidl('log (Comoving Volume) (Mpc^3)'),$
;   plot,logv[spec],nobj[spec],/ylog,xtit=textoidl('log (Comoving Volume) (h_{70}^{-3} Mpc^3)'),$
;   plot,logv[spec],nobj[spec],/ylog,xtit=textoidl('log Volume (!8h_{70}^{-3}!3 Mpc^3)'),$
      ytit='Number of Galaxies',psym=4,xrange=xr,yr=[1e3,2e6],$
      symsize=1.0, color=keycolor, xsty=1, ysty=1

;   if keyword_set(keynote) then col = im_color('goldenrod') else col = im_color('black')
    for i=0,n_elements(spec)-1 do begin
       oplot,[logv[spec[i]],logv[spec[i]]],[nobj[spec[i]],nobj[spec[i]]],$
         psym=symcat(14),symsize=symm[spec[i]],color=im_color(color[i])
    endfor
;   for i=0,n_elements(phot)-1 do begin
;      oplot,[logv[phot[i]],logv[phot[i]]],[nobj[phot[i]],nobj[phot[i]]],$
;        psym=symcat(17),symsize=symm[phot[i]]-0.4,color=im_color('royal blue',10)
;   endfor

; PRIMUS:
    ww = where(strmatch(survey,'*PRIMUS*',/fold))
    oplot,[logv[ww],logv[ww]],[nobj[ww],nobj[ww]],psym=symcat(15),$
      symsize=symm[ww]-0.8,color=im_color('orange')

    off=0.07
    csize = 1.4
    xyouts,logv[0]+0.02,nobj[0]+5e3,survey[0],charsize=csize, color=keycolor
;   xyouts,logv[1]+off,nobj[1]-1.4e3,survey[1],charsize=csize, color=keycolor
;   xyouts,logv[2]+off,nobj[2],survey[2],charsize=csize, color=keycolor
    xyouts,logv[3]+off,nobj[3]-7e2,survey[3],charsize=csize, color=keycolor

;   xyouts,logv[4]-0.3,nobj[4]-7e4,'COSMOS-30',charsize=csize, color=keycolor
;   xyouts,6.1+ alog10(0.7^3),nobj[4]+1e4,'COSMOS-30',charsize=csize, color=keycolor
;   xyouts,6.1+ alog10(0.7^3),nobj[4]-3e4,'phot-z',charsize=csize, color=keycolor
;   xyouts,logv[5]+off,nobj[5]-1e3,survey[5],charsize=csize, color=keycolor
    xyouts,logv[6]+off,nobj[6],survey[6],charsize=csize, color=keycolor
;   xyouts,logv[7]+off+0.01,nobj[7]-2e3,survey[7],charsize=csize, color=keycolor
;   xyouts,logv[8]+off,nobj[8],survey[8],charsize=csize, color=keycolor
    xyouts,logv[9]+off,nobj[9],survey[9],charsize=csize, color=keycolor
;   xyouts,logv[10]+off,nobj[10],survey[10],charsize=csize, color=keycolor

    xyouts,logv[11]-0.25,nobj[11]+4e4,'PRIMUS',charsize=2.2, color=keycolor

    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    

; --------------------------------------------------
; SDSS pie diagram
    psfile = talkpath+'pie_sdss.ps'
    im_plotconfig, 0, pos, psfile=psfile, width=7.0, $
      height=7.0, charsize=2, keynote=keynote

    zz = read_vagc_garching(/postlss)
    ww = where(zz.dec gt -2 and zz.dec lt 2)
    pie_plot, zz[ww].z, zz[ww].ra, psym=3, symsize=0.1, rlabel='Redshift', $
      rrange=[0,0.17], thaxes=90.0, axiscolor=keycolor, $
      color=keycolor, position=pos, rotate=270

    im_plotconfig, psfile=psfile, /psclose, keynote=keynote

stop    
    
; --------------------------------------------------
; pie diagram for all surveys, with and without primus
    zrange = [0,1.5]
    rot = 270.0
    thaxes = 90.0
    
    psfile = talkpath+'pie_allsurveys_noprimus.ps'
    im_plotconfig, 0, pos, psfile=psfile, width=7.0, $
      height=7.0, charsize=2, keynote=keynote
    pie_allsurveys, pos, keynote=keynote, keycolor=keycolor, $
      zrange=zrange, rot=rot, thaxes=thaxes
    im_plotconfig, psfile=psfile, /psclose, keynote=keynote
    
; now add primus
    psfile = talkpath+'pie_allsurveys.ps'
    im_plotconfig, 0, pos, psfile=psfile, width=7.0, $
      height=7.0, charsize=2, keynote=keynote
    pie_allsurveys, pos, keynote=keynote, keycolor=keycolor, $
      zrange=zrange, rot=rot, thaxes=thaxes
    
    field = ['cdfs','xmm','es1','deep2_02hr','deep2_23hr','cosmos','dls']
    for ii = 0, n_elements(field)-1 do begin
       zz = primus_read_zerod(field=field[ii])
       gal = where(strtrim(zz.zprimus_class,2) eq 'GALAXY')
       pie_plot, zz[gal].zprimus, zz[gal].ra, psym=3, symsize=0.2, $
         /over, rrange=zrange, color='orange', rotate=rot
    endfor
    im_plotconfig, psfile=psfile, /psclose, keynote=keynote

stop    
    
; --------------------------------------------------
; number and mass-density evolution in bins of mass
    zbins = [mf_zbins(/sdss),mf_zbins()]
    qq = [read_mf_vmax('sdss',/rho,/final,/log,/quiescent),$
      read_mf_vmax(/rho,/avgfield,/final,/log,/quiescent)]
    sf = [read_mf_vmax('sdss',/rho,/final,/log,/active),$
      read_mf_vmax(/rho,/avgfield,/final,/log,/active)]

    drory = [[get_mf_literature(/drory,/quiescent)],[get_mf_literature(/drory,/active)]]
    
    tot = [[qq],[sf]]
    psym = [14,16]
    psym_open = [4,9]
    color = ['coral','powder blue']
    symsize = [2.2,2.0]*1.3

    psfile = talkpath+'z_vs_numden_rhoden_bymass.ps'
    im_plotconfig, 3, pos, psfile=psfile, xmargin=[1.4,0.4], width=[3,3,3], $
      charsize=1.7, keynote=keynote

    xrange = [0.001,1.01]
    rhorange = [6.6,8.3]
    numrange = [-4.5,-2.0]

; 10-10.5
    plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='Redshift', $
      ytitle=textoidl('log (Number Density) (Mpc^{3})'), color=keycolor
    im_legend, '10^{10}-10^{10.5} M'+sunsymbol(), /right, /top, box=0, margin=0, textcolor=keycolor
;   im_legend, 'log (M/M'+sunsymbol()+')=10-10.5', /right, /top, box=0

    zcut = [0.65,0.8]
    for ii = 0, 1 do begin
;      gd = where(tot[*,ii].num_10_105 gt -999.0,ngd)
       gd = where(tot[*,ii].num_10_105 gt -999.0 and zbins.zbin lt zcut[ii],ngd)
       if (ngd ne 0) then begin
          oploterror, zbins[gd].zbin, tot[gd,ii].num_10_105, zbins[gd].zsigma*0, tot[gd,ii].numerr_cv_10_105, $
            psym=-symcat(psym[ii]), symsize=symsize[ii], color=im_color(color[ii]), errcolor=im_color(color[ii]), $
            errstyle=0, errthick=8, /nohat
       endif
       djs_oplot, drory[*,ii].z, drory[*,ii].num_10_105, psym=-symcat(psym_open[ii]), $
         color=im_color(color[ii]), line=1
    endfor 

; 10.5-11
    plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='Redshift', $
      ytitle='', ytickname=replicate(' ',10), color=keycolor
    im_legend, '10^{10.5}-10^{11} M'+sunsymbol(), /right, /top, box=0, margin=0, textcolor=keycolor
;   im_legend, 'log (M/M'+sunsymbol()+')=10.5-11', /right, /top, box=0

    zcut = [0.8,1.0]
    for ii = 0, 1 do begin
       splog, 'Hack!'
;      gd = where(tot[*,ii].num_105_11 gt -999.0,ngd)
       gd = where(tot[*,ii].num_10_105 gt -999.0 and zbins.zbin lt zcut[ii],ngd)
       if (ngd ne 0) then begin
          oploterror, zbins[gd].zbin, tot[gd,ii].num_105_11, zbins[gd].zsigma*0, tot[gd,ii].numerr_cv_105_11, $
            psym=-symcat(psym[ii]), symsize=symsize[ii], color=im_color(color[ii]), errcolor=im_color(color[ii]), $
            errstyle=0, errthick=8, /nohat
       endif
       djs_oplot, drory[*,ii].z, drory[*,ii].num_105_11, psym=-symcat(psym_open[ii]), $
         color=im_color(color[ii]), line=1
    endfor 

; 11-11.5
    plot, [0], [0], /nodata, /noerase, position=pos[*,2], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='Redshift', $
      ytickname=replicate(' ',10), ytitle='', color=keycolor
    im_legend, '10^{11}-10^{11.5} M'+sunsymbol(), /right, /top, box=0, margin=0, textcolor=keycolor
;   im_legend, 'log (M/M'+sunsymbol()+')=11-11.5', /right, /top, box=0

    zcut = [0.8,0.8]
    for ii = 0, 1 do begin
;      gd = where(tot[*,ii].num_11_115 gt -999.0,ngd)
       gd = where(tot[*,ii].num_11_115 gt -999.0  and zbins.zbin lt zcut[ii],ngd)
       if (ngd ne 0) then begin
          oploterror, zbins[gd].zbin, tot[gd,ii].num_11_115, zbins[gd].zsigma*0, tot[gd,ii].numerr_cv_11_115, $
            psym=-symcat(psym[ii]), symsize=symsize[ii], color=im_color(color[ii]), errcolor=im_color(color[ii]), $
            errstyle=0, errthick=8, /nohat
       endif
       djs_oplot, drory[*,ii].z, drory[*,ii].num_11_115, psym=-symcat(psym_open[ii]), $
         color=im_color(color[ii]), line=1
    endfor 

    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

stop    


; ---------------------------------------------------------------------------
; redshift example
    nzz = 11
    zz = range(0,1,nzz)
    doff = 1.1
    off = findgen(nzz)*doff
    xrange = [3600,1E4]
    
    ss = im_read_bc03(age=10.0,/silent,minwave=xrange[0],maxwave=xrange[1])
    
    psfile = talkpath+'redshift.ps'
    im_plotconfig, 8, pos, psfile=psfile, keynote=keynote
    loadct, 13, /silent
    
    plot, [0], [0], /nodata, xsty=9, ysty=5, position=pos, $
      xrange=xrange, yrange=[-doff,max(off)], $
      ytitle=ytitle1, color=keycolor, $
      xtickinterval=2000, xtitle=textoidl('Rest-Frame Wavelength (\AA)')
    for iz = 0, nzz-1 do begin
       ww = ss.wave*(1+zz[iz])
       ff = ss.flux/interpol(ss.flux,ss.wave,5500)+off[iz]-1
       inrange = where(ww ge xrange[0] and ww le xrange[1],npix)
       
       cc = 255*findgen(npix)/npix
       djs_oplot, ww[inrange], ff[inrange], psym=10, color=im_color(cc)
;      for ii = 0, npix-1 do oplot, [ww[inrange[ii]]], [ff[inrange[ii]]], $
;        psym=6, symsize=0.05, color=cc[ii]
    endfor
;   loadct, 0
    
    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

; ---------------------------------------------------------------------------
; sdss+primus - quiescent and SF stellar mass function in three 
; redshift bins 
    supergrid = 1
    xrange = [8.7,12.2]
    yrange = [-5,-1]

    zbins_sdss = mf_zbins(/sdss)
    zbins = mf_zbins(nzbins)
    
    maxis1 = range(xrange[0]+0.2,12.5,100)
    scolor = 'black' & sline = 5

    qqcolor = 'firebrick'  & qqpolycolor = 'coral'
    qqfitcolor = 'dark red' & qqline = 0

    sfcolor = 'navy' & sfpolycolor = 'powder blue'
    sffitcolor = 'navy'     & sfline = 5
    
    psfile = talkpath+'mf_qq_sf_evol.ps'
    im_plotconfig, 3, pos, psfile=psfile, keynote=keynote, $
      ymargin=[1.0,1.5]

; -------------------------
; SDSS
    plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle='', $
      ytitle=textoidl('log (Number Density) (Mpc^{-3} dex^{-1})'), $
      xtickinterval=1, xminor=4, color=keycolor, ytickinterval=1
    im_legend, '<Age> = '+strtrim(string(getage(zbins_sdss.zbin),format='(F12.1)'),2)+' Gyr', $
      /right, /top, box=0, margin=0, charsize=1.5, textcolor=keycolor

    qq = read_mf_vmax('sdss',noevol=noevol,/final,/log,/quiescent)
    sf = read_mf_vmax('sdss',noevol=noevol,/final,/log,/active)
    qqfit = read_mf_vmax('sdss',/bestfit,noevol=noevol,/final,/quiescent)
    sffit = read_mf_vmax('sdss',/bestfit,noevol=noevol,/final,/active)

;   djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=qqfit)), $
;     line=qqline, thick=5, color=im_color(qqfitcolor)
;   djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=sffit)), $
;     line=sfline, thick=5, color=im_color(sffitcolor)

; quiescent    
    good = where(qq.limit eq 1 and qq.number ge 3)
    phimass = qq.mass[good]
    polyfill, [phimass,reverse(phimass)],[qq.phi_lower_stat[good]-0.03,$
      reverse(qq.phi_upper_stat[good]+0.03)], $
      /data, color=im_color(qqpolycolor), noclip=0, /fill

; star-forming        
    good = where(sf.limit eq 1 and sf.number ge 3)
    phimass = sf.mass[good]
    polyfill, [phimass,reverse(phimass)],[sf.phi_lower_stat[good]-0.03,$
      reverse(sf.phi_upper_stat[good]+0.03)], $
      /data, color=im_color(sfpolycolor), noclip=0, /fill

; -------------------------
; PRIMUS
    qq = read_mf_vmax(noevol=noevol,/avgfield,/final,/log,/quiescent)
    sf = read_mf_vmax(noevol=noevol,/avgfield,/final,/log,/active)

; zbin2
    iz = 1
    plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, ytickname=replicate(' ',10), $
      xtitle=textoidl('log (Stellar Mass) (M'+sunsymbol()+')'), $
      ytitle='', xtickinterval=1, xminor=4, color=keycolor
    im_legend, '<Age> = '+strtrim(string(getage(zbins[iz].zbin),format='(F12.1)'),2)+' Gyr', $
      /right, /top, box=0, margin=0, charsize=1.5, textcolor=keycolor

    good = where(sf[iz].limit eq 1 and sf[iz].number ge 3)
    phimass = sf[iz].mass[good]
    polyfill, [phimass,reverse(phimass)],[sf[iz].phi_lower_stat[good],$
      reverse(sf[iz].phi_upper_stat[good])], $
      /data, color=im_color(sfpolycolor), noclip=0, /fill
    
    good = where(qq[iz].limit eq 1 and qq[iz].number ge 3)
    phimass = qq[iz].mass[good]
    polyfill, [phimass,reverse(phimass)],[qq[iz].phi_lower_stat[good],$
      reverse(qq[iz].phi_upper_stat[good])], $
      /data, color=im_color(qqpolycolor), noclip=0, /fill

; zbin5
    iz = 4
    plot, [0], [0], /nodata, /noerase, position=pos[*,2], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, ytickname=replicate(' ',10), $
      xtitle='', ytitle='', xtickinterval=1, xminor=4, color=keycolor
    im_legend, '<Age> = '+strtrim(string(getage(zbins[iz].zbin),format='(F12.1)'),2)+' Gyr', $
      /right, /top, box=0, margin=0, charsize=1.5, textcolor=keycolor

    good = where(sf[iz].limit eq 1 and sf[iz].number ge 3)
    phimass = sf[iz].mass[good]
    polyfill, [phimass,reverse(phimass)],[sf[iz].phi_lower_stat[good],$
      reverse(sf[iz].phi_upper_stat[good])], $
      /data, color=im_color(sfpolycolor), noclip=0, /fill
    
    good = where(qq[iz].limit eq 1 and qq[iz].number ge 3)
    phimass = qq[iz].mass[good]
    polyfill, [phimass,reverse(phimass)],[qq[iz].phi_lower_stat[good],$
      reverse(qq[iz].phi_upper_stat[good])], $
      /data, color=im_color(qqpolycolor), noclip=0, /fill

    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

stop    
    
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
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, charsize=1.7, $
      width=7.0, height=5.0

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

    zrange = [-0.05,1.2]
    xtitle1 = 'Redshift from High-Resolution Spectroscopy'
    ytitle1 = 'PRIMUS Redshift'

    zaxis = findgen((2.0-(-1.0))/0.01)*0.01+(-1.0)
    dzcut_03 = 0.03*(1.0+zaxis)
    dzcut_10 = 0.10*(1.0+zaxis)

; Q=4
    plot, [0], [0], /nodata, xsty=1, ysty=1, position=pos[*,0], $
      xrange=zrange, yrange=zrange, xtitle=xtitle1, ytitle=ytitle1, $
      color=keycolor
    best = where(zconf eq 4,nbest)
    for ii = 0, nbest-1 do begin
       this = where(field[best[ii]] eq strtrim(prefs.field,2))
       plots, zspec[best[ii]], zprimus[best[ii]], psym=symcat(prefs[this].psym,thick=6), $
         color=im_color(prefs[this].color,101), symsize=prefs[this].symsize
    endfor
    oplot, zaxis, zaxis, line=0, thick=5, color=keycolor
    oplot, zaxis, zaxis-dzcut_03, line=5, thick=5, color=keycolor
    oplot, zaxis, zaxis+dzcut_03, line=5, thick=5, color=keycolor
;   oplot, zaxis, zaxis-dzcut_10, line=1, thick=5, color=keycolor
;   oplot, zaxis, zaxis+dzcut_10, line=1, thick=5, color=keycolor
    
;    stats = primus_zvz_stats(zspec[best],zprimus[best])
;    label = [$
;      'F(>0.03)='+strtrim(string(100.0*stats.frac_bad_03,format='(F12.1)'),2)+'%',$
;      'F(>0.10)='+strtrim(string(100.0*stats.frac_bad_10,format='(F12.1)'),2)+'%',$
;      '\sigma_{\Delta'+'z/(1+z)}='+strtrim(string(stats.mad,format='(F12.4)'),2)]
;    legend, textoidl(label), /left, /top, box=0, charsize=1.4, $
;      margin=0, textcolor=keycolor
;;   xyouts, 0.43, 0.82, 'Q=4', align=0.5, /data, color=keycolor
    
; legend
    im_legend, prefs.nicefield, textcolor=keycolor, $
      psym=prefs.psym, color=strtrim(prefs.color,2), $
      /right, /bottom, box=0, margin=0, charsize=1.4, $
      symsize=prefs.symsize*2.3, symthick=6
    
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; SED-fitting example
    isedpath = mf_path(/isedfit)
    isedfit_sfhgrid_dir = mf_path(/montegrids)

    paramfile = isedpath+'cfhtls_xmm_supergrid01_isedfit.par'
    sfhgrid_paramfile = getenv('PRIMUS_DIR')+'/pro/science/mf/mf_sfhgrid.par'

    ndraw = isedfit_ndraw() ; number of random draws

; pick the galaxy
    pp = read_mf_ubersample('cfhtls_xmm')
    ra = 34.776369D & dec = -5.3947370D ; picked this one by eye
    spherematch, pp.ra, pp.dec, ra, dec, 1D/3600.0, m1, m2
    index = m1[0]
    galaxy = hogg_iau_name(ra,dec,'CFHTLS-XMM')


    if (n_elements(model) eq 0) then begin
       model = isedfit_restore(paramfile,isedfit,iopath=isedpath,$
         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,index=index)
    endif

    if (n_elements(mstar) eq 0L) then begin
       mstar = isedfit_reconstruct_posterior(paramfile,post=post,$
         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath,$
;        age=age,Z=Z,tau=tau,sfr0=sfr0,b100=b100,av=av,sfrage=sfrage,$
         index=index)
    endif

    filters = strtrim(get_mf_filters('cfhtls_xmm',nice_filte=nice_filters),2)
    filtinfo = im_filterspecs(filterlist=filters)

    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    ytitle1 = textoidl('AB Magnitude')

    yrange = [28.5,18.5]
    xrange1 = [1000.0,6E4]
    ticks1 = [1000,4000,10000,40000]
    
    psfile = talkpath+'isedfit_example.ps'
    im_plotconfig, 8, pos, psfile=psfile, ymargin=[1.0,1.1], keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=9, xtitle=xtitle1, ytitle=ytitle1, /xlog, $
      position=pos, color=keycolor, $
      xtickformat='(I0)', xticks=n_elements(ticks1)-1, xtickv=ticks1
    oplot, model.wave, model.flux, line=0, color=keycolor; color='grey'

; overplot the filter names
;   nice_filters = repstr(repstr(nice_filters,'ch1','[3.6]'),'ch2','[4.5]')
;   for ii = 0, n_elements(filters)-1 do xyouts, filtinfo[ii].weff, $
;     28.0, nice_filters[ii], /data, align=0.5, color=keycolor, $
;     charsize=1.4
    xyouts, mean(k_lambda_eff(filterlist=galex_filterlist())), 28.0, $
      'ultraviolet', /data, align=0.5, color=im_color('powder blue'), charsize=1.4
    xyouts, mean(k_lambda_eff(filterlist=sdss_filterlist())), 28.0, $
      'optical', /data, align=0.5, color=im_color('forest green'), charsize=1.4
    xyouts, mean(k_lambda_eff(filterlist=(irac_filterlist())[0:1])), 28.0, $
      'mid-infrared', /data, align=0.5, color=im_color('tan'), charsize=1.4
    
; overplot the observed and model photometry
    used = where((isedfit.maggies gt 0.0) and $ ; used in the fitting
      (isedfit.ivarmaggies gt 0.0),nused)
    notused = where((isedfit.maggies gt 0.0) and $ ; not used in the fitting
      (isedfit.ivarmaggies eq 0.0),nnotused)
    nodata = where((isedfit.maggies eq 0.0) and $ ; no measurement
      (isedfit.ivarmaggies eq 0.0),nnodata)

;   djs_oplot, filtinfo[used].weff, -2.5*alog10(isedfit[used].bestmaggies), $
;     psym=symcat(6,thick=6), symsize=2.5

; overplot galex, cfhtls, and irac
    galex = where(strmatch(filters,'*galex*',/fold))
    mab = maggies2mag(isedfit.maggies[galex],$
      ivar=isedfit.ivarmaggies[galex],magerr=mab_err)
    oploterror, filtinfo[galex].weff, mab, mab_err, psym=symcat(15), $
      symsize=3.0, color=im_color('powder blue',101), $
      errcolor=im_color('powder blue',101), errthick=!p.thick
    
    optical = where(strmatch(filters,'*capak*',/fold))
    mab = maggies2mag(isedfit.maggies[optical],$
      ivar=isedfit.ivarmaggies[optical],magerr=mab_err)
    oploterror, filtinfo[optical].weff, mab, mab_err, psym=symcat(15), $
      symsize=3.0, color=im_color('forest green',101), $
      errcolor=im_color('forest green',101), errthick=!p.thick
    
    irac = where(strmatch(filters,'*irac*',/fold))
    mab = maggies2mag(isedfit.maggies[irac],$
      ivar=isedfit.ivarmaggies[irac],magerr=mab_err)
    oploterror, filtinfo[irac].weff, mab, mab_err, psym=symcat(15), $
      symsize=3.0, color=im_color('tan',101), $
      errcolor=im_color('tan',101), errthick=!p.thick
    
; inset with P(M)     /noerase, 
    im_plothist, mstar, bin=0.04, /noplot, xx, yy
    im_plothist, mstar, bin=0.04, xsty=9, ysty=5, /noerase, yrange=[0,max(yy)*1.05], $
      position=[0.55,0.35,0.9,0.65], /fill, fcolor=im_color('grey60'), $
      ytitle='P(M)', xtitle='log (Stellar Mass)  (M'+sunsymbol()+')', color=keycolor, $
      ytickname=replicate(' ',10), charsize=1.5, xrange=isedfit.mass_50+4*isedfit.mass_err*[-1,1], $
      xtickinterval=0.2
;   oplot, isedfit.mass_50*[1,1], !y.crange, line=0, thick=6, color=djs_icolor('black')
;   oplot, isedfit.mass*[1,1], !y.crange, line=5, thick=6, color=djs_icolor('black')

    im_plotconfig, psfile=psfile, /psclose, keynote=keynote

stop    


    
    
; --------------------------------------------------
; spectrum of 
    ss = rd1dspec('ugca_116_drift_030_035.ms.fits',datapath=getenv('IM_PROJECTS_DIR')+'/atlas/atlas1d/')

    psfile = talkpath+'ugca116.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=4.5, charsize=2, keynote=keynote

    plot, [0], [0], /nodata, xrange=[3620,6880], yrange=[0,1.05], $
      xsty=9, ysty=5, xtitle=textoidl('Wavelength (\AA)'), $
      ytickname=replicate(' ',10), color=keycolor
    djs_oplot, ss.wave, ss.spec/max(ss.spec), psym=10, color=keycolor
    
    im_plotconfig, psfile=psfile, /psclose, keynote=keynote

stop    
    
stop
    
; --------------------------------------------------
; build a Montage of RC3 galaxies (must be run on offshore!)
;   hoggdir = '$IM_ARCHIVE_DIR/projects/hogg_rc3/labeled/'
    hoggdir = '$IM_ARCHIVE_DIR/projects/hogg_rc3/'
    pushd, hoggdir
    rc3 = file_search('*.jpg',count=nall)

    size = '150'
    nrow = '10'
    ngal = 100
;   these = shuffle_indx(nall,num_sub=ngal)

; 155, 381, 322 are bad
    these = [139,43,177,492,74,416,134,477,120,423,199,196,159,367,361,490,180,377,$
      443,405,229,378,460,474,233,480,322,88,499,0,105,438,205,189,440,505,235,156,$
      395,295,308,430,5,437,255,227,36,374,94,174,380,125,297,239,168,274,285,52,$
      394,132,369,343,75,415,269,493,65,366,266,410,126,265,325,136,311,301,242,89,$
      109,157,464,296,219,345,277,382,22,113,178,419,241,8,447,432,212,389,463,263,226,472]
    stop
    outfile = datapath+'rc3_montage.jpg'
    cmd = 'montage -bordercolor black -borderwidth 1 '+ $
      '-tile '+strtrim(nrow,2)+'x'+strtrim(nrow,2)+' -geometry +0+0 '+$
      '-quality 100 -resize '+size+'x'+size+' '+strjoin(rc3[these],' ')+$
      ' '+outfile
;   splog, cmd
    spawn, cmd, /sh

    popd

stop    
    
return
end
