pro talk_11sep_stanford, keynote=keynote, noevol=noevol
; jm11sep22ucsd - miscellaneous plots for my 2011 Sep talk at Stanford 

    mfpath = mf_path()
    mfspath = mf_path(/mfs)
    isedpath = mf_path(/isedfit)
    mzpath = mz_path()
    
    datapath = getenv('IM_RESEARCH_DIR')+'/talks/2011/11sep_stanford/'
    if keyword_set(keynote) then talkpath = datapath+'keynote/' else $
      talkpath = datapath

; --------------------------------------------------
; primus - 0.1(NUV-r) vs mass
    all = mrdfits(mfspath+'mfdata_all_supergrid01.fits.gz',1) ; fiducial supergrid
    zbins = mf_zbins(nzbins)
    
    psfile = talkpath+'mass_vs_nuvmr.ps'
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.6, $
      height=2.8*[1,1,1], width=3.5*[1,1], keynote=keynote

    if keyword_set(keynote) then keycolor = im_color('white',10) else $
      keycolor = im_color('black',10)

    xrange = [8.0,11.95]
    yrange = [0.01,7.5]
    xtitle1 = mfplot_masstitle()
    ytitle1 = textoidl('^{0.1}(NUV - r)')
;   maxis = range(xrange[0],xrange[1],50)
    maxis1 = range(xrange[0]+1.0,xrange[1]-0.08,50)
    maxis2 = range(xrange[0]+0.5,xrange[1]-0.08,50)
    maxis3 = range(xrange[0]+0.05,xrange[1]-0.08,50)

    for iz = 0, nzbins-1 do begin
       if odd(iz) then begin
          ytitle = '' & ytickname = replicate(' ',10)
       endif else begin
          ytitle = ytitle1 & delvarx, ytickname
       endelse
       if (iz lt 4) then xtickname = replicate(' ',10) else $
         delvarx, xtickname 

       indx = where((all.z ge zbins[iz].zlo) and (all.z lt zbins[iz].zup),nindx)
       mass = all[indx].mass
       sfrm = all[indx].sfrm100

       t2 = -alog10(2*getage(zbins[iz].zbin)) ; doubling time [Gyr^-1]
       t5 = -alog10(5*getage(zbins[iz].zbin)) ; tetratuple time [Gyr^-1]
       splog, zbins[iz].zbin, t2, t5
       low = where(sfrm lt t5,nlow)
       high = where(sfrm gt t2,nhigh)
       interm = where(sfrm gt t5 and sfrm lt t2,ninterm)

       mfplot_scatterplot, all[indx].mass, all[indx].nuvmr_01, $
         weight=all[indx].weight, noerase=(iz gt 0), /nogrey, $
         position=pos[*,iz], xrange=xrange, yrange=yrange, xsty=1, ysty=1, $
         ytitle=ytitle, xtickname=xtickname, ytickname=ytickname, $
         npix=16, ccolor=im_color('tan'), color=keycolor
;      mfplot_scatterplot, all[indx[high]].mass, all[indx[high]].nuvmr_01, $
;        weight=all[indx[high]].weight, noerase=(iz gt 0), /nogrey, $
;        position=pos[*,iz], xrange=xrange, yrange=yrange, xsty=1, ysty=1, $
;        ytitle=ytitle, xtickname=xtickname, ytickname=ytickname, $
;        npix=16, ccolor=im_color('sky blue'), outcolor=im_color('sky blue'), color=keycolor
;      mfplot_scatterplot, all[indx[low]].mass, all[indx[low]].nuvmr_01, $
;        weight=all[indx[low]].weight, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
;        xtitle='', ytitle='', /noerase, position=pos[*,iz], ccolor=im_color('tomato'), $
;        /nogrey, outcolor=im_color('tomato'), npix=16, color=keycolor
;      mfplot_scatterplot, all[indx[interm]].mass, all[indx[interm]].nuvmr_01, $
;        weight=all[indx[interm]].weight, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
;        xtitle='', ytitle='', /noerase, position=pos[*,iz], ccolor=im_color('spring green1'), $
;        /nogrey, outcolor=im_color('spring green1'), npix=16, color=keycolor

       junk = mf_select_nuvmr_quiescent(z=0.1,coeff=coeff,masspivot=masspivot)
       junk = mf_select_nuvmr_quiescent(z=zbins[iz].zbin,coeff=zcoeff)
       djs_oplot, maxis2, poly(maxis2-masspivot,coeff), line=0, color=keycolor, thick=6
       djs_oplot, maxis2, poly(maxis2-masspivot,zcoeff), line=5, color=keycolor, thick=6

       qq = mf_select_nuvmr_quiescent(all[indx].nuvmr_01,all[indx].mass,z=all[indx].z)
       
       im_legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[iz].zup,format='(F12.2)'),2), $
         /left, /top, box=0, margin=0, charsize=1.4, textcolor=keycolor

    endfor 
    xyouts, pos[0,5], pos[1,5]-0.07, align=0.5, /norm, xtitle1, color=keycolor

    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

stop    


; --------------------------------------------------
; number and mass-density evolution for all, quiescent, and
; star-forming galaxies in bins of mass
    all = mrdfits(mfpath+'rhoden_numden_vs_redshift_all.fits.gz',1)
    qq = mrdfits(mfpath+'rhoden_numden_vs_redshift_quiescent.fits.gz',1)
    sf = mrdfits(mfpath+'rhoden_numden_vs_redshift_active.fits.gz',1)

    tot = [[all],[qq],[sf]]
    psym = [15,14,16]
    if keyword_set(keynote) then $
      color = ['white','tomato','sky blue'] else $
        color = ['black','dark red','steel blue']
    symsize = [2.0,2.2,2.0]
    off = [0.0,+0.01,-0.01]

    if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'

    allerrors = 0
    if allerrors then $
      psfile = talkpath+'rhoden_numden_bymass_allerrors.ps' else $
        psfile = talkpath+'rhoden_numden_bymass.ps'
    im_plotconfig, 17, pos, psfile=psfile, xmargin=[1.2,0.3], ymargin=[0.7,0.85], $
      width=2.9*[1,1,1], height=[3.0,3.0], charsize=1.6, yspace=0.1, xspace=[0.1,0.1], keynote=keynote

    xrange = [0.001,1.01]
    rhorange = [6.6,8.3]
    numrange = [-4.3,-2.35]

; #########################
; 10-10.5
    plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='', xtickname=replicate(' ',10), $
      ytitle=mfplot_numdentitle(), color=im_color(keycolor)
    im_xyouts_title, title='log (M/M'+sunsymbol()+')=10-10.5', color=im_color(keycolor)
    im_legend, ['All','Quiescent','Star-Forming'], /left, /bottom, box=0, $
      charsize=1.2, psym=psym, symsize=symsize*0.9, color=color, textcolor=im_color(keycolor)

    for ii = 0, 2 do begin
       gd = where(tot[*,ii].num_10_105 gt -999.0,ngd)
       if (ngd ne 0) then begin
          djs_oploterr, tot[gd,ii].zbin+off[ii], tot[gd,ii].num_10_105, xerr=tot[gd,ii].zsigma*0, $
            yerr=tot[gd,ii].numerr_cv_10_105, $
            psym=symcat(psym[ii],color=im_color(keycolor)), symsize=symsize[ii], $
            color=im_color(color[ii]), errcolor=im_color(color[ii]), errstyle=0, errthick=6, /cap
          if allerrors then begin
; models             
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.9
                yhei = tot[gd[jj],ii].num_models_max_10_105-tot[gd[jj],ii].num_models_min_10_105
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].num_models_min_10_105
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=3
             endfor
;; priors
;             for jj = 0, ngd-1 do begin
;                xwid = tot[gd[jj],ii].zsigma*1.5
;                yhei = tot[gd[jj],ii].num_priors_max_10_105-tot[gd[jj],ii].num_priors_min_10_105
;                xcen = tot[gd[jj],ii].zbin+off[ii]
;                ycen = yhei/2.0+tot[gd[jj],ii].num_priors_min_10_105
;                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
;                  /data, color=im_color(color[ii]), thick=5, line=1
;             endfor
          endif
       endif
    endfor 

    plot, [0], [0], /nodata, /noeras, position=pos[*,3], xsty=1, ysty=1, $
      yrange=rhorange, xrange=xrange, xtitle='Redshift', $
      ytitle=mfplot_rhodentitle(), color=im_color(keycolor)
    for ii = 0, 2 do begin
       gd = where(tot[*,ii].rho_10_105 gt -999.0,ngd)
       if (ngd ne 0) then begin
          djs_oploterr, tot[gd,ii].zbin+off[ii], tot[gd,ii].rho_10_105, xerr=tot[gd,ii].zsigma*0, yerr=tot[gd,ii].rhoerr_cv_10_105, $
            psym=symcat(psym[ii],color=im_color(keycolor)), symsize=symsize[ii], color=im_color(color[ii]), $
            errcolor=im_color(color[ii]), errstyle=0, errthick=6, /cap
          if allerrors then begin
; models             
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.9
                yhei = tot[gd[jj],ii].rho_models_max_10_105-tot[gd[jj],ii].rho_models_min_10_105
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].rho_models_min_10_105
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=3
             endfor
;; priors
;             for jj = 0, ngd-1 do begin
;                xwid = tot[gd[jj],ii].zsigma*1.5
;                yhei = tot[gd[jj],ii].rho_priors_max_10_105-tot[gd[jj],ii].rho_priors_min_10_105
;                xcen = tot[gd[jj],ii].zbin+off[ii]
;                ycen = yhei/2.0+tot[gd[jj],ii].rho_priors_min_10_105
;                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
;                  /data, color=im_color(color[ii]), thick=5, line=1
;             endfor
          endif
       endif
    endfor 
    
; #########################
; 10.5-11
    plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='', xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), ytitle='', color=im_color(keycolor)
    im_xyouts_title, title='log (M/M'+sunsymbol()+')=10.5-11', color=im_color(keycolor)
;   djs_oplot, combo.z, combo.num_105_11, psym=6, color='green', symsize=2.0
;   djs_oplot, ilbert_qq.z, ilbert_qq.num_105_11, psym=6, color='orange', symsize=2.0
    for ii = 0, 2 do begin
       gd = where(tot[*,ii].num_105_11 gt -999.0,ngd)
       if (ngd ne 0) then begin
          djs_oploterr, tot[gd,ii].zbin+off[ii], tot[gd,ii].num_105_11, xerr=tot[gd,ii].zsigma*0, yerr=tot[gd,ii].numerr_cv_105_11, $
            psym=symcat(psym[ii],color=im_color(keycolor)), symsize=symsize[ii], color=im_color(color[ii]), $
            errcolor=im_color(color[ii]), errstyle=0, errthick=6, /cap
          if allerrors then begin
; models             
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.9
                yhei = tot[gd[jj],ii].num_models_max_105_11-tot[gd[jj],ii].num_models_min_105_11
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].num_models_min_105_11
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=3
             endfor
;; priors
;             for jj = 0, ngd-1 do begin
;                xwid = tot[gd[jj],ii].zsigma*1.5
;                yhei = tot[gd[jj],ii].num_priors_max_105_11-tot[gd[jj],ii].num_priors_min_105_11
;                xcen = tot[gd[jj],ii].zbin+off[ii]
;                ycen = yhei/2.0+tot[gd[jj],ii].num_priors_min_105_11
;                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
;                  /data, color=im_color(color[ii]), thick=5, line=1
;             endfor
          endif
       endif
    endfor 

    plot, [0], [0], /nodata, /noerase, position=pos[*,4], xsty=1, ysty=1, $
      yrange=rhorange, xrange=xrange, xtitle='Redshift', ytickname=replicate(' ',10), ytitle='', $
      color=im_color(keycolor)
;   djs_oplot, combo.z, combo.rho_105_11, psym=6, color='green', symsize=2.0
;   djs_oplot, ilbert_qq.z, ilbert_qq.rho_105_11, psym=6, color='orange', symsize=2.0
    for ii = 0, 2 do begin
       gd = where(tot[*,ii].rho_105_11 gt -999.0,ngd)
       if (ngd ne 0) then begin
          djs_oploterr, tot[gd,ii].zbin+off[ii], tot[gd,ii].rho_105_11, xerr=tot[gd,ii].zsigma*0, yerr=tot[gd,ii].rhoerr_cv_105_11, $
            psym=symcat(psym[ii],color=im_color(keycolor)), symsize=symsize[ii], color=im_color(color[ii]), $
            errcolor=im_color(color[ii]), errstyle=0, errthick=6, /cap
          if allerrors then begin
; models             
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.9
                yhei = tot[gd[jj],ii].rho_models_max_105_11-tot[gd[jj],ii].rho_models_min_105_11
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].rho_models_min_105_11
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=3
             endfor
;; priors
;             for jj = 0, ngd-1 do begin
;                xwid = tot[gd[jj],ii].zsigma*1.5
;                yhei = tot[gd[jj],ii].rho_priors_max_105_11-tot[gd[jj],ii].rho_priors_min_105_11
;                xcen = tot[gd[jj],ii].zbin+off[ii]
;                ycen = yhei/2.0+tot[gd[jj],ii].rho_priors_min_105_11
;                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
;                  /data, color=im_color(color[ii]), thick=5, line=1
;             endfor
          endif
       endif
    endfor 

; #########################
; 11-11.5
    plot, [0], [0], /nodata, /noerase, position=pos[*,2], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='', xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), ytitle='', color=im_color(keycolor)
    im_xyouts_title, title='log (M/M'+sunsymbol()+')=11-11.5', color=im_color(keycolor)
;   djs_oplot, combo.z, combo.num_11_115, psym=6, color='green', symsize=2.0
;   djs_oplot, ilbert_qq.z, ilbert_qq.num_11_115, psym=6, color='orange', symsize=2.0
    for ii = 0, 2 do begin
       gd = where(tot[*,ii].num_11_115 gt -999.0,ngd)
       if (ngd ne 0) then begin
          djs_oploterr, tot[gd,ii].zbin+off[ii], tot[gd,ii].num_11_115, xerr=tot[gd,ii].zsigma*0, yerr=tot[gd,ii].numerr_cv_11_115, $
            psym=symcat(psym[ii],color=im_color(keycolor)), symsize=symsize[ii], $
            color=im_color(color[ii]), errcolor=im_color(color[ii]), errstyle=0, errthick=6, /cap
          if allerrors then begin
; models             
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.9
                yhei = tot[gd[jj],ii].num_models_max_11_115-tot[gd[jj],ii].num_models_min_11_115
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].num_models_min_11_115
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=3
             endfor
;; priors
;             for jj = 0, ngd-1 do begin
;                xwid = tot[gd[jj],ii].zsigma*1.5
;                yhei = tot[gd[jj],ii].num_priors_max_11_115-tot[gd[jj],ii].num_priors_min_11_115
;                xcen = tot[gd[jj],ii].zbin+off[ii]
;                ycen = yhei/2.0+tot[gd[jj],ii].num_priors_min_11_115
;                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
;                  /data, color=im_color(color[ii]), thick=5, line=1
;             endfor
          endif
       endif
    endfor 

    plot, [0], [0], /nodata, /noerase, position=pos[*,5], xsty=1, ysty=1, $
      yrange=rhorange, xrange=xrange, xtitle='Redshift', ytickname=replicate(' ',10), ytitle='', $
      color=im_color(keycolor)
;   djs_oplot, combo.z, combo.rho_11_115, psym=6, color='green', symsize=2.0
;   djs_oplot, ilbert_qq.z, ilbert_qq.rho_11_115, psym=6, color='orange', symsize=2.0
    for ii = 0, 2 do begin
       gd = where(tot[*,ii].rho_11_115 gt -999.0,ngd)
       if (ngd ne 0) then begin
          djs_oploterr, tot[gd,ii].zbin+off[ii], tot[gd,ii].rho_11_115, xerr=tot[gd,ii].zsigma*0, yerr=tot[gd,ii].rhoerr_cv_11_115, $
            psym=symcat(psym[ii],color=im_color(keycolor)), symsize=symsize[ii], $
            color=im_color(color[ii]), errcolor=im_color(color[ii]), errstyle=0, errthick=6, /cap
          if allerrors then begin
; models             
             for jj = 0, ngd-1 do begin
                xwid = tot[gd[jj],ii].zsigma*1.9
                yhei = tot[gd[jj],ii].rho_models_max_11_115-tot[gd[jj],ii].rho_models_min_11_115
                xcen = tot[gd[jj],ii].zbin+off[ii]
                ycen = yhei/2.0+tot[gd[jj],ii].rho_models_min_11_115
                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
                  /data, color=im_color(color[ii]), thick=3
             endfor
;; priors
;             for jj = 0, ngd-1 do begin
;                xwid = tot[gd[jj],ii].zsigma*1.5
;                yhei = tot[gd[jj],ii].rho_priors_max_11_115-tot[gd[jj],ii].rho_priors_min_11_115
;                xcen = tot[gd[jj],ii].zbin+off[ii]
;                ycen = yhei/2.0+tot[gd[jj],ii].rho_priors_min_11_115
;                im_oplot_box, xwid, yhei, 0.0, xoffset=xcen, yoffset=ycen, $
;                  /data, color=im_color(color[ii]), thick=5, line=1
;             endfor
          endif
       endif
    endfor 

    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

stop    
    
; ---------------------------------------------------------------------------
; primus - 6x2 panel plot showing the field-averaged MFs for 'all', 
; 'quiescent', and 'star-forming galaxies'
    supergrid = 1
    xrange = [8.3,12.1]
;   xrange = [8.3,12.5]
    yrange = [-5.0,-1.1]
    showerr = 0
    zbins = mf_zbins(nzbins)
    maxis1 = range(8.8,12.5,100)

    nlimitshow = 3
    pcolor = 'black' & pline = 2
    scolor = 'black' & sline = 0
    
    psfile = talkpath+'mf_vmax_everything2.ps'
    im_plotconfig, 26, pos, psfile=psfile, charsize=1.4, $
      height=2.3*[1,1], keynote=keynote

    if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'

; 'all' galaxies
    mfdata = read_mf_vmax(noevol=noevol,/avgfield,supergrid=supergrid,/log)
    mffit = read_mf_vmax(/bestfit,noevol=noevol,/avgfield,supergrid=supergrid)
    mffit_sdss = read_mf_vmax('sdss',/bestfit,noevol=noevol,supergrid=supergrid)
    for iz = 0, nzbins-1 do begin
       if (iz gt 0) then begin
          ytitle = ''
          ytickname = replicate(' ',10)
       endif else begin
          ytitle = '' ; mfplot_phititle()
          delvarx, ytickname
       endelse
       plot, [0], [0], /nodata, noerase=(iz gt 0), position=pos[*,iz+nzbins*0], xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, xtickname=replicate(' ',10), ytickname=ytickname, xtitle='', $
         ytitle=ytitle, xtickinterval=1, xminor=4, color=im_color(keycolor)
;      if (iz eq 0) then legend, 'All', /left, /bottom, box=0, charsize=1.0, margin=0
       legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[iz].zup,format='(F12.2)'),2), $
         /right, /top, box=0, margin=0, charsize=1.1, textcolor=im_color(keycolor)

       good = where(mfdata[iz].limit eq 1)
       limit = where(mfdata[iz].limit eq 0 and mfdata[iz].number gt 0 and mfdata[iz].mass gt xrange[0]+0.05,nlimit)
       if nlimit gt 0 then limit = (reverse(limit))[0:nlimitshow<nlimit-1]
       if keyword_set(showerr) then begin
          oploterror, mfdata[iz].mass[good], mfdata[iz].phi[good], $
            mfdata[iz].phierr[good], psym=symcat(8,thick=3), symsize=0.75, $
            errthick=!p.thick
       endif else begin
          phimass = mfdata[iz].mass[good]
; Poisson+CV errors
          phimin = mfdata[iz].phi[good]-sqrt(mfdata[iz].phierr_lower[good]^2+mfdata[iz].phierr_cv[good]^2)
          phimax = mfdata[iz].phi[good]+sqrt(mfdata[iz].phierr_upper[good]^2+mfdata[iz].phierr_cv[good]^2)
          polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
            /data, color=im_color('dim grey'), noclip=0, /fill
; Poisson errors
          phimin = mfdata[iz].phi[good]-mfdata[iz].phierr_lower[good]
          phimax = mfdata[iz].phi[good]+mfdata[iz].phierr_upper[good]
          polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
            /data, color=im_color('light gray'), noclip=0, /fill

          if (nlimit gt 0) then begin
             oploterror, mfdata[iz].mass[limit], mfdata[iz].phi[limit], mfdata[iz].phierr_upper[limit], $
               psym=symcat(9,thick=3), symsize=1, /nohat, /hibar, $
               color=im_color('gray40'), errcolor=im_color('gray40')
             oploterror, mfdata[iz].mass[limit], mfdata[iz].phi[limit], mfdata[iz].phierr_lower[limit], $
               psym=3, symsize=1, /nohat, /lobar, color=im_color('gray40'), errcolor=im_color('gray40')
          endif
       endelse
       djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=mffit_sdss)), $
         line=sline, thick=4, color=im_color(scolor)
;      djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=mffit[iz])), $
;        line=pline, thick=4, color=im_color(pcolor)
    endfor 

; 'active' and 'quiscent' galaxies
    sf_mfdata = read_mf_vmax(noevol=noevol,/avgfield,supergrid=supergrid,/active,/log)
    sf_mffit = read_mf_vmax(/bestfit,noevol=noevol,/avgfield,supergrid=supergrid,/active)
    sf_mffit_sdss = read_mf_vmax('sdss',/bestfit,noevol=noevol,supergrid=supergrid,/active)
    qq_mfdata = read_mf_vmax(noevol=noevol,/avgfield,supergrid=supergrid,/quiescent,/log)
    qq_mffit = read_mf_vmax(/bestfit,noevol=noevol,/avgfield,supergrid=supergrid,/quiescent)
    qq_mffit_sdss = read_mf_vmax('sdss',/bestfit,noevol=noevol,supergrid=supergrid,/quiescent)
    for iz = 0, nzbins-1 do begin
       if (iz gt 0) then begin
          ytitle = ''
          ytickname = replicate(' ',10)
       endif else begin
;         ytitle = mfplot_phititle()
          delvarx, ytickname
       endelse
       plot, [0], [0], /nodata, /noerase, position=pos[*,iz+nzbins*1], xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, ytickname=ytickname, $
         xtitle='', ytitle=ytitle, xtickinterval=1, xminor=4, color=im_color(keycolor)
;      plot, [0], [0], /nodata, /noerase, position=pos[*,iz+nzbins*1], xsty=1, ysty=1, $
;        yrange=yrange, xrange=xrange, xtickname=replicate(' ',10), ytickname=ytickname, xtitle='', $
;        ytitle=ytitle, xtickinterval=1, xminor=4
;      if (iz eq 0) then legend, 'Star-Forming', /left, /bottom, box=0, charsize=1.0, margin=0

; 'active' galaxies
       good = where(sf_mfdata[iz].limit eq 1)
       limit = where(sf_mfdata[iz].limit eq 0 and sf_mfdata[iz].number gt 0 and sf_mfdata[iz].mass gt xrange[0]+0.05,nlimit)
       if nlimit gt 0 then limit = (reverse(limit))[0:nlimitshow<nlimit-1]
       if keyword_set(showerr) then begin
          oploterror, sf_mfdata[iz].mass[good], sf_mfdata[iz].phi[good], $
            sf_mfdata[iz].phierr[good], psym=symcat(8,thick=3), symsize=0.75, $
            errthick=!p.thick
       endif else begin
;         djs_oplot, sf_mfdata[iz].mass[good], sf_mfdata[iz].phi[good], $
;           psym=symcat(8,thick=3), symsize=0.75
          phimass = sf_mfdata[iz].mass[good]
; Poisson+CV errors
          phimin = sf_mfdata[iz].phi[good]-sqrt(sf_mfdata[iz].phierr_lower[good]^2+sf_mfdata[iz].phierr_cv[good]^2)
          phimax = sf_mfdata[iz].phi[good]+sqrt(sf_mfdata[iz].phierr_upper[good]^2+sf_mfdata[iz].phierr_cv[good]^2)
          polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
            /data, color=im_color('blue'), noclip=0, /fill
; Poisson errors
          phimin = sf_mfdata[iz].phi[good]-sf_mfdata[iz].phierr_lower[good]
          phimax = sf_mfdata[iz].phi[good]+sf_mfdata[iz].phierr_upper[good]
          polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
            /data, color=im_color('sky blue'), noclip=0, /fill

          if (nlimit gt 0) then begin
             oploterror, sf_mfdata[iz].mass[limit], sf_mfdata[iz].phi[limit], sf_mfdata[iz].phierr_upper[limit], $
               psym=symcat(9,thick=3), symsize=1, /nohat, /hibar, $
               color=im_color('sky blue'), errcolor=im_color('sky blue')
             oploterror, sf_mfdata[iz].mass[limit], sf_mfdata[iz].phi[limit], sf_mfdata[iz].phierr_lower[limit], $
               psym=3, symsize=1, /nohat, /lobar, color=im_color('sky blue'), errcolor=im_color('sky blue')
          endif
       endelse
       djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=sf_mffit_sdss)), $
         line=0, thick=4, color=im_color(scolor)
;      djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=sf_mffit[iz])), $
;        line=pline, thick=4, color=im_color(pcolor)
             
; 'quiescent' galaxies
       good = where(qq_mfdata[iz].limit eq 1)
       limit = where(qq_mfdata[iz].limit eq 0 and qq_mfdata[iz].number gt 0 and qq_mfdata[iz].mass gt xrange[0]+0.05,nlimit)
       if nlimit gt 0 then limit = (reverse(limit))[0:nlimitshow<nlimit-1]
       if keyword_set(showerr) then begin
          oploterror, qq_mfdata[iz].mass[good], qq_mfdata[iz].phi[good], $
            qq_mfdata[iz].phierr[good], psym=symcat(8,thick=3), symsize=0.75, $
            errthick=!p.thick
       endif else begin
;         djs_oplot, qq_mfdata[iz].mass[good], qq_mfdata[iz].phi[good], $
;           psym=symcat(8,thick=3), symsize=0.75
          phimass = qq_mfdata[iz].mass[good]
; Poisson+CV errors
          phimin = qq_mfdata[iz].phi[good]-sqrt(qq_mfdata[iz].phierr_lower[good]^2+qq_mfdata[iz].phierr_cv[good]^2)
          phimax = qq_mfdata[iz].phi[good]+sqrt(qq_mfdata[iz].phierr_upper[good]^2+qq_mfdata[iz].phierr_cv[good]^2)
          polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
            /data, color=im_color('dark red'), noclip=0, /fill
; Poisson errors
          phimin = qq_mfdata[iz].phi[good]-qq_mfdata[iz].phierr_lower[good]
          phimax = qq_mfdata[iz].phi[good]+qq_mfdata[iz].phierr_upper[good]
          polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
            /data, color=im_color('tomato'), noclip=0, /fill

          if (nlimit gt 0) then begin
             oploterror, qq_mfdata[iz].mass[limit], qq_mfdata[iz].phi[limit], qq_mfdata[iz].phierr_upper[limit], $
               psym=symcat(9,thick=3), symsize=1, /nohat, /hibar, $
               color=im_color('orange red'), errcolor=im_color('orange red')
             oploterror, qq_mfdata[iz].mass[limit], qq_mfdata[iz].phi[limit], qq_mfdata[iz].phierr_lower[limit], $
               psym=3, symsize=1, /nohat, /lobar, color=im_color('orange red'), errcolor=im_color('orange red')
          endif
       endelse
       djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=qq_mffit_sdss)), $
         line=5, thick=4, color=im_color(scolor)
;      djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=qq_mffit[iz])), $
;        line=pline, thick=4, color=im_color(pcolor)
    endfor

; x-title
    xyouts, pos[0,9], pos[1,9]-0.11, mfplot_masstitle(), align=0.5, /normal, color=im_color(keycolor)
    xyouts, pos[0,0]-0.05, pos[1,0], mfplot_phititle(), align=0.5, /normal, orientation=90, color=im_color(keycolor)
    
    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

stop    



; --------------------------------------------------
; SED-fitting example
    field = 'cfhtls_xmm'
    sfhgrid = 1
    supergrid = 1
    ndraw = isedfit_ndraw() ; number of random draws

    isedpath = mf_path(/isedfit)
    isedfit_sfhgrid_dir = mf_path(/montegrids)

; need all the chunks in memory!       
    chunkfile = file_search(isedpath+'sfhgrid01/charlot/'+field+$
      '_bc03_chab_chunk_????.fits.gz',count=nchunk)
    for jj = 0, nchunk-1 do begin
       modelgrid1 = mrdfits(chunkfile[jj],1)
       modelgrid1 = struct_trimtags(temporary(modelgrid1),except='modelmaggies')
       if (jj eq 0) then modelgrid = temporary(modelgrid1) else $
         modelgrid = [temporary(modelgrid),temporary(modelgrid1)]
    endfor
    nmodel = n_elements(modelgrid)
    nage = n_elements(modelgrid[0].age)
    nallmodel = nmodel*nage
    
    bigage    = reform(modelgrid.age,nallmodel)
    bigmass   = reform(modelgrid.mstar,nallmodel)
    bigtau    = reform(rebin(reform(modelgrid.tau,1,nmodel),nage,nmodel),nallmodel)
    bigZ      = reform(rebin(reform(modelgrid.Z,1,nmodel),nage,nmodel),nallmodel)
    bigav     = reform(rebin(reform(modelgrid.mu*modelgrid.av,1,nmodel),nage,nmodel),nallmodel)
    bignburst = reform(rebin(reform(modelgrid.nburst,1,nmodel),nage,nmodel),nallmodel)
    
;   splog, 'Reconstructing SFHs'
;   bigsfr = bigage*0D
;   bigsfr100 = bigage*0D       ; average over the previous 100 Myr
;   bigb100 = bigage*0D         ; birthrate parameter
;   bigmgal = bigage*0D         ; galaxy mass ignoring mass loss 
;   for imod = 0L, nmodel-1 do begin
;      tindx = lindgen(nage)+imod*nage
;      sfr = isedfit_reconstruct_sfh(modelgrid[imod],outage=bigage[tindx],$
;        sfr100=sfr100,b100=b100,mgalaxy=mgal)
;      bigsfr[tindx] = sfr          ; alog10(sfr)
;      bigsfr100[tindx] = sfr100    ; alog10(sfr100) 
;      bigb100[tindx] = b100
;      bigmgal[tindx] = mgal
;   endfor

; build the parameter file name    
    paramfile = isedpath+field+'_supergrid01_isedfit.par'
    params = read_isedfit_paramfile(paramfile,sfhgrid=sfhgrid)
    params.sfhgrid = string(sfhgrid,format='(I2.2)')
    params.imf = 'chab'
    params.synthmodels = 'bc03'
    params.redcurve = redcurve2string(1) ; charlot

    pp = read_mf_ubersample(field)
    ra = 34.776369D & dec = -5.3947370D ; picked this one by eye
    spherematch, pp.ra, pp.dec, ra, dec, 1D/3600.0, m1, m2
    index = m1[0]
    galaxy = hogg_iau_name(ra,dec,'CFHTLS-XMM')

; build the posterior distribution in mass    
    isedfit = read_mf_isedfit(field,supergrid=supergrid,rows=index,post=post)
    isedfit = isedfit[0] & post = post[0]

    logscale_err = post.scale_err/post.scale/alog(10)
    logscale = alog10(post.scale) + randomn(seed,ndraw)*logscale_err
    pofm = alog10(bigmass[post.draws])+logscale

    model = isedfit_restore(paramfile,junk,params=params,$
      iopath=isedpath,index=index,isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)
    filters = strtrim(get_mf_filters(field,nice_filte=nice_filters),2)
    filtinfo = im_filterspecs(filterlist=filters)

    xtitle1 = textoidl('Observed Wavelength (\AA)')
    xtitle2 = textoidl('Rest Wavelength (\AA)')
    ytitle1 = textoidl('m_{AB}')

    yrange = [28.5,18.5]
    xrange1 = [1000.0,6E4]
    xrange2 = xrange1/(1.0+isedfit.zobj)
    
    psfile = talkpath+'isedfit_example.ps'
    im_plotconfig, 8, pos, psfile=psfile, ymargin=[1.0,1.1], keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, /xlog, $
      xtickinterval=1000, position=pos, color=keycolor
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, /xlog, color=keycolor
    oplot, model.wave, model.flux, line=0, color=keycolor; color='grey'

; overplot the filter names
    nice_filters = repstr(repstr(nice_filters,'ch1','[3.6]'),'ch2','[4.5]')
    for ii = 0, n_elements(filters)-1 do xyouts, filtinfo[ii].weff, $
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
      'log (M/M'+sunsymbol()+')_{best} = '+strtrim(string(isedfit.mass,format='(F12.2)'),2),$
      'log (M/M'+sunsymbol()+')_{median} = '+strtrim(string(isedfit.mass_50,format='(F12.2)'),2)]
;   label = [$
;     galaxy,'z = '+string(isedfit.zobj,format='(F6.4)'),$
;     'log (M/M'+sunsymbol()+') = '+strtrim(string(isedfit.mass_avg,format='(F12.2)'),2)+$
;     '\pm'+strtrim(string(isedfit.mass_err,format='(F12.2)'),2)]

     legend, textoidl(label), /left, /top, box=0, charsize=1.5, $
       margin=0, textcolor=keycolor, spacing=2.2

; inset with P(M)     /noerase, 
    im_plothist, pofm, bin=0.02, /noplot, xx, yy
    im_plothist, pofm, bin=0.02, xsty=3, ysty=1, /noerase, yrange=[0,max(yy)*1.05], $
      position=[0.55,0.3,0.9,0.6], /fill, fcolor=im_color('goldenrod'), $
      ytitle='P(M)', xtitle='Stellar Mass', color=keycolor, $
      ytickname=replicate(' ',10), charsize=1.5, xrange=isedfit.mass_50+3.5*isedfit.mass_err*[-1,1]
    oplot, isedfit.mass_50*[1,1], !y.crange, line=0, thick=6, color=djs_icolor('black')
    oplot, isedfit.mass*[1,1], !y.crange, line=5, thick=6, color=djs_icolor('black')

    im_plotconfig, psfile=psfile, /psclose, keynote=keynote

stop    


; ---------------------------------------------------------------------------
; sdss - compare MFs derived using different prior assumptions 
    super = get_mf_supergrid(nsuper,superstring=superstring)

    psfile = talkpath+'mf_vmax_sdss_differentpriors.ps'
    im_plotconfig, 3, pos, psfile=psfile, charsize=1.4, keynote=keynote
    
    if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'

    xrange = [8.8,12.2]
    yrange = [-7.2,-1.3]
    maxis1 = range(xrange[0]+0.1,12.5,100)
    showerr = 0
    
    plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytitle=mfplot_phititle(), xtickinterval=1, color=im_color(keycolor)
    im_legend, 'All', /right, /top, box=0, charsize=1.4, margin=0, textcolor=keycolor
    stanford_render_priorplot, 'all', noevol=noevol, super=super, $
      xrange=xrange, yrange=yrange, /dolegend, textcolor=keycolor
    
    plot, [0], [0], /nodata, position=pos[*,1], /noerase, xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytickname=replicate(' ',10), xtickinterval=1, color=im_color(keycolor)
    im_legend, 'Quiescent', /right, /top, box=0, charsize=1.4, margin=0, textcolor=keycolor
    stanford_render_priorplot, 'quiescent', noevol=noevol, super=super, $
      xrange=xrange, yrange=yrange, textcolor=keycolor
    
    plot, [0], [0], /nodata, position=pos[*,2], /noerase, xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytickname=replicate(' ',10), xtickinterval=1, color=im_color(keycolor)
    im_legend, 'Star-Forming', /right, /top, box=0, charsize=1.4, margin=0, textcolor=keycolor
    stanford_render_priorplot, 'active', noevol=noevol, super=super, $
      xrange=xrange, yrange=yrange, textcolor=keycolor
    
    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

stop    
    
; ---------------------------------------------------------------------------
; sdss - compare all, quiescent, and active samples
    supergrid = 1
    psfile = talkpath+'mf_vmax_sdss.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.8, height=6.5, charsize=2, keynote=keynote

    if keyword_set(keynote) then keycolor = im_color('white',10) else $
      keycolor = im_color('black',10)

    xrange = [8,12.0]
    yrange = [-6.1,-0.8]
    maxis1 = range(xrange[0]+0.15,12.5,100)

    subsample = ['active','quiescent','all']
    if keyword_set(keynote) then begin
       fitcolor = ['medium blue','red','black']
       mfcolor = ['powder blue','tomato','snow']
    endif else begin
       fitcolor = ['medium blue','red','black']
       mfcolor = ['dodger blue','tomato','dark grey']
    endelse

    fitcolor_double = 'black'
    mfline_double = 1
    
    psymsize = [1.2,1.3,1.3]*1.0
    psymsize_limit = [0.8,0.9,0.8]*1.2
;   psymsize = [2.0,2.5,2.0]*0.7
;   mfpsym = [9,4,6]
    mfpsym = [9,4,15]
;   mfpsym = [16,15,14]
    mfpsym_limit = [9,6,4]
    mfline = [5,3,0]
    nlimitshow = 2

    plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytitle=mfplot_phititle(), color=keycolor
    im_legend, ['All','Quiescent','Star-Forming'], /left, /bottom, box=0, $
      color=reverse(mfcolor), psym=reverse(mfpsym), symsize=[1.7,2.1,1.9], $
      symthick=8, spacing=2.4, margin=0, charsize=1.7, textcolor=keycolor
    
; print the "quenching mass"
    qq = alog10(mf_schechter(maxis1,schechter=read_mf_vmax('sdss',/quie,/best)))
    aa = alog10(mf_schechter(maxis1,schechter=read_mf_vmax('sdss',/act,/best)) )
    splog, 'quenching mass ', interpol(maxis1,qq-aa,0.0)

    mffitall = maxis1*0.0 ; add the active and quiescent Schechter fits
    for jj = 0, n_elements(subsample)-1 do begin
       case jj of
          0: begin
             mfdata = read_mf_vmax('sdss',supergrid=supergrid,/active,/log)
             mffit = read_mf_vmax('sdss',supergrid=supergrid,/active,/bestfit)
             model = mf_schechter(maxis1,schechter=mffit)
             mffitall = model
          end
          1: begin
             mfdata = read_mf_vmax('sdss',supergrid=supergrid,/quiescent,/log)
             mffit = read_mf_vmax('sdss',supergrid=supergrid,/quiescent,/bestfit)
             model = mf_schechter(maxis1,schechter=mffit)
             mffitall += model
          end
          2: begin
             mfdata = read_mf_vmax('sdss',supergrid=supergrid,/log)
             mffit = read_mf_vmax('sdss',supergrid=supergrid,/bestfit)
;            mffit_double = read_mf_vmax('sdss',supergrid=supergrid,/bestfit,/double)
             model = mf_schechter(maxis1,schechter=mffit)
;            model_double = mf_schechter_plus(maxis1,schechter=mffit_double)
          end
       endcase
       
       good = where(mfdata.number gt 0 and mfdata.mass gt xrange[0]+0.08,nthese)
;      good = where(mfdata.limit eq 1,nthese)
       limit = where(mfdata.number gt 0 and mfdata.limit eq 0 and mfdata.mass gt xrange[0]+0.08,nlimit)

       oploterror, mfdata.mass[good], mfdata.phi[good], mfdata.phierr_upper[good], $
         psym=symcat(mfpsym[jj],thick=5), symsize=psymsize[jj], errthick=5, $
         color=im_color(mfcolor[jj]), errcolor=im_color(mfcolor[jj]), /hibar, /nohat
       oploterror, mfdata.mass[good], mfdata.phi[good], mfdata.phierr_lower[good], psym=3, $
         symsize=psymsize[jj], errthick=5, color=im_color(mfcolor[jj]), $
         errcolor=im_color(mfcolor[jj]), /lobar, /nohat

       if nlimit gt 0 then limit = (reverse(limit))[0:nlimitshow<nlimit-1]
;      if (nlimit gt 0) then begin
;         oploterror, mfdata.mass[limit], mfdata.phi[limit], mfdata.phierr_upper[limit], $
;           psym=symcat(mfpsym_limit[jj],thick=2), symsize=psymsize_limit[jj], errthick=5, $
;           color=im_color(mfcolor[jj]), errcolor=im_color(mfcolor[jj]), /hibar, /nohat
;         oploterror, mfdata.mass[limit], mfdata.phi[limit], mfdata.phierr_upper[limit], psym=3, $
;           symsize=psymsize[jj], errthick=5, color=im_color(mfcolor[jj]), $
;           errcolor=im_color(mfcolor[jj]), /lobar, /nohat
;      endif
       
;       if (nlimit gt 0) then begin
;          limit = reverse(these[0]-(lindgen(nlimit)+1))
;          plotsym, 2, thick=6
;          djs_oplot, mfdata.mass[limit], alog10(mfdata.phi[limit]), psym=8, $
;            color=im_color(mfcolor[jj],100+jj), symsize=psymsize[jj]
;       endif
       
        if subsample[jj] eq 'all' then begin
 ;         djs_oplot, maxis1, alog10(mffitall), line=mfline[jj], $
 ;           thick=7, color=im_color('dark green')
           djs_oplot, maxis1, alog10(model), line=mfline[jj], thick=7, color=im_color(fitcolor[jj])
 ;         djs_oplot, maxis1, alog10(model_double), line=mfline_double, thick=7, color=im_color(fitcolor_double)
        endif else begin
;          djs_oplot, maxis1, alog10(model), line=mfline[jj], thick=7, color=im_color(fitcolor[jj])
        endelse
    endfor
    mfoplot_lit, maxis1, color=im_color('black'), line=5, thick=8, /baldry, /log
;   mfoplot_lit, maxis1, color=keycolor, line=0, thick=8, /bell, /log

    im_plotconfig, /psclose, psfile=psfile, keynote=keynote

; ---------------------------------------------------------------------------
; make a legend for the NUV-r vs mass plot
    psfile = talkpath+'mass_vs_nuvmr_legend.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote

    djs_plot, [0], [0], /nodata, position=pos, xsty=5, ysty=5, $
      xrange=[0,1], yrange=[0,1]
    label = ['-1.2<log (\psi/M)','-1.2<log (\psi/M)<-1.6',$
      '-1.6>log (\psi/M)']
    im_legend, label, /left, /top, box=0, charsize=1.5, $
      textcolor=['sky blue','spring green1','tomato']
    
    im_plotconfig, /psclose, psfile=psfile, keynote=keynote



; ---------------------------------------------------------------------------    
; effect of star formation and mergers on the MF; for the "mergers"
; mass function see Baldry+04 (the appendix)
    minmass = 7.0
    maxmass = 12.5
    mass = range(minmass,maxmass,500)

; pure star formation for 5 Gyr
    time = 7D9
    logsfr = 0.67*mass - 6.19      ; from Noeske+07
    newmass_sf = alog10(10^mass+0.5*10^logsfr*time) ; account for recycling
    
    all = {phistar: 4.26D-3, logmstar: 10.648+0.25D, alpha: -0.46D, $
      phiplus: 0.58D-3, alphaplus: -1.58D}
    merge = {phistar: 1.6D-3, logmstar: 11.1D, alpha: -0.4D}

    xrange = [8.0,15.0]
    gg = where((mass gt xrange[0]+0.05) and (mass lt xrange[1]-0.05))
    gmass = mass[gg]
    g2 = where((newmass_sf gt xrange[0]+0.1) and (newmass_sf lt xrange[1]-0.1))
    gnewmass_sf = newmass_sf[g2]

    psfile = talkpath+'toymodel.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, charsize=2.0, $
      height=5.8

    if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'
    
    plot, gmass, alog10(mf_schechter_plus(gmass,schechter=all)), xsty=1, $
      ysty=1, xrange=xrange, yrange=[-7.1,1.5], xtitle='Mass (M'+sunsymbol()+')', $
      ytitle=textoidl('Number Density (Mpc^{-3} dex^{-1})'), $
      color=im_color(keycolor), thick=4, xtickinterval=1, line=0
    djs_oplot, gmass, alog10(mf_schechter(gmass,schechter=merge)), $
      color=im_color('orange',100), thick=12, line=3
    djs_oplot, gnewmass_sf, alog10(mf_schechter_plus(mass[g2],schechter=all)), $
      color=im_color('cyan',101), thick=12, line=5

; overplot the halo mass function; see email from Kyle Stewart and his
; eps2.pro routine
    readfast, datapath+'sigma_lcdm.dat', data, skip=1
    mass = reform(10^data[0,*])
    sigma = reform(data[1,*])
    nn = n_elements(sigma)-1 ; no endpoint
    mass2 = mass[1:nn-1] ; no endpoint
    
; calculate the analytic formula from Warren, et al 2006.
; dn/dlogM = -rho/M * dlogsigma/dlogM * ( A(sigma^-a + b)*exp(-c/sigma^2) )
;  rho = .3*rho_crit = .3*1.88h^2*10^.29 g/cm^3
; best fits for the fitting variables A,a,b,c are:
; A=.7234  a=1.625  b=.2538  c=1.1982
    const = -83.294D9 ; = (-.3*(1.88e-29))*((3.086e24)^3.0)*(1/(1.99e33)) to cancel Msun units and cm^3->Mpc^3

    dsigdm = ((shift(alog(sigma),1)-shift(alog(sigma),-1))/(shift(alog(mass),1)-shift(alog(mass),-1)))[1:nn-1]
    dndm = (const/mass2)*dsigdm*0.7234*(sigma[1:nn-1]^(-1.625)+.2538)*exp(-1.1982/(sigma[1:nn-1]^2))
    djs_oplot, alog10(mass2), alog10(dndm), line=0, thick=12, color=im_color(keycolor)

    im_legend, ['Star Formation','Mergers','Halo Mass Function'], $
      /right, /top, box=0, textcolor=['cyan','orange',keycolor], charsize=1.5, $
      line=[5,3,0], thick=10, pspacing=2.0, color=['cyan','orange',keycolor]
;   im_legend, ['Initial Mass Function','Mass-Dependent Star Formation','Simple Merger Model'], $
;     /left, /bottom, box=0, textcolor=['white','cyan','orange'], charsize=1.5, $
;     line=[0,5,3], thick=10, pspacing=2.0, color=['white','cyan','orange']

    
    
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; cosmic variance errors in bins of mass and redshift
    recalculate = 0
    if recalculate then begin
       zmin = 0.0 & zmax = 1.0 & zbin = 0.1
       nzedge = (zmax-zmin)/zbin+1
       zedge = range(zmin,zmax,nzedge)
       zz = zedge[0:nzedge-2]+zbin/2.0
       nzz = n_elements(zz)
       
       b090=0.074  & b190=2.58 & b290=1.039
       b0100=0.053 & b1100=3.07 & b2100=1.225
       b0110=0.173 & b1110=2.89 & b2110=1.438

       nsurvey = 3
       cv = {survey: '', s1: 0.0, s2: 0.0, zz: zz, $
         cv_dm: zz*0, cv_90: zz*0, cv_10: zz*0, cv_11: zz*0}
       cv = replicate(cv,nsurvey)

       cv[0].survey = 'cosmos' & cv[0].s1 = 84.0/60 & cv[0].s2 = 84.0/60
       cv[1].survey = 'goods'  & cv[1].s1 = 10.0/60 & cv[1].s2 = 16.0/60
       cv[2].survey = 'egs'    & cv[2].s1 = 10.0/60 & cv[2].s2 = 70.0/60

       for ss = 0, nsurvey-1 do begin
          for ii = 0, nzz-1 do begin
             cv[ss].cv_dm[ii] = quickcv(cv[ss].s1,cv[ss].s2,zz[ii],zbin,sig8=0.77) ; dark matter CV
             cv[ss].cv_90[ii] = cv[ss].cv_dm[ii]*(b090*(1+zz[ii])^b190+b290)
             cv[ss].cv_10[ii] = cv[ss].cv_dm[ii]*(b0100*(1+zz[ii])^b1100+b2100)
             cv[ss].cv_11[ii] = cv[ss].cv_dm[ii]*(b0110*(1+zz[ii])^b1110+b2110)
          endfor
       endfor
       im_mwrfits, cv, datapath+'cv_errors.fits', /nogzip
    endif else cv = mrdfits(datapath+'cv_errors.fits',1)
    nsurvey = n_elements(cv)
    
; test that we match moster+11
;   plot, cv[0].zz, cv[0].cv_dm, psym=6
;   bb = range(0,1,20)
;   djs_oplot, bb, 0.069/(0.234+bb^0.834)

    psfile = talkpath+'cv_errors.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, $
      xmargin=[1.5,0.4], width=6.6, height=5.2, charsize=2.4

    if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'

    plot, [0], [0], /nodata, xsty=1, ysty=1, position=pos, $
      xrange=[0.0,1.0], yrange=[0.0,0.8], color=im_color(keycolor), $
      xtitle='Redshift', ytitle=textoidl('\sigma_{V}')
    djs_xyouts, 0.2, 0.1, '\Delta'+'z=0.1', align=0.5, $
      charsize=2.2, color=im_color(keycolor);, charthick=4
    
    color1 = 'midnight blue'
    djs_oplot, cv[0].zz, cv[0].cv_10, line=0, color=im_color(color1), thick=7
    djs_oplot, cv[0].zz, cv[0].cv_11, line=0, color=im_color(color1), thick=12

    color2 = 'powder blue'
    djs_oplot, cv[1].zz, cv[1].cv_10, line=5, color=im_color(color2), thick=7
    djs_oplot, cv[1].zz, cv[1].cv_11, line=5, color=im_color(color2), thick=12

    djs_xyouts, 0.4, 0.54, '10^{11.25} M'+sunsymbol(), align=0, $
      charsize=1.6, color=im_color(keycolor), charthick=4
    djs_xyouts, 0.4, 0.40, '10^{10.25} M'+sunsymbol(), align=0, $
      charsize=1.6, color=im_color(keycolor), charthick=4

    djs_xyouts, 0.95, 0.24, '10^{11.25} M'+sunsymbol(), align=1, $
      charsize=1.6, color=im_color(keycolor), charthick=4
    djs_xyouts, 0.95, 0.08, '10^{10.25} M'+sunsymbol(), align=1, $
      charsize=1.6, color=im_color(keycolor), charthick=4

    area = string(cv[0:1].s2*cv[0:1].s1,format='(F4.2)')+' deg^{2}'
    im_legend, ['COSMOS '+area[0],'GOODS '+area[1]], /right, /top, box=0, $
      color=[color1,color2], line=[0,5], charsize=1.7, margin=0, thick=10, $
      textcolor=keycolor
    
;   scolor = ['midnight blue','powder blue','forest green']
;   mline = [3,0,5]
;   for ss = 0, nsurvey-1 do begin
;      djs_oplot, cv[ss].zz, cv[ss].cv_90, line=mline[0], color=im_color(scolor[ss]), thick=6
;      djs_oplot, cv[ss].zz, cv[ss].cv_10, line=mline[1], color=im_color(scolor[ss]), thick=5
;      djs_oplot, cv[ss].zz, cv[ss].cv_11, line=mline[2], color=im_color(scolor[ss]), thick=10
;   endfor
    
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose



; --------------------------------------------------
; effect of redshift errors on absolute magnitudes
    psfile = talkpath+'photoz_errors.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, $
      xmargin=[1.5,0.4], width=6.6, height=5.2, charsize=2.4

    if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'

    nzz = 100
    nzerr = 1000
    zz = range(0.05,1.5,nzz)
    zerr = range(0.0001,0.12,nzerr,/log) ; =dz/(1+z)
;   dzbin = median(zerr-shift(zerr,1))
    dzbin = median(alog10(zerr)-shift(alog10(zerr),1))

    dmerr = fltarr(nzerr)
    for ii = 0, nzerr-1 do begin
       zz1 = (zz + randomn(seed,nzz)*(zerr[ii]*(1+zz)))>0.01
       dmerr[ii] = stddev(lf_distmod(zz)-lf_distmod(zz1))
;      if finite(dmerr[ii]) eq 0 then stop
;      plot, zz, 100*(dmodulus(zz)-dmodulus(zz1)), psym=6, yr=[-30,30] & cc = get_kbrd(1)
    endfor
    zerr = 100*zerr ; %
    dmerr = 100*dmerr ; %
    
    coeff = linfit(zerr,dmerr,yfit=yfit,sigma=sig)
    splog, coeff
    emax = yfit+poly(zerr,[0.0,2.5])
    emin = yfit-poly(zerr,[0.0,2.5])
    
;   dfpsclose & plot, zerr, dmerr, psym=6
;   djs_oplot, zerr, yfit, line=0, color='red'
;   djs_oplot, zerr, emax, line=0, color='blue'
;   djs_oplot, zerr, emin, line=0, color='blue'
    
;   emax = yfit+1.5*sig
;   emin = (yfit-1.5*sig)>0

    xrange = minmax(zerr)
    yrange = [0,120]
    plot, [0], [0], /nodata, xsty=5, ysty=5, position=pos, $
      xrange=xrange, yrange=yrange, /xlog, color=im_color(keycolor)
    legend, '0 < z < 1.5', /left, /top, box=0, textcolor=im_color(keycolor), charsize=2

; primus
    polyfill, [0.2,0.8,0.8,0.2], [0,0,15,15], $
      /data, color=im_color('tomato'), noclip=0, /fill
    xyouts, 0.4, 18, 'PRIMUS', align=0.5, /data, color=im_color(keycolor), charsize=2

; specz
    polyfill, [min(zerr),0.07,0.07,min(zerr)], [0,0,5,5], $
      /data, color=im_color('dodger blue'), noclip=0, /fill
    xyouts, 0.027, 8, 'specz', align=0.5, /data, color=im_color(keycolor), charsize=2

; photoz
    polyfill, [2,12,12,2], [0,0,100,100], $
      /data, color=im_color('forest green'), noclip=0, /fill
    xyouts, 4.0, 105, 'photoz', align=0.5, /data, color=im_color(keycolor), charsize=2

; plot the errors    
    if keyword_set(keynote) then col = 'white' else col = 'grey'
    polyfill, [zerr,reverse(zerr)],[emin,reverse(emax)], $
      /data, color=im_color(col), noclip=0, /fill

; overplot the axes
    plot, [0], [0], /nodata, /noerase, xsty=1, ysty=1, position=pos, $
      xrange=xrange, yrange=yrange, xtitle='Relative Redshift Error (%)', $
      ytitle='Luminosity Error (%)', /xlog, color=im_color(keycolor)
    
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose
    
stop    




; --------------------------------------------------
; survey comparison - magnitude limit vs area
    psfile = talkpath+'maglimit_vs_area_comp.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, $
      xmargin=[1.5,0.2], width=6.8, height=5.5

    if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'

    readcol,datapath+'survey_info.dat',survey,nobj,zmin,$
      zmax,zmed,area,mag,depth,idepth,v,logv,$
      format='a,l,f,f,f,f,a,f,f,f', /silent, comment=';'
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

    if keyword_set(keynote) then col = im_color('goldenrod',101) else col = im_color('charcoal',10)
    for i=0,n_elements(spec)-1 do begin
       oplot,[logv[spec[i]],logv[spec[i]]],[nobj[spec[i]],nobj[spec[i]]],$
         psym=symcat(14),symsize=symm[spec[i]],color=col
    endfor
    for i=0,n_elements(phot)-1 do begin
       oplot,[logv[phot[i]],logv[phot[i]]],[nobj[phot[i]],nobj[phot[i]]],$
         psym=symcat(17),symsize=symm[phot[i]]-0.4,color=im_color('royal blue',10)
    endfor

; PRIMUS:
    oplot,[logv[11],logv[11]],[nobj[11],nobj[11]],psym=symcat(15),$
      symsize=symm[11]-0.8,color=im_color('firebrick',10)

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

    
    
    
    
    
    
return
end
