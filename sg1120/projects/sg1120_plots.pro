pro sg1120_plots, gcat1, rcat1, hstcat1, rim, hstim, hststamps=hststamps, ps=ps
; jm07apr26nyu - generate some preliminary plots

    hststamps = 1L ; NOTE!
    
    mosaicpath = sg1120_path(/mosaics)
    catpath = sg1120_path(/catalogs)
    analysis_path = sg1120_path(/analysis)

; read the g- and r-prime photometric catalogs, and the redshift
; catalog, and cross-match

    if (n_elements(rim) eq 0L) then rim = mrdfits(mosaicpath+'sg1120_rprime.fits',0,/silent)
    if (n_elements(gcat1) eq 0L) then gcat1 = rsex(catpath+'sg1120_gprime_matched.cat')
    if (n_elements(rcat1) eq 0L) then rcat1 = rsex(catpath+'sg1120_rprime_matched.cat')
    if (n_elements(hstcat1) eq 0L) then hstcat1 = rsex(catpath+'sg1120_hst_matched_merged.cat')
    zcat1 = rsex(catpath+'sg1120_07jun25.zcat.sex')

    spherematch, rcat1.alpha_j2000, rcat1.delta_j2000, zcat1.ra, zcat1.dec, $
      1.0/3600.0, catmatch, zcatmatch, distance12, maxmatch=1
    gcat = gcat1[catmatch]
    rcat = rcat1[catmatch]
    hstcat = hstcat1[catmatch]
    zcat = zcat1[zcatmatch]
    
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.7, height=6.7, $
      xmargin=[1.5,0.3], ymargin=[0.3,1.5], xpage=8.5, ypage=8.5, $
      position=pos, /normal
    
    if keyword_set(ps) then begin
       dfpsplot, analysis_path+'sg1120_plots.ps', xsize=8.5, ysize=8.5, /color
       postthick = 5.0
    endif else begin
       im_window, 0, xratio=0.45, /square
       postthick = 2.0
    endelse

; ---------------------------------------------------------------------------
; color-magnitude diagram of spectroscopically confirmed cluster members
; ---------------------------------------------------------------------------

    plotsym, 0, 1.2, /fill
    
    indx = where((rcat.mag_auto lt 90.0) and (gcat.mag_auto lt 90.0) and $
      (rcat.mag_aper lt 90.0) and (gcat.mag_aper lt 90.0) and $
      (zcat.z gt 0.32) and (zcat.z lt 0.42),nindx)
    x = rcat[indx].mag_auto
    y = gcat[indx].mag_aper-rcat[indx].mag_aper

;   xrange = minmax(x)
    xrange = [-10.5,-4.5]
    yrange = [0.0,2.5] ; minmax(y)
    
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, charsize=2.5, charthick=postthick, $
      xthick=postthick, ythick=postthick, xtitle='r_{auto}', ytitle='(g-r)_{1"}', $
      xrange=xrange, yrange=yrange, position=pos
    djs_oplot, x, y, ps=8
    legend, '0.32 < z < 0.42', /left, /top, box=0, charsize=2.0, charthick=postthick
    
; fit the red sequence

    grcut = 1.5
;   oplot, !x.crange, grcut*[1,1]
    
    red = where((y gt grcut),nred)
    pivot = -8.5
    xdata = x[red]-pivot
    ydata = y[red]
    yivar = ydata*0.0+1.0

;   coeff = robust_linefit(xdata,ydata,yfit,sig,coeff_err)
    aa = fltarr(2,nred)
    aa[0,*] = 1.0D & aa[1,*] = xdata
    hogg_iter_linfit, aa, ydata, yivar, coeff, nsigma=2.0, covar=covar
    coeff_err = coeff*0.0
    wh = where(covar[lindgen(2),lindgen(2)] ge 0.0,ct) ; from MPFIT, line 3060
    if (ct gt 0) then coeff_err[wh] = sqrt(covar[wh,wh])

    coeff[0] = coeff[0] - coeff[1]*pivot
    sig = djsig(ydata-poly(xdata,coeff))

    splog, coeff[0], coeff[1], sig
    splog, coeff_err[0], coeff_err[1]

    xfit = findgen((!x.crange[1]-!x.crange[0])/0.005)*0.005+!x.crange[0]
    yfit = poly(xfit,coeff)
    djs_oplot, xfit, yfit, line=0, thick=postthick
    djs_oplot, xfit, yfit-sig, line=2, thick=postthick
    djs_oplot, xfit, yfit+sig, line=2, thick=postthick

    legend, textoidl('\Delta(g-r)_{rms} = '+strtrim(string(sig,format='(F12.3)'),2)+' mag'), $
      /left, /bottom, box=0, charsize=2.0, charthick=postthick
;   legend, 'Slope = '+strtrim(string(coeff[1],format='(F12.3)'),2)+', '+$
;     'Scatter = '+strtrim(string(sig,format='(F12.3)'),2)+' mag', $
;     /left, /bottom, box=0, charsize=2.0, charthick=postthick
    
;   help, where(poly(coeff,x)
    
;   rej = where(yivar eq 0.0,nrej)
;   if (nrej ne 0L) then for irej = 0L, nrej-1L do plots, xdata[rej[irej]]+pivot, $
;     ydata[rej[irej]], ps=7, sym=3

    if (not keyword_set(ps)) then cc = get_kbrd(1)

; ---------------------------------------------------------------------------
; color-magnitude diagram of all galaxies
; ---------------------------------------------------------------------------

    plotsym, 0, 0.6, /fill
    
    indx = where((rcat1.mag_auto lt 90.0) and (gcat1.mag_auto lt 90.0) and $
      (rcat1.mag_aper lt 90.0) and (gcat1.mag_aper lt 90.0) and $
      (rcat1.class_star lt 0.8),nindx)
    x = rcat1[indx].mag_auto
    y = gcat1[indx].mag_aper-rcat1[indx].mag_aper

;   xrange = minmax(x)
;   yrange = [0.0,2.3] ; minmax(y)
    
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, charsize=2.5, charthick=postthick, $
      xthick=postthick, ythick=postthick, xtitle='r_{auto}', ytitle='(g-r)_{1"}', $
      xrange=xrange, yrange=yrange, position=pos
    djs_oplot, x, y, ps=8
    legend, 'All Galaxies', /left, /top, box=0, charsize=2.0, charthick=postthick

    xfit = findgen((!x.crange[1]-!x.crange[0])/0.005)*0.005+!x.crange[0]
    yfit = poly(xfit,coeff)
    djs_oplot, xfit, yfit, line=0, thick=postthick
    djs_oplot, xfit, yfit-sig, line=2, thick=postthick
    djs_oplot, xfit, yfit+sig, line=2, thick=postthick

    if keyword_set(ps) then dfpsclose else cc = get_kbrd(1)

; ---------------------------------------------------------------------------
; generate postage stamps of the red spectroscopic members
; ---------------------------------------------------------------------------

    r = rcat.mag_auto
    gr = gcat.mag_aper-rcat.mag_aper

    indx = where((rcat.mag_auto lt 90.0) and (gcat.mag_auto lt 90.0) and $
      (rcat.mag_aper lt 90.0) and (gcat.mag_aper lt 90.0) and $
      (zcat.z gt 0.32) and (zcat.z lt 0.42) and $
      (gr lt (poly(r,coeff)+sig)) and (gr gt (poly(r,coeff)-sig)),nindx)
;   indx = where((rcat.mag_auto lt 90.0) and (gcat.mag_auto lt 90.0) and $
;     (rcat.mag_aper lt 90.0) and (gcat.mag_aper lt 90.0) and $
;     (zcat.z gt 0.32) and (zcat.z lt 0.42),nindx)

    if keyword_set(hststamps) then begin
       pixscale = 0.035 
       stamp = 5.0/pixscale    ; [pixel]
    endif else begin
       pixscale = 0.188         ; [arcsec/pixel]
       stamp = 20.0/pixscale    ; [pixel]
    endelse
    
    nperpage = 16.0
    npage = ceil(nindx/nperpage)
    
;   ncols = 4.0 & nrows = 3.0
    ncols = sqrt(nperpage) & nrows = ncols
    xspace = 0.0 & yspace = 0.0
    xmargin = [0.2,0.2]-xspace*(ncols-1.0)*[1,1]/2.0
    ymargin = [0.2,0.2]-yspace*(nrows-1.0)*[1,1]/2.0

    psize = (8.5-total(xmargin)-total(xspace))/ncols
    width = replicate(psize,ncols)
    height = replicate(psize,nrows)

    xpage = total(width)+total(xmargin)+total(xspace)
    ypage = total(height)+total(ymargin)+total(yspace)

    psname = analysis_path+'sg1120_stamps.ps'
    if keyword_set(ps) then dfpsplot, psname, xsize=xpage, ysize=ypage, /color, bits=24

;   for jpage = 0L, 0L do begin
    for jpage = 0L, npage-1L do begin
  
       print, format='("Page ",I0,"/",I0,".",A10,$)', jpage+1, npage, string(13b)
          
       xspace1 = xspace & yspace1 = yspace
       xmargin1 = xmargin & ymargin1 = ymargin
       width1 = width & height1 = height
       xpage1 = xpage & ypage1 = ypage

;      set_plot, 'Z'
       arm_plotconfig, /landscape, nx=ncols, ny=nrows, xmargin=xmargin1, $
         ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
         height=height1, coord=pos, xpage=ypage1, ypage=xpage1, bw=0
 
       for kperpage = 0L, nperpage-1L do begin

          iobj = fix(jpage*nperpage+kperpage)
          if (iobj lt nindx) then begin

             if keyword_set(hststamps) then begin
                xcen = hstcat[indx[iobj]].x_image & ycen = hstcat[indx[iobj]].y_image
                image = (mrdfits(mosaicpath+'sg1120_hst.fits',0,/silent,range=ycen+stamp/2.0*[-1,1]))$
                  [xcen-stamp/2.0:xcen+stamp/2.0,*]
             endif else begin
                xcen = rcat[indx[iobj]].x_image & ycen = rcat[indx[iobj]].y_image
                image = rim[xcen-stamp/2.0:xcen+stamp/2.0,ycen-stamp/2.0:ycen+stamp/2.0]
             endelse

             imsize = size(image,/dimension)
             xsize = imsize[0] & xcen = xsize/2.0 & ysize = imsize[1] & ycen = ysize/2.0
             xaxis = (findgen(xsize)-xcen);*pixscale ; [arcsec]
             yaxis = (findgen(ysize)-ycen);*pixscale ; [arcsec]

             plotimage, asinhscl(image,negative=1,alpha=5.0,beta=16.0,omin=0,omax=250), /normal, $
;            plotimage, logscl(image,negative=1,exp=0.5,mean=1.0,omin=0,omax=255), /normal, $
               position=pos[*,kperpage], margin=0, imgxrange=minmax(xaxis), $
               imgyrange=minmax(yaxis), noerase=(kperpage gt 0L), $
               xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
               xsty=5, ysty=5;, /preserve
;            tvellipse, rcat[indx[iobj]].a_image, rcat[indx[iobj]].b_image, $
;              0.0, 0.0, rcat[indx[iobj]].theta_image, line=0, thick=1.0, /data

             legend, strtrim(string(zcat[indx[iobj]].z,format='(F12.4)'),2), $
               /right, /bottom, box=0, charsize=1.2, charthick=2.0, margin=0
;            legend, strtrim(string(zcat[indx[iobj]].oii,format='(F12.2)'),2), $
;              /left, /bottom, box=0, charsize=1.2, charthick=2.0, margin=0
;            legend, strtrim(string(gcat[indx[iobj]].mag_aper-rcat[indx[iobj]].mag_aper,format='(F8.2)'),2),$
;              /left, /bottom, box=0, charsize=1.5, charthick=2.0, margin=0
;            legend, strtrim(zcat[indx[iobj]].id,2), /left, /top, box=0, $
;              charsize=1.5, charthick=2.0, margin=0
             
          endif
             
       endfor

       if (not keyword_set(ps)) then cc = get_kbrd(1)

    endfor

    if keyword_set(ps) then dfpsclose
    cleanplot, /silent
    
stop    
    
return
end
    
