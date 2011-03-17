function zpt_mpfunc, x, p, color=color, mag=mag
; zeropoint function for MPFIT    
; X and COLOR are SDSS magnitudes
    model = x + p[0] + p[1]*color
return, model
end

function ndwfs_zpt_mpfit, sdssmag, mag, magerr, color=color
; fit for the zeropoint
    parinfo1 = replicate({value: 1.0, fixed: 0},2)
    coeff = mpfitfun('zpt_mpfunc',sdssmag,mag,magerr,$
      parinfo=parinfo1,functargs={color: color, mag:mag},$
      yfit=bfit,/quiet)
return, coeff
end

pro ndwfs_zpt_scatterplot, sdssmag, catmag, color, coeff, $
  pos, equation=equation, xrange1=xrange1, yrange1=yrange1, $
  xrange2=xrange2, yrange2=yrange2, xrange3=xrange3, $
  yrange3=yrange3, xtitle1=xtitle1, xtitle2=xtitle2, $
  xtitle3=xtitle3, ytitle1=ytitle1, ytitle2=ytitle2, $
  ytitle3=ytitle3
  ; scatterplot for the zeropoint code

    resid1 = catmag - (sdssmag + coeff[0])
    resid2 = catmag - (sdssmag + coeff[0] + coeff[1]*color)
    
    xaxis1 = findgen((xrange1[1]-xrange1[0])/0.01+1)*0.01+xrange1[0]
    xaxis2 = findgen((xrange2[1]-xrange2[0])/0.01+1)*0.01+xrange2[0]

    plotsym, 0, 0.3, /fill
; magnitude vs magnitude
    djs_plot, sdssmag, catmag, psym=8, xsty=1, ysty=1, $
      xtitle=xtitle1, ytitle=ytitle1, xrange=xrange1, $
      yrange=yrange1, position=pos[*,0], charsize=1.7
;   djs_oplot, xaxis1, xaxis1, line=5, thick=3.0
    djs_oplot, xaxis1, xaxis1-coeff[0], line=0, color='red'
    legend, textoidl('A1 = '+strtrim(string(coeff[0],format='(F12.3)'),2)), $
      /right, /bottom, box=0, charsize=1.7
    legend, textoidl(equation), /left, /top, box=0, charsize=1.5
; color residuals
    djs_plot, color, resid1, psym=8, xsty=1, ysty=1, $
      xtitle=xtitle2, ytitle=ytitle2, xrange=xrange2, $
      yrange=yrange2, position=pos[*,1], /noerase, charsize=1.7
    djs_oplot, xaxis2, xaxis2*0.0, line=5, thick=3.0
    djs_oplot, xaxis2, poly(xaxis2,[0.0,coeff[1]]), line=0, color='red'
    legend, textoidl('A2 = '+strtrim(string(coeff[1],format='(F12.4)'),2)), $
      /right, /bottom, box=0, charsize=1.7
; magnitude residuals
    sig = djsig(resid2,sigrej=4.5)
    djs_plot, sdssmag, resid2, psym=8, xsty=1, ysty=1, $
      xtitle=xtitle3, ytitle=ytitle3, xrange=xrange3, $
      yrange=yrange3, position=pos[*,2], /noerase, charsize=1.7
    djs_oplot, xaxis1, xaxis1*0.0, line=5, thick=3.0
    legend, textoidl('\sigma = '+strtrim(string(sig,format='(F12.3)'),2)), $
      /right, /top, box=0, charsize=1.5
    
return
end    

pro ndwfs_zeropoints
; jm09aug25ucsd - calibrate the NDWFS photometry
    
    common calibrate_ndwfs, ndwfs1, sdss1

    photodir = getenv('RESEARCHPATH')+'/data/ndwfs/'

; read and match the NDWFS and SDSS catalogs    
    if (n_elements(ndwfs1) eq 0) then begin
       photofile = file_search(photodir+'NDWFS_??_??.fits.gz',count=nfile)
       for ii = 0, nfile-1 do begin
          temp = mrdfits(photofile[ii],1)
          keep = where((temp.bw_flags eq 0) and (temp.r_flags eq 0) and $
            (temp.i_flags eq 0),nkeep)
          splog, nkeep
          temp = temp[keep]
          if (ii eq 0) then ndwfs1 = temporary(temp) else $
            ndwfs1 = [temporary(ndwfs1),temporary(temp)]
       endfor
    endif

    if (n_elements(sdss1) eq 0) then begin
       photofile = photodir+'NDWFS_SDSS_stars.fits.gz'
       sdss1 = mrdfits(photofile,1)
    endif

    sdss_filterlist = sdss_filterlist()
    bwrik_filterlist = ndwfs_filterlist()
    bwri_filterlist = bwrik_filterlist[0:2]

; spherematch, convert to maggies, and select a fiducial sample of
; objects with good multiband photometry
    spherematch, ndwfs1.i_alpha_j2000, ndwfs1.i_delta_j2000, $
      sdss1.ra, sdss1.dec, 1.0/3600.0, m1, m2

    sdss = ndwfs_zpt_sdsscat2mag(sdss1[m2],maggies=sdssmaggies,$
      ivarmaggies=sdssivarmaggies)
    ndwfs = zpt_ndwfscat2mag(ndwfs1[m1])

    these = where((ndwfs.bw gt 0.0) and (ndwfs.r gt 0.0) and $
      (ndwfs.i gt 0.0) and (ndwfs.bw_flags eq 0) and $
      (ndwfs.r_flags eq 0) and (ndwfs.i_flags eq 0) and $
      (sdss.u gt 0.0) and (sdss.g gt 0.0) and $
      (sdss.r gt 0.0) and (sdss.i gt 0.0) and (sdss.z gt 0.0) and $
      (sdss.r ge 18.0) and (sdss.r le 20.5) and $
      (sdss.g-sdss.r le 1.1))
    sdss = sdss[these]
    ndwfs = ndwfs[these]
    sdssmaggies = sdssmaggies[*,these]
    sdssivarmaggies = sdssivarmaggies[*,these]

; ---------------------------------------------------------------------------
; zeropoint test 1 - fit Kurucz models to all the stars and compare
; the synthesized vs observed BwRI photometry

    kall = fit_kurucz_models(sdssmaggies,sdssivarmaggies,$
      filterlist=sdss_filterlist)
    these = where(kall.kurucz_chi2min lt 3.0,nstar)
    kthese = kall[these]

;   im_plothist, sdss.g-sdss.r, bin=0.05
;   im_plothist, sdss[these].g-sdss[ww].r, /over, /fill, bin=0.05

    bwri = fltarr(3,nstar)
    for ii = 0L, nstar-1L do bwri[*,ii] = reform(k_project_filters($
      k_lambda_to_edges(kthese[ii].lambda),kthese[ii].spec,$
      filterlist=bwri_filterlist))

    psfile = photodir+'ndwfs_kurucz_zeropoints.ps'
    im_plotconfig, 6, pos1, psfile=psfile, xmargin=[1.4,0.3], $
      height=[4.5,3.0], width=6.8

    mrange = [[18.0,22.5],[17.0,21.5],[17.0,21.5]]
    rrange = 0.39*[-1,1]

    band = ['Bw','R','I']
    for ii = 0, 2 do begin
       magtag = tag_indx(ndwfs[0],band[ii])
       xx = ndwfs[these].(magtag)
       yy = reform(-2.5*alog10(bwri[ii,*]))
       djs_plot, xx, yy, position=pos1[*,0], $
         xsty=1, ysty=1, xrange=mrange[*,ii], yrange=mrange[*,ii], $
         psym=4, xtitle='', ytitle=band[ii]+' (SDSS synthesized, AB mag)', $
         xtickname=replicate(' ',10)
       djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
;      im_legend, band[ii]+'_{SDSS}-'+band[ii]+'_{NDWFS} = '+$
       im_legend, '<\Delta'+band[ii]+'> = '+$
         im_string_stats(yy-xx,sigrej=3.0), /left, /top, box=0
       djs_plot, xx, yy-xx, position=pos1[*,1], /noerase, xsty=1, ysty=1, $
         xrange=mrange[*,ii], yrange=rrange, psym=4, xtitle=band[ii]+' (NDWFS, AB mag)', $
         ytitle='Residuals (AB mag)'
       djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
    endfor
       
    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep

; show the spectra themselves    
    psfile = photodir+'ndwfs_kurucz_sedfits.ps'
    im_plotconfig, 8, pos1, psfile=psfile

    weff = k_lambda_eff(filterlist=sdss_filterlist())
    ndwfs_weff = k_lambda_eff(filterlist=bwri_filterlist)

    light = 2.99792458D18       ; speed of light [A/s]
    xrange1 = [3050.0,9500.0]
    xtitle1 = 'Wavelength (\AA)'
    ytitle1 = 'm_{AB}'
    
;   for ii = 0, 5 do begin
    for ii = 0, nstar-1 do begin
       wave = kthese[ii].lambda
       flux = kthese[ii].spec*wave^2/light
       flux = -2.5*alog10(flux>1D-50)-48.6

       ndwfs_mab = [ndwfs[these[ii]].bw,ndwfs[these[ii]].r,$
         ndwfs[these[ii]].i]
       mab = -2.5*alog10(kthese[ii].primus_maggies)
       bestmab = -2.5*alog10(kthese[ii].kurucz_maggies)
          
       get_element, wave, xrange1, xx
       yrange = fltarr(2)
       yrange[0] = (max(bestmab)>max(flux[xx[0]:xx[1]]))*1.02
       yrange[1] = (min(bestmab)<min(flux[xx[0]:xx[1]]))*0.95

       djs_plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
         xsty=1, ysty=1, xtitle=xtitle1, ytitle=ytitle1
       djs_oplot, wave, flux, line=0, color='grey'
       djs_oplot, weff, bestmab, psym=symcat(6,thick=6), $
         symsize=2.5

       mab = maggies2mag(kthese[ii].primus_maggies,$
         ivar=kthese[ii].primus_maggiesivar,magerr=mab_err)
       oploterror, weff, mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=djs_icolor('dark green'), $
         errcolor=djs_icolor('dark green'), errthick=!p.thick

       djs_oplot, ndwfs_weff, ndwfs_mab, psym=symcat(4,thick=6), $
         symsize=2.5, color='red'

       im_legend, [$
         '[Fe/H]='+string(kthese[ii].kurucz_feh,format='(F4.1)'),$
         'T_{eff}='+string(kthese[ii].kurucz_teff,format='(I0)'),$
         'log (g)='+string(kthese[ii].kurucz_g,format='(F3.1)'),$
         '\chi^{2}_{\nu}='+string(kthese[ii].kurucz_chi2min/4.0,format='(F3.1)')], $
         /left, /top, box=0, charsize=1.4

       im_legend, ['Kurucz ugriz','SDSS ugriz','NDWFS BwRI'], $
         psym=[6,16,4], color=['','dark green','red'], /right, $
         /bottom, box=0, charsize=1.4, symthick=3.0
    endfor
    
    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep

stop    
    
; ---------------------------------------------------------------------------    
    
; fit for the color terms and zeropoints
    cut = where(sdss.g-sdss.r le 1.2)
    bfit_ug1 = ndwfs_zpt_mpfit(sdss[cut].u,ndwfs[cut].bw,$
      ndwfs[cut].bwerr,color=sdss[cut].u-sdss[cut].g)
    bfit_ug2 = ndwfs_zpt_mpfit(sdss[cut].g,ndwfs[cut].bw,$
      ndwfs[cut].bwerr,color=sdss[cut].u-sdss[cut].g)
    bfit_gr = ndwfs_zpt_mpfit(sdss[cut].g,ndwfs[cut].bw,$
      ndwfs[cut].bwerr,color=sdss[cut].g-sdss[cut].r)
    rfit_gr = ndwfs_zpt_mpfit(sdss[cut].r,ndwfs[cut].r,$
      ndwfs[cut].rerr,color=sdss[cut].g-sdss[cut].r)
    rfit_ri = ndwfs_zpt_mpfit(sdss[cut].r,ndwfs[cut].r,$
      ndwfs[cut].rerr,color=sdss[cut].r-sdss[cut].i)
    ifit_ri = ndwfs_zpt_mpfit(sdss[cut].i,ndwfs[cut].i,$
      ndwfs[cut].ierr,color=sdss[cut].r-sdss[cut].i)
    ifit_iz = ndwfs_zpt_mpfit(sdss[cut].i,ndwfs[cut].i,$
      ndwfs[cut].ierr,color=sdss[cut].i-sdss[cut].z)

; finally make a QAplot
    psfile = photodir+'ndwfs_zeropoints.ps'
    im_plotconfig, 4, pos, yspace=[1.0,1.0], ymargin=[0.3,0.9], $
      height=[2.6,2.6,2.6], ypage=11.0, xmargin=[1.4,0.3], $
      xspace=0.0, width=6.8, xpage=8.5, psfile=psfile
; Bw/ug, u
    xrange1 = [18.5,22.5] & yrange1 = [18.5,22.5]
    xrange2 = [0.0,3.0]  & yrange2 = [-1,1]
    xrange3 = [18.5,22.5] & yrange3 = [-0.5,0.5]
    ndwfs_zpt_scatterplot, sdss.u, ndwfs.bw, sdss.u-sdss.g, bfit_ug1, pos, $
      equation='Bw_{NDWFS}=u_{SDSS}+A1+A2*(u-g)_{SDSS}', $
      xrange1=xrange1, yrange1=yrange1, xrange2=xrange2, $
      yrange2=yrange2, xrange3=xrange3, yrange3=yrange3, $
      xtitle1='u_{SDSS}', ytitle1='Bw_{NDWFS}', $
      xtitle2='(u - g)_{SDSS}', ytitle2='Bw_{NDWFS} - (u_{SDSS} + A1)', $
      xtitle3='u_{SDSS}', ytitle3='Residuals'
; Bw/ug, g
    xrange1 = [18.5,22.5] & yrange1 = [18.5,22.5]
    xrange2 = [0.0,3.0]  & yrange2 = [-1,1]
    xrange3 = [18.5,22.5] & yrange3 = [-0.5,0.5]
    ndwfs_zpt_scatterplot, sdss.g, ndwfs.bw, sdss.u-sdss.g, bfit_ug2, pos, $
      equation='Bw_{NDWFS}=g_{SDSS}+A1+A2*(u-g)_{SDSS}', $
      xrange1=xrange1, yrange1=yrange1, xrange2=xrange2, $
      yrange2=yrange2, xrange3=xrange3, yrange3=yrange3, $
      xtitle1='g_{SDSS}', ytitle1='Bw_{NDWFS}', $
      xtitle2='(u - g)_{SDSS}', ytitle2='Bw_{NDWFS} - (g_{SDSS} + A1)', $
      xtitle3='g_{SDSS}', ytitle3='Residuals'
; Bw/gr
    xrange1 = [18.5,22.5] & yrange1 = [18.5,22.5]
    xrange2 = [-0.3,1.3]  & yrange2 = [-1,1]
    xrange3 = [18.5,22.5] & yrange3 = [-0.5,0.5]
    ndwfs_zpt_scatterplot, sdss.g, ndwfs.bw, sdss.g-sdss.r, bfit_gr, pos, $
      equation='Bw_{NDWFS}=g_{SDSS}+A1+A2*(g-r)_{SDSS}', $
      xrange1=xrange1, yrange1=yrange1, xrange2=xrange2, $
      yrange2=yrange2, xrange3=xrange3, yrange3=yrange3, $
      xtitle1='g_{SDSS}', ytitle1='Bw_{NDWFS}', $
      xtitle2='(g - r)_{SDSS}', ytitle2='Bw_{NDWFS} - (g_{SDSS} + A1)', $
      xtitle3='g_{SDSS}', ytitle3='Residuals'
; R/gr
    xrange1 = [18.0,21.0] & yrange1 = [18.0,21.0]
    xrange2 = [-0.3,1.3]  & yrange2 = [-1,1]
    xrange3 = [18.0,21.0] & yrange3 = [-0.5,0.5]
    ndwfs_zpt_scatterplot, sdss.r, ndwfs.r, sdss.g-sdss.r, rfit_gr, pos, $
      equation='R_{NDWFS}=r_{SDSS}+A1+A2*(g-r)_{SDSS}', $
      xrange1=xrange1, yrange1=yrange1, xrange2=xrange2, $
      yrange2=yrange2, xrange3=xrange3, yrange3=yrange3, $
      xtitle1='r_{SDSS}', ytitle1='R_{NDWFS}', $
      xtitle2='(g - r)_{SDSS}', ytitle2='R_{NDWFS} - (r_{SDSS} + A1)', $
      xtitle3='r_{SDSS}', ytitle3='Residuals'
; R/ri
    xrange1 = [17.5,21.0] & yrange1 = [17.5,21.0]
    xrange2 = [-0.2,0.7]  & yrange2 = [-1,1]
    xrange3 = [17.5,21.0] & yrange3 = [-0.5,0.5]
    ndwfs_zpt_scatterplot, sdss.r, ndwfs.r, sdss.r-sdss.i, rfit_ri, pos, $
      equation='R_{NDWFS}=r_{SDSS}+A1+A2*(r-i)_{SDSS}', $
      xrange1=xrange1, yrange1=yrange1, xrange2=xrange2, $
      yrange2=yrange2, xrange3=xrange3, yrange3=yrange3, $
      xtitle1='r_{SDSS}', ytitle1='R_{NDWFS}', $
      xtitle2='(r - i)_{SDSS}', ytitle2='R_{NDWFS} - (r_{SDSS} + A1)', $
      xtitle3='r_{SDSS}', ytitle3='Residuals'
; I/ri
    xrange1 = [17.5,21.0] & yrange1 = [17.5,21.0]
    xrange2 = [-0.2,0.7]  & yrange2 = [-1,1]
    xrange3 = [17.5,21.0] & yrange3 = [-0.5,0.5]
    ndwfs_zpt_scatterplot, sdss.i, ndwfs.i, sdss.r-sdss.i, ifit_ri, pos, $
      equation='I_{NDWFS}=i_{SDSS}+A1+A2*(r-i)_{SDSS}', $
      xrange1=xrange1, yrange1=yrange1, xrange2=xrange2, $
      yrange2=yrange2, xrange3=xrange3, yrange3=yrange3, $
      xtitle1='i_{SDSS}', ytitle1='I_{NDWFS}', $
      xtitle2='(r - i)_{SDSS}', ytitle2='I_{NDWFS} - (i_{SDSS} + A1)', $
      xtitle3='i_{SDSS}', ytitle3='Residuals'
; I/iz
    xrange1 = [17.5,21.0] & yrange1 = [17.5,21.0]
    xrange2 = [-0.3,0.5]  & yrange2 = [-1,1]
    xrange3 = [17.5,21.0] & yrange3 = [-0.5,0.5]
    ndwfs_zpt_scatterplot, sdss.i, ndwfs.i, sdss.i-sdss.z, ifit_iz, pos, $
      equation='I_{NDWFS}=i_{SDSS}+A1+A2*(i-z)_{SDSS}', $
      xrange1=xrange1, yrange1=yrange1, xrange2=xrange2, $
      yrange2=yrange2, xrange3=xrange3, yrange3=yrange3, $
      xtitle1='i_{SDSS}', ytitle1='I_{NDWFS}', $
      xtitle2='(i - z)_{SDSS}', ytitle2='I_{NDWFS} - (i_{SDSS} + A1)', $
      xtitle3='i_{SDSS}', ytitle3='Residuals'

    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep

stop    
    
return
end
    
