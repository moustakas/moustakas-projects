;+
; NAME:
;   PLOT_COSMICIMF_MF
;
; PURPOSE:
;   Plot the mass function results from COSMICIMF_MF.
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 05, UCSD
;
; Copyright (C) 2010, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function get_mstar, data, mstar_err=mstar_err
    mstar = data.mstar
    if tag_exist(data,'mstar_cv_err') then $
      mstar_err = sqrt(data.mstar_err^2+data.mstar_cv_err^2) else $
      mstar_err = data.mstar_err
return, mstar    
end

function get_phistar, data, phistar_err=phistar_err
    phistar = data.phistar
    if tag_exist(data,'phistar_cv_err') then $
      phistar_err = sqrt(data.phistar_err^2+data.phistar_cv_err^2+data.phistar_lss_err^2) else $
        phistar_err = data.phistar_err
; take the log
    phistar_err = phistar_err/phistar/alog(10)
    phistar = alog10(phistar)
return, phistar
end

function get_rho, data, rho_err=rho_err
    rho = data.rho
    if tag_exist(data,'rho_cv_err') then $
      rho_err = sqrt(data.rho_err^2+data.rho_cv_err^2+data.rho_lss_err^2) else $
      rho_err = data.rho_err
; take the log        
    rho_err = rho_err/rho/alog(10)
    rho = alog10(rho)
return, rho
end

pro plot_mf, mf_fit, mf_data, zbins, pos=pos, title=title
; simple wrapper to make plotting the MFs of the various subsamples
; (e.g., all, quiescent, active) easier 

    maxis1 = im_array(8.1,12.5,0.02)
    maxis2 = im_array(8.1,12.5,0.1)
    xrange = [9.4,12.3] ; [8.0,12.5]
    yrange = [1E-6,4E-2]
    
    cole = mf_cole()
    bell = mf_bell()
    for ii = 0, n_elements(mf_fit)-1 do begin
       if odd(ii) then begin
          ytitle = ''
          ytickname = replicate(' ',10)
       endif else begin
          if (ii eq 2) then $
            ytitle = cosmicimf_phititle() else delvarx, ytitle
          delvarx, ytickname 
       endelse

       if (ii lt 4) then begin
          xtitle = ''
          xtickname = replicate(' ',10)
       endif else begin
          delvarx, xtickname 
       endelse
       
       djs_plot, [0], [0], /nodata, noerase=(ii gt 0), $
         position=pos[*,ii], xsty=1, ysty=1, /ylog, $
         yrange=yrange, xrange=xrange, xtickname=xtickname, $
         ytickname=ytickname, xtitle=xtitle, ytitle=ytitle, $
         xtickinterval=1

; generate a shaded histogram showing all the errors including cosmic
; variance 
       phistar_err = sqrt(mf_fit[ii].phistar_err^2+mf_fit[ii].phistar_cv_err^2+$
         mf_fit[ii].phistar_lss_err^2)
       mstar_err = sqrt(mf_fit[ii].mstar_err^2+mf_fit[ii].mstar_cv_err^2)
       alpha_err = sqrt(mf_fit[ii].alpha_err^2+mf_fit[ii].alpha_cv_err^2)
       
       phimax1 = mf_schechter(maxis2,mf_fit[ii].phistar+phistar_err/2,$
         mf_fit[ii].mstar,mf_fit[ii].alpha)
       phimax2 = mf_schechter(maxis2,mf_fit[ii].phistar,$
         mf_fit[ii].mstar+mstar_err/2,mf_fit[ii].alpha)
       phimax3 = mf_schechter(maxis2,mf_fit[ii].phistar,$
         mf_fit[ii].mstar,mf_fit[ii].alpha+alpha_err/2)
       phimax = phimax1>phimax2>phimax3

       phimin1 = mf_schechter(maxis2,mf_fit[ii].phistar-phistar_err/2,$
         mf_fit[ii].mstar,mf_fit[ii].alpha)
       phimin2 = mf_schechter(maxis2,mf_fit[ii].phistar,$
         mf_fit[ii].mstar-mstar_err/2,mf_fit[ii].alpha)
       phimin3 = mf_schechter(maxis2,mf_fit[ii].phistar,$
         mf_fit[ii].mstar,mf_fit[ii].alpha-alpha_err/2)
       phimin = phimin1<phimin2<phimin3
       polyfill, [maxis2,reverse(maxis2)],[phimin,reverse(phimax)], $
         /data, /fill, color=fsc_color('tan',100), $
         noclip=0
       
; plot the mf_data, distinguishing between objects that are above and
; below the 75% completeness limit; only plot objects that are within
; 1.5 dex of the completeness limit
       mlimit = mf_data[ii].mass[(where(mf_data[ii].fullbin gt 0))[0]]
       these = where((mf_data[ii].phi gt 0.0) and $
         (mf_data[ii].mass gt mlimit-1.0) and $
         (mf_data[ii].mass gt xrange[0]),nthese)
       fullbin = mf_data[ii].fullbin[these]
       mass = mf_data[ii].mass[these]
       phi = mf_data[ii].phi[these]
       phierr = mf_data[ii].phierr[these]
       phierr_cv = mf_data[ii].phierr_cv[these]
;      phierr_model = mf_data[ii].phierr_model[these]
       
       above = where(fullbin,comp=below)
       oploterror, mass[above], phi[above], phierr[above], $
         psym=symcat(9,thick=!p.thick), symsize=1.2
       oploterror, mass[below], phi[below], phierr[below], $
         psym=symcat(9,thick=!p.thick), symsize=1.2, color=djs_icolor('grey'), $
         errcolor=djs_icolor('grey')

; overplot Cole+ and Bell+
;      djs_oplot, maxis1, mf_schechter(maxis1,cole.phistar,$
;        cole.mstar,cole.alpha), line=0, color=fsc_color(cole.color2,100)
       djs_oplot, maxis1, mf_schechter(maxis1,bell.phistar,$
         cole.mstar,bell.alpha), line=5, color=fsc_color(bell.color2,101)

; overplot the Schechter fit and the SDSS results for reference
;      djs_oplot, maxis1, mf_schechter(maxis1,mf_fit[0]), $
;        line=5, color='blue'
;      djs_oplot, maxis1, mf_schechter(maxis1,mf_fit[ii]), $
;        line=0, color='red'
       legend, 'z='+strtrim(string(zbins[ii].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[ii].zup,format='(F12.2)'),2), $
         /right, /top, box=0, margin=0, charsize=1.5
    endfor
    xyouts, pos[0,5], pos[1,5]-0.07, align=0.5, /norm, cosmicimf_masstitle()
    if (n_elements(title) ne 0) then xyouts, pos[2,0], pos[3,0]+0.02, $
      align=0.5, /norm, textoidl(title)

return
end

pro plot_cosmicimf_mf, ps=ps, keynote=keynote
; jm10mar23ucsd - plot the stellar mass functions in AGES

    common mf_sdss, sdss
    
    cosmicimfpath = ages_path(/projects)+'cosmicimf/'
    paperpath = ages_path(/papers)+'cosmicimf/'
    if keyword_set(keynote) then paperpath = $
      getenv('RESEARCHPATH')+'/meetings/10apr_florida/keynote/'

    zbins = cosmicimf_zbins(nzbins)
    peg = read_cosmicimf_sample(/pegase)
    imf = ['Salpeter',strtrim(peg.imf,2)]
    nimf = n_elements(imf)

    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

    maxis1 = im_array(7.5,12.5,0.01)

; --------------------------------------------------
; plot the fiducial (Salpeter) mass functions in six redshift bins 
    psfile = paperpath+'mf_salp'+suffix
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.6, $
      height=2.9*[1,1,1], ymargin=[0.6,1.1], keynote=keynote
    mf_data = mrdfits(cosmicimfpath+'mf_data.fits.gz',1,/silent)
    mf_fit = mrdfits(cosmicimfpath+'mf_fit.fits.gz',1,/silent)
    plot_mf, mf_fit, mf_data, zbins, pos=pos

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.eps','.pdf'), /sh
       rmfile, psfile
    endif

stop    
    
; --------------------------------------------------
; fit and plot the SDSS mass function at z=0.1 using the SFRM sample
; (temporary hack to get through my Florida talk!) 
    h100 = 0.7
    binsize = 0.1
    parinfo = init_mf_parinfo(fixslope=0,quiescent=quiescent,$
      active=active,/double_schechter)
    if (n_elements(sdss) eq 0) then sdss = $
      mrdfits(sdss_path(/lowz)+'lowz_catalog.dr6.fits.gz',1)
    vmax = sdss.vmax/h100^3.0 ; h=1-->h=0.7
    bmass = reform(k_sdss_bell(sdss.absmag[0:4])) ; stellar mass
    bmass = alog10(bmass) - alog10(h100^2) + 0.14 ; h=1-->0.7; diet Salpeter-->Salpeter
    
    mf_vmax, bmass, 1.0/vmax/binsize, binsize=binsize, $
      histmin=histmin, histmax=histmax, binmass=binmass, $
      phi=phi, errphi=phierr, minmass=8.3, fullbin=fullbin, $
      number=number
    full = where(fullbin)
    mf_fit_schechter_plus, binmass[full], phi[full], $
      phierr[full], fit, parinfo=parinfo

    psfile = paperpath+'sdss_mf'+suffix
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, $
      xmargin=[1.3,0.4], width=6.8, height=6.8

    if keyword_set(keynote) then keycolor = djs_icolor('white')

    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      /ylog, yrange=[1E-5,3E-2], xrange=[8.3,12.0], $
      xtitle=cosmicimf_masstitle(), ytitle=cosmicimf_phititle()
    above = where(fullbin and (number gt 3))
;   above = where(fullbin,comp=below)
    oploterror, binmass[above], phi[above], phierr[above], $
      psym=symcat(6,thick=6), symsize=1.4
;   oploterror, binmass[below], phi[below], phierr[below], $
;     psym=symcat(9,thick=!p.thick), symsize=1.2, color=djs_icolor('grey'), $
;     errcolor=djs_icolor('grey')
    djs_oplot, maxis1, mf_schechter_plus(maxis1,fit), $
      color=fsc_color('tomato',101), line=0, thick=8
;   oplot_bgd08_mf, maxis1, params=params, /salp, color='red', line=5

    im_legend, 'Salpeter IMF (0.1-100 M_{\odot})', /left, /bottom, $
      box=0, margin=0, textcolor=keycolor
    
    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.eps','.pdf'), /sh
       rmfile, psfile
    endif

; --------------------------------------------------
; M*, phi*, rho* vs redshift
    psfile = paperpath+'rhostar_vs_redshift'+suffix
    im_plotconfig, 4, pos, psfile=psfile, height=3.0*[1,1,1], $
      charsize=1.6, keynote=keynote

    mf_fit = mrdfits(cosmicimfpath+'mf_fit.fits.gz',1,/silent) ; Salpeter

    zrange = [0.0,0.75]
    mstarrange = [10.9,11.35]
    phistarrange = [-3.1,-2.2]
    rhorange = [8.13,8.9]

    psym1 = 15 ; 6
    symsize1 = 3.2 ; 4
    errthick1 = 8
    errthick2 = 5
    color1 = 'black' ; 'firebrick'

; literature data
    bell = mf_bell()
    cole = mf_cole()
    borch = mf_borch()
    perez = mf_perez()
;   bundy = mf_bundy()
;   pozzetti07 = mf_07pozzetti()
    
; ###############
; M*
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=mstarrange, xrange=zrange, xtitle='', xtickname=replicate(' ',10), $
      ytitle=cosmicimf_masstitle(/mstar), ytickinterval=0.2
    mstar = get_mstar(mf_fit,mstar_err=mstar_err)
    oploterror, zbins.zbin, mstar, zbins.zbin_err*0.0, mstar_err, $
      psym=-symcat(psym1,thick=errthick1), symsize=symsize1, $
      errthick=errthick1, color=fsc_color(color1,150), $
      errcolor=fsc_color(color1,150), thick=errthick1

; Borch+06    
    if (n_elements(borch) ne 0) then begin
       mstar = get_mstar(borch,mstar_err=mstar_err)
       oploterror, borch.z+0.015, mstar, borch.zerr, mstar_err, $
         symsize=borch.symsize, psym=symcat(borch.psym,thick=errthick1), errthick=errthick2, $
         color=fsc_color(borch.color1,101), errcolor=fsc_color(borch.color1,101)
    endif
; Perez+08
    mstar = get_mstar(perez,mstar_err=mstar_err)
    oploterror, perez.z-0.015, mstar, perez.zerr, mstar_err, $
      symsize=perez.symsize, psym=symcat(perez.psym,thick=errthick1), errthick=errthick2, $
      color=fsc_color(perez.color1,101), errcolor=fsc_color(perez.color1,101)
; Bundy+06    
    if (n_elements(bundy) ne 0) then begin
       mstar = get_mstar(bundy,mstar_err=mstar_err)
       oploterror, bundy.z, mstar, bundy.zerr, mstar_err, $
         symsize=bundy.symsize, psym=symcat(bundy.psym,thick=errthick1), errthick=errthick2, $
         color=fsc_color(bundy.color1,101), errcolor=fsc_color(bundy.color1,101)
    endif
; Pozzetti+07; M* at z~0.23 is not constrained, so don't plot it!  
    if (n_elements(pozzetti07) ne 0) then begin
       indx = lindgen(n_elements(mstar)-1)+1
       mstar = get_mstar(pozzetti07[indx],mstar_err=mstar_err)
       oploterror, pozzetti07[indx].z, mstar, pozzetti07[indx].zerr, mstar_err[indx], $
         symsize=pozzetti07.symsize, psym=symcat(pozzetti07.psym,thick=errthick1), errthick=errthick2, $
         color=fsc_color(pozzetti07.color1,101), errcolor=fsc_color(pozzetti07.color1,101)
    endif
; Bell+03
    mstar = get_mstar(bell,mstar_err=mstar_err)
    oploterror, bell.z-0.015, mstar, bell.zerr, mstar_err, $
      symsize=bell.symsize, psym=symcat(bell.psym,thick=errthick1), errthick=errthick2, $
      color=fsc_color(bell.color1,101), errcolor=fsc_color(bell.color1,101)
; Cole+01
    mstar = get_mstar(cole,mstar_err=mstar_err)
    oploterror, cole.z, mstar, cole.zerr, mstar_err, $
      symsize=cole.symsize, psym=symcat(cole.psym,thick=errthick1), errthick=errthick2, $
      color=fsc_color(cole.color1,101), errcolor=fsc_color(cole.color1,101)

; ###############
; phi* 
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      yrange=phistarrange, xrange=zrange, xtitle='', xtickname=replicate(' ',10), $
      ytitle=cosmicimf_phititle(/phistar,/log), ytickinterval=0.2
    phistar = get_phistar(mf_fit,phistar_err=phistar_err)
    oploterror, zbins.zbin, phistar, zbins.zbin_err*0.0, phistar_err, $
      psym=-symcat(psym1,thick=errthick1), symsize=symsize1, $
      errthick=errthick1, color=fsc_color(color1,150), $
      errcolor=fsc_color(color1,150), thick=errthick1

; Borch+06    
    if (n_elements(borch) ne 0) then begin
       phistar = get_phistar(borch,phistar_err=phistar_err)
       oploterror, borch.z+0.015, phistar, borch.zerr, phistar_err, $
         symsize=borch.symsize, psym=symcat(borch.psym,thick=errthick1), errthick=errthick2, $
         color=fsc_color(borch.color1,101), errcolor=fsc_color(borch.color1,101)
    endif
; Perez+08
    phistar = get_phistar(perez,phistar_err=phistar_err)
    oploterror, perez.z-0.015, phistar, perez.zerr, phistar_err, $
      symsize=perez.symsize, psym=symcat(perez.psym,thick=errthick1), errthick=errthick2, $
      color=fsc_color(perez.color1,101), errcolor=fsc_color(perez.color1,101)
; Bundy+06    
    if (n_elements(bundy) ne 0) then begin
       phistar = get_phistar(bundy,phistar_err=phistar_err)
       oploterror, bundy.z, phistar, bundy.zerr, phistar_err, $
         symsize=bundy.symsize, psym=symcat(bundy.psym,thick=errthick1), errthick=errthick2, $
         color=fsc_color(bundy.color1,101), errcolor=fsc_color(bundy.color1,101)
    endif
; Pozzetti07
    if (n_elements(pozzetti07) ne 0) then begin
       phistar = get_phistar(pozzetti07,phistar_err=phistar_err)
       oploterror, pozzetti07.z, phistar, pozzetti07.zerr, phistar_err, $
         symsize=pozzetti07.symsize, psym=symcat(pozzetti07.psym,thick=errthick1), errthick=errthick2, $
         color=fsc_color(pozzetti07.color1,101), errcolor=fsc_color(pozzetti07.color1,101)
    endif
; Bell+03
    phistar = get_phistar(bell,phistar_err=phistar_err)
    oploterror, bell.z-0.015, phistar, bell.zerr, phistar_err, $
      symsize=bell.symsize, psym=symcat(bell.psym,thick=errthick1), errthick=errthick2, $
      color=fsc_color(bell.color1,101), errcolor=fsc_color(bell.color1,101)
; Cole+01
    phistar = get_phistar(cole,phistar_err=phistar_err)
    oploterror, cole.z, phistar, cole.zerr, phistar_err, $
      symsize=cole.symsize, psym=symcat(cole.psym,thick=errthick1), errthick=errthick2, $
      color=fsc_color(cole.color1,101), errcolor=fsc_color(cole.color1,101)

; ###############
; rho*    
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,2], xsty=1, ysty=1, $
      yrange=rhorange, xrange=zrange, xtitle='Redshift', ytitle=cosmicimf_rhotitle(), $
      ytickinterval=0.2
    rho = get_rho(mf_fit,rho_err=rho_err)
    oploterror, zbins.zbin, rho, zbins.zbin_err*0.0, rho_err, $
      psym=-symcat(psym1,thick=errthick1), symsize=symsize1, $
      errthick=errthick1, color=fsc_color(color1,150), $
      errcolor=fsc_color(color1,150), thick=errthick1

;; fit a line to the z>0.2 points plus Cole+Bell and overplot
;    zz = [z_cole,z_bell03,zbins[1:nzbins-1].zbin]
;    rho = [rhostar_cole,rhostar_bell03,rhostar[1:nzbins-1]]
;    rho_err = [rhostarerr_cole,rhostarerr_bell03,rhostar_err[1:nzbins-1]]
;
;    pivot = 0.0 ; 0.4
;    coeff = im_linefit(zz-pivot,rho,yerr=rho_err,chi2=chi2,coeff_err=coeff_err)
;    zaxis = im_array(0.02,0.7,0.02)
;
;    rhostarfit = poly(zaxis-pivot,coeff)
;    rhostarmax = poly(zaxis-pivot,coeff+[0.0,coeff_err[1]])
;    rhostarmin = poly(zaxis-pivot,coeff-[0.0,coeff_err[1]])
;    polyfill, [zaxis,reverse(zaxis)],[rhostarmin,reverse(rhostarmax)], $
;      /data, /fill, color=fsc_color('tan',100), noclip=0
;;   djs_oplot, zaxis, rhostarfit, line=0, color='grey', thick=3
;;   djs_oplot, zaxis, rhostarmax, line=0;, color='grey', thick=3
;;   djs_oplot, zaxis, rhostarmin, line=0;, color='grey', thick=3
    
; Borch+06    
    rho = get_rho(borch,rho_err=rho_err)
    oploterror, borch.z+0.015, rho, borch.zerr, rho_err, $
      symsize=borch.symsize, psym=symcat(borch.psym,thick=errthick1), errthick=errthick2, $
      color=fsc_color(borch.color1,101), errcolor=fsc_color(borch.color1,101)
; Perez+08
    rho = get_rho(perez,rho_err=rho_err)
    oploterror, perez.z-0.015, rho, perez.zerr, rho_err, $
      symsize=perez.symsize, psym=symcat(perez.psym,thick=errthick1), errthick=errthick2, $
      color=fsc_color(perez.color1,101), errcolor=fsc_color(perez.color1,101)
; Bundy+06    
    if (n_elements(bundy) ne 0) then begin
       rho = get_rho(bundy,rho_err=rho_err)
       oploterror, bundy.z, rho, bundy.zerr, rho_err, $
         symsize=bundy.symsize, psym=symcat(bundy.psym,thick=errthick1), errthick=errthick2, $
         color=fsc_color(bundy.color1,101), errcolor=fsc_color(bundy.color1,101)
    endif
; Pozzetti07
    if (n_elements(pozzetti07) ne 0) then begin
       rho = get_rho(pozzetti07,rho_err=rho_err)
       oploterror, pozzetti07.z, rho, pozzetti07.zerr, rho_err, $
         symsize=pozzetti07.symsize, psym=symcat(pozzetti07.psym,thick=errthick1), errthick=errthick2, $
         color=fsc_color(pozzetti07.color1,101), errcolor=fsc_color(pozzetti07.color1,101)
    endif
; Bell+03    
    rho = get_rho(bell,rho_err=rho_err)
    oploterror, bell.z-0.015, rho, bell.zerr, rho_err, $
      symsize=bell.symsize, psym=symcat(bell.psym,thick=errthick1), errthick=errthick2, $
      color=fsc_color(bell.color1,101), errcolor=fsc_color(bell.color1,101)
; Cole+01
    rho = get_rho(cole,rho_err=rho_err)
    oploterror, cole.z, rho, cole.zerr, rho_err, $
      symsize=cole.symsize, psym=symcat(cole.psym,thick=errthick1), errthick=errthick2, $
      color=fsc_color(cole.color1,101), errcolor=fsc_color(cole.color1,101)

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.eps','.pdf'), /sh
       rmfile, psfile
    endif

stop    
    
; --------------------------------------------------
; mass functions in six redshift bins - test plot (this does not work
; as an EPS file!)
    psfile = paperpath+'mf_allimfs.ps'
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.6, $
      height=2.9*[1,1,1], ymargin=[0.6,1.1], keynote=keynote
;   for ii = 0, nimf-1 do begin
    for ii = 0, 0 do begin
       mf_data = mrdfits(cosmicimfpath+'mf_data.fits.gz',ii+1,/silent)
       mf_fit = mrdfits(cosmicimfpath+'mf_fit.fits.gz',ii+1,/silent)
       plot_mf, mf_fit, mf_data, zbins, pos=pos, $
         title=repstr(imf[ii],'cosmicimf_','\alpha_{2}=')
    endfor
    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.eps','.pdf'), /sh
       rmfile, psfile
    endif

return
end
    


;; --------------------------------------------------
;; rhostar versus redshift
;    psfile = paperpath+'rhostar_vs_redshift'+suffix
;    im_plotconfig, 0, pos, psfile=psfile
;
;    xtitle = 'Redshift'
;    ytitle = cosmicimf_rhotitle()
;    xrange = [-0.02,0.75]
;    yrange = [8.3,9.2]
;
;    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;      yrange=yrange, xrange=xrange, xtitle=xtitle, ytitle=ytitle
;    for ii = 0, nimf-1 do begin
;       mf_data = mrdfits(cosmicimfpath+'mf_data.fits.gz',ii+1,/silent)
;       oploterror, zbins.zbin, alog10(mf_fit.rho), mf_fit.rho_err/mf_fit.rho/alog(10), $
;         psym=-symcat(16)
;    endfor
;    im_plotconfig, /psclose
;

;; --------------------------------------------------
;; plot some local mass functions for comparison
;    massaxis = im_array(7.0,12.0,0.02)
;    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=[9,12], yrange=[-4.8,-1]
;; Borch+06
;    djs_oplot, massaxis, alog10(mf_schechter(10^massaxis,35D-4,10^10.91,-1.1)), color='grey' ; local
;    djs_oplot, massaxis, alog10(mf_schechter(10^massaxis,19D-4,10^11.03,-1.1)), color='orange' ; z~0.3
;    djs_oplot, massaxis, alog10(mf_schechter(10^massaxis,18D-4,10^11.02,-1.1)), color='blue'   ; z~0.5
;    djs_oplot, massaxis, alog10(mf_schechter(10^massaxis,16D-4,10^11.09,-1.1)), color='red'    ; z~0.7
;    djs_oplot, massaxis, alog10(mf_schechter(10^massaxis,12D-4,10^11.08,-1.1)), color='green'  ; z~0.9
;    
;    print, alog10(total(0.02*10^massaxis*mf_schechter(10^massaxis,35D-4,10^10.91,-1.1)))

