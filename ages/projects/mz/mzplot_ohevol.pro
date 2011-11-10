function sfrd_frac, sfr, inzaxis, outzaxis, plotzaxis
; jm09apr17nyu - compute the cumulative fraction of stellar mass (or
; metals) formed (synthesized) versus time/redshift

    lookback = getage(0.0)-getage(outzaxis) ; lookback time
    ntime = n_elements(lookback)
    
    sfr = interpol(sfr,inzaxis,outzaxis)
    sfrd_total = im_integral(lookback,sfr*1D9,$ ; [M_sun]
      min(lookback),max(lookback)) 
    sfrd_time = im_integral(lookback,sfr*1D9,$  ; [M_sun]
      lookback,replicate(max(lookback),ntime))
    sfrd_frac = sfrd_time/sfrd_total

; extrapolate to z<0.1    
;   plotzaxis = im_array(0.0,0.6,0.01)
    plotzaxis = [0.0,outzaxis]
    notzero = where(sfrd_frac gt 0.0)
    sfrd_frac = interpol(1.0*(sfrd_frac[notzero]),$
      outzaxis[notzero],plotzaxis)
    
return, sfrd_frac
end

function zformpeg, peg, zform=pegzform, peg_zaxis=peg_zaxis
; jm09apr17nyu - attach a formation redshift to a Pegase model

    ppeg = peg
    rev = reverse(sort(ppeg.age)) ; reverse the time array!
    ppeg = ppeg[rev]

    pegtform = getage(pegzform)
    peg_zaxis = getredshift(ppeg.age/1D3+pegtform)
    peg_lookback = getage(0.0)-getage(peg_zaxis)

    peg_good = where(peg_zaxis gt 0.0)
    peg_zaxis = peg_zaxis[peg_good]
    peg_lookback = peg_lookback[peg_good]
    ppeg = ppeg[peg_good]
    
return, ppeg
end
    
pro mzplot_ohevol, ps=ps
; jm09mar27nyu - plots of the evolution of (O/H)
; jm10nov02ucsd - major update

    litpath = getenv('IM_PAPERS_DIR')+'/literature/data/'

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

    masscut = 9.5D
    zbins = mz_zbins(nz)

; --------------------------------------------------
; Figure 17 - plot the evolution of the MZ relation on a single panel
    calib = 't04'
    mzlocal = mrdfits(mzpath+'mzlocal_sdss_ews_'+calib+'.fits.gz',1)
    mzevol = mrdfits(mzpath+'mzevol_'+calib+'.fits.gz',1,/silent)

    psfile = pspath+'mzevol_model'+suffix
    im_plotconfig, 0, pos, psfile=psfile, height=4.7, width=6.5, $
      xmargin=[1.6,0.4]

    xrange = [9.6,11.2]
    if calib eq 't04' then yrange = [8.75,9.12]
    if calib eq 'kk04' then yrange = [8.9,9.15]
    xtitle1 = mzplot_masstitle()
    ytitle1 = mzplot_ohtitle(kk04=calib eq 'kk04',t04=calib eq 't04',m91=calib eq 'm91')
    maxis = range(xrange[0]+0.02,xrange[1]-0.1,50)

;   zval = zbins.zbin
    zval = mzavg.medz[1:nz]
    zlabel = 'z='+string(zval,format='(F4.2)')
    zline = [0,3,5,1,4,2]
    zcolor = ['black','firebrick','dodger blue','tan','forest green','orange']
    below_zpsym = [ 9, 5, 6,45, 4,11]
    above_zpsym = [16,17,15,46,14,18]
    zsymsize = [1.0,1.2,1.0,1.4,1.0,1.1]*1.8

    zthese = lindgen(nz)
;   zthese = [0,2,4]
    nzthese = n_elements(zthese)
    
    djs_plot, [0], [0], /nodata, ysty=1, xsty=1, position=pos, $
      xrange=xrange, yrange=yrange, xtitle=xtitle1, ytitle=ytitle1
    im_legend, zlabel[zthese], /right, /bottom, box=0, line=zline[zthese], $
      color=zcolor[zthese], pspacing=2, thick=6, margin=0, $
      charsize=1.4, psym=-above_zpsym, symsize=zsymsize*0.7

; plot the data
    plotdata = 1
    if plotdata then begin
;      djs_oplot, mzlocal.bin_mass, mzlocal.bin_oh, psym=symcat(16), symsize=0.7
;      djs_oplot, mzlocal.bin_mass, mzlocal.bin_oh_mean, psym=symcat(16), symsize=0.7, color='purple'

       good = where(mzevol.sdss_ohmean_bymass gt -900 and mzevol.sdss_medmass_bymass gt xrange[0]+0.05)
;      oploterror, mzevol.sdss_medmass_bymass[good], mzevol.sdss_ohmean_bymass[good], $
;        mzevol.sdss_ohmean_bymass_err[good], psym=symcat(6,thick=5), symsize=zsymsize[0], color=im_color(zcolor[0]), $
;        errcolor=im_color(zcolor[0])

       for iz = 0, nzthese-1 do begin
          above = where(mzevol.complete_bymass[*,zthese[iz]] eq 1 and mzevol.ohmean_bymass[*,zthese[iz]] gt -900 and $
            mzevol.medmass_bymass[*,zthese[iz]] gt xrange[0]+0.05,nabove)
          niceprint, mzevol.ngal_bymass[above,zthese[iz]], mzevol.medmass_bymass[above,zthese[iz]], $
            mzevol.ohmean_bymass[above,zthese[iz]] & print
          if (nabove ne 0) then begin
             oploterror, mzevol.medmass_bymass[above,zthese[iz]], mzevol.ohmean_bymass[above,zthese[iz]], $
               mzevol.ohmean_bymass_err[above,zthese[iz]], psym=symcat(above_zpsym[zthese[iz]],thick=5), $
               symsize=zsymsize[zthese[iz]], color=im_color(zcolor[zthese[iz]]), $
               errcolor=im_color(zcolor[zthese[iz]])
          endif

          below = where(mzevol.complete_bymass[*,zthese[iz]] eq 0 and mzevol.ohmean_bymass[*,zthese[iz]] gt -900,nbelow)
          if (nbelow ne 0) then begin
             oploterror, mzevol.medmass_bymass[below,zthese[iz]], mzevol.ohmean_bymass[below,zthese[iz]], $
               mzevol.ohmean_bymass_err[below,zthese[iz]], psym=symcat(below_zpsym[zthese[iz]],thick=5), $
               symsize=zsymsize[zthese[iz]], color=im_color(zcolor[zthese[iz]]), $
               errcolor=im_color(zcolor[zthese[iz]])
          endif
       endfor 
    endif
       
; and now plot the model
    ww = where(strtrim(mzavg.calib,2) eq calib)
    for iz = 0, nzthese-1 do begin
;      djs_oplot, maxis, mzevol_func(maxis,mzavg.mzevol_coeffs[*,ww],$
;        z=zval[zthese[iz]],qz0=mzavg.qz0), line=zline[zthese[iz]], color=im_color(zcolor[zthese[iz]]), $
;        thick=2
       djs_oplot, maxis, mzevol_func(maxis,mzavg.mzevol_coeffs_r0zero[*,ww],$
         z=zval[zthese[iz]],qz0=mzavg.qz0), line=zline[zthese[iz]], color=im_color(zcolor[zthese[iz]]), $
         thick=10
    endfor
    
    im_plotconfig, /psclose, psfile=psfile

stop    
    
; --------------------------------------------------
; Figure 15 - metallicity evolution rate (slope) vs stellar mass
;   slopemass = massbins[0:nmassbins-2].massbin
    masscut = 9.5D
    massbins = mz_massbins(/rev,masscut=masscut,masskeep=masskeep)
    nmassbins = n_elements(massbins)

    calib = strtrim(mzavg.calib,2)
    calibcolor = ['firebrick','forest green','dodger blue']
    calibpsym = [16,17,15]
;   calibpsym = [9,5,6]
    calibsymsize = [1.5,1.8,1.5]*0.8
    calibline = [0,3,5]
    off = [+0.02,-0.01,+0.01]

    psfile = pspath+'mass_vs_dlogohdz'+suffix
    im_plotconfig, 0, pos, psfile=psfile, height=4.7, width=6.5, $
      xmargin=[1.6,0.4]

    xrange = [9.6,11.2]
    yrange = [-0.44,-0.05]
;   yrange = [-0.4,-0.1]
    xtitle1 = mzplot_masstitle()
    ytitle1 = textoidl('d[log(O/H)] / dz (dex z^{-1})')
    maxis = range(xrange[0]+0.1,xrange[1]-0.1,50)

    djs_plot, [0], [0], /nodata, ysty=1, xsty=1, position=pos, $
      xrange=xrange, yrange=yrange, xtitle=xtitle1, ytitle=ytitle1, $
      xtickinterval=0.5
    legend, ['z=0.05-0.75'], /left, /top, box=0

    im_legend, [strupcase(calib),'Average'], /right, /bottom, box=0, margin=0, $
      charsize=1.5, pspacing=1.9, line=[calibline,0], thick=8, $
      psym=-[calibpsym,15], color=[calibcolor,'black'], symsize=[calibsymsize*1.0,1.5], $
      symthick=4.0

;   legend, textoidl('d[log(O/H)]/dz = '+string(mzavg.dlogohdz_coeff[0],format='(F6.3)')+'+'+$
;     string(mzavg.dlogohdz_coeff[1],format='(F5.3)')+' log (M/10^{10.5} M_{'+sunsymbol()+'})'), $
;     /right, /bottom, box=0, charsize=1.4

;; overplot the 1-sigma model
;    rand = mrandomn(seed,mzavg.dlogohdz_covar,1000)
;    rand[*,0] = rand[*,0] + mzavg.dlogohdz_coeff[0]
;    rand[*,1] = rand[*,1] + mzavg.dlogohdz_coeff[1]
;    
;    ell = covar2ellipse(mzavg.dlogohdz_covar,nsigma=1.0) ; get models within 1-sigma
;    indx = get_ellipse_indices(rand[*,0],rand[*,1],$
;      major=ell.major,minor=ell.minor,angle=ell.angle, $
;      xcenter=mzavg.dlogohdz_coeff[0],ycenter=mzavg.dlogohdz_coeff[1])
;    nindx = n_elements(indx)
;
;    allgrad = fltarr(n_elements(maxis),nindx)
;    for jj = 0, nindx-1 do allgrad[*,jj] = poly(maxis-mzavg.dlogohdz_normmass,rand[indx[jj],*])
;    mingrad = min(allgrad,dim=2)
;    maxgrad = max(allgrad,dim=2)
;    polyfill, [maxis,reverse(maxis)],[mingrad,reverse(maxgrad)], $
;      /data, /fill, color=im_color('tan'), noclip=0

;   maxmodel = poly(maxis-mzavg.dlogohdz_normmass,mzavg.dlogohdz_coeff+mzavg.dlogohdz_coeff_err/2)
;   minmodel = poly(maxis-mzavg.dlogohdz_normmass,mzavg.dlogohdz_coeff-mzavg.dlogohdz_coeff_err/2)
;   polyfill, [maxis,reverse(maxis)], [minmodel,reverse(maxmodel)], $
;     /data, /fill, color=im_color('tan',101), noclip=0

;   djs_oplot, maxis, poly(maxis-mzavg.dlogohdz_normmass,mzavg.dlogohdz_coeff), line=0, $
;     color=im_color('firebrick'), thick=7
    
; the list of calibrations here has to match fit_mzlzevol!!    
    for ii = 0, n_elements(calib)-1 do begin
       oploterror, mzavg.dlogohdz_medmass_all[*,ii]+off[ii], mzavg.coeffs_bymass_all[1,*,ii], $
         mzavg.coeffs_bymass_err_all[1,*,ii], color=im_color(calibcolor[ii]), $
         psym=symcat(calibpsym[ii]), symsize=calibsymsize[ii], thick=5, errthick=5
       djs_oplot, maxis, poly(maxis-mzavg.dlogohdz_normmass,mzavg.dlogohdz_coeff_all[*,ii]), $
         line=calibline[ii], color=im_color(calibcolor[ii]), thick=4
    endfor
    oploterror, mzavg.dlogohdz_medmass, mzavg.coeffs_bymass[1,*], $
      mzavg.coeffs_bymass_err[1,*], psym=symcat(15,thick=8), symsize=3, $
      errthick=6
    djs_oplot, maxis, poly(maxis-mzavg.dlogohdz_normmass,mzavg.dlogohdz_coeff), $
      line=0, thick=8
    print, mzavg.dlogohdz_coeff, mzavg.dlogohdz_coeff_err

    im_plotconfig, /psclose, psfile=psfile

; --------------------------------------------------
; Figure 14 - metallicity vs redshift in bins of mass - all three calibrations
    massbins = mz_massbins(/rev,masscut=masscut,masskeep=masskeep)
    nmassbins = n_elements(massbins)
    struct_print, massbins

    psfile = pspath+'z_vs_oh_bymass'+suffix
    im_plotconfig, 19, pos, psfile=psfile, height=2.2*[1,1,1,1], $
      charsize=1.6

    calib = ['kk04','t04','m91']
    calibcolor = ['firebrick','forest green','dodger blue']
    open_calibpsym = [9,5,6]
    fill_calibpsym = [16,17,15]

    calibsymsize = [1.5,1.8,1.5]*1.1
    calibline = [0,3,5]
    ncalib = n_elements(calib)

    xrange = [0.0,0.72]
    yrange = [8.5,9.27]
    xtitle1 = 'Redshift'
    ytitle1 = mzplot_ohtitle(/mean)

    zaxis = range(xrange[0]+0.03,xrange[1]-0.03,100)
    
    for mm = 0, n_elements(massbins)-1 do begin
       if (mm lt 6) then begin
          xtitle = ''
          xtickname = replicate(' ',10)
       endif else begin
          xtitle = xtitle1
          delvarx, xtickname
       endelse
       if odd(mm) then begin
          ytitle = ''
          ytickname = replicate(' ',10)
       endif else begin
          ytitle = ytitle1
          delvarx, ytickname
       endelse
       
       djs_plot, [0], [0], /nodata, noerase=(mm gt 0), ysty=1, xsty=1, xrange=xrange, $
         yrange=yrange, xtitle=xtitle, ytitle='', position=pos[*,mm], $
         xtickname=xtickname, ytickname=ytickname, yminor=2, xtickinterval=0.2
       im_legend, massbins[mm].label, /right, /top, box=0, margin=0, $
         charsize=1.3
       if (mm eq 0) then begin
          im_legend, strupcase(calib), /left, /bottom, box=0, margin=0, $
            charsize=1.2, pspacing=1.9, line=calibline, thick=8, $
            psym=-fill_calibpsym, color=calibcolor, symsize=calibsymsize*0.8, $
            symthick=4.0
       endif
;; overplot the model assuming no evolution in M*
;       for ii = 0, ncalib-1 do begin
;          mzevol = mrdfits(mzpath+'mzevol_'+calib[ii]+'.fits.gz',1,/silent)
;          minmodel = mzevol_func(zaxis*0+massbins[mm].massbin,$
;            mzevol.coeffs_r0zero-mzevol.coeffs_r0zero_err/2.0,z=zaxis,qz0=mzevol.qz0)
;          maxmodel = mzevol_func(zaxis*0+massbins[mm].massbin,$
;            mzevol.coeffs_r0zero+mzevol.coeffs_r0zero_err/2.0,z=zaxis,qz0=mzevol.qz0)
;          polyfill, [zaxis,reverse(zaxis)], [minmodel,reverse(maxmodel)], $
;            /data, /fill, color=im_color(polycolor[ii],101), noclip=0
;;         polyfill, [zaxis,reverse(zaxis)], [minmodel,reverse(maxmodel)], $
;;           /data, /line_fill, spacing=0.05, orientation=135, $
;;           color=im_color(polycolor[ii],101), noclip=0
;        endfor
; and now the data
       for ii = 0, ncalib-1 do begin
          mzevol = mrdfits(mzpath+'mzevol_'+calib[ii]+'.fits.gz',1,/silent)
;         good = where(mzevol.ohmean_bymass[masskeep[mm],*] gt -900.0)
;         djs_oplot, zbins[good].zbin, mzevol.ohmean[masskeep[mm],good], $
;           psym=symcat(calibpsym[ii],thick=5), symsize=calibsymsize[ii], $
;           color=im_color(calibcolor[ii],101)


          above = where(mzevol.ohmean_bymass[masskeep[mm],*] gt -900.0 and $
            mzevol.complete_bymass[masskeep[mm],*] eq 1,nabove)
          if (nabove ne 0) then begin
             oploterror, [mzevol.sdss_medz_bymass[masskeep[mm]],reform(mzevol.medz_bymass[masskeep[mm],above])], $
               [mzevol.sdss_ohmean_bymass[masskeep[mm]],reform(mzevol.ohmean_bymass[masskeep[mm],above])], $
               [mzevol.sdss_ohmean_bymass_err[masskeep[mm]],reform(mzevol.ohmean_bymass_err[masskeep[mm],above])], $
               psym=symcat(fill_calibpsym[ii],thick=5), symsize=calibsymsize[ii], $
               color=im_color(calibcolor[ii],101), errcolor=im_color(calibcolor[ii],101)
          endif

          below = where(mzevol.ohmean_bymass[masskeep[mm],*] gt -900.0 and $
            mzevol.complete_bymass[masskeep[mm],*] eq 0,nbelow)
          if (nbelow ne 0) then begin
             oploterror, mzevol.medz_bymass[masskeep[mm],below], mzevol.ohmean_bymass[masskeep[mm],below], $
               mzevol.ohmean_bymass_err[masskeep[mm],below], $
               psym=symcat(open_calibpsym[ii],thick=5), symsize=calibsymsize[ii], $
               color=im_color(calibcolor[ii],101), errcolor=im_color(calibcolor[ii],101)
          endif

;; AGES          
;          oploterror, mzevol.medz_bymass[masskeep[mm],good], mzevol.ohmean_bymass[masskeep[mm],good], $
;            mzevol.ohmean_bymass_err[masskeep[mm],good], $
;            psym=symcat(calibpsym[ii],thick=5), symsize=calibsymsize[ii], $
;            color=im_color(calibcolor[ii],101), errcolor=im_color(calibcolor[ii],101)
;; SDSS
;          oploterror, mzevol.sdss_medz_bymass[masskeep[mm]], mzevol.sdss_ohmean_bymass[masskeep[mm]], $
;            mzevol.sdss_ohmean_bymass_err[masskeep[mm]], $
;            psym=symcat(sdss_calibpsym[ii],thick=5), symsize=calibsymsize[ii], $
;            color=im_color(calibcolor[ii],101), errcolor=im_color(calibcolor[ii],101)          
          
          djs_oplot, zaxis, poly(zaxis-mzevol.qz0,mzevol.coeffs_bymass[*,masskeep[mm]]), $
            line=calibline[ii], color=im_color(calibcolor[ii],101)
;         djs_oplot, zaxis, mzevol_func(zaxis*0+massbins[masskeep[mm]].massbin,mzevol.coeffs,$
;           z=zaxis,qz0=mzevol.qz0), line=calibline[ii], color=im_color(calibcolor[ii],101)
       endfor
       xyouts, pos[0,0]-0.075, pos[1,0], ytitle1, align=0.5, orientation=90, /norm
       xyouts, pos[0,4]-0.075, pos[1,4], ytitle1, align=0.5, orientation=90, /norm
    endfor

    im_plotconfig, /psclose, psfile=psfile

stop    
    
;; --------------------------------------------------
;; QB (luminosity evolution rate) vs redshift for various bins of
;; stellar mass
;    mzlocal = mrdfits(mzpath+'mzlocal_sdss_ews_t04.fits.gz',1)
;    psfile = pspath+'mzevol_model'+suffix
;    im_plotconfig, 0, pos, psfile=psfile, height=4.7, width=6.5, $
;      xmargin=[1.6,0.4]
;
;    xrange = [9.4,11.5]
;    yrange = [8.55,9.15]
;    xtitle1 = mzplot_masstitle()
;    ytitle1 = mzplot_ohtitle()+'_{T04}'
;    maxis = range(9.5,xrange[1]-0.05,100)
;
;    zval = [0.1,0.3,0.5,0.7,0.9]
;    zlabel = 'z='+string(zval,format='(F3.1)')
;    zline = [0,3,5,4,2]*0
;    zline2 = zval*0+5
;    zcolor = ['black','firebrick','dodger blue','forest green','orange']
;    
;    djs_plot, [0], [0], /nodata, ysty=1, xsty=1, position=pos, $
;      xrange=xrange, yrange=yrange, xtitle=xtitle1, ytitle=ytitle1
;    im_legend, zlabel, /right, /bottom, box=0, line=zline, $
;      color=zcolor, pspacing=1.8, thick=6, margin=0, $
;      charsize=1.6
;    for ii = 0, n_elements(zval)-1 do begin
;       ohmodel = mz_closedbox(maxis,mzlocal.coeff) + (zval[ii]-mzavg.qz0)*$
;         poly(maxis-mzavg.dlogohdz_normmass,mzavg.dlogohdz_coeff)
;       keep = where((ohmodel gt yrange[0]+0.03))
;       djs_oplot, maxis[keep], ohmodel[keep], line=zline[ii], $
;         color=im_color(zcolor[ii],101), thick=7
;; P,R model
;       djs_oplot, maxis, mzevol_func(maxis,mzavg.mlfit_coeffs,z=zval[ii],qz0=mzavg.qz0), $
;         line=zline2[ii], color=im_color(zcolor[ii],101), thick=7
;    endfor
;    
;    im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)

stop    
    
stop
stop
stop
    
; --------------------------------------------------
; evolution of (O/H)* and c0 with redshift
    ztitle1 = 'Redshift'
    ohtitle1 = textoidl('12+log(O/H) at 10^{'+string(mz_pivotmass(),format='(F4.1)')+'} M_{'+sunsymbol()+'}')
;   ohtitle1 = textoidl('12+log(O/H)^{*}')
    ohtitle2 = textoidl('12+log(O/H) at M_{B}='+string(lz_pivotmag(),format='(F5.1)'))
    zrange = [-0.02,0.8]
    ohrange = [8.49,9.2]

    calibcolor = ['orange','firebrick','forest green']
    calibline = [0,5,3]
    calibpsym = [15,17,16]
    calibsymsize = [2.5,2.8,2.5]

    sdsspsym = 6
    sdsssymsize = 3.0
    sdsscolor = 'dark orchid'
    
    psfile = pspath+'z_vs_ohstar'+suffix
    im_plotconfig, 6, pos, psfile=psfile, charsize=1.8, yspace=0.9, $
      height=4.0*[1,1], xmargin=[1.4,0.4], width=6.7

; MZ evolution with r0=0.0 dex/z
    djs_plot, [0], [0], /nodata, ysty=1, xsty=1, xrange=zrange, $
      yrange=ohrange, xtitle=ztitle1, ytitle=ohtitle1, position=pos[*,0]
    im_legend, strupcase(mzlocal.calib), /left, /bottom, box=0, $
      charsize=1.6, pspacing=1.9, line=calibline, thick=8, $
      psym=-calibpsym, color=calibcolor, symsize=calibsymsize*0.7
    
    r0indx = 0
    for ii = 0, ncalib-1 do begin
       djs_oplot, zaxis1, mzevol_func(replicate(mz_pivotmass(),nzaxis1),$
         mzevol[ii].coeffs[*,r0indx],z=zaxis1,qz0=qz0), $
         line=calibline[ii], color=im_color(calibcolor[ii],101), thick=8
       djs_oplot, zbins.zbin, mzevol[ii].ohstar[*,r0indx], psym=symcat(calibpsym[ii],thick=8), $
         symsize=calibsymsize[ii], color=im_color(calibcolor[ii],101)
;      djs_oplot, zbins.zbin, mzevol[ii].ohstar[*,2], psym=symcat(calibpsym[ii],thick=8), $
;        symsize=calibsymsize[ii], color=im_color(calibcolor[ii],101)
;      plots, 0.07, mz_brokenpl(mz_pivotmass(),mzlocal[ii].coeff_bin), $
;        psym=symcat(sdsspsym,thick=10), symsize=sdsssymsize, color=im_color(sdsscolor,101)
    endfor

; LZ evolution with q0=1.6 dex/z
    djs_plot, [0], [0], /nodata, /noerase, ysty=1, xsty=1, xrange=zrange, $
      yrange=ohrange, xtitle=ztitle1, ytitle=ohtitle2, position=pos[*,1]

    q0indx = 2
    for ii = 0, ncalib-1 do begin
       djs_oplot, zaxis1, lzevol_func(replicate(lz_pivotmag(),nzaxis1),$
         lzevol[ii].coeffs[*,q0indx],z=zaxis1,qz0=qz0,pivotmag=lz_pivotmag()), $
         line=calibline[ii], color=im_color(calibcolor[ii],101), thick=8
       djs_oplot, zbins.zbin, lzevol[ii].ohstar[*,q0indx], psym=symcat(calibpsym[ii],thick=8), $
         symsize=calibsymsize[ii], color=im_color(calibcolor[ii],101)
;      plots, 0.07, poly(0.0,lzlocal[ii].coeff), $
;        psym=symcat(sdsspsym,thick=10), symsize=sdsssymsize, color=im_color(sdsscolor,101)
    endfor

    im_plotconfig, /psclose

stop
stop
stop    

; --------------------------------------------------
; evolution of Delta-log (O/H)* with redshift
    ztitle1 = 'Redshift'
    dohtitle1 = textoidl('<\Delta'+'log(O/H)> from !8M-Z!6')
    dohtitle2 = textoidl('<\Delta'+'log(O/H)> from !8L-Z!6')

    zrange = [-0.02,0.8]
    dohrange = [-0.49,0.09]

    r0 = [0.0,0.5]
    q0 = [0.0,2.0]
;   r0 = [0.0,0.25,0.5]
;   q0 = [0.0,1.2,1.6,2.0]
    
    r0color = ['orange','firebrick','tan']
    r0line = [0,5] & r0psym = [15,17] & r0symsize = [2.5,2.8]
    q0color = ['forest green','dodger blue']
    q0line = [0,5] & q0psym = [15,17] & q0symsize = [2.5,2.8]
    
    psfile = pspath+'z_vs_dlogoh'+suffix
    im_plotconfig, 6, pos, psfile=psfile, charsize=1.8, yspace=0.9, $
      height=4.0*[1,1], xmargin=[1.4,0.4], width=6.7

; dlogoh from MZ relation for various values of r0
    djs_plot, [0], [0], /nodata, ysty=1, xsty=1, xrange=zrange, $
      yrange=dohrange, xtitle=ztitle1, ytitle=dohtitle1, position=pos[*,0]
;   im_legend, ['R='+string(r0,format='(F4.1)')], /left, $
    im_legend, ['R='+string(r0,format='(F3.1)')+' dex z^{-1}'], /left, $
      /bottom, box=0, charsize=1.6, pspacing=1.9, line=r0line, thick=8, $
      psym=-r0psym, color=r0color, symsize=r0symsize*0.7
    
    for ii = 0, n_elements(r0)-1 do begin
       get_element, mzevol[0].coeffs[5,*], r0[ii], r0indx
       oploterror, zbins.zbin, mzevol[0].dlogoh_avg[*,r0indx], mzevol[0].dlogoh_avg_err[*,r0indx], $
         psym=symcat(r0psym[ii],thick=8), symsize=r0symsize[ii], errthick=8, $
         color=im_color(r0color[ii],101), errcolor=im_color(r0color[ii],101)
       djs_oplot, zaxis1, poly(zaxis1-qz0,[0.0,mzevol[0].pavg[r0indx]]), line=r0line[ii], $
         color=im_color(r0color[ii],101), thick=8
    endfor
       
; dlogoh from LZ relation for various values of q0
    djs_plot, [0], [0], /nodata, /noerase, ysty=1, xsty=1, xrange=zrange, $
      yrange=dohrange, xtitle=ztitle1, ytitle=dohtitle2, position=pos[*,1]
;   im_legend, ['Q='+string(q0,format='(F3.1)')], /left, $
    im_legend, ['Q_{B}='+string(q0,format='(F3.1)')+' mag z^{-1}'], /left, $
      /bottom, box=0, charsize=1.6, pspacing=1.9, line=q0line, thick=8, $
      psym=-q0psym, color=q0color, symsize=q0symsize*0.7

    for ii = 0, n_elements(q0)-1 do begin
       get_element, lzevol[0].coeffs[3,*], q0[ii], q0indx
       oploterror, zbins.zbin, lzevol[0].dlogoh_avg[*,q0indx], lzevol[0].dlogoh_avg_err[*,q0indx], $
         psym=symcat(q0psym[ii],thick=8), symsize=q0symsize[ii], errthick=8, $
         color=im_color(q0color[ii],101), errcolor=im_color(q0color[ii],101)
       djs_oplot, zaxis1, poly(zaxis1-qz0,[0.0,lzevol[0].savg[q0indx]]), line=q0line[ii], $
         color=im_color(q0color[ii],101), thick=8
    endfor

    im_plotconfig, /psclose

stop    
    
    timelabel1 = [1.0,3.0,5.0] ; [Gyr]
    timelabel2 = [1.0,3.0,5.0,7.0,10.0,12.0] ; [Gyr]
    nzaxis1 = 100
    zaxis1 = range(-0.05,0.8,nzaxis1)
    zz = im_array(0.0,1.0,0.02)

; representative Pegase models
    splog, 'Reading the Pegase models'
    pegpath = getenv('PEGASE_HR_SFHGRID_DIR')+'/measure/'
    peg1 = mrdfits(pegpath+'salp_tau_001.0Gyr.info.fits',1,silent=0)
    peg2 = mrdfits(pegpath+'salp_tau_001.0Gyr.info.fits',1,silent=0)
    peg3 = mrdfits(pegpath+'salp_tau_003.0Gyr.info.fits',1,silent=0)
    ipeg3 = mrdfits(pegpath+'salp_kennlaw_0.10_infall_003.0Gyr.info.fits',1,silent=0)
    pegconst = mrdfits(pegpath+'salp_tau_999.0Gyr.info.fits',1,silent=0)

    zf = 5.0
    zf2_peg1 = zformpeg(peg1,zform=zf,peg_zaxis=zf2_peg_zaxis)
    zf2_peg2 = zformpeg(peg2,zform=zf)
    zf2_peg3 = zformpeg(peg3,zform=zf)
    zf2_ipeg3 = zformpeg(ipeg3,zform=zf)
    zf2_pegconst = zformpeg(pegconst,zform=zf)

    zf = 1.5
    zf1_peg1 = zformpeg(peg1,zform=zf,peg_zaxis=zf1_peg_zaxis)
    zf1_peg2 = zformpeg(peg2,zform=zf)
    zf1_peg3 = zformpeg(peg3,zform=zf)
    zf1_ipeg3 = zformpeg(ipeg3,zform=zf)
    zf1_pegconst = zformpeg(pegconst,zform=zf)
    
; read the output from FIT_MZLZEVOL
    ohevol = mrdfits(mzpath+'ohevol.fits.gz',1)

    glzevol_ages = mrdfits(mzpath+'lzevol_g.fits.gz',1)
    glzevol_noevol = mrdfits(mzpath+'lzevol_g.fits.gz',2)
    glzevol_levol = mrdfits(mzpath+'lzevol_g.fits.gz',3)

    mzevol_ages = mrdfits(mzpath+'mzevol.fits.gz',1)
    mzevol_noevol = mrdfits(mzpath+'mzevol.fits.gz',2)
    mzevol_levol = mrdfits(mzpath+'mzevol.fits.gz',3)

; ---------------------------------------------------------------------------
; redshift versus SFR density

    psfile = pspath+'z_vs_sfrd'+suffix
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, height=5.5, $
      ymargin=[0.9,1.1], xmargin=[1.5,0.2], thick=6
    
    h = rsex(litpath+'04hopkins.sex')
    uv = where(strmatch(h.indicator,'*UV*',/fold))
    ha = where(strmatch(h.indicator,'*Ha*',/fold) or $
      strmatch(h.indicator,'*Hb*',/fold) or $
      strmatch(h.indicator,'*OII*',/fold))
    ir = where(strmatch(h.indicator,'*IR*',/fold))
    rad = where(strmatch(h.indicator,'*RADIO*',/fold) or $
      strmatch(h.indicator,'*xray*',/fold))

    uvcolor = 'royal blue'
    hacolor = 'forest green'
    ircolor = 'firebrick'
    radcolor = 'orange'

    uvpsym = 16
    hapsym = 14
    irpsym = 15
    radpsym = 17

    uvsymsize = 1.5
    hasymsize = 2.0
    irsymsize = 1.5
    radsymsize = 1.5
        
    xrange1 = [-0.05,0.92]
    yrange1 = [-2.3,-0.2]
    xtitle1 = 'log (1 + z)'
    ytitle1 = 'log (d\rho_{*}/dt) (M'+sunsymbol()+' yr^{-1} Mpc^{-3} h_{70}^{3})'

    xrange2 = [0.0,0.9]
    yrange2 = [0.7,1.04]
    xtitle2 = 'Redshift'
    ytitle2 = 'F(<z)/F(z=0.1)'
    
    yrange3 = alog10(yrange2)
    ytitle3 = textoidl('\Delta'+'log(O/H)')
;   ytitle3 = textoidl('<\Delta'+'log(O/H)>')

    djs_plot, [0], [0], /nodata, xtitle=xtitle1, ytitle=ytitle1, $
      xsty=9, ysty=1, xrange=xrange1, yrange=yrange1, position=pos[*,0]
    axis, /xaxis, xthick=postthick1, xsty=1, xrange=xrange, xtitle='Lookback Time (Gyr)', $
      xtickv=alog10(1.0+getredshift(getage(0.0)-timelabel2)), $
      xticks=n_elements(timelabel2)-1L, xtickname=string(timelabel2,format='(I0)')
; UV    
    oploterror, alog10(1+h[uv].z), h[uv].sfrd, h[uv].zerr_lo/(1+h[uv].z)/alog(10.0), $
      h[uv].sfrderr_lo, psym=symcat(uvpsym), color=im_color(uvcolor,10), $
      errcolor=im_color(uvcolor,10), errthick=!p.thick, /lobar, symsize=uvsymsize
    oploterror, alog10(1+h[uv].z), h[uv].sfrd, h[uv].zerr_hi/(1+h[uv].z)/alog(10.0), $
      h[uv].sfrderr_hi, psym=symcat(uvpsym), color=im_color(uvcolor,10), $
      errcolor=im_color(uvcolor,10), errthick=!p.thick, /hibar, symsize=uvsymsize
; Ha
    oploterror, alog10(1+h[ha].z), h[ha].sfrd, h[ha].zerr_lo/(1+h[ha].z)/alog(10.0), $
      h[ha].sfrderr_lo, psym=symcat(hapsym), color=im_color(hacolor,10), $
      errcolor=im_color(hacolor,10), errthick=!p.thick, /lobar, symsize=hasymsize
    oploterror, alog10(1+h[ha].z), h[ha].sfrd, h[ha].zerr_hi/(1+h[ha].z)/alog(10.0), $
      h[ha].sfrderr_hi, psym=symcat(hapsym), color=im_color(hacolor,10), $
      errcolor=im_color(hacolor,10), errthick=!p.thick, /hibar, symsize=hasymsize
; IR
    oploterror, alog10(1+h[ir].z), h[ir].sfrd, h[ir].zerr_lo/(1+h[ir].z)/alog(10.0), $
      h[ir].sfrderr_lo, psym=symcat(irpsym), color=im_color(ircolor,10), $
      errcolor=im_color(ircolor,10), errthick=!p.thick, /lobar, symsize=irsymsize
    oploterror, alog10(1+h[ir].z), h[ir].sfrd, h[ir].zerr_hi/(1+h[ir].z)/alog(10.0), $
      h[ir].sfrderr_hi, psym=symcat(irpsym), color=im_color(ircolor,10), $
      errcolor=im_color(ircolor,10), errthick=!p.thick, /hibar, symsize=irsymsize
; Radio
    oploterror, alog10(1+h[rad].z), h[rad].sfrd, h[rad].zerr_lo/(1+h[rad].z)/alog(10.0), $
      h[rad].sfrderr_lo, psym=symcat(radpsym), color=im_color(radcolor,10), $
      errcolor=im_color(radcolor,10), errthick=!p.thick, /lobar, symsize=radsymsize
    oploterror, alog10(1+h[rad].z), h[rad].sfrd, h[rad].zerr_hi/(1+h[rad].z)/alog(10.0), $
      h[rad].sfrderr_hi, psym=symcat(radpsym), color=im_color(radcolor,10), $
      errcolor=im_color(radcolor,10), errthick=!p.thick, /hibar, symsize=radsymsize

; legend
    label = ['UV','H\alpha/H\beta/[OII]','IR','Radio/Xray']
    im_legend, label, /left, /top, box=0, psym=[uvpsym,hapsym,irpsym,radpsym], $
      symsize=[uvsymsize,hasymsize,irsymsize,radsymsize], $
      color=[uvcolor,hacolor,ircolor,radcolor], $
      spacing=1.7, charsize=1.5
    
; overplot the Pegase models here
    
; coefficients from plot_redshift_vs_sfr_mass_density    
    sfrd_zaxis = im_array(0.1,6.0,0.01)
    sfrd_lookback = getage(0.0)-getage(sfrd_zaxis)  ; lookback time
    plot_sfrd_zaxis = [0.0,sfrd_zaxis]

;   sfrd = 10.0^h.sfrd
;   logsfrderr = sqrt((total([[h.sfrderr_lo],[h.sfrderr_hi]],2)/2.0)^2+0.1^2)
;   sfrderr = logsfrderr*sfrd*alog(10.0)*0.0+0.1 ; no weighting!
;   coeff = mpfitexpr('(P[0]+P[1]*X)/(1.0+(X/P[2])^P[3])',h.z,$
;     sfrd,sfrderr,fltarr(4)+1.0,/quiet)
;   print, coeff

;   coeff = [0.0170,0.13,3.3,5.3]
    coeff = [0.010464392D,0.11303816D,3.4534319D,3.7109016D]
    sfrd = (coeff[0]+coeff[1]*sfrd_zaxis)/(1.0+(sfrd_zaxis/coeff[2])^coeff[3])
    sfrd_up = sfrd*10^0.1
    sfrd_lo = sfrd/10^0.1
    
    plot_sfrd = (coeff[0]+coeff[1]*plot_sfrd_zaxis)/(1.0+(plot_sfrd_zaxis/coeff[2])^coeff[3])
    djs_oplot, alog10(1+plot_sfrd_zaxis), alog10(plot_sfrd), line=0
    djs_oplot, alog10(1+plot_sfrd_zaxis), alog10(plot_sfrd)+0.1, line=5
    djs_oplot, alog10(1+plot_sfrd_zaxis), alog10(plot_sfrd)-0.1, line=5

; now compute and overplot the fraction of stars formed versus
; redshift
    sfrd_total = im_integral(sfrd_lookback,sfrd*1D9,$              ; M_sun/Mpc^3
      min(sfrd_lookback),max(sfrd_lookback))                       
    sfrd_time = im_integral(sfrd_lookback,sfrd*1D9,sfrd_lookback,$ ; M_sun/Mpc^3
      replicate(max(sfrd_lookback),n_elements(sfrd_lookback)))
    sfrd_frac = sfrd_time/sfrd_total

; upper    
    sfrd_total_up = im_integral(sfrd_lookback,sfrd_up*1D9,$              ; M_sun/Mpc^3
      min(sfrd_lookback),max(sfrd_lookback))                       
    sfrd_time_up = im_integral(sfrd_lookback,sfrd_up*1D9,sfrd_lookback,$ ; M_sun/Mpc^3
      replicate(max(sfrd_lookback),n_elements(sfrd_lookback)))
    sfrd_frac_up = sfrd_time_up/sfrd_total_up
; lower
    sfrd_total_lo = im_integral(sfrd_lookback,sfrd_lo*1D9,$              ; M_sun/Mpc^3
      min(sfrd_lookback),max(sfrd_lookback))                       
    sfrd_time_lo = im_integral(sfrd_lookback,sfrd_lo*1D9,sfrd_lookback,$ ; M_sun/Mpc^3
      replicate(max(sfrd_lookback),n_elements(sfrd_lookback)))
    sfrd_frac_lo = sfrd_time_lo/sfrd_total_lo

; extrapolate to z<0.1    
    notzero = where(sfrd_frac gt 0.0)
    plot_sfrd_frac1 = interpol((sfrd_frac[notzero]),$
      sfrd_zaxis[notzero],plot_sfrd_zaxis,/quad)

; inset figure    
    djs_plot, [0], [0], /noerase, xtitle=xtitle2, ytitle=ytitle2, $
      xsty=1, ysty=9, xrange=xrange2, yrange=yrange2, charsize=1.3, $
      pos=[0.52,0.24,0.85,0.48], ytickinterval=0.1
    axis, /yaxis, ysty=1, yrange=yrange3, ytitle=ytitle3, charsize=1.3
    djs_oplot, zz, zz*0.0+1.0, line=1, thick=3
    djs_oplot, plot_sfrd_zaxis, plot_sfrd_frac1, line=0
    
    im_plotconfig, /psclose

; ------------------------------------------------------------
; luminosity and mass evolution for blue galaxies (literature
; compilation)

    h100 = mz_h100()
    Bvega2ab = (k_vega2ab(filterlist='bessell_B.par',/kurucz,/silent))[0]
    B2r01 = -0.6429 ; ^{0.1}r = B-0.6429-1.0845*[(B-V)-0.5870] [AB, Blanton & Roweis]
    B2g01 = +0.0759 ; ^{0.1}g = B+0.0759+0.0620*[(B-V)-0.5870] [AB, Blanton & Roweis]

    psfile = pspath+'z_vs_mr_mass_lit'+suffix
    im_plotconfig, 6, pos, psfile=psfile, charsize=1.8, $
      xmargin=[1.4,0.2], ymargin=[1.0,1.1]

    xtitle = 'Redshift'
    ytitle1 = textoidl('M_{0.1r}^{*}')
    ytitle2 = textoidl('log (M^{*}/M'+sunsymbol()+')')

    xrange = [-0.02,0.82]
    yrange1 = [-20.2,-22.4]
    yrange2 = [10.6,11.29]

    w06color   = 'orange'    & w06psym = 16 & w06symsize = 2.4
    f06color   = 'navy'      & f06psym = 14 & f06symsize = 2.8
    e07color   = 'firebrick' & e07psym = 15 & e07symsize = 2.0
    i09color   = 'firebrick' & i09psym = 15 & i09symsize = 2.0
    b06color   = 'navy'      & b06psym = 14 & b06symsize = 2.8

; luminosity evolution plot    
    djs_plot, [0], [0], /nodata, ysty=1, xsty=9, xrange=xrange, $
      yrange=yrange1, xtitle='', ytitle=ytitle1, position=pos[*,0], $
      xtickname=replicate(' ',10)
    axis, /xaxis, xsty=1, xrange=xrange, xtickv=getredshift(getage(0.0)-timelabel1), $
      xtitle='Lookback Time (Gyr)', xticks=n_elements(timelabel1)-1L, $
      xtickname=string(timelabel1,format='(I0)')

; overplot the line for 1.5 mag/z of luminosity evolution    
    djs_oplot, zaxis1, poly(zaxis1-0.1,[-20.95,-1.5]), line=0

; -------------------------
; Eisenstein et al. 2009 [AGES]; Omega_0=0.3, Omega_lamba=0.7, h=1.0;
; AB; alpha fixed at -1.10 for all blue galaxies; evolving color cut
; based on the (u0.1-r0.1) color: A=u0.1-r0.1+0.08(M_r0.1+20);  

    z_e07    = [0.1,0.2,0.3,0.4,0.5,0.65]
    zerr_e07 = [0.05,0.05,0.05,0.05,0.05,0.1]

    mstar_e07 = [-19.95,-20.39,-20.43,-20.57,-20.86,-20.83] + 5.0*alog10(h100) ; h=1-->h=0.7
    mstarerr_e07 = [0.12,0.07,0.04,0.05,0.11,0.08]

    oploterror, z_e07, mstar_e07, zerr_e07, mstarerr_e07, symsize=e07symsize, $
      psym=symcat(e07psym), color=im_color(e07color,1E8), errcolor=im_color(e07color,1E8)

; -------------------------
; Willmer et al. 2006 [DEEP2]; "minimal" model adopted, as
; recommended; Omega_0=0.3, Omega_lamba=0.7, h=0.7; Vega; alpha=-1.30
; fixed for all blue galaxies, based on a fit to all galaxies at
; 0.2<z<0.6 in COMBO-17 (Faber et al. 2005); color cut: (U-B) =
; -0.032(M_B+21.52)+0.454-0.25; drop the last bin

    z_w06 = [0.3,0.5,0.7,0.9,1.1];,1.3]
    zerr_w06 = [0.1,0.1,0.1,0.1,0.1];,0.1]

    mstar_w06 = [-20.36,-20.72,-21.15,-21.21,-21.38] + Bvega2ab + B2r01 ; Vega-->AB; B-->r0.1 ; ,-21.86]
    mstarerr_up_w06 = [0.13,0.05,0.07,0.01,0.04];,0.07]
    mstarerr_lo_w06 = [0.11,0.07,0.07,0.03,0.05];,0.08]
    mstarerr_w06 = mstar_w06*0.0
    for i = 0, n_elements(mstar_w06)-1L do mstarerr_w06[i] = mean([mstarerr_up_w06[i],mstarerr_lo_w06[i]])
    
    oploterror, z_w06, mstar_w06, zerr_w06, mstarerr_up_w06, /hi, symsize=w06symsize, $
      psym=symcat(w06psym), color=im_color(w06color,1E8), errcolor=im_color(w06color,1E8)
    oploterror, z_w06, mstar_w06, zerr_w06, mstarerr_lo_w06, /lo, symsize=w06symsize, $
      psym=symcat(w06psym), color=im_color(w06color,1E8), errcolor=im_color(w06color,1E8)

; -------------------------
; Faber et al. 2007 [COMBO17]; Omega_0=0.3, Omega_lamba=0.7, h=0.7;
; Vega; alpha=-1.3 fixed for all blue galaxies;

    z_f06 = [0.3,0.5,0.7,0.9,1.1]
    zerr_f06 = [0.1,0.1,0.1,0.1,0.1]

    mstar_f06 = [-20.74,-21.10,-21.30,-21.10,-21.25] + Bvega2ab + B2r01 ; Vega-->AB; B-->r0.1
    mstarerr_f06 = [0.20,0.15,0.16,0.17,0.18]

    oploterror, z_f06, mstar_f06, zerr_f06, mstarerr_f06, symsize=f06symsize, $
      psym=symcat(f06psym), color=im_color(f06color,1E8), errcolor=im_color(f06color,1E8)
    
; legend
    label = ['Willmer+06','Faber+07','Eisenstein+09']
    im_legend, label, /left, /top, box=0, psym=[w06psym,f06psym,e07psym], $
      symsize=[w06symsize,f06symsize,e07symsize], color=[w06color,f06color,e07color], $
      spacing=2.0, charsize=1.4

; mass evolution plot    
    djs_plot, [0], [0], /nodata, /noerase, ysty=1, xsty=1, xrange=xrange, $
      yrange=yrange2, xtitle=xtitle, ytitle=ytitle2, position=pos[*,1], $
      yminor=3

; overplot the line for 1.5 mag/z of luminosity evolution    
    djs_oplot, zaxis1, poly(zaxis1-0.1,[10.87,0.25]), line=0

; -------------------------
; Ilbert+09 [COSMOS]; Omega_0=0.3, Omega_lamba=0.7, h=0.7; Chabrier
; IMF; see ILBERT_MASS_FUNCTION for details

    z_i09 = [0.3,0.5,0.7,0.9,1.1,1.35,1.75]
    zerr_i09 = [0.1,0.1,0.1,0.1,0.1,0.15,0.25]

; intermediate-activity star-forming galaxies    
    mstar_i09 = [11.01,10.95,10.89,10.78,10.73,10.76,10.86]
    mstarerr_i09 = [0.03,0.03,0.02,0.02,0.02,0.02,0.02]
; high-activity star-forming galaxies    
;   mstar_i09 = [10.26,10.30,10.39,10.42,10.42,10.49,10.86]
;   mstarerr_i09 = [0.05,0.03,0.02,0.02,0.02,0.01,0.5]
; all star-forming galaxies    
    mstar_i09 = [11.01,10.96,10.92,10.83,10.80,10.81,10.94]
    mstarerr_i09 = mstar_i09*0.0+0.1

    oploterror, z_i09, mstar_i09, zerr_i09, mstarerr_i09, symsize=i09symsize, $
      psym=symcat(i09psym), color=im_color(i09color,1E8), errcolor=im_color(i09color,1E8)

; -------------------------
; Borch+06 [COMBO-17]; Omega_0=0.3, Omega_lamba=0.7, h=0.7; Kroupa+93
; IMF (same as Chabrier) 
    
    z_b06 = [0.1,0.3,0.5,0.7,0.9]
    zerr_b06 = [0.1,0.1,0.1,0.1,0.1]

    mstar_b06 = [10.80,10.86,10.93,11.07,11.00]
    mstarerr_b06 = [0.10,0.10,0.19,0.12,0.07,0.25]

    oploterror, z_b06, mstar_b06, zerr_b06, mstarerr_b06, symsize=b06symsize, $
      psym=symcat(b06psym), color=im_color(b06color,1E8), errcolor=im_color(b06color,1E8)

; fit to both samples
    xx = [z_b06,z_i09] & yy = [mstar_b06,mstar_i09] & yyerr = [mstarerr_b06,mstarerr_i09]
    ww = where(xx lt 0.8)
;   print, linfit(xx[ww],yy[ww],measure_errors=yyerr[ww])
    
; legend
    label = ['Borch+06','Ilbert+09']
    im_legend, label, /left, /top, box=0, psym=[b06psym,i09psym], $
      symsize=[b06symsize,i09symsize], color=[b06color,i09color], $
      spacing=2.0, charsize=1.4

    im_plotconfig, /psclose

; ------------------------------------------------------------
; redshift vs fraction of stars/metals formed
    
    psfile = pspath+'z_vs_ohfraction'+suffix
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, width=6.5, $
      height=6.5, xmargin=[1.6,0.4], ymargin=[1.0,1.0], thick=6.0

    xrange = [0.0,0.77]
    yrange = [0.6,1.049]
;   yrange = [0.74,1.0499]
;   yrange = [73.0,104.0] ; per cent
;   yrange = [60.0,109.0] ; per cent
    xtitle = 'Redshift'
    ytitle = 'Fraction of Metals Synthesized by z=0.1'
    
; coefficients from plot_redshift_vs_sfr_mass_density    
    sfrd_zaxis = findgen((6.0-0.1)/0.01+1)*0.01+0.1 ; for integral
    sfrd_lookback = getage(0.0)-getage(sfrd_zaxis)  ; lookback time
    
    coeff = [0.010464392D,0.11303816D,3.4534319D,3.7109016D]
    sfrd = (coeff[0]+coeff[1]*sfrd_zaxis)/(1.0+(sfrd_zaxis/coeff[2])^coeff[3])
;   plot, alog10(1+sfrd_zaxis), alog10(sfrd), xsty=3, ysty=3, yr=[-2.3,-0.3]
    
; what fraction of the stars (and therefore metals) were formed by
; various redshifts; integrate the model; convert the SFRD to
; M_sun/Gyr; the factor of 1/64 (1.6%) is the IMF-weighted yield based
; on the Woosley & Weaver (1995) models, assuming a Salpeter IMF
; 0.1-125 M_sun (Conti et al. 2003; Pettini et al. 2006); note that
; Madau et al. (1996) used 1/42 (2.4%) for a Salpeter 0.1-100 IMF;
; also note that the exact value is ~independent of the IMF; and in
; detail, of course, the factor drops out because we are doing
; fractions here
    sfrd_frac = sfrd_frac(sfrd,sfrd_zaxis,sfrd_zaxis,plot_zaxis)

; interpolate the Pegase models and then integrate; note that I'm not
; dividing by the co-moving volume, but I could
    zf2_peg1_sfrd_frac = sfrd_frac(zf2_peg1.sfr,zf2_peg_zaxis,sfrd_zaxis)
    zf1_peg1_sfrd_frac = sfrd_frac(zf1_peg1.sfr,zf1_peg_zaxis,sfrd_zaxis)

    zf2_peg3_sfrd_frac = sfrd_frac(zf2_peg3.sfr,zf2_peg_zaxis,sfrd_zaxis)
    zf1_peg3_sfrd_frac = sfrd_frac(zf1_peg3.sfr,zf1_peg_zaxis,sfrd_zaxis)

    zf2_pegconst_sfrd_frac = sfrd_frac(zf2_pegconst.sfr,zf2_peg_zaxis,sfrd_zaxis)
    zf1_pegconst_sfrd_frac = sfrd_frac(zf1_pegconst.sfr,zf1_peg_zaxis,sfrd_zaxis)

; now make the plot    
    djs_plot, [0], [0], /nodata, ysty=1, xsty=9, xrange=xrange, $
      yrange=yrange, xtitle=xtitle, ytitle=ytitle, position=pos
    axis, /xaxis, xsty=1, xrange=xrange, xtickv=getredshift(getage(0.0)-timelabel1), $
      xtitle='Lookback Time (Gyr)', xticks=n_elements(timelabel1)-1L, $
      xtickname=string(timelabel1,format='(I0)')
    djs_oplot, plot_zaxis, sfrd_frac, line=0, thick=12

; overplot the Pegase models
; tau=1
;   djs_oplot, plot_zaxis, zf2_peg1_sfrd_frac, line=5, $
;     color=im_color('dodger blue',148), thick=8
    djs_oplot, plot_zaxis, zf1_peg1_sfrd_frac, line=5, $
      color=im_color('dodger blue',148), thick=8

; tau=3    
;   djs_oplot, plot_zaxis, zf2_peg3_sfrd_frac, line=3, $
;     color=im_color('orange',149), thick=8
    djs_oplot, plot_zaxis, zf1_peg3_sfrd_frac, line=3, $
      color=im_color('orange',149), thick=8

; tau=inf    
;   djs_oplot, plot_zaxis, zf2_pegconst_sfrd_frac, line=1, $
;     color=im_color('purple',150), thick=10
    djs_oplot, plot_zaxis, zf1_pegconst_sfrd_frac, line=1, $
      color=im_color('purple',150), thick=10

; and finally the data    
    color1 = 'firebrick' & psym1 = 6 & symsize1 = 3.0 & line1 = 3
    color2 = 'royal blue' & psym2 = 9 & symsize2 = 3.0 & line2 = 5
    color3 = 'forest green' & psym3 = 4 & symsize3 = 3.5 & line3 = 0
    color4 = 'orange'       & psym3 = 7 & symsize3 = 3.5 & line3 = 4

    zerr = [0.05,0.05,0.05,0.05,0.05,0.1] ; no error bar on the first point
    z1 = ohevol.z
    z2 = ohevol.z-0.015
    z3 = ohevol.z+0.015
    z4 = ohevol.z+0.015

    oh1 = 1.0*10^ohevol.ldlogoh_noevol_cor
    oh2 = 1.0*10^ohevol.ldlogoh_levol_cor
    oh3 = 1.0*10^ohevol.mdlogoh_noevol_cor
    oh4 = 1.0*10^ohevol.mdlogoh_levol_cor

    oh1_err = ohevol.ldlogoh_noevol_cor_err*oh1*alog(10.0)
    oh2_err = ohevol.ldlogoh_levol_cor_err*oh2*alog(10.0)
    oh3_err = ohevol.mdlogoh_noevol_cor_err*oh3*alog(10.0)
    oh4_err = ohevol.mdlogoh_levol_cor_err*oh4*alog(10.0)

;   oh1_err[0] = 0.0
;   oh2_err[0] = 0.0
;   oh3_err[0] = 0.0
    
; plot the mean values
    zz = z1
    moh = zz*0.0
    moh_err = zz*0.0
    weighted = 0
    for ii = 0, n_elements(zz)-1L do begin
       if weighted then begin
          moh[ii] = im_weighted_mean([oh2[ii],oh3[ii]],$
            [oh2_err[ii],oh3_err[ii]],wsigma=jj)
          moh_err[ii] = jj
       endif else begin
          oo = [oh1[ii],oh2[ii],oh3[ii],oh4[ii]]
;         oo = [oh2[ii],oh3[ii]]
          moh[ii] = mean(oo)
          moh_err[ii] = stddev(oo)
       endelse
    endfor

;   for ii = 0, n_elements(zz)-1L do im_oplot_box, 1.9*zerr[ii], $
;     2*moh_err[ii], 0.0, xoffset=zz[ii], yoffset=moh[ii]
    oploterror, zz, moh, zerr, moh_err, psym=-symcat(6,thick=10), $
      symsize=5.0, color=djs_icolor('dark red'), errcolor=djs_icolor('dark red'), $
      thick=8.0, errthick=10

; now plot the individual points    
    plotpoints = 0
    if plotpoints then begin
       oploterror, z1, oh1, zerr, oh1_err, psym=-symcat(psym1,thick=8), $
         symsize=symsize1, line=line1, color=im_color(color1,91), $
         errcolor=im_color(color1,91)
       oploterror, z2, oh2, zerr, oh2_err, psym=-symcat(psym2,thick=8), $
         symsize=symsize2, color=im_color(color2,92), errcolor=im_color(color2,92), line=line2
       oploterror, z3, oh3, zerr, oh3_err, psym=-symcat(psym3,thick=8), $
         symsize=symsize3, color=im_color(color3,93), errcolor=im_color(color3,93), line=line3
       oploterror, z4, oh4, zerr, oh4_err, psym=-symcat(psym4,thick=8), $
         symsize=symsize4, color=im_color(color4,94), errcolor=im_color(color4,94), line=line4
    endif

; label the Pegase curves
    xyouts, 0.52, 1.01, textoidl('\tau=1 Gyr'), align=0.5, $
      charsize=2.0, color=im_color('dodger blue',10)
    xyouts, 0.52, 0.96, textoidl('\tau=3 Gyr'), align=0.5, $
      charsize=2.0, color=im_color('orange',11), $
      orientation=-20
    xyouts, 0.3, 0.83, textoidl('\tau=!7y!6'), align=0.5, $
      charsize=2.2, color=im_color('purple',12), $
      orientation=-60
;   xyouts, 0.6, 0.79, textoidl('!MI!N'+'d\rho_{*}(t)'), $
;     align=0.5, charsize=2.0, color=im_color('black',13)
;   xyouts, 0.6, 0.77, textoidl('!MI!N'+'dt\prime'+'(d\rho_{*}/dt\prime)'), $
    xyouts, 0.605, 0.78, textoidl('!MI!N'+'dt'+'(d\rho_{*}/dt)'), $
      align=0.5, charsize=1.9, color=im_color('black',13)
    
; legend    
    label = ['\tau=1 Gyr','\tau=3 Gyr','\tau=const']
    cc = ['dodger blue','orange','purple']
    line = [5,3,1]
    spacing = 2.2

;   label = ['!MI!N'+'dt'+'\rho_{*}(t)','\tau=1 Gyr','\tau=3 Gyr','\tau=const']
;   cc = ['','dodger blue','orange','purple']
;   line = [0,5,3,1]
;   spacing = 2.0
;;  pp = [0.01,!y.crange[0]+1.5]
    
;   im_legend, label, /left, /bottom, box=0, pspacing=1.3, color=cc, $
;     textcolor=cc, line=line, spacing=spacing, margin=0, charsize=1.7, $
;     position=pp, thick=10

    im_plotconfig, /psclose

; ------------------------------------------------------------
; delta-log(O/H) with redshift

    psfile = pspath+'dlogohevol'+suffix
    im_plotconfig, 6, pos, psfile=psfile, charsize=2.0, $
      height=4.0*[1,1], ymargin=[0.8,1.0], xmargin=[1.4,0.1], thick=10

    ztitle1 = 'Redshift'
    dohtitle1 = textoidl('<\Delta'+'log(O/H)> from !8L-Z!6')
    dohtitle2 = textoidl('<\Delta'+'log(O/H)> from !8M-Z!6')
    zrange = [0.0,0.75]
    dohrange = [-0.499999,0.09]
;   dohrange = [-0.58,0.13]

    color1 = 'black'        & line1 = 5 & psym1 = 6  & symsize1 = 3.0
    color2 = 'red'          & line2 = 0 & psym2 = 15 & symsize2 = 3.0
    color3 = 'forest green' & line3 = 5 & psym3 = 4  & symsize3 = 4.0
    color4 = 'royal blue'   & line4 = 0 & psym4 = 14 & symsize4 = 4.0

    color5 = 'black'     & line5 = 5 & psym5 = 6 & symsize5 = 3.0
    color6 = 'navy'      & line6 = 0 & psym6 = 9 & symsize6 = 3.5
    color7 = 'firebrick' & line7 = 0 & psym7 = 16 & symsize7 = 3.5

; LZ    
    djs_plot, [0], [0], /nodata, ysty=1, xsty=9, position=pos[*,0], $
      xrange=zrange, yrange=dohrange, xtitle='', ytitle=dohtitle1, $
      xtickname=replicate(' ',10)
    axis, /xaxis, xsty=1, xrange=xrange, xtickv=getredshift(getage(0.0)-timelabel1), $
      xtitle='Lookback Time (Gyr)', xticks=n_elements(timelabel1)-1L, $
      xtickname=string(timelabel1,format='(I0)')
    djs_oplot, zz, zz*0.0, line=1, thick=3
;   cc = poly_fit([0.1,0.65],[0.0,-alog10(2.0)],1) ; factor of 2 to z=0.65
;   djs_oplot, zz, poly(zz,cc), line=1, thick=3
;   cc = poly_fit([0.1,0.65],[0.0,-alog10(4.0)],1) ; factor of 4
;   djs_oplot, zz, poly(zz,cc), line=1, thick=3
; Q=0
    oploterror, ohevol.z, ohevol.ldlogoh_noevol, ohevol.ldlogoh_noevol_err, $
      psym=-symcat(psym1,thick=8), color=im_color(color1,91), errcolor=im_color(color1,91), $
      symsize=symsize1, line=line1
    oploterror, ohevol.z, ohevol.ldlogoh_noevol_cor, ohevol.ldlogoh_noevol_cor_err, $
      psym=-symcat(psym2), color=im_color(color2,92), errcolor=im_color(color2,92), $
      symsize=symsize2, line=line2
; Q=1.5
    oploterror, ohevol.z, ohevol.ldlogoh_levol, ohevol.ldlogoh_levol_err, $
      psym=-symcat(psym3,thick=8), color=im_color(color3,93), errcolor=im_color(color3,93), $
      symsize=symsize3, line=line3
    oploterror, ohevol.z, ohevol.ldlogoh_levol_cor, ohevol.ldlogoh_levol_cor_err, $
      psym=-symcat(psym4), color=im_color(color4,94), errcolor=im_color(color4,94), $
      symsize=symsize4, line=line4

    im_legend, ['Q=0, Observed','Q=0, Corrected','Q=1.5, Observed','Q=1.5, Corrected'], $
      /left, /bottom, box=0, charsize=1.8, color=[color1,color2,color3,color4], $
      linestyle=[line1,line2,line3,line4], psym=[psym1,psym2,psym3,psym4], $
      symthick=7, symsize=0.6*[symsize1,symsize2,symsize3,symsize4]

; MZ
    djs_plot, [0], [0], /nodata, /noerase, ysty=1, xsty=1, position=pos[*,1], $
      xrange=zrange, yrange=dohrange, xtitle=ztitle1, ytitle=dohtitle2
    djs_oplot, zz, zz*0.0, line=1, thick=3

    oploterror, ohevol.z, ohevol.mdlogoh, ohevol.mdlogoh_err, $
      psym=-symcat(psym5,thick=8), color=im_color(color5,91), errcolor=im_color(color5,91), $
      symsize=symsize5, line=line5
; Q=1.5
    oploterror, ohevol.z, ohevol.mdlogoh_levol_cor, ohevol.mdlogoh_levol_cor_err, $
      psym=-symcat(psym7), color=im_color(color7,92), errcolor=im_color(color7,92), $
      symsize=symsize7, line=line7
; Q=0
    oploterror, ohevol.z, ohevol.mdlogoh_noevol_cor, ohevol.mdlogoh_noevol_cor_err, $
      psym=-symcat(psym6,thick=8), color=im_color(color6,92), errcolor=im_color(color6,92), $
      symsize=symsize6, line=line6

    im_legend, ['Observed','Q=0, Corrected','Q=1.5, Corrected'], $
      /left, /bottom, box=0, charsize=1.8, color=[color5,color6,color7], $
      linestyle=[line5,line6,line7], psym=[psym5,psym6,psym7], $
      symthick=7, symsize=0.7*[symsize5,symsize6,symsize7]

    im_plotconfig, /psclose

; ------------------------------------------------------------
; evolution of characteristic (O/H) with redshift

    psfile = pspath+'ohevol'+suffix
    im_plotconfig, 6, pos, psfile=psfile, charsize=1.8, $
      height=4.0*[1,1], ymargin=[0.8,1.0], xmargin=[1.3,0.2], thick=10

    ztitle1 = 'Redshift'
    ohtitle1 = textoidl('12+log(O/H) at M_{0.1g}=-20.5')
    ohtitle2 = textoidl('12+log(O/H) at 10^{10.5} M'+sunsymbol())
    zrange = [0.0,0.75]
    ohrange = [8.62,9.33]
;   ohrange = [8.599,9.32]

    color1 = 'black'        & line1 = 5 & psym1 = 6  & symsize1 = 3.0
    color2 = 'red'          & line2 = 0 & psym2 = 15 & symsize2 = 3.0
    color3 = 'forest green' & line3 = 5 & psym3 = 4  & symsize3 = 4.0
    color4 = 'royal blue'   & line4 = 0 & psym4 = 14 & symsize4 = 4.0

    color5 = 'black'     & line5 = 0 & psym5 = 6 & symsize5 = 3.0
    color6 = 'navy'      & line6 = 5 & psym6 = 9 & symsize6 = 3.5
    color7 = 'firebrick' & line7 = 5 & psym7 = 16 & symsize7 = 3.5
    
; LZ
    djs_plot, [0], [0], /nodata, ysty=1, xsty=9, position=pos[*,0], $
      xrange=zrange, yrange=ohrange, xtitle='', ytitle=ohtitle1, $
      xtickname=replicate(' ',10)
    axis, /xaxis, xsty=1, xrange=xrange, xtickv=getredshift(getage(0.0)-timelabel1), $
      xtitle='Lookback Time (Gyr)', xticks=n_elements(timelabel1)-1L, $
      xtickname=string(timelabel1,format='(I0)')
; AGES
    oploterror, glzevol_ages.z, glzevol_ages.levol_mean, glzevol_ages.levol_mean_err, $
      psym=-symcat(psym2), color=im_color(color2,92), errcolor=im_color(color2,92), $
      symsize=symsize2, line=line2
    oploterror, glzevol_ages.z, glzevol_ages.noevol_mean, glzevol_ages.noevol_mean_err, $
      psym=-symcat(psym1,thick=8), color=im_color(color1,91), errcolor=im_color(color1,91), $
      symsize=symsize1, line=line1
; SDSS2AGES/levol/noevol
    oploterror, glzevol_levol.z, glzevol_levol.levol_mean, glzevol_levol.levol_mean_err, $
      psym=-symcat(psym4), color=im_color(color4,94), errcolor=im_color(color4,94), $
      symsize=symsize4, line=line4
    oploterror, glzevol_noevol.z, glzevol_noevol.noevol_mean, glzevol_noevol.noevol_mean_err, $
      psym=-symcat(psym3,thick=8), color=im_color(color3,93), errcolor=im_color(color3,93), $
      symsize=symsize3, line=line3

    im_legend, ['AGES, Q=0','AGES, Q=1.5','Mock AGES, Q=0','Mock AGES, Q=1.5'], $
      /left, /bottom, box=0, charsize=1.7, color=[color1,color2,color3,color4], $
      linestyle=[line1,line2,line3,line4], psym=[psym1,psym2,psym3,psym4], $
      symthick=7, symsize=0.6*[symsize1,symsize2,symsize3,symsize4]

; MZ
    djs_plot, [0], [0], /nodata, /noerase, ysty=1, xsty=1, position=pos[*,1], $
      xrange=zrange, yrange=ohrange, xtitle=ztitle1, ytitle=ohtitle2
; AGES
    oploterror, mzevol_ages.z, mzevol_ages.mean, mzevol_ages.mean_err, $
      psym=-symcat(psym5,thick=8), color=im_color(color5,91), errcolor=im_color(color5,92), $
      symsize=symsize5, line=line5
; SDSS2AGES/levol
    oploterror, mzevol_levol.z, mzevol_levol.mean, mzevol_levol.mean_err, $
      psym=-symcat(psym7), color=im_color(color7,94), errcolor=im_color(color7,94), $
      symsize=symsize7, line=line7
; SDSS2AGES/noevol
    oploterror, mzevol_noevol.z, mzevol_noevol.mean, mzevol_noevol.mean_err, $
      psym=-symcat(psym6,thick=8), color=im_color(color6,93), errcolor=im_color(color6,93), $
      symsize=symsize6, line=line6
   
    im_legend, ['AGES','Mock AGES, Q=0','Mock AGES, Q=1.5'], $
      /left, /bottom, box=0, charsize=1.7, color=[color5,color6,color7], $
      linestyle=[line5,line6,line7], psym=[psym5,psym6,psym7], $
      symthick=7, symsize=0.6*[symsize5,symsize6,symsize7]

    im_plotconfig, /psclose

return
end
    

;;
;;; ---------------------------------------------------------------------------
;;; redshift versus SFR density
;;
;;    psfile = pspath+'z_vs_sfrd'+suffix
;;    im_plotconfig, 6, pos, psfile=psfile, charsize=1.7, height=[4.0,2.0], $
;;      ymargin=[0.9,1.1], xmargin=[1.3,0.2], thick=6
;;    
;;    h = rsex(litpath+'04hopkins.sex')
;;    uv = where(strmatch(h.indicator,'*UV*',/fold))
;;    ha = where(strmatch(h.indicator,'*Ha*',/fold) or $
;;      strmatch(h.indicator,'*Hb*',/fold) or $
;;      strmatch(h.indicator,'*OII*',/fold))
;;    ir = where(strmatch(h.indicator,'*IR*',/fold))
;;    rad = where(strmatch(h.indicator,'*RADIO*',/fold) or $
;;      strmatch(h.indicator,'*xray*',/fold))
;;
;;    uvcolor = 'royal blue'
;;    hacolor = 'forest green'
;;    ircolor = 'firebrick'
;;    radcolor = 'orange'
;;
;;    uvpsym = 16
;;    hapsym = 14
;;    irpsym = 12
;;    radpsym = 18
;;    
;;    xrange = [-0.05,0.92]
;;    yrange1 = [-2.3,-0.2]
;;    yrange2 = [0.6,1]
;;    xtitle = 'log (1 + z)'
;;    ytitle1 = 'log \rho_{SFR} (M'+sunsymbol()+' yr^{-1} Mpc^{-3} h_{70}^{3})'
;;    ytitle2 = 'Fraction of Stars Formed by z=0.1'
;;    
;;    djs_plot, [0], [0], /nodata, xtitle='', ytitle=ytitle, $
;;      xsty=9, ysty=1, xrange=xrange, yrange=yrange1, position=pos[*,0], $
;;      xtickname=replicate(' ',10)
;;    axis, /xaxis, xthick=postthick1, xsty=1, xrange=xrange, xtitle='Lookback Time (Gyr)', $
;;      xtickv=alog10(1.0+getredshift(getage(0.0)-timelabel2)), $
;;      xticks=n_elements(timelabel2)-1L, xtickname=string(timelabel2,format='(I0)')
;;
;;; UV    
;;    oploterror, alog10(1+h[uv].z), h[uv].sfrd, h[uv].zerr_lo/(1+h[uv].z)/alog(10.0), $
;;      h[uv].sfrderr_lo, psym=symcat(uvpsym), color=im_color(uvcolor,10), $
;;      errcolor=im_color(uvcolor,10), errthick=!p.thick, /lobar
;;    oploterror, alog10(1+h[uv].z), h[uv].sfrd, h[uv].zerr_hi/(1+h[uv].z)/alog(10.0), $
;;      h[uv].sfrderr_hi, psym=symcat(uvpsym), color=im_color(uvcolor,10), $
;;      errcolor=im_color(uvcolor,10), errthick=!p.thick, /hibar
;;; Ha
;;    oploterror, alog10(1+h[ha].z), h[ha].sfrd, h[ha].zerr_lo/(1+h[ha].z)/alog(10.0), $
;;      h[ha].sfrderr_lo, psym=symcat(hapsym), color=im_color(hacolor,10), $
;;      errcolor=im_color(hacolor,10), errthick=!p.thick, /lobar
;;    oploterror, alog10(1+h[ha].z), h[ha].sfrd, h[ha].zerr_hi/(1+h[ha].z)/alog(10.0), $
;;      h[ha].sfrderr_hi, psym=symcat(hapsym), color=im_color(hacolor,10), $
;;      errcolor=im_color(hacolor,10), errthick=!p.thick, /hibar
;;; IR
;;    oploterror, alog10(1+h[ir].z), h[ir].sfrd, h[ir].zerr_lo/(1+h[ir].z)/alog(10.0), $
;;      h[ir].sfrderr_lo, psym=symcat(irpsym), color=im_color(ircolor,10), $
;;      errcolor=im_color(ircolor,10), errthick=!p.thick, /lobar
;;    oploterror, alog10(1+h[ir].z), h[ir].sfrd, h[ir].zerr_hi/(1+h[ir].z)/alog(10.0), $
;;      h[ir].sfrderr_hi, psym=symcat(irpsym), color=im_color(ircolor,10), $
;;      errcolor=im_color(ircolor,10), errthick=!p.thick, /hibar
;;; Radio
;;    oploterror, alog10(1+h[rad].z), h[rad].sfrd, h[rad].zerr_lo/(1+h[rad].z)/alog(10.0), $
;;      h[rad].sfrderr_lo, psym=symcat(radpsym), color=im_color(radcolor,10), $
;;      errcolor=im_color(radcolor,10), errthick=!p.thick, /lobar
;;    oploterror, alog10(1+h[rad].z), h[rad].sfrd, h[rad].zerr_hi/(1+h[rad].z)/alog(10.0), $
;;      h[rad].sfrderr_hi, psym=symcat(radpsym), color=im_color(radcolor,10), $
;;      errcolor=im_color(radcolor,10), errthick=!p.thick, /hibar
;;    
;;; coefficients from plot_redshift_vs_sfr_mass_density    
;;    sfrd_zaxis = im_array(0.0,6.0,0.05)
;;    coeff = [0.010464392D,0.11303816D,3.4534319D,3.7109016D]
;;    sfrd = (coeff[0]+coeff[1]*sfrd_zaxis)/(1.0+(sfrd_zaxis/coeff[2])^coeff[3])
;;    djs_oplot, alog10(1+sfrd_zaxis), alog10(sfrd), line=0
;;
;;; now compute and overplot the fraction of stars formed versus
;;; redshift
;;    sfrd_lookback = getage(0.0)-getage(sfrd_zaxis)  ; lookback time
;;    sfrd_total = im_integral(sfrd_lookback,sfrd*1D9,$              ; M_sun/Mpc^3
;;      min(sfrd_lookback),max(sfrd_lookback))                       
;;    sfrd_time = im_integral(sfrd_lookback,sfrd*1D9,sfrd_lookback,$ ; M_sun/Mpc^3
;;      replicate(max(sfrd_lookback),n_elements(sfrd_lookback)))
;;    sfrd_frac = sfrd_time/sfrd_total
;;
;;; extrapolate to z<0.1    
;;    plot_sfrd_zaxis = [0.0,sfrd_zaxis]
;;    notzero = where(sfrd_frac gt 0.0)
;;    plot_sfrd_frac1 = interpol((sfrd_frac[notzero]),$
;;      sfrd_zaxis[notzero],plot_sfrd_zaxis)
;;
;;    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle=ytitle2, $
;;      xsty=1, ysty=1, xrange=xrange, yrange=yrange2, position=pos[*,1]
;;    djs_oplot, plot_sfrd_zaxis, plot_sfrd_frac1, line=0, thick=12
;;    
;;    im_plotconfig, /psclose
;;

;;
;;; --------------------------------------------------
;;; average dlog-O/H vs redshift in bins of mass
;;    mzavg = mrdfits(mzpath+'mzevol_avg.fits.gz',1,/silent)
;;    
;;    psfile = pspath+'z_vs_dlogoh_avg'+suffix
;;    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.7, $
;;      xmargin=[1.4,0.4]
;;
;;    xrange = [-0.02,0.75]
;;    yrange = [-0.4,0.15]
;;    xtitle1 = 'Redshift'
;;    ytitle1 = textoidl('<\Delta'+'log(O/H)>')
;;
;;    zaxis = range(0.0,xrange[1]-0.03,50)
;;
;;    djs_plot, [0], [0], /nodata, ysty=1, xsty=1, position=pos, $
;;      xrange=xrange, yrange=yrange, xtitle=xtitle1, ytitle=ytitle1
;;    djs_oplot, !x.crange, [0,0], line=5, thick=3
;;    im_legend, massbins[0:nmassbins-2].label, /left, /bottom, box=0, $
;;      charsize=1.3, psym=massbins[0:nmassbins-2].psym, color=massbins[0:nmassbins-2].color
;;    
;;    for mm = 0, nmassbins-2 do begin
;;       good = where(mzavg.dlogoh[mm] gt -900.0)
;;       oploterror, zbins[good].zbin, mzavg[good].dlogoh[mm], mzavg[good].dlogoh_err[mm], $
;;         psym=symcat(massbins[mm].psym,thick=5), symsize=massbins[mm].symsize, $
;;         color=im_color(massbins[mm].color,101), errcolor=im_color(massbins[mm].color,101)
;;    endfor
;;    
;;    im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)
;;

;;
;;
;;; --------------------------------------------------
;;; old Figure 15 - average d[log(O/H)] vs redshift in bins of mass 
;;    psfile = pspath+'z_vs_dlogoh_avg'+suffix
;;    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.7, $
;;      xmargin=[1.4,0.4]
;;
;;    xrange = [-0.02,0.7]
;;    yrange = [-0.15,0.03]
;;    xtitle1 = 'Redshift'
;;;   ytitle1 = textoidl('log (O/H)_{z} - log (O/H)_{z=0.1}')
;;    ytitle1 = textoidl('\Delta'+'<12 + log (O/H)>')
;;
;;    zaxis = range(0.0,xrange[1]-0.02,50)
;;
;;    djs_plot, [0], [0], /nodata, ysty=1, xsty=1, position=pos, $
;;      xrange=xrange, yrange=yrange, xtitle=xtitle1, ytitle=ytitle1, $
;;      ytickname=['-0.15','-0.10','-0.05','0.00']
;;    djs_oplot, zaxis, zaxis*0, line=0, thick=3, color='grey'
;;;   maxmodel = poly(zaxis-mzavg.qz0,mzavg.coeffs+mzavg.coeffs_err/2)
;;;   minmodel = poly(zaxis-mzavg.qz0,mzavg.coeffs-mzavg.coeffs_err/2)
;;;   polyfill, [zaxis,reverse(zaxis)], [minmodel,reverse(maxmodel)], $
;;;     /data, /fill, color=im_color('tan',101), noclip=0
;;;   oploterror, mzavg.zbin, mzavg.dlogoh, mzavg.dlogoh_err, $
;;;     psym=symcat(15,thick=5), symsize=2.0, color=im_color('firebrick',101), $
;;;     errcolor=im_color('firebrick',101)
;;
;;;   xyouts, 0.15, -0.10, textoidl('log (M/M_{'+sunsymbol()+'})'), align=0.5, $
;;;     /data, charsize=1.4
;;;   im_legend, ['>11','10.5-11','10-10.5','9.5-10'], /left, /bottom, box=0, $
;;;     charsize=1.3, psym=-massbins[0:nmassbins-2].psym, $
;;;     color=massbins[0:nmassbins-2].color, line=massbins[0:nmassbins-2].line, $
;;;     thick=6, pspacing=1.8
;;
;;;   im_legend, massbins.label, /left, /bottom, box=0, charsize=1.3, $
;;;     psym=-massbins.psym, color=massbins.color, line=massbins.line, $
;;;     thick=6, pspacing=1.8, spacing=1.9, margin=0
;;
;;;   for mm = 1, 1 do begin
;;    for mm = 0, nmassbins-1 do begin
;;; AGES
;;       good = where(mzavg.dlogoh_bymass[masskeep[mm],*] gt -900.0,ngood)
;;       oploterror, mzavg.medz_bymass[masskeep[mm],good], $
;;         mzavg.dlogoh_bymass[masskeep[mm],good], mzavg.dlogoh_bymass_err[masskeep[mm],good], $
;;         psym=symcat(massbins[mm].psym,thick=7), symsize=massbins[mm].symsize, $
;;         color=im_color(massbins[mm].color), errcolor=im_color(massbins[mm].color), $
;;         errthick=!p.thick
;;
;;; SDSS
;;       oploterror, mzavg.sdss_medz_bymass[masskeep[mm]], $
;;         mzavg.sdss_dlogoh_bymass[masskeep[mm]], mzavg.sdss_dlogoh_bymass_err[masskeep[mm],good], $
;;         psym=symcat(massbins[mm].psym,thick=7), symsize=massbins[mm].symsize, $
;;         color=im_color(massbins[mm].color), errcolor=im_color(massbins[mm].color), $
;;         errthick=!p.thick
;;
;;; the model fit
;;       djs_oplot, zaxis, poly(zaxis-mzavg.qz0,mzavg.coeffs_bymass[*,masskeep[mm]])-$
;;         poly(0.0,mzavg.coeffs_bymass[*,masskeep[mm]]), thick=7, $
;;         color=im_color(massbins[mm].color), line=massbins[mm].line
;;    endfor
;;
;;;   djs_oplot, zaxis, poly(zaxis-mzavg.qz0,mzavg.coeffs), line=0
;;;   djs_oplot, zaxis, poly(zaxis-mzavg.qz0,mzavg.coeffs+mzavg.coeffs_err/2), line=0
;;;   djs_oplot, zaxis, poly(zaxis-mzavg.qz0,mzavg.coeffs-mzavg.coeffs_err/2), line=0
;;    
;;    im_plotconfig, /psclose, psfile=psfile
;;
