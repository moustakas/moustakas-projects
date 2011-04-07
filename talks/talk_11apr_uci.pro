pro talk_11apr_uci, keynote=keynote, noevol=noevol
; jm11apr01ucsd - miscellaneous plots for my 2011-Apr UC Irvine talk  

    common com_talk, sdss_parent, sdss_ised
    
    mfpath = mf_path(zerodversion=zerodver)
    isedpath = mf_path(/isedfit)
    
    rootpath = getenv('IM_RESEARCH_DIR')+'/talks/2011/11apr_uci/'
    datapath = rootpath+'data/'
    if keyword_set(keynote) then talkpath = rootpath+'keynote/' else $
      talkpath = rootpath+'figures/'

    noevol = 1 ; NOTE!
    if keyword_set(noevol) then esuffix = 'noevol' else esuffix = 'evol'

; read the supergrid parameter file    
    supergrid_paramfile = mf_supergrid_parfile()
    super = yanny_readone(supergrid_paramfile)
    superstring = string(super.supergrid,format='(I2.2)')
    nsuper = n_elements(super)

; --------------------------------------------------
; SED-fitting example
    field = 'cfhtls_xmm'
    sfhgrid = 3
    supergrid = 4

    isedpath = mf_path(/isedfit)
    sfhgrid_basedir = mf_path(/montegrids)

    paramfile = isedpath+field+'_isedfit.par'
    params = read_isedfit_paramfile(paramfile,sfhgrid=sfhgrid)
    params.sfhgrid = string(sfhgrid,format='(I2.2)')
    params.imf = 'salp'
    params.synthmodels = 'bc03'
    params.redcurve = redcurve2string(1) ; charlot

    ii = read_mf_ubersample(field)
    spherematch, ii.ra, ii.dec, 15D*hms2dec('02:19:12.1'), hms2dec('-04:01:08.4'),$
      5.0/3600.0, m1, m2 & help, m1
    galaxy = 'J021912.1-040108.4'    
    index = m1[0]
    isedfit = read_mf_isedfit(field,supergrid=supergrid,rows=index,post=post)
    isedfit = isedfit[0] & post = post[0]

    model = isedfit_restore(paramfile,junk,params=params,$
      iopath=isedpath,index=index,sfhgrid_basedir=sfhgrid_basedir)
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

; inset with P(M)     
    im_plothist, post.mass, bin=0.005, xsty=3, ysty=1, yrange=[0,120], $
      position=[0.55,0.35,0.9,0.55], /noerase, /fill, fcolor='dark green', $
      ytitle='P(M)', xtitle='Stellar Mass', color=keycolor, $
      ytickname=replicate(' ',10), charsize=1.5, xrange=[10.65,11.04]
    oplot, isedfit.mass_50*[1,1], !y.crange, line=0, thick=6, color=djs_icolor('black')
    oplot, isedfit.mass*[1,1], !y.crange, line=5, thick=6, color=djs_icolor('black')

    im_plotconfig, psfile=psfile, /psclose, keynote=keynote

stop    

; --------------------------------------------------
; Omega* vs redshift + SDSS + literature    
    H0 = (100.0*mf_h100())*3.1556926D7*1D6/3.086D19 ; [Myr^-1]
    bigG = 4.4983D-21 ; [Mpc^3/Myr^2/M_sun]
    rhocrit = 3.0*H0^2.0/(8.0*!pi*bigG) ; [M_sun/Mpc^3, h=0.7]
;   rhocrit = 1.35D11 ; [M_sun/Mpc^3, h=0.7, from Pettini 2006]

    zbins = mf_zbins(nzbins)
    field = get_mf_fields(nice=nice)
    nfield = n_elements(field)
    prefs = mfplot_prefs(field)
    prefs.color = ['purple','dark green','firebrick','dodger blue','orange'];,'orange','navy','blue']

;   yrange = [-3.2,-1.7]
    yrange = [-2.9,-1.8]
    
    psfile = talkpath+'z_vs_omega.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.4,0.3], $
      width=6.8, height=5.4, charsize=2, keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')
    
    plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      yrange=yrange, xrange=[0,1.03], xtitle='Redshift', $
      ytitle=mfplot_omegatitle(), ytickinterval=0.5, color=keycolor
    
; wilkins+; Chabrier --> Salpeter
    wilk = rsex(getenv('IM_PAPERS_DIR')+'/literature/data/08wilkins.sex')
    oploterror, total([[wilk.zmin],[wilk.zmax]],2)/2.0, $
      alog10(wilk.omega)+0.25, wilk.omega_err/wilk.omega/alog(10), psym=symcat(9,thick=6), $
      color=fsc_color('tan',101), errcolor=fsc_color('tan',101), $
      symsize=1.0

; primus - average and then field-by-field
    allomega = fltarr(nzbins,nfield)
    for ii = 0, n_elements(field)-1 do begin
       rho = get_mf_param(field[ii],/rho,param_err=rho_err,noevol=noevol,supergrid=8);/avgsuper)
       omega = alog10(rho/rhocrit)
       omega_err = rho_err/rho/alog(10)
       oploterror, zbins.zbin, omega, omega_err, psym=symcat(prefs[ii].psym,thick=8), $
         symsize=2.2, errthick=8, color=fsc_color(prefs[ii].color,101), $
         errcolor=fsc_color(prefs[ii].color,101)
       allomega[*,ii] = omega
    endfor
    splog, 'HACK!!!!'
;   rho = get_mf_param(/rho,param_err=rho_err,noevol=noevol,/avgfield,supergrid=8);,/avgsuper)
;   omega = alog10(rho/rhocrit)
;   omega_err = rho_err/rho/alog(10)
    omega = fltarr(nzbins)
    omega_err = fltarr(nzbins)
    for ii = 0, nzbins-1 do omega[ii] = djs_mean(allomega[ii,*])
    for ii = 0, nzbins-1 do omega_err[ii] = djsig(allomega[ii,*])
    oploterror, zbins.zbin, omega, omega_err, psym=symcat(6,thick=10), $
      symsize=2.7, errthick=8
    niceprint, zbins.zbin, omega_err

; sdss    
    sdsszbins = mf_zbins(/sdss)
    rho = get_mf_param('sdss',/rho,param_err=rho_err,noevol=noevol,/avgsuper)
    sdssomega = alog10(rho/rhocrit)
    sdssomega_err = rho_err/rho/alog(10)
    plots, sdsszbins.zbin, sdssomega, psym=symcat(36,thick=8), $
      symsize=3.0, color=fsc_color('navy',101)
;   oploterror, sdsszbins.zbin, sdssomega, sdssomega_err, psym=symcat(6,thick=8), $
;     symsize=3.9, errthick=8, color=fsc_color('navy',101)

    im_legend, ['SDSS',nice], psym=[36,prefs.psym], color=['navy',strtrim(prefs.color,2)], $
      /left, /bottom, box=0, margin=0, charsize=1.4, $
      symsize=[1.2,prefs.symsize]*1.8, textcolor=keycolor
    
    
    im_plotconfig, psfile=psfile, /psclose, keynote=keynote

stop    
    
; -------------------------
; primus - loop on each redshift interval and intercompare the MFs
; from the various fields, for each subsample
    supergrid = 8
    fields = get_mf_fields()
    psfile = talkpath+'binned_mf_all_'+esuffix+'.ps'
    plot_binned_mf, noevol=noevol, psfile=psfile, field=fields, supergrid=supergrid, keynote=keynote ; hack!

stop    
    
; --------------------------------------------------
; paper plot: redshift success rate in bins of magnitude and color

; read the completeness map    
    compfile = mfpath+'mf_completeness_'+zerodver+'.fits.gz'
    splog, 'Reading '+compfile
    completeness = mrdfits(compfile,1)

    levels = [0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.99]
    cannotation = string(levels,format='(F4.2)')
    
    nmock = 500 ; 1000 ; 2000L
    mincol = 50 ; 100 
    ndiv = 5
    ticknames = string(range(0.0,1.0,ndiv+1),format='(F5.2)')

    loadct, 16 ; 15
    psfile = talkpath+'completeness_mag_vs_color.ps'
    im_plotconfig, 17, pos, psfile=psfile, charsize=1.4, $
      xmargin=[0.8D,0.3D], ymargin=[0.2D,0.9D], xspace=0.8D*[1,1], $
      yspace=1D, width=2.7D*[1,1,1], height=2.8D*[1,1], keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')
    
    field = get_mf_fields(/nodeep2,/nodls,nice_field=nice_field)
    nfield = n_elements(field)
    for ii = 0, nfield-1 do begin
       params = get_mf_completeness_params(field[ii])
       magaxis = mf_grid_axis(params,0)
       coloraxis = mf_grid_axis(params,1)
       mag2daxis = magaxis#(coloraxis*0.0+1.0)
       color2daxis = transpose(coloraxis#(magaxis*0.0+1.0))
       
       magrange = [params.hmin[0],params.hmax[0]] ;  [18.0,23.5]
       colorrange = [params.hmin[1],params.hmax[1]] ; [-0.5,2.0]
       mag = randomu(seed,nmock,nmock)*(magrange[1]-magrange[0])+magrange[0]
       color = randomu(seed,nmock,nmock)*(colorrange[1]-colorrange[0])+colorrange[0]
       
       grid = completeness[where(field[ii] eq strtrim(completeness.field,2))]
       zsuccess = get_mf_zsuccess(field[ii],mag=mag,color1=color,color2=color*0-999.0,/noextrap)

       good = where(zsuccess gt 0.0)
       mf_hogg_scatterplot, mag[good], color[good], weight=zsuccess[good], position=pos[*,ii], $
         noerase=(ii gt 0), xsty=1, ysty=1, xrange=magrange, yrange=colorrange, xnpix=params.nbins[0], $
         ynpix=params.nbins[1], xtitle=textoidl(params.nice_targ_filter+' (AB mag)'), $
         ytitle=textoidl(params.nice_color_filter[0]), exp=2.0, $
         darkest=mincol, levels=levels, /nocontour, color=keycolor, keynote=keynote
       contour, grid.fracobj_color1, mag2daxis, color2daxis, levels=levels, $
         /overplot, c_annotation=cannotation

       legend, nice_field[ii], position=[pos[0,ii]-0.01,pos[3,ii]], $
         /norm, box=0, charsize=1.4, charthick=3.5, textcolor=keycolor
    endfor       

    off1 = 0.03 & hei = 0.3 & wid = 0.08 & off2 = 0.14
    cpos = [pos[0,ii]+off1,pos[1,ii],pos[0,ii]+wid,pos[3,ii]]
;   off1 = 0.03 & hei = 0.3 & wid = 0.05 & off2 = 0.13
;   cpos = [pos[2,1]+off1,(pos[1,1]-pos[3,3])/2+pos[3,3]-hei,$
;     pos[2,1]+off1+wid,(pos[1,1]-pos[3,3])/2+pos[3,3]+hei]
    primus_fsc_colorbar, position=cpos, range=[0.0,1.0], /vertical, $
      charsize=1.4, minor=4, division=ndiv, ticknames=ticknames, $
      /invert, bottom=mincol, /right, color=keycolor
    xyouts, cpos[0]+off2, (cpos[3]-cpos[1])/2.0+cpos[1], 'Redshift Success', $
      orientation=270, align=0.5, charsize=1.8, /norm, color=keycolor
    
    im_plotconfig, psfile=psfile, /psclose, keynote=keynote
    loadct, 0

stop
    
; -------------------------
; sdss - compare MFs derived using different prior assumptions
    psfile = talkpath+'binned_mf_sdss_differentpriors.ps'
    im_plotconfig, 3, pos, psfile=psfile, charsize=1.4, keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')
    
    xrange = [8.8,12.6]
    yrange = [-7.2,-1.3]
    maxis1 = range(xrange[0]+0.1,12.5,100)

    showerr = 0
    
; left - all galaxies    
    plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytitle=mfplot_phititle(), color=keycolor
    legend, 'All', /right, /top, box=0, charsize=1.2, margin=0, textcolor=keycolor
    render_priorplot, 'all', esuffix=esuffix, super=super, $
      xrange=xrange, yrange=yrange, /dolegend, textcolor=keycolor
    
; middle - quiescent galaxies    
    plot, [0], [0], /nodata, position=pos[*,1], /noerase, xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytickname=replicate(' ',10), color=keycolor
    legend, 'Quiescent', /right, /top, box=0, charsize=1.2, margin=0, textcolor=keycolor
    render_priorplot, 'quiescent', esuffix=esuffix, super=super, $
      xrange=xrange, yrange=yrange, textcolor=keycolor

; right - active galaxies    
    plot, [0], [0], /nodata, position=pos[*,2], /noerase, xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytickname=replicate(' ',10), color=keycolor
    legend, 'Star-Forming', /right, /top, box=0, charsize=1.2, margin=0, textcolor=keycolor
    render_priorplot, 'active', esuffix=esuffix, super=super, $
      xrange=xrange, yrange=yrange, textcolor=keycolor

    im_plotconfig, psfile=psfile, /psclose, keynote=keynote

stop    
    
; -------------------------
; sdss - compare the supergridavg MFs for all, quiescent, and active 
    psfile = talkpath+'binned_mf_sdss.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.8, height=6.5, keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')
    
    xrange = [8.8,12.4]
    yrange = [-7.4,-1.1]
    maxis1 = range(xrange[0]+0.1,12.5,100)

    fitcolor = ['black','red','navy blue']
    mfcolor = ['white','tomato','dodger blue']
    symsize = [2.0,2.0,2.5]
    mfpsym = [6,9,4]
    mfline = [0,5,3]

    plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytitle=mfplot_phititle(), color=keycolor
    legend, textoidl('Salpeter IMF (0.1-100 M_{'+sunsymbol()+'})'), $
      /right, /top, box=0, charsize=1.7, textcolor=keycolor
    im_legend, ['All','Quiescent','Star-Forming'], /left, /bottom, box=0, $
      color=mfcolor, psym=mfpsym, symsize=symsize, symthick=8, $
      textcolor=keycolor

    usepolyfill = 0
    
    subsample = ['all','quiescent','active']
    for jj = 0, n_elements(subsample)-1 do begin
;      mfdata = read_binned_mf('sdss',quiescent=jj eq 1,active=jj eq 2,supergrid=8)
;      mffit = read_binned_mf('sdss',quiescent=jj eq 1,active=jj eq 2,supergrid=8,/bestfit)
       mfdata = read_binned_mf('sdss',quiescent=jj eq 1,active=jj eq 2,/avgsuper)
       mffit = read_binned_mf('sdss',quiescent=jj eq 1,active=jj eq 2,/avgsuper,/bestfit)
       print, mffit.rho_tot, mffit.alpha
;      if (jj eq 0) then mffit_double = read_binned_mf('sdss',/avgsuper,/double_bestfit)

       good = where(mfdata.phierr gt 0.0)
       mass = mfdata.mass[good]
       phi = mfdata.phi[good]
       phierr = mfdata.phierr[good]

       if usepolyfill then begin
          phimin = alog10(phi-phierr)
          phimax = alog10(phi+phierr)
          polyfill, [mass,reverse(mass)],[phimin,reverse(phimax)], $
            /data, fill=(jj ne 1), color=fsc_color(mfcolor[jj],100), noclip=0

;         polyfill, [mfdata.mass,reverse(mfdata.mass)],[phimin,reverse(phimax)], $
;           /data, color=fsc_color(mfcolor[jj],100), noclip=0, $
;           /line_fill, orientation=45, spacing=0.01, linestyle=1
;         polyfill, [mfdata.mass,reverse(mfdata.mass)],[phimin,reverse(phimax)], $
;           /data, color=fsc_color(mfcolor[jj],100), noclip=0, $
;           /line_fill, orientation=135, spacing=0.01, linestyle=1
       endif else begin
          oploterror, mass, alog10(phi), phierr/phi/alog(10), $
            psym=symcat(mfpsym[jj],thick=8), symsize=1.6, errthick=!p.thick, $
            color=fsc_color(mfcolor[jj],100+jj), $
            errcolor=fsc_color(mfcolor[jj],100+jj)
       endelse
          
;      djs_oplot, maxis1, alog10(mf_schechter(maxis1,mffit)), $
;        line=mfline[jj], thick=10, color=fsc_color(fitcolor[jj],101)
;      if (jj eq 0) then djs_oplot, maxis1, alog10(mf_schechter_plus($
;        maxis1,mffit_double)), line=1, thick=8, color=fsc_color(fitcolor[jj],101)
    endfor

    im_plotconfig, psfile=psfile, /psclose, keynote=keynote

stop    





; --------------------------------------------------
; prior grid distributions
    priorfile = mf_path(/monte)+'sfhgrid03/bc03/charlot/salp_montegrid.fits.gz'
    prior = mrdfits(priorfile,1)

    psfile = talkpath+'priors.ps'
    im_plotconfig, 15, pos, psfile=psfile, charsize=1.6, charthick=3.5, keynote=keynote, $
      xmargin=[0.3,0.25], xspace=[0.15,0.15], yspace=[0.8,0.8], width=[2.5,2.5,2.5], $
      height=2.0*[1,1,1]

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    taubin = 0.4
    agebin = 0.2
    zbin = 0.002
    dustbin = 0.03
    mubin = 0.02
    nbin = 0.5
    tbbin = 0.3
    dtbin = 0.005
    fbbin = 0.05
    
    im_plothist, prior.tau, bin=taubin, /noplot, xtau, ytau
    im_plothist, float(prior.age), bin=agebin, /noplot, xage, yage
    im_plothist, prior.Z, bin=zbin, /noplot, xz, yz
    im_plothist, prior.ebv, bin=dustbin, /noplot, xdust, ydust
    im_plothist, prior.mu, bin=mubin, /noplot, xmu, ymu
    im_plothist, float(prior.nburst), bin=nbin, /noplot, xnb, ynb

    isb = where(prior.tburst gt -1)
    im_plothist, float((prior.tburst)[isb]), bin=tbbin, /noplot, xtb, ytb
    im_plothist, float((prior.dtburst)[isb]), bin=dtbin, /noplot, xdtb, ydtb
    im_plothist, alog10(float((prior.fburst)[isb])), bin=fbbin, /noplot, xfb, yfb

    im_plothist, prior.tau, bin=taubin, position=pos[*,0], yrange=[0,max(ytau)*1.05], $
      xsty=1, ysty=1, xtitle='\tau (Gyr)', ytickname=replicate(' ',10), $
      /fill, fcolor=fsc_color('powder blue',101), color=keycolor
    im_plothist, float(prior.age), bin=agebin, /noerase, position=pos[*,1], yrange=[0,max(yage)*1.1], $
      xsty=1, ysty=1, xtitle='t_{form} (Gyr)', ytickname=replicate(' ',10), $
      /fill, fcolor=fsc_color('powder blue',101), color=keycolor
    im_plothist, prior.Z, bin=zbin, /noerase, position=pos[*,2], yrange=[0,max(yz)*1.1], $
      xsty=1, ysty=1, xtitle='Metallicity', ytickname=replicate(' ',10), $
      xtickinterval=0.015, /fill, fcolor=fsc_color('powder blue',101), color=keycolor

    im_plothist, prior.ebv, bin=dustbin, /noerase, position=pos[*,3], yrange=[0,max(ydust)*1.1], $
      xsty=1, ysty=1, xtitle='E(B-V)', ytickname=replicate(' ',10), $
      /fill, fcolor=fsc_color('powder blue',101), color=keycolor
    im_plothist, prior.mu, bin=mubin, /noerase, position=pos[*,4], yrange=[0,max(ymu)*1.1], $
      xsty=1, ysty=1, xtitle='\mu', ytickname=replicate(' ',10), $
      /fill, fcolor=fsc_color('powder blue',101), color=keycolor
    
    im_plothist, float(prior.nburst), bin=nbin, /noerase, position=pos[*,5], yrange=[0,max(ynb)*1.1], $
      xsty=1, ysty=1, xtitle='N_{burst}', ytickname=replicate(' ',10), $
      /fill, fcolor=fsc_color('powder blue',101), color=keycolor
    im_plothist, float((prior.tburst)[isb]), bin=tbbin, /noerase, position=pos[*,6], yrange=[0,max(ytb)*1.1], $
      xsty=1, ysty=1, xtitle='t_{burst}', ytickname=replicate(' ',10), $
      /fill, fcolor=fsc_color('powder blue',101), color=keycolor
    im_plothist, float((prior.dtburst)[isb]), bin=dtbin, /noerase, position=pos[*,7], yrange=[0,max(ydtb)*1.1], $
      xsty=1, ysty=1, xtitle='\Delta'+'t_{burst}', ytickname=replicate(' ',10), $
      /fill, fcolor=fsc_color('powder blue',101), xtickinterval=0.1, color=keycolor
    im_plothist, alog10(float((prior.fburst)[isb])), bin=fbbin, /noerase, position=pos[*,8], yrange=[0,max(yfb)*1.1], $
      xsty=1, ysty=1, xtitle='log (F_{burst})', ytickname=replicate(' ',10), $
      /fill, fcolor=fsc_color('powder blue',101), color=keycolor, xtickinterval=1.0
    
    
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; SDSS - SFR/M vs M coded by UV-optical color
    if (n_elements(sdss_parent) eq 0L) then sdss_parent = read_mf_parent('sdss')
    if (n_elements(sdss_ised) eq 0L) then sdss_ised = read_mf_isedfit('sdss',/parent,supergrid=4)
    mz = sdss_parent.k_ugrizjhk_absmag_05[4]
    umz = sdss_parent.k_ugrizjhk_absmag_05[0]-sdss_parent.k_ugrizjhk_absmag_05[4]
    mass = sdss_ised.mass_50
    sfrm = (alog10(sdss_ised.sfr100_avg)+9-mass)>(-4) ; [Gyr^-1]
    weight = sdss_parent.final_weight

    psfile = talkpath+'mass_vs_sfrm_sdss.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, $
      width=6.7, xmargin=[1.4,0.4], height=6.3, keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    xrange = [8.6,12.3]
    yrange = [-4.5,+0.7]

    xtitle = mfplot_masstitle()
    ytitle = mfplot_sfrmtitle()
    maxis = range(xrange[0]+0.1,xrange[1]-0.1,50)

    qq = mf_select_quiescent(umz,mz,z=sdss_parent.zprimus)
    aa = mf_select_quiescent(umz,mz,z=sdss_parent.zprimus,/active)

    mfplot_scatterplot, mass[aa], sfrm[aa], weight=weight[aa], $
      position=pos, xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, /sdss, ccolor=fsc_color('sky blue',101), $
      /nogrey, outcolor=fsc_color('sky blue',101), color=keycolor, $
      cthick=8
    mfplot_scatterplot, mass[qq], sfrm[qq], weight=weight[qq], $
      position=pos, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
      xtitle='', ytitle='', /sdss, ccolor=fsc_color('tomato'), $
      /overplot, /nogrey, outcolor=fsc_color('tomato'), color=keycolor, $
      cthick=8

    im_legend, ['Star-Forming','Quiescent'], /right, /top, box=0, $
      charsize=1.9, textcolor=['sky blue','tomato']
    
;   legend, ['Salim+07'], /right, /top, box=0, charsize=1.4, line=0, $
;     thick=6, pspacing=1.9, color=keycolor
;   djs_oplot, maxis, poly(maxis-10+0.25,[-9.83,-0.35])+9, line=0, thick=6 ; Salim+07, eq 11; Chabrier-->Salpeter

    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    
    
; --------------------------------------------------
; SDSS - 0.5(u-z) vs 0.5z
    zbins = mf_zbins(/sdss)
    if (n_elements(sdss_parent) eq 0L) then sdss_parent = read_mf_parent('sdss')
    mz = sdss_parent.k_ugrizjhk_absmag_05[4]
    umz = sdss_parent.k_ugrizjhk_absmag_05[0]-mz
    weight = sdss_parent.final_weight
    
    psfile = talkpath+'mz_vs_umz_sdss.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    xrange = [-18.3,-23.8]
    yrange = [0.3,7.2]
    xtitle = mfplot_mztitle()
    ytitle = textoidl('^{0.5}(u - z)')
    mzaxis = range(xrange[0]-0.3,xrange[1]+0.3,50)

    mfplot_scatterplot, mz, umz, weight=weight, position=pos, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, /sdss, color=keycolor;, $
;     ccolor=keycolor ; fsc_color('tan',101)
; overplot the axes in black to make the tick marks visible 
    djs_plot, [0], [0], /nodata, /noerase, position=pos, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    xyouts, -19.2, 6.5, 'Quiescent', align=0.5, charsize=1.8
    xyouts, -19.5, 0.8, 'Star-Forming', align=0.5, charsize=1.8
;   xyouts, -22.8, 0.8, 'Star-Forming', align=0.5, charsize=1.6

    junk = mf_select_quiescent(z=0.1,coeff=coeff,$
      mzpivot=mzpivot,wyder_coeff=wyder_coeff)
    djs_oplot, mzaxis, poly(mzaxis-mzpivot,wyder_coeff), $
      thick=6, color=fsc_color('firebrick',101) ; from Wyder+07
    djs_oplot, mzaxis, poly(mzaxis-mzpivot,coeff), line=5, $
      color='dark green', thick=6

; show a reddening vector
    uweff = k_lambda_eff(filterlist='sdss_u0.par',band_shift=0.5)
    zweff = k_lambda_eff(filterlist='sdss_z0.par',band_shift=0.5)

    ebv = 0.3
    mz_true = -23.3
    umz_true = 0.9

    mz_red = mz_true + 0.4*ebv*k_lambda(zweff,/char)
    umz_red = umz_true + 0.4*ebv*(k_lambda(uweff,/char)-k_lambda(zweff,/char))

    arrow, mz_true, umz_true, mz_red, umz_red, /data, $
      hsize=-0.25, hthick=8, thick=8 ;, color=djs_icolor('blue')
    xyouts, (mz_red-mz_true)/2.0+mz_true+0.3, (umz_red-umz_true)/2.0+umz_true-0.22, $
      orientation=-38, /data, 'E(B-V) = '+string(ebv,format='(F3.1)'), $
      charsize=1.3, align=0.5   ;, color=djs_icolor('blue')

    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

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
    
    all = {phistar: 4.26D-3, mstar: 10.648+0.25D, alpha: -0.46D, $
      phiplus: 0.58D-3, alphaplus: -1.58D}
    merge = {phistar: 1.6D-3, mstar: 11.1D, alpha: -0.4D}

    xrange = [9,12.4]
    gg = where((mass gt xrange[0]+0.05) and (mass lt xrange[1]-0.05))
    gmass = mass[gg]
    g2 = where((newmass_sf gt xrange[0]+0.1) and (newmass_sf lt xrange[1]-0.1))
    gnewmass_sf = newmass_sf[g2]

    psfile = talkpath+'toymodel.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, charsize=2.0, $
      height=5.8

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')
    
    plot, gmass, alog10(mf_schechter_plus(gmass,all)), xsty=3, $
      ysty=3, xrange=[9.0,12.5], yrange=[-8,-1], xtitle='Stellar Mass', $
      ytitle=textoidl('Number Density (Mpc^{-3} dex^{-1})'), $
      color=keycolor, thick=10, xtickinterval=1, line=0
    djs_oplot, gmass, alog10(mf_schechter(gmass,merge)), $
      color=fsc_color('orange',100), thick=12, line=3
    djs_oplot, gnewmass_sf, alog10(mf_schechter_plus(mass[g2],all)), $
      color=fsc_color('cyan',101), thick=12, line=5

    im_legend, ['Initial Mass Function','Mass-Dependent Star Formation','Simple Merger Model'], $
      /left, /bottom, box=0, textcolor=['white','cyan','orange'], charsize=1.5, $
      line=[0,5,3], thick=10, pspacing=2.0, color=['white','cyan','orange']
    
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

return
end


;; ---------------------------------------------------------------------------    
;; example of how the mass function might change with star formation
;; and mergers
;
;; incomplete piece of code    
;    
;    minmass = 8.8
;    maxmass = 12.5
;    binsize = 0.01
;    nbin = ceil((maxmass-minmass)/binsize)
;    maxis = range(minmass,maxmass,nbin)
;
;    schechter = {phistar: 7D, mstar: 10.648D, alpha: -0.46D, $
;      phiplus: 1D, alphaplus: -1.58D}
;    phimodel = mf_schechter_plus(maxis,schechter)
;;   schechter = {phistar: 4.26D-3, mstar: 10.648D, alpha: -0.46D, $
;;     phiplus: 0.58D-3, alphaplus: -1.58D}
;    phi = mf_schechter_plus_random(minmass,maxmass,nbin,$
;      schechter=schechter,mass=mass)
;    srt = sort(mass) & mass = mass[srt] & phi = phi[srt]
;    phi = phi/min(phi)
;    phierr = phi*0.0
;    for ii = 0, nbin-1 do phierr[ii] = 0.2*sqrt(djs_mean(im_poisson_limits(phi[ii])))
;;   phi = (phi+randomn(seed,nbin)*phierr)>0
;
;; pure star formation for 5 Gyr
;    time = 7D9
;    logsfr = 0.67*mass - 6.19      ; from Noeske+07
;    newmass_sf = alog10(10^mass+0.5*10^logsfr*time) ; account for recycling
;
;    djs_plot, mass, alog10(phi), xsty=3, ysty=3, xrange=[9,12.0], yrange=[0,3], psym=6
;    djs_oplot, maxis, alog10(phimodel*total(phi)/total(phimodel)/sqrt(2)), color='green'
;    djs_oplot, newmass_sf, alog10(phi), color='blue'
;    
;stop    
;    
;; now simulate mergers
;    all = lindgen(ngal)
;    nleft = ngal
;    for ii = 0, ngal/2-1 do begin ; let 50% merge
;       indx = all[long(randomu(seed,2)*nleft)]
;       remove, indx, all
;       nleft = n_elements(all)
;       mergemass = alog10(total(10^mass[indx]))
;       mergephi = interpolate(phi,findex(mass,mergemass))*2
;       if (ii eq 0) then begin
;          newmass_merge = mergemass
;          newphi_merge = mergephi
;       endif else begin
;          newmass_merge = [newmass_merge,mergemass]
;          newphi_merge = [newphi_merge,mergephi]
;       endelse
;    endfor
;;   newmass_merge = 
;
;stop    
;    
;    
;    djs_plot, mass, alog10(phi), xsty=3, ysty=3, xrange=[9,12.0]
;    djs_oplot, maxis, alog10(mf_schechter_plus(maxis,schechter)), color='green'
;    djs_oplot, newmass_sf, alog10(phi), color='blue'
;
