pro pie_allsurveys, pos, keynote=keynote, keycolor=keycolor, $
  zrange=zrange, rot=rot, thaxes=thaxes

; SDSS    
    zz = read_vagc_garching(/postlss)
    ww = where(zz.dec gt -2 and zz.dec lt 2)
    pie_plot, zz[ww].z, zz[ww].ra, psym=3, symsize=0.1, rlabel='Redshift', $
      rrange=zrange, thaxes=thaxes, /hours, axiscolor=keycolor, $
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


    
; ---------------------------------------------------------------------------
; sdss - compare all, quiescent, and active samples
    psfile = talkpath+'mf_sdss.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.8, height=6.0, charsize=2.0, keynote=keynote

    xrange = [8.8,12.2]
    yrange = [-7,-1.5]
    maxis1 = range(xrange[0]+0.15,12.5,100)

    subsample = ['active','quiescent']
    mfcolor = ['powder blue','tomato']

    psymsize = [0.9,1.1,0.9]*1.2
    mfpsym = [16,14,15]

    plot, [0], [0], /nodata, position=pos, xsty=9, ysty=9, $
      yrange=yrange, xrange=xrange, xtitle='log (Stellar Mass) (M'+sunsymbol()+')', $
      ytitle=textoidl('log (Number Density) (Mpc^{-3} dex^{-1})'), xtickinterval=1, color=keycolor
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
; SDSS pie diagram
    psfile = talkpath+'pie_sdss.ps'
    im_plotconfig, 0, pos, psfile=psfile, width=7.0, $
      height=7.0, charsize=2, keynote=keynote

    zz = read_vagc_garching(/postlss)
    ww = where(zz.dec gt -2 and zz.dec lt 2)
    pie_plot, zz[ww].z, zz[ww].ra, psym=3, symsize=0.1, rlabel='Redshift', $
      rrange=[0,0.17], thaxes=90.0, /hours, axiscolor=keycolor, $
      color=keycolor, position=pos, rotate=270

    im_plotconfig, psfile=psfile, /psclose, keynote=keynote

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
