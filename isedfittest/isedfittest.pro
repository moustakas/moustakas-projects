pro test_oplot_reddening, ubvri=ubvri

    ebv = 0.5

    if keyword_set(ubvri) then begin
       bmv_true = 1.2
       umb_true = 0.1
       uweff = k_lambda_eff(filterlist='bessell_U.par',band_shift=0.0)
       bweff = k_lambda_eff(filterlist='bessell_B.par',band_shift=0.0)
       vweff = k_lambda_eff(filterlist='bessell_R.par',band_shift=0.0)
       
       umb_red = umb_true + 0.4*ebv*(k_lambda(uweff,/calz)-k_lambda(bweff,/calz))
       bmv_red = bmv_true + 0.4*ebv*(k_lambda(bweff,/calz)-k_lambda(vweff,/calz))
       
       arrow, bmv_true, umb_true, bmv_red, umb_red, /data, $
         hsize=-0.25, hthick=8, thick=8
       xyouts, bmv_red, umb_red+0.15, /data, 'E(B-V) = '+$
         string(ebv,format='(F3.1)'), charsize=1.3, align=0.5

    endif else begin

       gmr_true = 1.2
       umg_true = 0.1
       uweff = k_lambda_eff(filterlist='sdss_u0.par',band_shift=0.1)
       gweff = k_lambda_eff(filterlist='sdss_g0.par',band_shift=0.1)
       rweff = k_lambda_eff(filterlist='sdss_r0.par',band_shift=0.1)

       umg_red = umg_true + 0.4*ebv*(k_lambda(uweff,/calz)-k_lambda(gweff,/calz))
       gmr_red = gmr_true + 0.4*ebv*(k_lambda(gweff,/calz)-k_lambda(rweff,/calz))

       arrow, gmr_true, umg_true, gmr_red, umg_red, /data, $
         hsize=-0.25, hthick=8, thick=8
       xyouts, gmr_red, umg_red+0.15, /data, 'E(B-V) = '+$
         string(ebv,format='(F3.1)'), charsize=1.3, align=0.5

    endelse
   
return
end

function test_get_modelcolors, col, ubvri=ubvri

    nage = n_elements(col[0].age)
    nmodel = n_elements(col)
    bigtau = reform(rebin(reform(col.tau,1,nmodel),nage,nmodel),nage*nmodel)
    bigebv = reform(rebin(reform(col.ebv,1,nmodel),nage,nmodel),nage*nmodel)
    bigZ = reform(rebin(reform(col.Z,1,nmodel),nage,nmodel),nage*nmodel)
    bigage = reform(col.age,nage*nmodel)

    if keyword_set(ubvri) then begin
       ub = reform(reform(col.abmag[0,*]-col.abmag[1,*]),nage*nmodel)
       bv = reform(reform(col.abmag[1,*]-col.abmag[2,*]),nage*nmodel)
       vr = reform(reform(col.abmag[2,*]-col.abmag[3,*]),nage*nmodel)
       ri = reform(reform(col.abmag[3,*]-col.abmag[4,*]),nage*nmodel)
       modelcolors = {ub: ub, bv: bv, vr: vr, ri: ri, $
         tau: bigtau, ebv: bigebv, age: bigage, z: bigz}
    endif else begin
       ug = reform(reform(col.abmag[0,*]-col.abmag[1,*]),nage*nmodel)
       gr = reform(reform(col.abmag[1,*]-col.abmag[2,*]),nage*nmodel)
       ri = reform(reform(col.abmag[2,*]-col.abmag[3,*]),nage*nmodel)
       iz = reform(reform(col.abmag[3,*]-col.abmag[4,*]),nage*nmodel)
       modelcolors = {ug: ug, gr: gr, ri: ri, iz: iz, $
         tau: bigtau, ebv: bigebv, age: bigage, z: bigz}
    endelse
       
return, modelcolors
end    

pro isedfittest, build_parent=build_parent, models=models, isedfit=isedfit, $
  get_colors=get_colors, kcorrect=kcorrect, qaplot=qaplot, clobber=clobber, $
  compare_with_kcorrect=compare_with_kcorrect

    iopath = getenv('IM_ARCHIVE_DIR')+'/projects/isedfittest/'
    paramfile = iopath+'test_isedfit.par'
    params = read_isedfit_paramfile(paramfile)

; build a test sample    
    if keyword_set(build_parent) then begin
       vagc = mz_get_vagc(sample=sample,letter=letter,poststr=poststr)
       post = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/postlss)
       kauff = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/mpamass)

       keep = where((post.z gt 0.05) and (post.z lt 0.1) and (kauff.mass_median gt 0.0),nkeep)
       ngal = 5000
       these = random_indices(nkeep,ngal)

       keep = keep[these]
       post = post[keep]
       kauff = kauff[keep]
       
       vagcpath = getenv('VAGC_REDUX')+'/'
       sdssphot = mrdfits(vagcpath+'object_sdss_imaging.fits.gz',$
         row=post.object_position,1)
       galex = mrdfits(vagcpath+'object_galex_gr6.fits.gz',$
         row=post.object_position,1)
;      if (n_elements(sdsstwomass) eq 0L) then sdsstwomass = $
;        mrdfits(vagcpath+'object_twomass.fits.gz',$
;        row=phot.object_position,1)

       sdss_to_maggies, maggies, ivarmaggies, calib=sdssphot, flux='cmodel'
       im_galex_to_maggies, galex, gmaggies, givarmaggies
       
       parent = struct_addtags(sdssphot,struct_trimtags(post,$
         select=['z','absm','object_position']))
       parent = struct_addtags(temporary(parent),replicate({maggies: fltarr(7), $
         ivarmaggies: fltarr(7)},ngal))
       parent.maggies = [gmaggies,maggies]
       parent.ivarmaggies = [givarmaggies,ivarmaggies]

       im_mwrfits, parent, 'parent.fits', /clobber
       im_mwrfits, kauff, 'kauffmann.fits', /clobber
    endif

    parent = mrdfits(iopath+'parent.fits.gz',1)
    kcorr = mrdfits(iopath+'kcorr.fits.gz',1)
    kauff = mrdfits(iopath+'kauffmann.fits.gz',1)
    ngal = n_elements(parent)

;   kk = sdss_kcorrect(parent.z,nmgy=parent.maggies*1D9,$
;     ivar=parent.ivarmaggies/1D18,band_shift=0.1,chi2=chi2,$
;     mass=mass,absmag=absmag,rmaggies=rmaggies)
    
; --------------------------------------------------
; compare the data and the models in color-color space
    if keyword_set(get_colors) then begin

; test with an SSP       
       test = mrdfits(getenv('ISEDFIT_SFHGRID_DIR')+'/ssp/bc03/bc03_chab_Z0.05.fits.gz',1)
       test.flux = test.flux/rebin(reform(test.mstar,1,220),6900,220)
       testmaggies = fltarr(5,n_elements(test.age))
       for jj = 0, n_elements(test.age)-1 do testmaggies[*,jj] = $
         reform(k_project_filters(k_lambda_to_edges(test.wave),$
         test.flux[*,jj],filterlist=sdss_filterlist(),band_shift=0.1))
       ug_test = -2.5*alog10(testmaggies[0,*]/testmaggies[1,*])
       gr_test = -2.5*alog10(testmaggies[1,*]/testmaggies[2,*])

; define the grids       
;      sfhgrid = 1
;      title = 'Dust, Z=0.004-0.05, No bursts'

       sfhgrid = [1,2,3,4]
       title = [$
         'Dust, Z=0.004-0.05, No bursts',$
         'Dust, Z=0.02, No bursts',$
         'No dust, Z=0.004-0.05, No bursts',$
         'Dust, Z=0.004-0.05, Bursts']

;       sfhgrid = [1,2,4,5,6,7]
;       title = [$
;         'Dust, Z=0.004-0.05, No bursts',$
;         'Dust, Z=0.02, No bursts',$
;;        'No dust, Z=0.004-0.05, No bursts',$
;         'Dust, Z=0.004-0.05, Bursts',$
;         'Dust, Z=0.004-0.05, Weaker Bursts',$
;         'Dust, Z=0.004-0.05, Truncated Bursts',$
;         'More Dust, Z=0.004-0.05, No bursts']

; pack in the observed colors
;      kk = kcorr
       kk = kcorr[where(kcorr.k_chi2 lt 10.0)]
;      data = {$
;        ub: reform(kk.k_ubvri_absmag_00[0]-kk.k_ubvri_absmag_00[1]), $
;        bv: reform(kk.k_ubvri_absmag_00[1]-kk.k_ubvri_absmag_00[2]), $
;        vr: reform(kk.k_ubvri_absmag_00[2]-kk.k_ubvri_absmag_00[3]), $
;        ri: reform(kk.k_ubvri_absmag_00[3]-kk.k_ubvri_absmag_00[4])}

       data = {$
         ug: reform(kk.k_ugriz_absmag_01[0]-kk.k_ugriz_absmag_01[1]), $
         gr: reform(kk.k_ugriz_absmag_01[1]-kk.k_ugriz_absmag_01[2]), $
         ri: reform(kk.k_ugriz_absmag_01[2]-kk.k_ugriz_absmag_01[3]), $
         iz: reform(kk.k_ugriz_absmag_01[3]-kk.k_ugriz_absmag_01[4])}

; for the plot
;      xcolor = ['bv','vr','ri']
;      ycolor = ['ub','bv','vr']
;      xrange = [[-0.5,2.0],[-0.5,1.1],[-0.5,1.0]]
;      yrange = [[-0.3,3.1],[-0.5,2.0],[-0.5,1.0]]
;      xtitle = '^{0.0}('+['B - V','V - R','R - I']+')'
;      ytitle = '^{0.0}('+['U - B','B - V','V - R']+')'

      xcolor = ['gr','ri','iz']
      ycolor = ['ug','gr','ri']
      xrange = [[-0.2,1.7],[-0.1,1.0],[-0.1,0.9]]
      yrange = [[-0.3,2.6],[-0.2,1.7],[-0.1,1.0]]
      xtitle = '^{0.1}('+['g - r','r - i','i - z']+')'
      ytitle = '^{0.1}('+['u - g','g - r','r - i']+')'
       
       for ii = 0, n_elements(sfhgrid)-1 do begin
          sfhg = string(sfhgrid[ii],format='(I2.2)')
          suffix = 'ugriz'
;         suffix = 'ubvri'
          colorfile = iopath+'chab_sfhgrid'+sfhg+'_colors_'+suffix+'.fits.gz'
          if (file_test(colorfile) eq 0) then begin
;            colors = get_sfhgrid_colors(sfhgrid[ii],$
;              filterlist=bessell_filterlist(),band_shift=0.0)
             colors = get_sfhgrid_colors(sfhgrid[ii],$
               filterlist=sdss_filterlist(),band_shift=0.1)
             im_mwrfits, colors, repstr(colorfile,'.gz',''), /clobber
          endif else colors = mrdfits(colorfile,1)

          model = test_get_modelcolors(colors)

; make the plot for every color-color combination
          psfile = repstr(colorfile,'.fits.gz','.ps')
          im_plotconfig, 0, pos, psfile=psfile, ymargin=[0.6,1.1], $
            xmargin=[1.3,0.4], width=6.8, height=6.8

          for jj = 0, n_elements(xcolor)-1 do begin
             dxtag = tag_indx(data,xcolor[jj]) ; data
             dytag = tag_indx(data,ycolor[jj])
             mxtag = tag_indx(model,xcolor[jj]) ; model
             mytag = tag_indx(model,ycolor[jj])
             
             hogg_scatterplot, data.(dxtag), data.(dytag), position=pos, xsty=1, ysty=1, $
               xtitle=textoidl(xtitle[jj]), ytitle=textoidl(ytitle[jj]), $
               xrange=xrange[*,jj], yrange=yrange[*,jj], levels=[0.1,0.25,0.5,0.75,0.9,0.95], $
               /outliers, outpsym=6, outcolor=djs_icolor('black'), /internal, $
               title='('+string(sfhgrid[ii],format='(I2.2)')+') '+title[ii]
;            test_oplot_reddening

             hogg_scatterplot, model.(mxtag), model.(mytag), position=pos, /overplot, xsty=9, ysty=9, $
               xrange=xrange[*,jj], yrange=yrange[*,jj], levels=[0.1,0.25,0.5,0.75,0.9,0.95], $
               ccolor=djs_icolor('red'), /nogrey, /outliers, outcolor=djs_icolor('red'), /internal
             ww = where(model.age gt 12.0 and model.tau lt 1 and model.ebv lt 0.1,nww)
;            ww = where(model.ebv gt 0.7,nww)
             if (nww ne 0L) then djs_oplot, (model.(mxtag))[ww], $
               (model.(mytag))[ww], psym=symcat(16), sym=0.8, color='blue'
;            if (jj eq 0) then begin
;               rc3 = read_rc3()
;               djs_oplot, rc3.b_vt, rc3.u_bt, psym=6, color='blue'
;            endif
             if (jj eq 0) then djs_oplot, gr_test, ug_test, line=0, thick=6
          endfor
          im_plotconfig, psfile=psfile, /psclose, /gzip
;         spawn, 'rsync -auv '+psfile+'.gz ~/', /sh
       endfor 
    endif

; --------------------------------------------------
; build the models
    if keyword_set(models) then isedfit_models, $
      paramfile, iopath=iopath, clobber=clobber

; --------------------------------------------------
; do the fitting!  
;    jj = mrdfits(iopath+'test1_bc03_chab_sfhgrid03.fits.gz',1)
;;   jj = mrdfits(iopath+'test1_bc03_chab_calzetti_sfhgrid02.fits.gz',1)
;    kk = mrdfits(iopath+'kcorr.fits.gz',1)
;    index = where(abs(jj.mass_50-kk.k_mass) gt 0.3)
    if keyword_set(isedfit) then begin
;      index = [690,1202,1790,2122,2866,3048,3502,4097,4147,4243,4766,4824]
;      index = lindgen(50)
;      outprefix = 'test2'
       isedfit, paramfile, parent.maggies, parent.ivarmaggies, parent.z, $
         result, iopath=iopath, outprefix=outprefix, galchunksize=5000, $
         clobber=clobber, index=index
    endif 

;   mz_kcorrect_qaplot, kk[index], psfile='junk.ps', in_filterlist=sdss_filterlist()

; --------------------------------------------------
    if keyword_set(kcorrect) then begin
       kcorr = test_isedfit_kcorrect(parent.z,parent.maggies,$
         parent.ivarmaggies,filterlist=[galex_filterlist(),sdss_filterlist()])
       im_mwrfits, kcorr, 'kcorr.fits', /clobber
    endif

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
       jj = mrdfits(iopath+'test1_bc03_chab_charlot_sfhgrid08.fits.gz',1)
;      kk = mrdfits(iopath+'kcorr.fits.gz',1)
       index = random_indices(n_elements(jj),50)
;      index = where(abs(jj.mass_50-kk.k_mass) gt 0.3)
       isedfit_qaplot, paramfile, isedfit, iopath=iopath, galaxy=galaxy, $
         index=index, clobber=clobber, outprefix=outprefix
    endif

; --------------------------------------------------
    if keyword_set(compare_with_kcorrect) then begin
       residrange = 0.9*[-1,1]

;      sfhgrid = 5
;      title = 'Dust, Z=0.004-0.05, Weaker bursts'
       sfhgrid = 10
       title = 'Dust, Z=0.004-0.05, Bursts'

;      sfhgrid = [1,2,3,4]
;      title = [$
;        'Dust, Z=0.004-0.05, No bursts',$
;        'Dust, Z=0.02, No bursts',$
;        'No dust, Z=0.004-0.05, No bursts',$
;        'Dust, Z=0.004-0.05, Bursts']
      
;       sfhgrid = [1,2,4,5,6,7]
;       title = [$
;         'Dust, Z=0.004-0.05, No bursts',$
;         'Dust, Z=0.02, No bursts',$
;;        'No dust, Z=0.004-0.05, No bursts',$
;         'Dust, Z=0.004-0.05, Bursts',$
;         'Dust, Z=0.004-0.05, Weaker Bursts',$
;         'Dust, Z=0.004-0.05, Truncated Bursts',$
;         'More Dust, Z=0.004-0.05, No bursts']
      
       psfile = iopath+'qaplot_ised_vs_kcorr.ps'
       im_plotconfig, 10, pos, psfile=psfile, charsize=1.6, $
         yspace=[1.0,1.0], xmargin=[1.2,0.4], ymargin=[0.6,1.1];, $
;        height=[3.0,3.0]

       kauffmass = kauff.mass_median-0.04 ; Kroupa-->Chabrier
       kcorrmass = kcorr.k_mass ; Chabrier

; ##################################################
; on page 1 compare the K-correct and Kauffmann stellar masses
; mass
       resid = kauffmass-kcorrmass
       hogg_scatterplot, kauffmass, resid, position=pos[*,0], xsty=1, ysty=1, $
         xrange=[9,11.8], yrange=residrange, xtitle=textoidl('log (M/M_{\odot})'), $
         ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
         levels=[0.1,0.25,0.5,0.75,0.9], /internal
       djs_oplot, !x.crange, [0,0], line=0, color='red'
       im_legend, '\Delta = '+im_string_stats(resid,ndec=3), $
         /left, /top, box=0, charsize=1.4
; color
       hogg_scatterplot, kcorr.k_ugriz_absmag_01[0]-kcorr.k_ugriz_absmag_01[1], resid, $
         /noerase, position=pos[*,1], xsty=1, ysty=1, $
         xrange=[0.2,2.5], yrange=residrange, xtitle=textoidl('^{0.1}(u - g)'), $
         ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
         levels=[0.1,0.25,0.5,0.75,0.9], /internal
       djs_oplot, !x.crange, [0,0], line=0, color='red'
       xyouts, 0.03, 0.5, 'Mass Residuals [Kauffmann-Kcorrect] (dex)', $
         /normal, orientation=90, align=0.5
       !p.multi = 0
       
; now compare K-correct and Kauffmann separately with K-correct       
       for ii = 0, n_elements(sfhgrid)-1 do begin
          sfhg = 'sfhgrid'+string(sfhgrid[ii],format='(I2.2)')
          if (sfhgrid[ii] eq 3) then red = '' else red = '_charlot'
;         if (sfhgrid[ii] eq 3) then red = '' else red = '_calzetti'
          sfhgridfile = iopath+'test3_bc03_chab'+red+'_'+sfhg+'.fits.gz'
;         sfhgridfile = iopath+'test1_bc03_chab'+red+'_'+sfhg+'.fits.gz'
          splog, 'Reading '+sfhgridfile
          ised = mrdfits(sfhgridfile,1)
          isedmass = ised.mass_50 ; Chabrier
; ##################################################
; comparison with K-correct          
; mass
          resid = isedmass-kcorrmass
          hogg_scatterplot, isedmass, resid, position=pos[*,0], xsty=1, ysty=1, $
            xrange=[9,11.8], yrange=residrange, xtitle=textoidl('log (M/M_{\odot})'), $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
          im_legend, '\Delta = '+im_string_stats(resid,ndec=3), $
            /left, /top, box=0, charsize=1.4
; color
          hogg_scatterplot, kcorr.k_ugriz_absmag_01[0]-kcorr.k_ugriz_absmag_01[1], resid, $
            /noerase, position=pos[*,1], xsty=1, ysty=1, $
            xrange=[0.2,2.5], yrange=residrange, xtitle=textoidl('^{0.1}(u - g)'), $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; AGE
          hogg_scatterplot, alog10(ised.age_50), resid, /noerase, position=pos[*,2], xsty=1, ysty=1, $
            xrange=alog10([1,20]), yrange=residrange, xtitle='log (Age) (Gyr)', $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; tau
          hogg_scatterplot, alog10(ised.tau_50>0.4), resid, /noerase, position=pos[*,3], xsty=3, ysty=1, $
            xrange=alog10([0.3,20]), yrange=residrange, xtitle=textoidl('\tau (Gyr)'), $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; dust
          hogg_scatterplot, alog10(ised.ebv_50>0.01), resid, /noerase, position=pos[*,4], xsty=3, ysty=1, $
            xrange=alog10([0.01,1]), yrange=residrange, xtitle='E(B-V) (mag)', $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; metallicity
          hogg_scatterplot, ised.Z_50, resid, /noerase, position=pos[*,5], xsty=3, ysty=1, $
            xrange=[0.003,0.052], yrange=residrange, xtitle='Metallicity', $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; titles
          xyouts, 0.03, 0.5, 'Mass Residuals [iSEDfit-Kcorrect] (dex)', $
            /normal, orientation=90, align=0.5
          xyouts, 0.5, 0.96, '('+string(sfhgrid[ii],format='(I2.2)')+') '+$
            title[ii], /normal, align=0.5
; ##################################################
; comparison with Kauffmann+
; mass
          resid = isedmass-kauffmass
          hogg_scatterplot, isedmass, resid, position=pos[*,0], xsty=1, ysty=1, $
            xrange=[9,11.8], yrange=residrange, xtitle=textoidl('log (M/M_{\odot})'), $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
          im_legend, '\Delta = '+im_string_stats(resid,ndec=3), $
            /left, /top, box=0, charsize=1.4
; color
          hogg_scatterplot, kcorr.k_ugriz_absmag_01[0]-kcorr.k_ugriz_absmag_01[1], resid, $
            /noerase, position=pos[*,1], xsty=1, ysty=1, $
            xrange=[0.2,2.7], yrange=residrange, xtitle=textoidl('^{0.1}(u - g)'), $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; AGE
          hogg_scatterplot, alog10(ised.age_50), resid, /noerase, position=pos[*,2], xsty=1, ysty=1, $
            xrange=alog10([1,20]), yrange=residrange, xtitle='log (Age) (Gyr)', $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; tau
          hogg_scatterplot, alog10(ised.tau_50>0.4), resid, /noerase, position=pos[*,3], xsty=3, ysty=1, $
            xrange=alog10([0.3,20]), yrange=residrange, xtitle=textoidl('\tau (Gyr)'), $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; dust
          hogg_scatterplot, alog10(ised.ebv_50>0.01), resid, /noerase, position=pos[*,4], xsty=3, ysty=1, $
            xrange=alog10([0.01,1]), yrange=residrange, xtitle='E(B-V) (mag)', $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; metallicity
          hogg_scatterplot, ised.Z_50, resid, /noerase, position=pos[*,5], xsty=3, ysty=1, $
            xrange=[0.003,0.052], yrange=residrange, xtitle='Metallicity', $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; titles
          xyouts, 0.03, 0.5, 'Mass Residuals [iSEDfit-Kauffmann] (dex)', $
            /normal, orientation=90, align=0.5
          xyouts, 0.5, 0.96, '('+string(sfhgrid[ii],format='(I2.2)')+') '+$
            title[ii], /normal, align=0.5
       endfor
       im_plotconfig, psfile=psfile, /psclose, /gzip ; /pdf
;      spawn, 'rsync -auv '+psfile+'.gz ~/', /sh
    endif
    stop
return
end
