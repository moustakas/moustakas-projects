pro bootes_zeropoints
; jm10jan05ucsd - calibrate the BOOTES/NDWFS photometry
    
    common calibrate_bootes, bootes1, sdss1

    ndwfspath = getenv('RESEARCHPATH')+'/data/ndwfs/'
    bootespath = getenv('RESEARCHPATH')+'/data/bootes/'

    psfile = bootespath+'bootes_kurucz_zeropoints.ps'

    band = ['Bw','R','I']
    nband = n_elements(band)
    mrange = [[19.0,22.2],[18.5,21.2],[18.2,21.0]]
    rrange = 0.39*[-1,1]
    
; read and match the BOOTES and SDSS catalogs    
    if (n_elements(sdss1) eq 0) or (n_elements(bootes1) eq 0) then begin
       sdssfile = ndwfspath+'NDWFS_SDSS_stars.fits.gz'
       splog, 'Reading '+sdssfile
       allsdss = mrdfits(sdssfile,1)
       
       bootesfile = bootespath+'bootes_I.fits.gz'
       splog, 'Reading '+bootesfile
       allbootes = mrdfits(bootesfile,1)

; match twice, calculating the median offset in the coordinate systems       
       m1 = im_spherematch(allbootes,allsdss,match2=m2,$
         ratagname1='alpha_j2000',dectagname1='delta_j2000',$
         raoff=raoff,decoff=decoff)
       splog, 3600.0*raoff, 3600.0*decoff
       
       sdss1 = allsdss[m2]
       tags = tag_names(allbootes[0])
       bootes1 = im_struct_trimtags(allbootes[m1],$
         select=tags,newtags='I_'+tags)
       
; now read the Bw- and R-band catalogs and merge the three optical
; catalogs 
       bwband = mrdfits(bootespath+'bootes_Bw.fits.gz',1,rows=m1)
       tags = tag_names(bwband[0])
       bootes1 = struct_addtags(im_struct_trimtags(bwband,$
         select=tags,newtags='Bw_'+tags),temporary(bootes1))

       rband = mrdfits(bootespath+'bootes_R.fits.gz',1,rows=m1)
       tags = tag_names(rband[0])
       bootes1 = struct_addtags(im_struct_trimtags(rband,$
         select=tags,newtags='R_'+tags),temporary(bootes1))
    endif

; convert to maggies and select a fiducial sample of objects with good
; multiband photometry
    cut1 = where($
      (bootes1.i_flag_duplicate eq 0) and $
      (bootes1.i_flag_subfield eq 1) and $
      (bootes1.i_segflags_aper_04 eq 0) and $
      (total(sdss1.psfflux gt 0.0,1) eq 5.0))
    sdss_to_maggies, sdssmaggies, sdssivarmaggies, calib=sdss1[cut1], flux='psf'
    rmag = -2.5*alog10(sdssmaggies[2,*])
    cut2 = where((rmag gt 19.0) and (rmag lt 20.5),ncut2)

; fit Kurucz models to all the stars and compare the synthesized vs
; observed photometry
    splog, 'Fitting '+string(ncut2,format='(I0)')+' objects...'
    kall = fit_kurucz_models(sdssmaggies[*,cut2],sdssivarmaggies[*,cut2],$
      filterlist=sdss_filterlist())
    cut3 = where(kall.kurucz_chi2min lt 3.0,nstar)
    splog, 'N = ', nstar

    sdssmaggies = sdssmaggies[*,cut2[cut3]]
    sdssivarmaggies = sdssivarmaggies[*,cut2[cut3]]
    
    kthese = kall[cut3]
    sdss = sdss1[cut1[cut2[cut3]]]
    bootes = bootes1[cut1[cut2[cut3]]]

; just keep BwRI    
    bootes_to_maggies, bootes, maggies, ivarmaggies, $
      filterlist=synthmaggies_filterlist, /psf, /nozpoffset
    maggies = maggies[0:2,*]
    ivarmaggies = ivarmaggies[0:2,*]
    synthmaggies_filterlist = strtrim(synthmaggies_filterlist[0:2],2)

; synthesize magnitudes
    synthmaggies = fltarr(nband,nstar)
    for ii = 0L, nstar-1L do synthmaggies[*,ii] = reform(k_project_filters($
      k_lambda_to_edges(kthese[ii].lambda),kthese[ii].spec,$
      filterlist=synthmaggies_filterlist))

; make the plot    
    im_plotconfig, 6, pos1, psfile=psfile, xmargin=[1.4,0.3], $
      height=[4.5,3.0], width=6.8
    for ii = 0, nband-1 do begin
       good = where((maggies[ii,*] gt 0.0),ngood)
       xx = reform(-2.5*alog10(maggies[ii,good]))
       yy = reform(-2.5*alog10(synthmaggies[ii,good]))
       djs_plot, xx, yy, position=pos1[*,0], xsty=1, ysty=1, $
         xrange=mrange[*,ii], yrange=mrange[*,ii], $
         xtitle='', ytitle=band[ii]+' (SDSS synthesized, AB mag)', $
         xtickname=replicate(' ',10), psym=6, symsize=0.1
       djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
       im_legend, '<\Delta'+band[ii]+'> = '+$
         im_string_stats(yy-xx,sigrej=3.0,ndecimal=3), /left, /top, box=0
       djs_plot, xx, yy-xx, position=pos1[*,1], /noerase, xsty=1, ysty=1, $
         xrange=mrange[*,ii], yrange=rrange, psym=6, symsize=0.1, $
         xtitle=band[ii]+' (NDWFS/BOOTES, AB mag)', $
         ytitle='Residuals (AB mag)'
       djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
    endfor
    im_plotconfig, psfile=psfile, /psclose, /gzip
    spawn, 'rsync -auv '+psfile+'.gz ~/', /sh

; show the spectra themselves    
    psfile = bootespath+'bootes_kurucz_sedfits.ps'
    im_plotconfig, 8, pos1, psfile=psfile

    sdss_weff = k_lambda_eff(filterlist=sdss_filterlist())
    bootes_weff = k_lambda_eff(filterlist=synthmaggies_filterlist)

    light = 2.99792458D18       ; speed of light [A/s]
    xrange1 = [3050.0,9500.0]
    xtitle1 = 'Wavelength (\AA)'
    ytitle1 = 'm_{AB}'
    
    bootes_mab = maggies2mag(maggies,magerr=bootes_mab_err,$
      ivarmaggies=ivarmaggies)
    bootes_kurucz_mab = -2.5*alog10(synthmaggies)

    sdss_mab = maggies2mag(sdssmaggies,magerr=sdss_mab_err,$
      ivarmaggies=sdssivarmaggies)
    sdss_kurucz_mab = -2.5*alog10(kthese.kurucz_maggies)

    for ii = 0, 99 do begin
;   for ii = 0, nstar-1 do begin
       wave = kthese[ii].lambda
       flux = kthese[ii].spec*wave^2/light
       flux = -2.5*alog10(flux>1D-50)-48.6

       get_element, wave, xrange1, xx
       yrange = fltarr(2)
       yrange[0] = (max(sdss_kurucz_mab[*,ii])>max(flux[xx[0]:xx[1]]))*1.02
       yrange[1] = (min(sdss_kurucz_mab[*,ii])<min(flux[xx[0]:xx[1]]))*0.98

       djs_plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
         xsty=1, ysty=1, xtitle=xtitle1, ytitle=ytitle1
       djs_oplot, wave, flux, line=0, color='grey'

       djs_oplot, sdss_weff, sdss_kurucz_mab[*,ii], psym=symcat(6,thick=6), $
         symsize=2.5, color=''
       djs_oplot, bootes_weff, bootes_kurucz_mab[*,ii], psym=symcat(6,thick=6), $
         symsize=2.5, color='blue'

       oploterror, sdss_weff, sdss_mab[*,ii], sdss_mab_err[*,ii], psym=symcat(16), $
         symsize=2.0, color=djs_icolor('dark green'), $
         errcolor=djs_icolor('dark green'), errthick=!p.thick
       oploterror, bootes_weff, bootes_mab[*,ii], bootes_mab_err[*,ii], psym=symcat(16), $
         symsize=2.0, color=djs_icolor('red'), $
         errcolor=djs_icolor('red'), errthick=!p.thick

       im_legend, [$
         '[Fe/H]='+string(kthese[ii].kurucz_feh,format='(F4.1)'),$
         'T_{eff}='+string(kthese[ii].kurucz_teff,format='(I0)'),$
         'log (g)='+string(kthese[ii].kurucz_g,format='(F3.1)'),$
         '\chi^{2}_{\nu}='+string(kthese[ii].kurucz_chi2min/4.0,format='(F3.1)')], $
         /left, /top, box=0, charsize=1.4

       im_legend, ['Kurucz ugriz','SDSS ugriz','Kurucz BwRI','BOOTES BwRI'], $
         psym=[6,16,6,16], color=['','dark green','blue','red'], /right, $
         /bottom, box=0, charsize=1.4, symthick=3.0
    endfor
    
    im_plotconfig, psfile=psfile, /psclose, /gzip
    spawn, 'rsync -auv '+psfile+'.gz ~/', /sh

stop    

return
end
