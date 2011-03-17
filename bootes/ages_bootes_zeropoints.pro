pro ages_bootes_zeropoints
; jm10jan06ucsd - check the relative and absolute zeropoints of the
;   BOOTES photometry by using the SDSS photometry of the AGES galaxies

    common ages_bootes_zeropoints, kcorr1, sdss1, bootes1

    agespath = ages_path(/analysis)
    bootespath = getenv('RESEARCHPATH')+'/data/bootes/'

    psfile = bootespath+'bootes_ages_zeropoints.ps'
    kcorr_qafile = bootespath+'bootes_ages_zeropoints_sedfits.ps'

    vname = 'default.nolines'

    band = ['Bw','R','I']
    nband = n_elements(band)
    mrange = [[18.0,23.2],[17.0,21.3],[17.0,21.0]]
    rrange = 0.89*[-1,1]
    
; read Eisenstein's k-corrections (for the redshifts) and the
; SDSS/AGES catalog and define a clean sample; then match to the
; BOOTES catalogs
    if (n_elements(kcorr1) eq 0) or (n_elements(sdss1) eq 0) or $
      (n_elements(bootes1) eq 0) then begin
       splog, 'Reading '+agespath+'ages_bootes.fits.gz'
       bootes1 = mrdfits(agespath+'ages_bootes.fits.gz',1)
       
       splog, 'Reading '+agespath+'catalog.kcorr.v3.fits.gz'
       kcorr1 = mrdfits(agespath+'catalog.kcorr.v3.fits.gz',1)
       
       splog, 'Reading '+agespath+'ages.sdss.phot.dr72.fits.gz'
       sdss1 = mrdfits(agespath+'ages.sdss.phot.dr72.fits.gz',1)
    endif    

    ages_id = lindgen(n_elements(bootes1))
    
; convert to maggies and select a fiducial sample of objects with good
; multiband photometry
    cut1 = where($
      (bootes1.i_mag_auto gt 0.0) and (bootes1.i_mag_auto lt 90.0) and $
      (bootes1.bw_mag_aper_04 gt 0.0) and (bootes1.bw_mag_aper_04 lt 90.0) and $
      (bootes1.r_mag_aper_04 gt 0.0) and (bootes1.r_mag_aper_04 lt 90.0) and $
      (bootes1.i_flag_duplicate eq 0) and $
      (bootes1.i_flag_subfield eq 1) and $
      (bootes1.i_segflags_aper_04 eq 0) and $
      (bootes1.bw_segflags_aper_04 eq 0) and $
      (bootes1.r_segflags_aper_04 eq 0) and $
      (bootes1.i_imaflags_aper_04 gt 5) and $
      (bootes1.bw_imaflags_aper_04 gt 5) and $
      (bootes1.r_imaflags_aper_04 gt 5) and $
      (kcorr1.z gt 0.05) and (kcorr1.z lt 0.3) and $
      (total(sdss1.petroflux gt 0.0,1) eq 5.0))
    sdss_to_maggies, sdssmaggies, sdssivarmaggies, $
      calib=sdss1[cut1], flux='petro'
    rmag = -2.5*alog10(sdssmaggies[2,*])
    cut2 = where((rmag gt 18.0) and (rmag lt 20.5),ncut2)
    sdssmaggies = sdssmaggies[*,cut2]
    sdssivarmaggies = sdssivarmaggies[*,cut2]
    
    ages_id = ages_id[cut1[cut2]]
    bootes = bootes1[cut1[cut2]]
    sdss = sdss1[cut1[cut2]]
    kcorr = kcorr1[cut1[cut2]]
    ngal = n_elements(bootes)
    splog, 'N = ', ngal

; fit K-correct models and compare the synthesized vs observed
; photometry
    synthmaggies_filterlist = (bootes_filterlist())[0:2]

    splog, 'Fitting '+string(ncut2,format='(I0)')+' objects...'
    kall = im_kcorrect(kcorr.z,sdssmaggies,sdssivarmaggies,$
      sdss_filterlist(),synthmaggies_filterlist,chi2=chi2,$
      coeffs=coeffs,bestmaggies=bestmaggies,mass=mass,$
      band_shift=0.0,/silent,synth_outmaggies_obs=synthmaggies,$
      vname=vname)
    
; just keep BwRI    
    bootes_to_maggies, bootes, maggies, ivarmaggies, /nozpoffset
    maggies = maggies[0:2,*]
    ivarmaggies = ivarmaggies[0:2,*]

; force the median offset in the I-band to be the same in order to
; take out differences in how total magnitudes are measured (e.g.,
; Petrosian vs Kron)...
    
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
         xtitle=band[ii]+' (NDWFS/Bootes, AB mag)', $
         ytitle='Residuals (AB mag)'
       djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
    endfor
    im_plotconfig, psfile=psfile, /psclose, /gzip
    spawn, 'rsync -auv '+psfile+'.gz ~/', /sh

; make a K-correct QAplot    
    str = replicate({galaxy: '', z: 0.0, chi2: 0.0, mass: 0.0, $
      coeffs: fltarr(5), in_filterlist: strarr(5), maggies: fltarr(5), $
      ivarmaggies: fltarr(5), bestmaggies: fltarr(5), $
      out_filterlist: strarr(3), synth_outmaggies_obs: fltarr(3), $
      outmaggies_obs: fltarr(3), outivarmaggies_obs: fltarr(3)},ngal)
    str.galaxy = 'AGES/'+string(ages_id,format='(I5.5)')
    str.z = kcorr.z
    str.chi2 = chi2
    str.mass = alog10(mass)
    str.coeffs = coeffs
    str.maggies = sdssmaggies
    str.ivarmaggies = sdssivarmaggies
    str.bestmaggies = bestmaggies
    str.in_filterlist = sdss_filterlist()
    str.out_filterlist = synthmaggies_filterlist

    str.outmaggies_obs = maggies
    str.outivarmaggies_obs = ivarmaggies
    str.synth_outmaggies_obs = synthmaggies

    kcorrect_qaplot, str[0:100], vname=vname, psfile=kcorr_qafile, /clobber
    spawn, 'rsync -auv '+kcorr_qafile+'.gz ~/', /sh
    
stop    

return
end
