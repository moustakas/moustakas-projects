function get_emlines, ewoii, cflux=cflux, cwave=cwave, z=z, $
  newratio=newratio
    light = im_light(/kms)

    coeff0 = alog10(3600D)       ; starting rest wavelength [log-10 Angstrom]
    coeff1 = 1D-4                ; [dispersion in log-10 Angstroms]
    velsz = (10D^coeff1-1)*light ; [km/s]

    maxwave = 6700D
    nwavepix = long((alog10(maxwave)-coeff0)/coeff1+1)
    logwave = coeff0 + dindgen(nwavepix)*coeff1

; specify the lines we care about; for [OII] use a 1:0.71 flux ratio
; (f_3726 = 0.71 f_3728.8) following Weiner, and for the other lines
; adopt Mostek's values
    restwave = [3726.032D,3728.815D,4861.33D,4958.91D,5006.84D,6562.8D]
    lineratio = [0.415,0.585,0.4577,0.094,0.2821,2.002]
    if keyword_set(newratio) then lineratio = [0.415,0.585,0.9,0.04,0.12,1.75]
    linewave = alog10(restwave)
    nline = n_elements(linewave)

; specify the rest-frame emission-line EW (this will eventually come
; from a Monte Carlo simulation); then compute the total (rest-frame)
; [OII] flux
    oiiwave = djs_mean([3726.032,3728.815])
    oiiflux = ewoii*cflux ; [erg/s/cm2] 

; add the emission lines (see Schlegel's linebackfit for the
; conversions to and from log-10 A)
    zshift = 0.0
    sigma = 100.0/light/alog(10) ; line-width [log-10 Angstrom]
    oiiamplitude = oiiflux/alog(10)/oiiwave/(sqrt(2.0*!pi)*sigma)
    lineflux = logwave*0
    for jj = 0, nline-1 do begin
; emission-line amplitude          
       amplitude = lineratio[jj]*oiiamplitude
       linewave1 = linewave[jj]+alog10(1+zshift)
       lineflux += amplitude*exp(-0.5*(logwave-linewave1)^2/sigma^2)
    endfor

    lineflux = interpolate(lineflux/(1.0+z),findex((10^logwave)*(1+z),cwave),missing=0.0)
    
return, lineflux
end

pro talk_13jul_desi
; jm13jul14siena - make some plots for the talk

    path = getenv('IM_RESEARCH_DIR')+'/talks/2013/13jul_desi/'

; ---------------------------------------------------------------------------
; plot the EW([OII]) distribution from DEEP2

    zz = read_deep2(/kcorr)
    ss = read_deep2(/ispec)
    ww = where(zz.zbest gt 0.7 and zz.zbest lt 1.4 and $
      ss.oii_3727_1[0]/ss.oii_3727_1[1] gt 5.0 and $
      ss.oii_3727_2[0]/ss.oii_3727_2[1] gt 5.0 and $
      ss.oii_3727_1_ew[0] gt 0.0 and $
      ss.oii_3727_2_ew[0] gt 0.0,ngal)
    zz = zz[ww]
    ss = ss[ww]

    umg = zz.ugriz_absmag_01[0]-zz.ugriz_absmag_01[1]
    ewoii = alog10(ss.oii_3727_1_ew[0] + ss.oii_3727_2_ew[0])

    psfile = path+'ewoii_hist.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, charthick=4, $
      width=5.5, margin=[1.6,0.2], charsize=2.0, xpage=7.3
    djs_plot, [0], [0], /nodata, position=pos, xsty=5, ysty=5, $
      xrange=[0.5,3.0], yrange=[0,700]
    im_plothist, ewoii, /fill, bin=0.03, fcolor=im_color('tan'), /overplot
    djs_plot, [0], [0], /nodata, /noerase, position=pos, xsty=1, ysty=1, $
      xrange=[0.5,3.0], yrange=[0,700], xtitle='log EW([O II]) (\AA)', $
      ytitle='Number of Galaxies'
;   im_legend, 'DEEP2 at z~1', /left, /top, margin=0
;   im_legend, 'DEEP2 at z~1', /left, /top, margin=0
    im_plotconfig, psfile=psfile, /psclose, /pdf

; now show a couple example model fits
    light = 2.99792458D18       ; speed of light [A/s]
    filt = deep2_filterlist()
    weff = k_lambda_eff(filterlist=filt)

    get_element, 10^ewoii, [30.0,220.0], ww
;   ww = (where(ewoii gt 1.8))[0:1]
;   ww = [500,850]
    splog, 10^ewoii[ww], zz[ww].zbest
    mag = maggies2mag(zz[ww].maggies,ivarmaggies=zz[ww].ivarmaggies,magerr=magerr)

    vname = 'default.nolines'
    k_load_vmatrix, vmatrix, lambda, vname=vname
    wave = k_lambda_to_centers(lambda)

    psfile = path+'example_sedfits.ps'
    im_plotconfig, 6, pos, psfile=psfile, ymargin=[1.1,1.1], charsize=2.0

; first object    
    cwave = wave*(1+zz[ww[0]].zbest)
    cflux = vmatrix#zz[ww[0]].coeffs/(1+zz[ww[0]].zbest)
    lineflux = get_emlines(10^ewoii[ww[0]],cflux=zz[ww[0]].cflux_3727,$
      cwave=cwave,z=zz[ww[0]].zbest)
    cflux = cflux + lineflux
    cflux = cflux*cwave^2/light           ; observed frame
    cflux = -2.5*alog10(cflux>1D-50)-48.6 ; [AB mag]

    djs_plot, [0], [0], /nodata, position=pos[*,0], $
      xsty=1, ysty=1, xrange=[0.15,1.25], yrange=[25.5,20], $
      ytitle='AB mag', xtickname=replicate(' ',10)
    djs_oplot, cwave/1D4, cflux, line=0, color='grey'
    oploterror, weff/1D4, mag[*,0], magerr[*,0], psym=symcat(15), $
      symsize=2.0, color=im_color('firebrick')
    im_legend, ['EW([O II]) = 30 \AA'], /left, /top, box=0, margin=0, charsize=1.8
    
; second object
    cwave = wave*(1+zz[ww[1]].zbest)
    cflux = vmatrix#zz[ww[1]].coeffs/(1+zz[ww[1]].zbest)
    lineflux = get_emlines(10^ewoii[ww[1]],cflux=zz[ww[1]].cflux_3727,$
      cwave=cwave,z=zz[ww[1]].zbest,/newratio)
    cflux = cflux + lineflux
    cflux = cflux*cwave^2/light           ; observed frame
    cflux = -2.5*alog10(cflux>1D-50)-48.6 ; [AB mag]

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], $
      xsty=1, ysty=1, xrange=[0.15,1.25], yrange=[25,18], $
      ytitle='AB mag', xtitle='Observed-Frame Wavelength (\mu'+'m)'
    djs_oplot, cwave/1D4, cflux, line=0, color='grey'
    oploterror, weff/1D4, mag[*,1], magerr[*,1], psym=symcat(15), $
      symsize=2.0, color=im_color('firebrick')
    im_legend, ['EW([O II]) = 220 \AA'], /left, /top, box=0, margin=0, charsize=1.8
    
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
    
stop    
    
; ---------------------------------------------------------------------------
; see ISEDFIT_CALIBRATE_EMLINES
    
    datapath = getenv('IM_PROJECTS_DIR')+'/isedfit/'
;   path = getenv('IM_PROJECTS_DIR')+'/desi/'

    allatlas = mrdfits(datapath+'integrated_atlas_speclinefit_v3.1_nodust.fits.gz',1)
;   allatlas = mrdfits(datapath+'integrated_atlas_speclinefit_v3.1.fits.gz',1)
    allsdss = read_vagc_garching(/ispecline)
    allhiiinfo = read_hiiregions(linefit=allhii)

    snrcut = 5.0
    
; remove AGN from the galaxy samples and also require S/N>5 on all
; emission lines of interest 
    sclass = iclassification(allsdss,snrcut=snrcut,ratios=sratios)
    keep = where($
      strmatch(sratios.final_class,'*AGN*') eq 0 and $
;     strmatch(sratios.final_class,'*SF*') and $
      allsdss.oiii_5007[0]/allsdss.oiii_5007[1] gt snrcut and $
      allsdss.oii_3727[0]/allsdss.oii_3727[1] gt snrcut and $
      allsdss.h_beta[0]/allsdss.h_beta[1] gt snrcut and $
      allsdss.h_alpha[0]/allsdss.h_alpha[1] gt snrcut and $
      allsdss.nii_6584[0]/allsdss.nii_6584[1] gt snrcut and $
      allsdss.sii_6716[0]/allsdss.sii_6716[1] gt snrcut and $
      allsdss.sii_6731[0]/allsdss.sii_6731[1] gt snrcut,nsdss)
    splog, 'SDSS ', nsdss
    sdss = allsdss[keep]
    
    aclass = iclassification(allatlas,snrcut=snrcut,ratios=aratios)
    keep = where($
      strmatch(aratios.final_class,'*AGN*') eq 0 and $
;     strmatch(aratios.final_class,'*SF*') and $
      allatlas.oiii_5007[0]/allatlas.oiii_5007[1] gt snrcut and $
      allatlas.oii_3727[0]/allatlas.oii_3727[1] gt snrcut and $
      allatlas.h_beta[0]/allatlas.h_beta[1] gt snrcut and $
      allatlas.h_alpha[0]/allatlas.h_alpha[1] gt snrcut and $
      allatlas.nii_6584[0]/allatlas.nii_6584[1] gt snrcut and $
      allatlas.sii_6716[0]/allatlas.sii_6716[1] gt snrcut and $
      allatlas.sii_6731[0]/allatlas.sii_6731[1] gt snrcut,natlas)
    splog, 'ATLAS ', natlas
    atlas = allatlas[keep]

    keep = where($
      allhii.oiii_5007[0]/allhii.oiii_5007[1] gt snrcut and $
      allhii.oii_3727[0]/allhii.oii_3727[1] gt snrcut and $
      allhii.h_beta[0]/allhii.h_beta[1] gt snrcut and $
      allhii.h_alpha[0]/allhii.h_alpha[1] gt snrcut and $
      allhii.nii_6584[0]/allhii.nii_6584[1] gt snrcut and $
      allhii.sii_6716[0]/allhii.sii_6716[1] gt snrcut and $
      allhii.sii_6731[0]/allhii.sii_6731[1] gt snrcut,nhii)
    splog, 'HII ', nhii
    hii = allhii[keep]
    hiiinfo = allhiiinfo[keep]
    
; build the line-ratios
    lineratio, sdss, 'OIII_5007', 'H_BETA', 'OII_3727', $
      'H_BETA', soiiihb, soiiihb_err, soiihb, soiihb_err, snrcut=snrcut
    lineratio, atlas, 'OIII_5007', 'H_BETA', 'OII_3727', $
      'H_BETA', aoiiihb, aoiiihb_err, aoiihb, aoiihb_err, snrcut=snrcut
    lineratio, hii, 'OIII_5007', 'H_BETA', 'OII_3727', $
      'H_BETA', hoiiihb, hoiiihb_err, hoiihb, hoiihb_err, snrcut=snrcut

    lineratio, sdss, 'OIII_5007', 'H_BETA', 'NII_6584', $
      'H_ALPHA', soiiihb, soiiihb_err, sniiha, sniiha_err, snrcut=snrcut
    lineratio, atlas, 'OIII_5007', 'H_BETA', 'NII_6584', $
      'H_ALPHA', aoiiihb, aoiiihb_err, aniiha, aniiha_err, snrcut=snrcut
    lineratio, hii, 'OIII_5007', 'H_BETA', 'NII_6584', $
      'H_ALPHA', hoiiihb, hoiiihb_err, hniiha, hniiha_err, snrcut=snrcut

    lineratio, sdss, 'OIII_5007', 'H_BETA', ['SII_6716','SII_6731'], $
      'H_ALPHA', soiiihb, soiiihb_err, ssiiha, ssiiha_err, snrcut=snrcut
    lineratio, atlas, 'OIII_5007', 'H_BETA', ['SII_6716','SII_6731'], $
      'H_ALPHA', aoiiihb, aoiiihb_err, asiiha, asiiha_err, snrcut=snrcut
    lineratio, hii, 'OIII_5007', 'H_BETA', ['SII_6716','SII_6731'], $
      'H_ALPHA', hoiiihb, hoiiihb_err, hsiiha, hsiiha_err, snrcut=snrcut

; fit the median relations, giving equal weight to every sample
    xx = [soiiihb,aoiiihb,hoiiihb]
    weight = [soiiihb*0.0+1.0/nsdss,aoiiihb*0.0+1.0/natlas,hoiiihb*0.0+1.0/nhii]

    binsz = 0.1
    oiihbmed = im_medxbin(xx,[soiihb,aoiihb,hoiihb],binsz,weight=weight)
    niihamed = im_medxbin(xx,[sniiha,aniiha,hniiha],binsz,weight=weight)
    siihamed = im_medxbin(xx,[ssiiha,asiiha,hsiiha],binsz,weight=weight)

    oiihbcoeff = poly_fit(oiihbmed.medx,oiihbmed.medy,3)
    niihacoeff = poly_fit(niihamed.medx,niihamed.medy,3)
    siihacoeff = poly_fit(siihamed.medx,siihamed.medy,3)

    splog, '[OII]/Hb coeff ', oiihbcoeff
    splog, '[NII]/Ha coeff ', niihacoeff
    splog, '[SII]/Ha coeff ', siihacoeff
    
; make the plots
    nogrey = 0
    oiiihbrange = [-1.0,1.0]
    oiiihbaxis = range(oiiihbrange[0]+0.1,oiiihbrange[1]-0.1,200)
    
    psfile = path+'oiiihb_oiihb.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.7, charsize=1.9, height=5.0
    
; [OII]/Hb vs [OIII]/Hb
    im_hogg_scatterplot, soiiihb, soiihb, position=pos[*,0], nogrey=nogrey, $
      xsty=1, ysty=1, xrange=oiiihbrange, yrange=[-1,1.2], $
      ytitle=textoidl('log ([OII] \lambda3727 / H\beta)'), $
      xtitle=textoidl('log ([O III] \lambda5007 / H\beta)')
    djs_oplot, hoiiihb, hoiihb, psym=symcat(6,thick=4), color='orange', symsize=0.5
    djs_oplot, aoiiihb, aoiihb, psym=symcat(16), color='blue'
;   oploterror, oiihbmed.medx, oiihbmed.medy, oiihbmed.sigy, $
;     psym=symcat(9,thick=8), symsize=3.0, errthick=6.0

;   djs_oplot, oiihbmed.medx, oiihbmed.medy, line=0, thick=10
;   djs_oplot, oiihbmed.medx, oiihbmed.medy+0.2, line=5, thick=10
;   djs_oplot, oiihbmed.medx, oiihbmed.medy-0.2, line=5, thick=10

;   djs_oplot, oiihbmed.medx, oiihbmed.quant75, line=5, thick=10
;   djs_oplot, oiihbmed.medx, oiihbmed.quant25, line=5, thick=10
    djs_oplot, oiiihbaxis, poly(oiiihbaxis,oiihbcoeff), line=0, thick=10
    djs_oplot, oiiihbaxis, poly(oiiihbaxis,oiihbcoeff)+0.2, line=5, thick=10
    djs_oplot, oiiihbaxis, poly(oiiihbaxis,oiihbcoeff)-0.2, line=5, thick=10
    
    im_legend, ['SDSS Galaxies','Integrated Spectra','HII Regions'], $
      psym=[15,16,6], color=['grey','blue','orange'], /left, /top, $
      box=0, margin=0, symthick=10, charsize=1.8
    im_legend, 'Median Fit', /right, /bottom, box=0, line=0, $
      thick=10, pspacing=1.9, charsize=2.1, margin=0
    
    im_plotconfig, /psclose, psfile=psfile, /pdf
    
stop    
    
return
end
    
