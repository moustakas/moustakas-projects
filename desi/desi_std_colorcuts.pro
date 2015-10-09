pro desi_std_colorcuts, qaplots=qaplots
; jm15oct01siena - develop a set of grz color-cuts for F-standard stars for DESI 

    common com_star, star_info, star_restflux, star_restwave, star_ugriz, star_grz_decals

    light = im_light(/km)
    
    std_rbright = 15.0
    std_rfaint = 19.0

    outdir = getenv('DESI_ROOT')+'/qaplots/'

    version = desi_star_templates_version()
    templatepath = getenv('DESI_ROOT')+'/spectro/templates/'+$
      'star_templates/'+version+'/'
    
    if n_elements(star_restflux) eq 0L then begin
       restfile = templatepath+'star_templates_'+version+'.fits'
       
       star_restflux = mrdfits(restfile,0,resthdr)
       star_restwave = 10D^make_wave(resthdr)
       star_info = mrdfits(restfile,1)

; synthesize ugriz photometry
       star_ugriz = fltarr(5,n_elements(star_info))
       for ii = 0, n_elements(star_info)-1 do star_ugriz[*,ii] = $
         reform(k_project_filters(k_lambda_to_edges(star_restwave),$
         star_restflux[*,ii],filterlist=sdss_filterlist()))
       star_ugriz = -2.5*alog10(star_ugriz)

; synthesize DECaLS grz photometry
       star_grz_decals = fltarr(3,n_elements(star_info))
       for ii = 0, n_elements(star_info)-1 do star_grz_decals[*,ii] = $
         reform(k_project_filters(k_lambda_to_edges(star_restwave),$
         star_restflux[*,ii],filterlist=decals_filterlist(/grz)))
       star_grz_decals = -2.5*alog10(star_grz_decals)
    endif
    ntemplate = n_elements(star_info)

; define a distance in color space; see equation 4 in Dawson+13; just 6 stars
; satisfy this criterion!  but they roughly match the distribution of physical
; properties of the BOSS standards
    ug = reform(star_ugriz[0,*]-star_ugriz[1,*])
    gr = reform(star_ugriz[1,*]-star_ugriz[2,*])
    ri = reform(star_ugriz[2,*]-star_ugriz[3,*])
    iz = reform(star_ugriz[3,*]-star_ugriz[4,*])
    
    mdist = sqrt((ug-0.82)^2 + (gr-0.30)^2 + (ri-0.09)^2 + (iz-0.02)^2)
    keep = where(mdist lt 0.08,nkeep)

; use the cuts proposed by Schlegel
    gr_decals = reform(star_grz_decals[0,*]-star_grz_decals[1,*])
    rz_decals = reform(star_grz_decals[1,*]-star_grz_decals[2,*])
    print, median(gr_decals[keep]), minmax(gr_decals[keep]), $
      max(gr_decals[keep])-min(gr_decals[keep])
    print, median(rz_decals[keep]), minmax(rz_decals[keep]), $
      max(rz_decals[keep])-min(rz_decals[keep])

    ddist = sqrt((gr_decals-0.32)^2 + (rz_decals-0.13)^2)
    these = where(ddist lt 0.06,nthese)
    print, these

    psfile = outdir+'std_grz_colorcuts.eps'
    im_plotconfig, 0, pos1, psfile=psfile, height=5.0, $
      xmargin=[1.3,0.4], width=6.8
    djs_plot, [0], [0], /nodata, position=pos1, $
      xsty=1, ysty=1, xrange=[-0.5,0.7], yrange=[-0.5,1.1], $
      xtitle='r - z', ytitle='g - r'

    djs_oplot, star_grz_decals[1,*]-star_grz_decals[2,*], $
      star_grz_decals[0,*]-star_grz_decals[1,*], psym=symcat(16), $
      symsize=0.7, color=cgcolor('grey')
    djs_oplot, star_grz_decals[1,these]-star_grz_decals[2,these], $
      star_grz_decals[0,these]-star_grz_decals[1,these], psym=symcat(16), $
      symsize=0.9, color=cgcolor('firebrick')
    djs_oplot, star_grz_decals[1,keep]-star_grz_decals[2,keep], $
      star_grz_decals[0,keep]-star_grz_decals[1,keep], psym=symcat(15), $
      symsize=0.7, color=cgcolor('dodger blue')
    im_legend, '[(g-r)-0.32]^2 + [(r-z)-0.13]^2 < 0.06', $
      /right, /bottom, box=0, charsize=1.5, margin=0
    im_legend, ['All Stars','ugriz-selected F-stars','grz-selected F-stars (proposed)'], $
      /left, /top, box=0, color=['grey','dodger blue','firebrick'], $
      psym=[16,16,15], charsize=1.5

;   im_legend, '[(u-g)-0.82]^2 + [(g-r)-0.30]^2 + !c!c '+$
;     '[(r-i)-0.09]^2 + [(i-z)-0.02]^2 < 0.08', $
;     /right, /bottom, box=0, charsize=1.5, margin=2
    im_plotconfig, psfile=psfile, /psclose, /png


stop
    
    boss = mrdfits(outdir+'boss-stds.fits',1)
    splog, median(boss.teff), median(boss.logg), median(boss.feh)
    struct_print, star_info[keep]
       
; assign r-band magnitudes uniformly
    rmag = randomu(seed,nstd)*(std_rbright-std_rfaint)+std_rfaint
    these = long(randomu(seed,nstd)*nkeep)
    
; build the output truth table
    outinfo = {$
      objtype:  'STD',$
      sdss_r:     0.0,$
      sdss_ug:    0.0,$
      sdss_gr:    0.0,$
      sdss_ri:    0.0,$
      sdss_iz:    0.0,$
      templateid: 0L}
    outinfo = replicate(outinfo,nstd)
    outinfo.sdss_r = rmag
    outinfo.sdss_ug = ug[keep[these]]
    outinfo.sdss_gr = gr[keep[these]]
    outinfo.sdss_ri = ri[keep[these]]
    outinfo.sdss_iz = iz[keep[these]]
    outinfo.templateid = keep[these]
    
; normalize each spectrum to the desired r-magnitude and write out
    restwave = k_lambda_to_edges(star_restwave)
    npix = n_elements(star_restwave)
    outflux = fltarr(npix,nstd)
    for ii = 0L, nstd-1 do begin
       robj = (reform(k_project_filters(restwave,$
         star_restflux[*,outinfo[ii].templateid],$
         filterlist='decam_r.par')))[0]
       
       outflux1 = reform(star_restflux[*,outinfo[ii].templateid])*$
         10.0^(-0.4*outinfo[ii].sdss_r)/robj
       outflux[*,ii] = outflux1
       
;         djs_plot, star_restwave, outflux1, xsty=1, ysty=1 ;, /xlog
;         struct_print, star_info[outinfo[ii].templateid], /no_head
;         cc = get_kbrd(1)
    endfor

; write out
    outfile = outdir+'std_templates.fits'

    cdelt1 = 1.4341626783D-05              ; [pixel size in log-10 A]
    velpixsize = cdelt1*light*alog(10)     ; [km/s]
    
    mkhdr, outhdr, outflux, /extend
    sxdelpar, outhdr, 'DATE'
    sxdelpar, outhdr, 'COMMENT'
;   sxaddpar, outhdr, 'VERSION', version, ' template version number'
    sxaddpar, outhdr, 'OBJTYPE', 'STD', ' object type'
    sxaddpar, outhdr, 'DISPAXIS', 1, ' dispersion axis'
    sxaddpar, outhdr, 'CTYPE1', 'WAVE-WAV', ' meaning of wavelength array'
    sxaddpar, outhdr, 'CUNIT1', 'Angstrom', ' units of wavelength array'
    sxaddpar, outhdr, 'CRPIX1', 1, ' reference pixel number'
    sxaddpar, outhdr, 'CRVAL1', min(alog10(star_restwave)), ' reference log10(Angstrom)'
    sxaddpar, outhdr, 'CDELT1', cdelt1, ' delta log10(Angstrom)'
    sxaddpar, outhdr, 'LOGLAM', 1, ' log10 spaced wavelengths?'
    sxaddpar, outhdr, 'WAVEMIN', min(star_restwave), ' minimum wavelength (Angstrom)'
    sxaddpar, outhdr, 'WAVEMAX', max(star_restwave), ' maximum wavelength (Angstrom)'
    sxaddpar, outhdr, 'WAVEUNIT', 'Angstrom', ' wavelength units'
    sxaddpar, outhdr, 'AIRORVAC', 'vac', ' wavelengths in vacuum (vac) or air'
    sxaddpar, outhdr, 'VELSCALE', velpixsize, ' pixel size in km/s'
    sxaddpar, outhdr, 'BUNIT', 'erg/s/cm2/A', ' spectrum flux units'
    sxaddpar, outhdr, 'FLUXUNIT', 'erg/s/cm2/A', ' spectrum flux units'
    
; we need to do some MWRFITS jujitsu to get the metadata header right 
    units = [$
      'OBJTYPE,,object type',$
      'SDSS_R,,r-band magnitude',$
      'SDSS_UG,,u-g color',$
      'SDSS_GR,,g-r color',$
      'SDSS_RI,,r-i color',$
      'SDSS_IZ,,i-z color',$
      'TEMPLATEID,,unique template ID number (0-indexed)']
    newheader = [$
      'EXTNAME'+" = '"+string('METADATA','(A-8)')+"'           /"]
    
    mwrfits, 0, outfile, /create
    mwrfits, outinfo, outfile, /silent
    metahdr =  im_update_header(outfile,units,newheader)
    
    im_mwrfits, outflux, outfile, outhdr, /clobber, /nogzip
    im_mwrfits, outinfo, outfile, metahdr, /append, /nogzip

; --------------------------------------------------
; make some QAplots
    if keyword_set(qaplots) then begin

       std_tempinfo = mrdfits(outdir+'std_templates.fits',1)
       boss = mrdfits(outdir+'boss-stds.fits',1)

       psfile = outdir+'qa_std.ps'
       im_plotconfig, 0, pos1, psfile=psfile, height=4.0, $
         ymargin=[0.4,6.6], xmargin=[1.2,0.3]
       djs_plot, [0], [0], /nodata, position=pos1, $
         xsty=1, ysty=1, xrange=[-1,3], yrange=[-0.7,2], $
         xtitle='r - z', ytitle='g - r'
       djs_oplot, star_ugriz[2,*]-star_ugriz[4,*], $
         star_ugriz[1,*]-star_ugriz[2,*], psym=symcat(16), $
         symsize=0.7, color=cgcolor('grey')
       djs_oplot, std_tempinfo.sdss_ri+std_tempinfo.sdss_iz, $
         std_tempinfo.sdss_gr, psym=symcat(9), $
         color=cgcolor('dodger blue'), symsize=0.7

       im_legend, '[(u-g)-0.82]^2 + [(g-r)-0.30]^2 + !c!c '+$
         '[(r-i)-0.09]^2 + [(i-z)-0.02]^2 < 0.08', $
         /right, /bottom, box=0, charsize=1.5, margin=2

; magnitude histogram
       im_plotconfig, 12, pos2, ymargin=[6.0,1.1], height=4.5, $
         xspace=1.0, xmargin=[1.2,0.3]
       djs_plot, [0], [0], /nodata, /noerase, position=pos2[*,0], $
         xsty=1, ysty=1, yrange=[0,1.05], $
         xrange=[std_rbright-0.4,std_rfaint+0.2], $
         ytitle='Relative Number', xtitle='r (AB mag)'
       im_plothist, std_tempinfo.sdss_r, bin=0.2, /overplot, /peak

       djs_plot, [0], [0], /nodata, /noerase, position=pos2[*,1], $
         xsty=1, ysty=1, yrange=[4.5,8], xrange=[-3.0,0.5], $
         ytitle='T_{eff} (1000 K)', xtitle='[Fe/H]', $
         xtickinterval=1.0, ytickinterval=1.0
       djs_oplot, boss.feh, boss.teff/1000, psym=6, symsize=1.0
       djs_oplot, star_info[std_tempinfo.templateid].feh, $
         star_info[std_tempinfo.templateid].teff/1000.0, $
         psym=symcat(16), color=cgcolor('dodger blue'), symsize=1.0
       im_legend, ['BOSS Standards','STD Templates'], /left, /top, box=0, $
         psym=[6,16], color=['black','dodger blue'], charsize=1.4, $
         margin=0

       im_plotconfig, psfile=psfile, /psclose, /pdf
       
stop

    endif

stop

return
end
