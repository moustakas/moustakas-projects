 pro ediscs_mass_vs_redshift, ediscs, postscript=postscript, write=write
; jm07sep16nyu - based on AGES_MB_VS_REDSHIFT

    outpath = ediscs_path(/projects)+'mz/'
    
    lsun = 3.826D33             ; [erg/s]
    light = 2.99792458D5        ; speed of light [km/s]
    mbolsun = 4.74              ; bolometric
    pc = 3.086D18               ; [cm]
    distance = 10.0*pc          ; fiducial distance [cm]
    fluxarea = 4.0*!dpi*distance^2
    
    I22 = 22.0
    I15 = 15.0
    
    red, h100=0.7, omega_0=0.3, omega_lambda=0.7 ; cosmological parameters

    if (n_elements(ediscs) eq 0L) then ediscs = read_ediscs()

; read the model    
    
    modelspath = getenv('SFHGRID_DIR')+'/basemodels/'
    csffits = 'burst_m62_salp_20.0Gyr.fits.gz'
    sspfits = 'tau_m62_salp_00.0Gyr.fits.gz'

    csf = mrdfits(modelspath+csffits,1,/silent)
    ssp = mrdfits(modelspath+sspfits,1,/silent)

    restwave = csf.wave
    restwave_edges = k_lambda_to_edges(restwave)

; construct the redshift and age grid; set the zeropoint to be 12 Gyr
; at z=0 (implying a formation redshift of z=4.111

;   zform = [3.9,4.0,4.1,4.2,4.3,4.4]
;   print, interpol(zform,getage(0.0)-getage(zform),12.0)
    
    minz = 0.01 & maxz = 1.05 & dz = 0.01
    z = findgen((maxz-minz)/dz+1)*dz+minz
    nz = n_elements(z)

    dlum = dluminosity(z,/Mpc) ; luminosity distance [Mpc]

    age = getage(z)-getage(5.0) ; z_form = 5
    get_element, csf.age, age*1E9, ageindx

; define the filter functions

    filterlist = ['bessell_B','dwfs_I']+'.par'
    info = im_filterspecs(filterlist=filterlist,/verbose)
    Binfo = info[0] & Iinfo = info[1]

    lbsun_lsun = 10.0^(-0.4*(Binfo.solarmags-mbolsun))
    lbsun = lsun*10.0^(-0.4*(Binfo.solarmags-mbolsun))

    result = replicate({z: 0.0, csf_MB_I22: 0.0, ssp_MB_I22: 0.0, $
      csf_MB_I15: 0.0, ssp_MB_I15: 0.0, csf_mass_I22: 0.0,  ssp_mass_I22: 0.0},nz)
    result.z = z
    
    for iz = 0L, nz-1L do begin

       csf_mass = csf.m_[ageindx[iz]]
       ssp_mass = ssp.m_[ageindx[iz]]
       
       csf_flux     = lsun*csf.flux[*,ageindx[iz]]/(1.0+z[iz])/fluxarea ; [erg/s/cm2/A/M_sun]
       csf_restflux = lsun*csf.flux[*,ageindx[iz]]/fluxarea
       ssp_flux     = lsun*ssp.flux[*,ageindx[iz]]/(1.0+z[iz])/fluxarea
       ssp_restflux = lsun*ssp.flux[*,ageindx[iz]]/fluxarea

       wave = restwave*(1.0+z[iz])
       wave_edges = k_lambda_to_edges(wave)

       csf_Brest_synth = im_filtermag(restwave_edges,csf_restflux,filterlist='bessell_B.par')
       csf_Iobs_synth  = im_filtermag(wave_edges,csf_flux,filterlist='dwfs_I.par')
       ssp_Brest_synth = im_filtermag(restwave_edges,ssp_restflux,filterlist='bessell_B.par')
       ssp_Iobs_synth  = im_filtermag(wave_edges,ssp_flux,filterlist='dwfs_I.par')

       m_csf_Brest_synth = -2.5*alog10(csf_Brest_synth) - Binfo.vega2ab ; [Vega]
       m_csf_Iobs_synth  = -2.5*alog10(csf_Iobs_synth) - Iinfo.vega2ab  ; [Vega]
       m_ssp_Brest_synth = -2.5*alog10(ssp_Brest_synth) - Binfo.vega2ab ; [Vega]
       m_ssp_Iobs_synth  = -2.5*alog10(ssp_Iobs_synth) - Iinfo.vega2ab  ; [Vega]

       result[iz].csf_MB_I22 = I22 + (m_csf_Brest_synth - m_csf_Iobs_synth) - 5*alog10(dlum[iz]) - 25
       result[iz].ssp_MB_I22 = I22 + (m_ssp_Brest_synth - m_ssp_Iobs_synth) - 5*alog10(dlum[iz]) - 25
       result[iz].csf_MB_I15 = I15 + (m_csf_Brest_synth - m_csf_Iobs_synth) - 5*alog10(dlum[iz]) - 25
       result[iz].ssp_MB_I15 = I15 + (m_ssp_Brest_synth - m_ssp_Iobs_synth) - 5*alog10(dlum[iz]) - 25

       result[iz].ssp_mass_I22 = ssp.m__lb[ageindx[iz]]*10^(-0.4*(result[iz].ssp_MB_I22-Binfo.solarmags))
       result[iz].csf_mass_I22 = csf.m__lb[ageindx[iz]]*10^(-0.4*(result[iz].csf_MB_I22-Binfo.solarmags))

    endfor

    if keyword_set(write) then begin
       splog, 'Writing '+outpath+'ediscs_mb_vs_redshift.fits.'
       mwrfits, result, outpath+'ediscs_mb_vs_redshift.fits', /create
       spawn, ['gzip -f '+outpath+'ediscs_mb_vs_redshift.fits'], /sh
    endif
    
; make a plot    

    if keyword_set(postscript) then begin
       postthick1 = 5.0
       dfpsplot, outpath+'ediscs_mass_vs_redshift.ps', /color, /square
    endif else begin
       im_window, 0, /square
       im_window, 2, /square
       postthick1 = 2.0
    endelse

    plotsym, 0, 0.7, /fill

    if (not keyword_set(postscript)) then wset, 2
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=postthick1, $
      ythick=postthick1, charthick=postthick1, charsize=1.8, xrange=[0.3,0.9], $
      yrange=[8.5,12.5], xtitle='Redshift', ytitle='log (M / M_{\odot})'
    djs_oplot, z, alog10(result.csf_mass_I22), line=2, thick=postthick1
    djs_oplot, z, alog10(result.ssp_mass_I22), line=0, thick=postthick1
    indx = where((ediscs.kcorr_mass gt -900.0) and ediscs.memberflag,nindx)
    djs_oplot, ediscs[indx].z, alog10(ediscs[indx].kcorr_mass), ps=8

    legend, textoidl(['\psi = const.','SSP']), /right, /bottom, $
      box=0, charthick=postthick1, charsize=1.7, line=[2,0], thick=postthick1
    legend, textoidl('I_{Vega} < 22'), /right, /top, $
      box=0, charthick=postthick1, charsize=1.7, thick=postthick1
      
    if (not keyword_set(postscript)) then wset, 0
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=postthick1, $
      ythick=postthick1, charthick=postthick1, charsize=1.8, xrange=[0.3,0.9], $
      yrange=[-15,-24], xtitle='Redshift', ytitle='M_{B}'
    djs_oplot, z, result.csf_MB_I22, line=2, thick=postthick1
    djs_oplot, z, result.ssp_MB_I22, line=0, thick=postthick1
;   djs_oplot, z, result.ssp_MB_I15, line=1, color='red'
;   djs_oplot, z, result.csf_MB_I15, line=2, color='red'
    indx = where((ediscs.m_b gt -900.0) and ediscs.memberflag,nindx)
    djs_oplot, ediscs[indx].z, ediscs[indx].m_b, ps=8

    legend, textoidl(['\psi = const.','SSP']), /right, /bottom, $
      box=0, charthick=postthick1, charsize=1.7, line=[2,0], thick=postthick1
    legend, textoidl('I_{Vega} < 22'), /right, /top, $
      box=0, charthick=postthick1, charsize=1.7, thick=postthick1
      
    if keyword_set(postscript) then dfpsclose

stop    
    
return
end
