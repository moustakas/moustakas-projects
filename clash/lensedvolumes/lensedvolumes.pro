function lensedvolumes_getsed, tau, z, zform=zform
; build an SED

    common com_getsed, igmgrid, ssp

    if n_elements(igmgrid) eq 0 then igmgrid = $
      mrdfits(getenv('IMPRO_DIR')+'/etc/igmtau_grid.fits.gz',1)
    if n_elements(ssp) eq 0 then ssp = im_read_fsps(met='Z0.0039',/flam)

    sed = {restwave: ssp.wave, restflux: ssp.wave*0.0, wave: ssp.wave*(1.0+z), flux: ssp.wave*0}
    
    h100 = 0.7D
    pc10 = 3.085678D19          ; fiducial distance [10 pc in cm]
    light = 2.99792458D18       ; speed of light [A/s]
    if n_elements(zform) eq 0 then zform = 15.0
    
    galage = getage(z)-getage(zform)
    
    sfh = {tau: tau, delayed: 0, nburst: 0, tautrunc: 0, maxage: 13.5D}
    csp1 = isedfit_convolve_sfh(ssp,info=sfh,time=galage)

    dlum = pc10*10D^(lf_distmod(z,omega0=0.3D,$ ; [cm]
      omegal0=0.7D)/5D)/h100
    sed.restflux = csp1*(pc10/dlum)^2D
    sed.flux = sed.restflux/(1.0+z)
    
    windx = findex(igmgrid.wave,sed.wave)
    zindx = findex(igmgrid.zgrid,z)
    igm = interpolate(igmgrid.igm,windx,zindx,/grid,missing=1.0)
    sed.flux = sed.flux*igm
    
return, sed
end

pro lensedvolumes, starforming=starforming
; jm12aug13siena - compute some limiting absolute magnitudes

    path = clash_path(/lensedvolumes)
    filt = ['galex_fuv','bessell_B','bessell_V']+'.par'
    nfilt = n_elements(filt)

; redshift parameters    
    zz = range(1.0,12.0,111)
;   zz = 7.0
    nzz = n_elements(zz)

    if keyword_set(starforming) then begin
       tau = 100D
       suffix = 'starforming'
    endif else begin
       tau = 0.01D
       suffix = 'passive'
    endelse
    
; set up the survey parameters    
    nsurvey = 12
    cat = replicate({survey: '', depth: 0.0, filt: '', absmag: fltarr(nfilt,nzz)},nsurvey)
    cat.survey = ['TENIS','UltraVISTA','CANDELS-Wide','CANDELS-Deep',$
      'HUDF09','HUDF09-Parallels','WFC3-ERS','HIPPIES','BoRG09','BoRG12',$
      'NICMOS+','MOIRCS/ISAAC']
    cat.depth = [25.3,24.0,27.1,27.7,29.9,29.3,27.2,26.7,26.3,26.6,26.5,25.5]
    cat.filt = ['twomass_Ks','twomass_H','clash_wfc3_f160w','clash_wfc3_f160w',$
      'clash_wfc3_f160w','clash_wfc3_f160w','clash_wfc3_f160w','clash_wfc3_f160w',$
      'clash_wfc3_f160w','clash_wfc3_f160w','clash_wfc3_f160w','twomass_H']+'.par'

    for iz = 0, nzz-1 do begin
; normalize the SED to the 5-sigma depth
       sed = lensedvolumes_getsed(tau,zz[iz]) ; SSP
       norm = reform(k_project_filters(k_lambda_to_edges(sed.wave),$
         sed.flux,filterlist=cat.filt))/10^(-0.4*cat.depth)
       for ii = 0, nsurvey-1 do begin
          cat[ii].absmag[*,iz] = -2.5*alog10(k_project_filters(k_lambda_to_edges(sed.restwave),$
            sed.restflux/norm[ii],filterlist=filt))-48.6
;         plot, sed.restwave, sed.restflux/norm[ii], xr=[1000,4D4], /xlog
;         djs_oplot, sed.wave, sed.flux/norm[ii], color='orange'
       endfor
    endfor

; write out
    im_mwrfits, cat, path+'lensedvolumes_absmag_'+suffix+'.fits', /clob
    
return
end

