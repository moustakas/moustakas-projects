pro sfrm_sdss_mf, quiescent=quiescent, active=active, $
  double_schechter=double_schechter, debug=debug
; jm10feb15ucsd - build the local stellar mass functions using the SDSS

    sfrmpath = ages_path(/projects)+'sfrm/'

    suffix = 'all'
    if keyword_set(quiescent) then suffix = 'quiescent'
    if keyword_set(active) then suffix = 'active'
    if keyword_set(double_schechter) then $
      dsuffix = '_double' else dsuffix = ''
    
    fitfile = sfrmpath+'sdss_mf_fit_'+suffix+dsuffix+'.fits'
    datafile = sfrmpath+'sdss_mf_data_'+suffix+'.fits'

    minmass = 9.5
    maxis1 = im_array(minmass,12.0,0.01)

; read the parent catalog
    parent = read_sfrm_sample(/sdss)
    binsize = sfrm_binsize(histmin=histmin,histmax=histmax)
    parinfo = init_mf_parinfo(fixslope=0,quiescent=quiescent,$
      active=active,double_schechter=double_schechter)

    init_mf_results, 1, mf_fit=mf_fit, mf_data=mf_data, $
      double_schechter=double_schechter

    nuvmr = parent.k_galex_absmag_01[1]-parent.k_ugriz_absmag_01[2]
;   rmj = parent.ugriz_absmag[2]-(parent.ubvrijhk_absmag[5]+jv2ab)
    qq = select_quiescent(nuvmr,active=aa)
    if keyword_set(quiescent) then parent = parent[qq]
    if keyword_set(active) then parent = parent[aa]
    ngal = n_elements(parent)

; build the Vmax-weighted mass function and then fit it
    mf_vmax, parent.k_mass, 1.0/parent.vmax_evol/binsize, binsize=binsize, $
      histmin=histmin, histmax=histmax, binmass=binmass, $
      phi=phi, errphi=phierr, minmass=minmass, fullbin=fullbin, $
      number=number

    full = where(fullbin)
    if keyword_set(double_schechter) then begin
       mf_fit_schechter_plus, binmass[full], phi[full], $
         phierr[full], fit, parinfo=parinfo
    endif else begin
       mf_fit_schechter, binmass[full], phi[full], $
         phierr[full], fit, parinfo=parinfo
    endelse
    
; store the results       
    fill_sfrm_results, 0, ngal, mf_fit, mf_data, fit, $
      binmass, phi, phierr, fullbin, number
    mf_fit.minmass = minmass
    
; QAplot
    if keyword_set(debug) then begin
       ploterror, binmass, phi, phierr, ps=10, xsty=1, ysty=1, $
         xrange=minmax(maxis1)+[-0.5,+0.0], yrange=[1E-6,0.1], /ylog
       if keyword_set(double_schechter) then $
         djs_oplot, maxis1, mf_schechter_plus(maxis1,fit), color='blue' else $
           djs_oplot, maxis1, mf_schechter(maxis1,fit), color='blue'
       oplot_bgd08_mf, maxis1, params=params, log=0, color='red', line=5
    endif

; write out
    im_mwrfits, mf_fit, fitfile, /clobber
    im_mwrfits, mf_data, datafile, /clobber
    
return
end
    
