pro sdss_isedfit_test, sdss
; jm06mar13uofa - test ISEDFIT on a random sample of SDSS galaxies

    datapath = ages_path(/mass)
    if (n_elements(sdss) eq 0L) then sdss = mrdfits(datapath+'sdss_isedfit_test.fits.gz',1,/silent)
    
; ---------------------------------------------------------------------------
; write out the random sub-sample of SDSS galaxies
; ---------------------------------------------------------------------------

;   sdss = mrdfits(sdss_path()+'sdss_ancillary_data_dr4.fits.gz',1,/silent)
;   indx = uniq_indx(1000,n_elements(sdss),/sort)
;   mwrfits, sdss[indx], datapath+'sdss_isedfit_test.fits', /create
;   spawn, ['gzip -f '+datapath+'sdss_isedfit_test.fits'], /sh
    
; ---------------------------------------------------------------------------
; photometric bands of interest
; ---------------------------------------------------------------------------

    modelsprefix = 'sdss_isedfit_ugriz'
    chi2prefix = 'sdss_isedfit_ugriz'

    filterlist = ['sdss_u0','sdss_g0','sdss_r0','sdss_i0','sdss_z0']+'.par'
    filtinfo = im_filterspecs(filterlist=filterlist,/verbose)
    nfilter = n_elements(filterlist)
    
; ---------------------------------------------------------------------------
; generate the model magnitudes for the ugriz photometry
; ---------------------------------------------------------------------------

    tau = [0.0,0.1,0.2,0.3,0.5,0.8,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.5,8.0,10.0,12.5,16.0,20.0]
    isedfit_models, filterlist, datapath=datapath, modelsprefix=modelsprefix, $
      minredshift=0.01D, maxredshift=0.35D, dredshift=0.01D, Zstellar=[2L,4L,5L], $
      tau=tau, minage=0.05D, maxage=13.5D, /allages, debug=0, write=1

;   tau = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,8.0,10.0,13.0,16.0,20.0]
;   isedfit_models, filterlist, datapath=datapath, modelsprefix=modelsprefix, $
;     minredshift=0.01D, maxredshift=0.35D, Zstellar=[2L,4L], tau=tau, $
;     minage=0.05D, maxage=13.5D, /allages, write=1, debug=0

;   isedfit_models, filterlist, datapath=datapath, modelsprefix=modelsprefix, $
;     minredshift=0.01D, maxredshift=0.35D, minage=0.1D, maxage=13.0D, $
;     dlogage=0.05D, Zstellar=[2L,4L], tau=tau, write=0, debug=1

; ---------------------------------------------------------------------------
; convert the photometry in each band to MAGGIES and do the fit!
; ---------------------------------------------------------------------------

    bands = 'sdss_'+['u','g','r','i','z']
    maggies = convert_to_maggies(sdss,bands,filtinfo.vega2ab*0.0,maggies_invvar=maggies_invvar)

    fitindx = lindgen(n_elements(sdss))
;   fitindx = lindgen(10)+100
;   fitindx = lindgen(100)+499
    thesefilters = filterlist

    isedfit_chi2,maggies[*,fitindx],maggies_invvar[*,fitindx],sdss[fitindx].z_obj,chi2,chi2_info,filterlist=thesefilters,$
      modelsprefix=modelsprefix,chi2prefix=chi2prefix,galaxy=sdss[fitindx].galaxy,maxold=0,/write

    isedfit, result, datapath=datapath, chi2prefix=chi2prefix, maxold=0, /write
    
    isedfit_qaplot, datapath=datapath, chi2prefix=chi2prefix, maxold=0, /postscript

    isedfit_measure, datapath=datapath, restfilterfile=restfilterfile, chi2prefix=chi2prefix, maxold=0

stop

    isedfit_qaplot, datapath=datapath, chi2prefix=chi2prefix;, /postscript

; read the results and generate a QA plot

;   result = mrdfits(datapath+'sdss_isedfit_ugriz_mass_maxold.fits.gz',1,/silent)
    result = mrdfits(datapath+'sdss_isedfit_ugriz_mass.fits.gz',1,/silent)
    good = where(result.mass_chi2min gt 0.0) & result.mass_chi2min = alog10(result.mass_chi2min)

    dfpsplot, datapath+'kauffmann_vs_isedfit.ps', /square, /color
    good = where((result.chi2min lt 10.0) and (sdss.kauffmann_mass gt -900.0),ngood)
    stats = im_stats(sdss[good].kauffmann_mass-result[good].mass_chi2min)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    djs_plot, result[good].mass_chi2min, sdss[good].kauffmann_mass, ps=4, syms=0.5, xsty=3, ysty=3, $
      xthick=5.0, ythick=5.0, charsize=1.8, charthick=5.0, xrange=[8,13.5], yrange=[8,13.5], $
      xtitle='log (M/M'+sunsymbol()+') [isedfit]', ytitle='log (M/M'+sunsymbol()+') [Kauffmann]'
    djs_oplot, !x.crange, !y.crange, line=0, thick=3.0, color='red'
    legend, textoidl('\Delta(log M) = '+xstr), /left, /top, box=0, charsize=1.5, charthick=5.0
; ---
    good = where((result.chi2min lt 10.0) and (sdss.kcorr_mass gt -900.0),ngood)
    stats = im_stats(sdss[good].kcorr_mass-result[good].mass_chi2min)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    djs_plot, result[good].mass_chi2min, sdss[good].kcorr_mass, ps=4, syms=0.5, xsty=3, ysty=3, $
      xthick=5.0, ythick=5.0, charsize=1.8, charthick=5.0, xrange=[8,13.5], yrange=[8,13.5], $
      xtitle='log (M/M'+sunsymbol()+') [isedfit]', ytitle='log (M/M'+sunsymbol()+') [k-correct]'
    djs_oplot, !x.crange, !y.crange, line=0, thick=3.0, color='red'
    legend, textoidl('\Delta(log M) = '+xstr), /left, /top, box=0, charsize=1.5, charthick=5.0
; ---
    good = where((sdss.kcorr_mass gt -900.0) and (sdss.kauffmann_mass gt -900.0),ngood)
    stats = im_stats(sdss[good].kauffmann_mass-sdss[good].kcorr_mass)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    djs_plot, sdss[good].kcorr_mass, sdss[good].kauffmann_mass, ps=4, syms=0.5, xsty=3, ysty=3, $
      xthick=5.0, ythick=5.0, charsize=1.8, charthick=5.0, xrange=[8,13.5], yrange=[8,13.5], $
      xtitle='log (M/M'+sunsymbol()+') [k-correct]', ytitle='log (M/M'+sunsymbol()+') [Kauffmann]'
    djs_oplot, !x.crange, !y.crange, line=0, thick=3.0, color='red'
    legend, textoidl('\Delta(log M) = '+xstr), /left, /top, box=0, charsize=1.5, charthick=5.0
; ---
    good = where((sdss.tremonti_bc_mass gt -900.0) and (sdss.kauffmann_mass gt -900.0),ngood)
    stats = im_stats(sdss[good].kauffmann_mass-sdss[good].tremonti_bc_mass)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    djs_plot, sdss[good].tremonti_bc_mass, sdss[good].kauffmann_mass, ps=4, syms=0.5, xsty=3, ysty=3, $
      xthick=5.0, ythick=5.0, charsize=1.8, charthick=5.0, xrange=[8,13.5], yrange=[8,13.5], $
      xtitle='log (M/M'+sunsymbol()+') [Tremonti BC03]', ytitle='log (M/M'+sunsymbol()+') [Kauffmann]'
    djs_oplot, !x.crange, !y.crange, line=0, thick=3.0, color='red'
    legend, textoidl('\Delta(log M) = '+xstr), /left, /top, box=0, charsize=1.5, charthick=5.0
    dfpsclose
    
return
end
    
