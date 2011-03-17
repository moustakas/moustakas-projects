;+
; NAME:
;       AGES_STELLAR_MASS
;
; PURPOSE:
;       Estimate the stellar masses of the AGES galaxies using SED
;       fitting. 
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Feb 10, U of A
;       jm05aug01uofa - updated for release 2.0
;-

pro ages_stellar_mass, ages, write=write, debug=debug

    if (n_elements(ages) eq 0L) then ages = ages_read_ancillary()
    ngalaxy = n_elements(ages)

    datapath = ages_path(/mass)

    modelsprefix = 'ages_BwRIKJKs_ugriz'
    chi2prefix = 'test_ages_BwRIKJKs'

    filterlist = ['dwfs_Bw','dwfs_R','dwfs_I','dwfs_onis_K','flamex_flamingos_J',$
      'flamex_flamingos_Ks','sdss_'+['u0','g0','r0','i0','z0'] ]+'.par'
    filtinfo = im_filterspecs(filterlist=filterlist,/verbose)
    filtinfo[where(filtinfo.sdssflag)].vega2ab = 0.0
    nfilter = n_elements(filterlist)
    
; ---------------------------------------------------------------------------
; generate the model magnitudes for the BwRIKJKs photometry
; ---------------------------------------------------------------------------

;   tau = [0.0,0.1,0.2,0.3,0.5,0.8,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.5,8.0,10.0,12.5,16.0,20.0]
;   isedfit_models, filterlist, datapath=datapath, modelsprefix=modelsprefix, $
;     minredshift=0.01D, maxredshift=1.0D, dredshift=0.01D, Zstellar=[2L,4L], $
;     tau=tau, minage=0.05D, maxage=13.5D, /allages, debug=0, write=1

; OLD
;   tau = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,8.0,10.0,13.0,16.0,20.0]
;   isedfit_models, filterlist, datapath=datapath, modelsprefix=modelsprefix, $
;     minredshift=0.01D, maxredshift=1.0D, minage=0.1D, maxage=13.0D, $
;     dlogage=0.05D, Zstellar=[2L,4L], tau=tau, /write, debug=debug
 
; ---------------------------------------------------------------------------
; convert the photometry in each band to MAGGIES and do the fit!
; ---------------------------------------------------------------------------

    bands = 'phot_'+['bw','r','i','k','flamj','flamk','sdss_u','sdss_g','sdss_r','sdss_i','sdss_z']
    allmaggies = convert_to_maggies(ages,bands,filtinfo.vega2ab,maggies_invvar=allmaggies_invvar)

    thesefilters = filterlist[0:5]
    maggies = allmaggies[0:5,*]
    maggies_invvar = allmaggies_invvar[0:5,*]
    zobj = ages.z
    galaxy = ages.ages_id

;   fitindx = where(ages.main_flag)
    fitindx = (where(ages.main_flag))[5000:5010]
;   fitindx = lindgen(10)+1100L
;   fitindx = lindgen(ngalaxy)

; compute and minimize chi2, then generate a QA plot

    isedfit_chi2,maggies[*,fitindx],maggies_invvar[*,fitindx],zobj[fitindx],chi2,chi2_info,filterlist=thesefilters,$
      modelsprefix=modelsprefix,chi2prefix=chi2prefix,galaxy=galaxy[fitindx],maxold=0,/write

stop
    
    isedfit, result, chi2prefix=chi2prefix, maxold=0, /write
    
    isedfit_qaplot, datapath=datapath, chi2prefix=chi2prefix, maxold=0, /postscript

    isedfit_measure, datapath=datapath, restfilterfile=restfilterfile, chi2prefix=chi2prefix, maxold=0

;   window, 0 & plotsym, 0, 1, /fill
;   plot, result.mass_chi2min, result.mass_median, ps=8, xr=[8,12.5], yr=[8,12.5], xsty=3, ysty=3
;   plot, result.mass_chi2min, result.mass_mean, ps=8, xr=[8,12.5], yr=[8,12.5], xsty=3, ysty=3
;   plot, result.mass_median, result.mass_mean, ps=8, xr=[8,12.5], yr=[8,12.5], xsty=3, ysty=3
;   oplot, !x.crange, !y.crange
;   plothist, result.ebv_chi2min, bin=0.1
    
return
end

;   for j=0,n_elements(result)-1 do for i=0,n_tags(result)-1 do if size(result[j].(i),/type) ne 7L then if total(finite(result[j].(i)) eq 0B) gt 0.0 then print, j, i
