;+
; NAME:
;       NFGS_STELLAR_MASS
;
; PURPOSE:
;       Estimate stellar masses using ISEDFIT.
;
; COMMENTS:
;       Use the emission-line corrected magnitudes.  
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 Mar 31, U of A
;-

pro nfgs_stellar_mass, nfgs, write=write

    datapath = nfgs_path(/analysis)
    if (n_elements(nfgs) eq 0L) then nfgs = read_nfgs()
    ngalaxy = n_elements(nfgs)
    galaxy = strtrim(nfgs.galaxy,2)

    red, h100=0.7, omega_0=0.3, omega_lambda=0.7
    H0 = 100.0*redh100() & omega0 = redomega0() & omegal = redomegal()
    q0 = omega0/2.0 - omegal
    light = 2.99792458D5        ; speed of light [km/s]

    zobj = nfgs.z_obj
;   zobj = H0*nfgs.distance/light

; ---------------------------------------------------------------------------
; photometric bands of interest
; ---------------------------------------------------------------------------

    modelsprefix = 'nfgs_UBVRJHKs'
    filterlist = ['bessell_'+['U','B','V','R'],'twomass_'+['J','H','Ks']]+'.par'
    filtinfo = im_filterspecs(filterlist=filterlist,/verbose)
    nfilter = n_elements(filterlist)

    niceprint, nfgs.nfgs_id, nfgs.galaxy, nfgs.nfgs_galaxy, nfgs.u, nfgs.b, nfgs.v, $
      nfgs.r, nfgs.twomass_j, nfgs.twomass_j, nfgs.twomass_ks
    
; ---------------------------------------------------------------------------
; generate the model magnitudes for the UBVRJHKs photometry
; ---------------------------------------------------------------------------

;   tau = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,8.0,10.0,13.0,16.0,20.0]
;   isedfit_models, filterlist, datapath=datapath, modelsprefix=modelsprefix, $
;     minredshift=0.01D, maxredshift=0.05D, minage=0.1D, maxage=13.0D, $
;     dlogage=0.05D, Zstellar=[2L,4L], tau=tau, /write, debug=debug

; ---------------------------------------------------------------------------
; convert the photometry in each band to MAGGIES and do the fit!
; ---------------------------------------------------------------------------

    bands = ['u','b','v','r','twomass_'+['j','h','ks']]
    maggies = convert_to_maggies(nfgs,bands,filtinfo.vega2ab,maggies_invvar=maggies_invvar)

    thesefilters = filterlist
;   thesefilters = filterlist[0:3]

    massprefix = 'nfgs_UBVRJHKs'

;   fitindx = lindgen(ngalaxy)
    fitindx = lindgen(3)+76

    isedfit,maggies[*,fitindx],maggies_invvar[*,fitindx],zobj[fitindx],result,result_info,filterlist=thesefilters,$
      datapath=datapath,modelsprefix=modelsprefix,massprefix=massprefix,galaxy=galaxy[fitindx],/debug;,/write

    oisedfit,maggies[*,fitindx],maggies_invvar[*,fitindx],zobj[fitindx],oresult,oresult_info,filterlist=thesefilters,$
      datapath=datapath,modelsprefix=modelsprefix,massprefix=massprefix,galaxy=galaxy[fitindx];,/debug;,/write

; test with the models themselves    

    models1 = mrdfits(datapath+result_info.modelsprefix+'.fits.gz',1,/silent)
    redshift = result_info.redshift & nredshift = n_elements(redshift)
    models = 1.45E11*float(reform(models1.modelmaggies,nfilter,n_elements(models1)/nredshift,nredshift))
    indx = lindgen(6)+250 & redindx = 0L & ngal = n_elements(indx)
    ebv = transpose(rebin(reform(10^(-0.4*k_lambda(filtinfo.weff,/charlot)*0.0),1,nfilter),ngal,nfilter))
    mmaggies_err = 0.1*models[*,indx,redindx]
    mmaggies = models[*,indx,redindx] + randomn(seed,nfilter,ngal)*mmaggies_err*1.0

    isedfit,mmaggies*ebv,1.0/(mmaggies_err*ebv)^2,replicate(redshift[redindx],ngal),mresult,mresult_info,$
      filterlist=thesefilters,datapath=datapath,modelsprefix=modelsprefix,massprefix='test_'+massprefix,/debug ;,/write

stop    
    
    isedfit,maggies[*,fitindx],maggies_invvar[*,fitindx],zobj[fitindx],result,filterlist=thesefilters,$
      modelsprefix=modelsprefix,massprefix=massprefix,galaxy=galaxy[fitindx],/maxold,/write

stop    

    isedfit_qaplot, result, result_info, datapath=datapath, massprefix=massprefix, /postscript

stop    
    
    isedfit_measure, measure, massprefix=massprefix, /write
    
return
end
