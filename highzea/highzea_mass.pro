pro highzea_mass

    path = highzea_path(/mass)
    photo = mrdfits(path+'hizea_sdss_photo.fit',1,/silent)
    
    datapath = path
    modelsprefix = 'highzea'
    chi2prefix = 'highzea'

    filterlist = ['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par','sdss_z0.par']
    filtinfo = im_filterspecs(filterlist=filterlist,/verbose)
    
;   tau = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,8.0,10.0,13.0,16.0,20.0]
;   isedfit_models, filterlist, datapath=datapath, modelsprefix=modelsprefix, $
;     minredshift=0.5D, maxredshift=0.8D, dredshift=0.05D, Zstellar=[2L,4L], $
;     tau=tau, minage=0.05D, maxage=8.5D, /allages, debug=0, write=1

    thesefilters = filterlist
    maggies = 10.0^(-0.4*photo.petrocounts)
    maggies_invvar = 1.0/(0.4*alog(10.0)*maggies*photo.petrocountserr)^2.0
    zobj = photo.zfinal
    galaxy = 'EA '+string(lindgen(n_elements(photo)),format='(I2.2)')

    fitindx = lindgen(n_elements(photo))

; compute and minimize chi2, then generate a QA plot

    isedfit_chi2,maggies[*,fitindx],maggies_invvar[*,fitindx],zobj[fitindx],chi2,chi2_info,filterlist=thesefilters,$
      modelsprefix=modelsprefix,chi2prefix=chi2prefix,galaxy=galaxy[fitindx],maxold=0,/write

    isedfit, result, chi2prefix=chi2prefix, maxold=0, /write
    
    isedfit_qaplot, datapath=datapath, chi2prefix=chi2prefix, maxold=0, /postscript

; maximally old    
    
    isedfit_chi2,maggies[*,fitindx],maggies_invvar[*,fitindx],zobj[fitindx],chi2,chi2_info,filterlist=thesefilters,$
      modelsprefix=modelsprefix,chi2prefix=chi2prefix,galaxy=galaxy[fitindx],maxold=1,/write

    isedfit, result, chi2prefix=chi2prefix, maxold=1, /write
    
    isedfit_qaplot, datapath=datapath, chi2prefix=chi2prefix, maxold=1, /postscript

return
end
    
