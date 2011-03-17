;+
; NAME:
;       AGES_ISEDFIT_QAPLOT
;
; PURPOSE:
;       Generate QA plots of the results of fitting the AGES
;       photometry with ISEDFIT. 
;
; CALLING SEQUENCE:
;       ages_isedfit_qaplot, 
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Feb 12, U of A
;       jm05aug04uofa - updated for release 2.0
;-

pro ages_isedfit_qaplot, result, filterlist=filterlist, prefix=prefix, datapath=datapath, $
  spec1dpath=spec1dpath, make_png=make_png

    if (n_elements(datapath) eq 0L) then datapath = ages_path(/mass)
    if (n_elements(spec1dpath) eq 0L) then spec1dpath = ages_path(/spec1d)

    webpath = ages_path(/web)+'mass/'

    if (n_elements(prefix) eq 0L) then prefix = 'ages_BwRIKugriz'
    allprefix = 'ages_BwRIKugriz'

; define some constants, the model SED path names, and the
; cosmological model 
    
    lsun = 3.826D33             ; [erg/s]
    light = 2.99792458D5        ; speed of light [km/s]
    
    rootpath = getenv('CATALOGS_DIR')+'/sfhgrid/' ; parameter I/O path
    modelspath = rootpath+'basemodels/'           ; base models output path

    red, h100=0.7, omega_0=0.3, omega_lambda=0.7 ; cosmological parameters
    
; restore the output from ISEDFIT_MODELS; if requested, only use a
; subset of the filters

    modelfile_info = allprefix+'_models_info.fits.gz'
    if (file_test(datapath+modelfile_info,/regular) eq 0L) then begin
       splog, 'ISEDFIT output '+datapath+modelfile_info+' not found.'
       return
    endif

    splog, 'Restoring '+datapath+modelfile_info+'.'
    isedfit_info = mrdfits(datapath+modelfile_info,1,/silent)

    if (n_elements(filterlist) eq 0L) then filterlist = isedfit_info.filters

    match, strtrim(isedfit_info.filters,2), filterlist, filterindx, indx
    if (n_elements(indx) ne n_elements(filterlist)) then begin
       splog, 'FILTERLIST must be a subset of the filters used in ISEDFIT_MODELS.'
       return
    endif else begin

; sort the filter list       
       
       splog, 'Selecting the following filters: ' & niceprint, '    '+strupcase(filterlist)
       srt = sort(filterindx)
       filterindx = filterindx[srt]

; subscript the information structure to the selected filters

       isedfit_info = struct_addtags({filters: isedfit_info.filters[filterindx]},$
         struct_trimtags(isedfit_info,except='FILTERS'))
       nfilter = n_elements(isedfit_info.filters)

    endelse

    metallicity = isedfit_info.metallicity

    tau      = isedfit_info.tau      & ntau      = n_elements(tau)
    tage     = isedfit_info.tage     & ntage     = n_elements(tage)
    ebv      = isedfit_info.ebv      & nebv      = n_elements(ebv)
    Z        = isedfit_info.Z        & nZ        = n_elements(Z)
    redshift = isedfit_info.redshift & nredshift = n_elements(redshift)

; restore the output from ISEDFIT for the requested set of filters

    isedfitfile = prefix+'_isedfit_mass.fits.gz'
    if (file_test(datapath+isedfitfile,/regular) eq 0L) then begin
       splog, 'ISEDFIT output '+datapath+isedfitfile+' not found.'
       return
    endif

    splog, 'Restoring '+datapath+isedfitfile+'.'
    result = mrdfits(datapath+isedfitfile,1,/silent)
    ngalaxy = n_elements(result)

    galaxy = strtrim(result.galaxy,2)
    passfield = repstr(strmid(galaxy,5),'_','/')
    pngfile = galaxy+'_mass.png'
    specfile = galaxy+'.fits' ; Hectospec spectrum
    
; recover the filter specifications
    
    filters = strtrim(isedfit_info.filters,2)
    nfilter = n_elements(filters)

    filtinfo = im_filterspecs(filterlist=filters)
    k_load_filters, filters, filter_nlambda, filter_lambda, filter_pass
    
    dlum = dluminosity(result.zobj,/cm) ; luminosity distance
    fluxarea = 4.0*!dpi*dlum*dlum
    constant = 3.631D-20*light*1D13/filtinfo.weff^2

; restore the AGES information structure and match with RESULT 

    ages = ages_read_info()
;   match = cmset_op(strtrim(ages.galaxy,2),'and',galaxy,/index)
;   agesmatch = ages[match]
    agesmatch = ages
    ndwfsgalaxy = agesmatch.ndwfs_galaxy
    
; generate QA plots for each galaxy

    if (not keyword_set(make_png)) then begin
       im_window, 0, xratio=0.8, yratio=0.6
    endif

;   for igalaxy = 12, 15 do begin
    for igalaxy = 0L, ngalaxy-1L do begin

; PNG output - black background          

       if keyword_set(make_png) then begin
          print, format='("Object ",I0,"/",I0,".",A1,$)', igalaxy+1, ngalaxy, string(13b)
          set_plot, 'Z'
          polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')
       endif

       zobj = result[igalaxy].zobj
       
       ssp = im_read_bc03(isedfile=result[igalaxy].isedfit_sspfile,isedpath=modelspath,$
         age=result[igalaxy].isedfit_tage,minwave=isedfit_info.minwave,bc=sspextras,$
         maxwave=isedfit_info.maxwave,/silent)

       ssprestwave = ssp.wave
       ssprestflux = lsun*(10^result[igalaxy].isedfit_mass)*ssp.flux/sspextras.m_ ; rest-frame SSP [erg/s/cm2/A]

       klambda = k_lambda(ssprestwave,/charlot)
       ssprestflux_red = ssprestflux * 10^(-0.4*result[igalaxy].isedfit_ebv*klambda)

       sspwave = ssp.wave*(1+zobj)
       sspflux = ssprestflux_red/fluxarea[igalaxy]/(1+zobj)

       modelmaggies_flambda = reform(result[igalaxy].isedfit_modelmaggies)*10^(result[igalaxy].isedfit_mass)*constant
       maggies_flambda = reform(result[igalaxy].isedfit_maggies)*constant
       maggies_sigma_flambda = reform(result[igalaxy].isedfit_maggies_sigma)*constant

; make the plot

       factor = 1D17
       
;      xrange = [1000.0,max(sspwave/(1+zobj))]
       xrange = minmax(sspwave/(1+zobj)) ; minmax(sspwave)
       yrange = factor*minmax(sspflux*(1+zobj))*[1.0,1.2] ; factor*minmax(sspflux)*[1.0,1.2]
       
       djs_plot, [0], [0], /nodata, noerase=keyword_set(make_png), xsty=3, ysty=3, $
         charsize=1.5, /xlog, charthick=2.0, xthick=2.0, ythick=2.0, $
         xtitle='Rest Wavelength (\AA)', ytitle='Flux (10^{-17} '+flam_units()+')', $
         yrange=yrange, xrange=xrange
       djs_oplot, sspwave/(1+zobj), factor*sspflux*(1+zobj), ps=10, thick=1.0, color='grey'

; overplot the observed and model photometry

       plotsym, 0, 1.1, /fill
       good = where(maggies_flambda gt 0.0,ngood)
       if (ngood ne 0L) then $
         oploterror, filtinfo[good].weff/(1+zobj), factor*maggies_flambda[good]*(1+zobj), $
         filtinfo[good].fwhm/(1+zobj), factor*maggies_sigma_flambda[good]*(1+zobj), ps=8, $
         color=djs_icolor('red'), errcolor=djs_icolor('red'), errthick=2.0

       plotsym, 8, 1.1, /fill
       djs_oplot, filtinfo.weff/(1+zobj), factor*modelmaggies_flambda*(1+zobj), ps=8, color='green'

; overplot the filters

;      for ifilter = 0L, nfilter-1L do begin
;         filtw = filter_lambda[0L:filter_nlambda[0]-1L,ifilter]
;         filtf = filter_pass[0L:filter_nlambda[0]-1L,ifilter]
;         filtf = filtf/max(filtf)
;         djs_oplot, filtw/(1+zobj), filtf*(!y.crange[1]-!y.crange[0])*0.10+!y.crange[0], line=0, color='cyan'
;      endfor
       
;; overplot the AGES spectrum, scaled to the R-band photometry
;
;       scale = 10^(-0.4*(ages[igalaxy].R-ages[igalaxy].synth_R))
;       
;       specdata  = read_ages_specfit(agesgalaxy[igalaxy],plate[igalaxy],/silent) ; this will need to be modified!
;       restwave  = reform(specdata[*,0])
;       restflux  = reform(specdata[*,1])
;       restcflux = reform(specdata[*,2])
;       resteflux = reform(specdata[*,3])>0.0
;       
;       wave  = restwave*(1+zobj)
;       flux  = scale*restflux/(1+zobj)
;       cflux = scale*restcflux/(1+zobj)
;       eflux = scale*resteflux/(1+zobj)
;
;       djs_oplot, wave, smooth(flux,10), thick=1.0, color='blue'

; generate the legend       

       label = [$
         'log M = '+strtrim(string(result[igalaxy].isedfit_mass,format='(F8.2)'),2)+' M'+sunsymbol(),$
         'Z/Z'+sunsymbol()+' = '+string(result[igalaxy].isedfit_metallicity/0.02,format='(G0.0)'), $
;        'Z/Z'+sunsymbol()+' = '+string(result[igalaxy].isedfit_metallicity/0.02,format='(F5.3)'), $
         '\tau = '+string(result[igalaxy].isedfit_tau,format='(F4.1)')+' Gyr',$
         't = '+string(result[igalaxy].isedfit_tage,format='(F5.2)')+' Gyr',$
         'E(B-V)= '+string(result[igalaxy].isedfit_ebv,format='(F3.1)'),$
         '\chi^{2}_{\nu} = '+strtrim(string(result[igalaxy].isedfit_chi2min,format='(F12.2)'),2)]
       legend, textoidl(label), /left, /top, box=0, charsize=1.1, charthick=2.0, spacing=1.5

       label = [repstr(ndwfsgalaxy[igalaxy],'_',' '),passfield[igalaxy],'z = '+string(zobj,format='(F6.4)')]
;      label = [repstr(ndwfsgalaxy[igalaxy],'_',' '),'z = '+string(zobj,format='(F6.4)')]
;      label = [repstr(result[igalaxy].galaxy,'_',' '),'z = '+string(zobj,format='(F6.4)')]
       legend, label, /right, /top, box=0, charsize=1.1, charthick=2.0, spacing=1.5

       if keyword_set(make_png) then begin
          img = tvrd()          
          tvlct, r, g, b, /get
          write_png, webpath+pngfile[igalaxy], img, r, g, b
          set_plot, 'X'
       endif else begin
;         wait, 0.75
          cc = get_kbrd(1)
       endelse

    endfor

;      contour, reform(chi2[0].chi2surface[result[igalaxy].tauindx,*,*,result[igalaxy].metallicityindx]), $
;        tage, ebv, thick=2.0, charsize=1.8, charthick=2.0, xtitle='Age [Gyr]', ytitle='E(B-V)'
;      oplot, result[igalaxy].tage*[1,1], !y.crange, thick=2.0, line=2.0
;      oplot, !x.crange, result[igalaxy].ebv*[1,1], thick=2.0, line=2.0

return
end
    
