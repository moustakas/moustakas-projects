pro atlas1d_demos, paper=paper
; jm03dec1uofa
; demonstrate various features of the template continuum fitting

    atlaspath = atlas_path(/atlas1d)
    specfitpath = atlas_path(/specfit)

    debug = 1
    pspath = ''
    if keyword_set(paper) then begin
       pspath = atlas_path(/papers)+'atlas/FIG_ATLAS/'
       debug = 0
    endif 
    
    atlas = read_integrated(/silent)

; ---------------------------------------------------------------------------    
; Balmer absorption
; ---------------------------------------------------------------------------    

;   psname_abs = 'demo_balmerabs.ps'
    psname = 'absorption_'+['hd','hg','hb','ha']+'.ps'
    
    galaxy = 'ngc2366'
    gatlas = atlas[speclinefit_locate(atlas,galaxy)]
    linefit = struct_trimtags(gatlas,select=[strtrim(gatlas[0].linename,2)+'*'])
    spec = read_atlas_specfit(galaxy,/silent)

    wave = spec[*,0]/(1+gatlas.z_obj) ; rest wavelength
    flux = spec[*,2]                  ; continuum fit
    ferr = flux*0.0+1.0               ; error spectrum
    eflux = spec[*,3]                 ; emission-line spectrum

    balmerwaves = [4101.7,4340.46,4861.33,6562.8]
    balmersigma = sqrt((balmerwaves*70.0/2.99E5)^2.0+(8.5/2.35)^2.0)
    nbalmer = n_elements(balmerwaves)

    for i = 0L, nbalmer-1L do begin
    
       im_openclose, pspath+psname[i], postscript=paper, /color
       balmerabs = ibalmerabs(wave,flux,ferr=ferr,balmerwaves=balmerwaves[i],$
         balmersigma=balmersigma[i],gsigma=gsigma,debug=debug,postscript=paper)
;      djs_oplot, wave, 1D15*(eflux+balmerabs[i].babs_continuum), ps=10, color='brown', thick=5.0
       im_openclose, postscript=paper, /close

    endfor
    
; ---------------------------------------------------------------------------    
; local continuum
; ---------------------------------------------------------------------------    

;   psname = 'demo_localcontinuum.ps'
    psname = 'local_continuum_'+['oii','oiii','nii','sii']+'.ps'
    
    galaxy = 'ngc4713' ; 'ngc4736'
    gatlas = atlas[speclinefit_locate(atlas,galaxy)]
    linefit = struct_trimtags(gatlas,select=[strtrim(gatlas[0].linename,2)+'*'])
    spec = read_atlas_specfit(galaxy,/silent)

    wave = spec[*,0]     ; rest wavelength
    flux = spec[*,1]     ; galaxy spectrum
    specfit = spec[*,3]  ; emission-line fit
    ferr = flux*0.0+1.0  ; error spectrum

    linename = ['OII_3727','OIII_5007','NII_6584','SII_6716']
    linewave = [gatlas.oii_3727_wave,gatlas.oiii_5007_wave,gatlas.nii_6584_wave,gatlas.sii_6716_wave]
    linesigmakms = [gatlas.oii_3727_sigma_total,gatlas.oiii_5007_sigma_total,$
      gatlas.nii_6584_sigma_total,gatlas.sii_6716_sigma_total]
    print, linesigmakms
    nline = n_elements(linename)

    for i = 0L, nline-1L do begin

       im_openclose, pspath+psname[i], postscript=paper, /color
       ilocal = ilocal_continuum(wave,flux,specfit,linename=linename[i],$
         linewave=linewave[i],linesigmakms=linesigmakms[i],glosigma=glosigma,$
         ghisigma=ghisigma,debug=debug,postscript=paper)
       im_openclose, postscript=paper, /close

    endfor
       
; ---------------------------------------------------------------------------    
; absorption-line fitting with and without dust reddening
; ---------------------------------------------------------------------------    

;   fnodust = ispeclinefit('ngc_3351_drift_360_200.ms.fits.gz',ndust=0,$
;     specres=8.5,/zcrosscor,specfit=snodust)
;   fdust = ispeclinefit('ngc_3351_drift_360_200.ms.fits.gz',$
;     specres=8.5,/zcrosscor,dustmodel=2,specfit=sdust)

stop

return
end    
