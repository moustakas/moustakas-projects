pro model_grids, filters=filters, templates=templates, noforest=noforest
;+
; NAME:
;	MODEL_GRIDS
;
; PURPOSE:
;	Generate grids of observations as a function of filters and
;	redshift for all the available SED models and write them as
;	IDL save sets.  The save sets will act as lookup tables when
;	the grids need to be compared to simulated galaxy observations
;	or to photometric catalogs.
;
; CALLING SEQUENCE:
;	model_grids, [filters=], [templates=], noforest=noforest
;
; INPUTS:
;
; OPTIONAL INPUTS:
;	filters   - filters in which to calculate model fluxes
;                   (default to all the filters in the library)
;	templates - SEDs to "observe" (default to sirtf.templates)
;	
; KEYWORD PARAMETERS:
;	noforest  - do not apply the Ly-forest correction
;	
; OUTPUTS:
;	An IDL save set is written in the lib subdirectory, in
;	addition to an informational structure summarizing the
;	parameters of the observation (redshift grid size, filters
;	used, etc.)
;
; COMMON BLOCKS:
;	sirtf_simulations
;
; COMMENTS:
;	Depending on the redshift binning and the number of filters,
;	this routine may take as long as thirty minutes.  
;
;	Does not include internal reddening.  Should we put
;	interpolation of colors elsewhere?  Add logarithmic redshift
;	binning.  
;	
;	Write this as a damn C code!
;
; PROCEDURES USED:
;	FILTER_MATCH(), BAND_FLUX(), CMSAVE, STRN()
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 May 7, U of A
;	jm01aug1uofa, documented
;	jm01aug15uofa, added interpolation between colors and IGM
;	extinction
;-

    tstart = systime(/seconds)

    common sirtf_simulations

    if not keyword_set(templates) then templates = sirtf.templates ; SEDs
    
; all filters [default].  check for out of range filters

    if not keyword_set(filters) then filters = sirtf.bandcube.bandnames 

    obands = filter_match(filters,sirtf.bandcube.bandnames)
    good = where(sirtf.bandcube[obands].lambda_eff*1E-4 lt max(*sirtf.sedcube[0].lambda),ngood,comp=bad,ncomp=nbad)

    if nbad ne 0L then splog, 'The following filters were out of range: ', filters[bad]
    obands = obands[good]
    nbands = n_elements(obands)
    
    dz = sirtf.redshift.dz
    zmax = max(*sirtf.redshift.zarray)
    zmin = dz
    zarray = (findgen((zmax-zmin+dz)/dz))*dz+zmin
    nz = n_elements(zarray)

    nseds = n_elements(sirtf.sedcube)
    kcorrection = dblarr(nseds,nz,nbands) ; for computing colors
    colorflux = dblarr(nseds,nz,nbands)   ; for computing colors
    filterflux = dblarr(nseds,nz,nbands)  ; for computing photometric redshifts

; read in the lyman-forest depletion corrections

    path = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='lib')
    if not keyword_set(noforest) then begin
    
       print, 'Reading the Lyman forest depletion table.'
       readfast, path+'lyforest.dat', lyfdat, skip=2, nlines=nlyf
       d_a = interpol(lyfdat[1,*],lyfdat[0,*],zarray) ; interpolate onto the defined redshift grid
       d_b = interpol(lyfdat[2,*],lyfdat[0,*],zarray)

    endif
       
    light = 2.99792458D14                  ; speed of light [micron/s]
    dlum = dluminosity(zarray)*3.085678D16 ; luminosity distance (m)
    
    for k = 0L, nbands-1L do begin   ; loop on each filter

       wband = *sirtf.bandcube[obands[k]].wband*1E-4  ; response function wavelength [micron]
       rband = *sirtf.bandcube[obands[k]].rband       ; response function
       rband = rband/max(rband)

       if sirtf.bandcube[obands[k]].surface eq float(0) then surface = 1.0 else $
         surface = sirtf.bandcube[obands[k]].surface*1E-4
       lambda_eff = sirtf.bandcube[obands[k]].lambda_eff*1E-4 ; effective wavelength [micron]
       flambda_2_fnu = lambda_eff*lambda_eff/light            ; [W/micron -> W/Hz]
       
       for i = 0L, nseds-1L do begin ; loop on each model SED

          wave = reform(*sirtf.sedcube[i].lambda)     ; SED wavelength vector [micron]
          mlum_nu = reform(*sirtf.sedcube[i].mlum_nu) ; SED L_nu or f_nu
          mlum_mu = reform(*sirtf.sedcube[i].mlum_mu) ; SED L_mu or f_mu
          
; apply the Lyman forest corrections to the appropriate wavelengths
             
          if not keyword_set(noforest) then begin

             lya = where((wave lt 1215.67*1E-4) and (wave gt 1025.72*1E-4),nlya)
             lyb = where((wave lt 1025.72*1E-4) and (wave gt 911.8*1E-4),nlyb)
             lylimit = where(wave lt 911.8*1E-4,nlylimit)
                
          endif

          lumrest = interpol(mlum_nu,wave,sirtf.bandcube[obands[k]].lambda_eff*1E-4)
          restflux = band_flux(wave,mlum_mu,wband,rband,z=0D) ; flux at z=0
          
          for j = 0L, nz-1L do begin ; loop on redshift

             if not keyword_set(noforest) then begin

                if nlya ne 0L then mlum_mu[lya] = mlum_mu[lya] * d_a[j]
                if nlyb ne 0L then mlum_mu[lyb] = mlum_mu[lyb] * d_b[j]
                if nlylimit ne 0L then mlum_mu[lylimit] = 0.0
                
             endif
             
             zflux = band_flux(wave,mlum_mu,wband,rband,z=zarray[j]) ; flux at z

             filterflux[i,j,k] = flambda_2_fnu * zflux / surface

             if restflux gt float(0) then kcorrection[i,j,k] = zflux / restflux; / (1+zarray[j])
             colorflux[i,j,k] = lumrest * kcorrection[i,j,k]

             print, format='("SED ",i0,"/",i0,", filter ",i0,"/",i0,".",a1,$)', $
               i+1, nseds, k, nbands, string(13b)

          endfor ; loop on redshift
          
       endfor    ; loop on SED type

    endfor       ; loop on filter

    modelgrid = {date       : strmid(systime(),4,20), $
                 zarray     : zarray,  $
                 filters    : filters, $
                 filterflux : filterflux, $
                 colorflux  : colorflux, $
                 kcorrection: kcorrection}

    infoname = path+strlowcase(templates)+'.idlsave'
    cmsave, filename=infoname, modelgrid  
    
    print & print, 'Elapsed time: '+strn((systime(/seconds)-tstart)/60)+' minutes.'
    
return
end
