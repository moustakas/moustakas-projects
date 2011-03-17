;+
; NAME:
;       KENN92_DISTANCES
;
; PURPOSE:
;       Compile distances for the Kennicutt (1992) galaxies.
;
; CALLING SEQUENCE:
;       kenn92_distances, distdata, /write
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       write - write out
;
; OUTPUTS:
;       distdata - distances data structure
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
;       J. Moustakas, 2005 Aug 02, U of A
;-

pro kenn92_distances, distdata, write=write, postscript=postscript

    if keyword_set(write) then postscript = 1L
    
    light = 2.99792458D5 ; speed of light [km/s]
    
    analysis_path = kenn92_path(/analysis)

    red, h100=0.7, omega_0=0.3, omega_lambda=0.7
    H0 = redh100()*100.0
    omega0 = redomega0()
    omega_lambda = redomegal()
    
    q0 = omega0/2.0 - omega_lambda

    kenn92 = mrdfits(analysis_path+'kenn92_ned.fits.gz',1,/silent)
    ngalaxy = n_elements(kenn92)

    ra = 15.0*im_hms2dec(kenn92.ra)
    dec = im_hms2dec(kenn92.dec)

; initialize the output data structure

    distdata = {$
      galaxy:              ' ', $
;     nedgalaxy:           ' ', $
      ra:                  ' ', $
      dec:                 ' ', $
      cz:                -999.0,$ ; redshift
      litdist:           -999.0,$ ; literature distance
      litdist_err:       -999.0,$ ; error
      litdist_method:       '-',$ ; direct distance method
      litdist_ref:          '-',$ ; reference
      modeldist:        -999.0, $ ; model distance
      z_cosmic:         -999.0, $ ; cosmic (Hubble-flow) redshift based on the proper distance
      distance:         -999.0, $ ; final distance
      distance_err:     -999.0, $ ; error
      distance_ref:        ' ', $ ; reference
      distance_texref:     ' ', $ ; latex reference
      distance_method:     ' '}   ; distance method
    distdata = replicate(distdata,ngalaxy)

    distdata.galaxy = strtrim(kenn92.galaxy,2)
;   distdata.nedgalaxy = strtrim(kenn92.ned_galaxy,2)
    distdata.ra = kenn92.ra
    distdata.dec = kenn92.dec

    goodz = where(kenn92.z gt -900)
    distdata[goodz].cz = kenn92[goodz].z*light

; ---------------------------------------------------------------------------
; compute model distances for all objects
; ---------------------------------------------------------------------------

    czgood = where(distdata.cz gt -900.0,nczgood,comp=czbad,ncomp=nczbad)
    mould = mould_distance(distdata[czgood].ra,distdata[czgood].dec,$
      distdata[czgood].cz,object=distdata[czgood].galaxy,/proper)
    
    distdata[czgood].modeldist = mould.distance

; ---------------------------------------------------------------------------
; assign final distances and distance errors (see ATLAS_DISTANCES) 
; ---------------------------------------------------------------------------

    litdist = where(distdata.litdist gt -900.0,ndist,comp=modeldist,ncomp=nmodeldist)
    if (ndist ne 0L) then begin
       distdata[litdist].distance        = distdata[litdist].litdist
       distdata[litdist].distance_err    = distdata[litdist].litdist_err
       distdata[litdist].distance_ref    = distdata[litdist].litdist_ref
       distdata[litdist].distance_method = distdata[litdist].litdist_method
    endif
    
    alpha = 0.006 ; [1/Mpc]

    if (nmodeldist ne 0L) then begin
       distdata[modeldist].distance        = distdata[modeldist].modeldist
       distdata[modeldist].distance_err    = distdata[modeldist].modeldist*0.13*exp(-alpha*distdata[modeldist].modeldist)
       distdata[modeldist].distance_ref    = 'Mould et al. 2000'
       distdata[modeldist].distance_texref = 'mould00'
       distdata[modeldist].distance_method = 'Infall Model'
    endif

    struct_print, struct_trimtags(distdata,select=['galaxy','ra','dec','cz','distance',$
      'distance_err','distance_ref','distance_method'])

    if keyword_set(write) then begin
       outfile = 'kenn92_distances.fits'
       splog, 'Writing '+analysis_path+outfile
       mwrfits, distdata, analysis_path+outfile, /create
       spawn, ['gzip -f '+analysis_path+outfile], /sh
    endif

return
end    
