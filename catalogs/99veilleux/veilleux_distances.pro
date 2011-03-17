;+
; NAME:
;       VEILLEUX_DISTANCES
;
; PURPOSE:
;       Compile distances for the VEILLEUX galaxies.
;
; CALLING SEQUENCE:
;       veilleux_distances, distdata, /write
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
;       J. Moustakas, 2006 Feb 09, U of A
;-

pro veilleux_distances, distdata, write=write

    light = 2.99792458D5 ; speed of light [km/s]

    veilleux_path = getenv('CATALOGS_DIR')+'/99veilleux/'

    veilleux95 = mrdfits(veilleux_path+'veilleux95_ned.fits.gz',1,/silent)
    veilleux99 = mrdfits(veilleux_path+'veilleux99_ned.fits.gz',1,/silent)
    veilleux = struct_append(veilleux95,veilleux99)
    ngalaxy = n_elements(veilleux)

    ra = 15.0*im_hms2dec(veilleux.ra)
    dec = im_hms2dec(veilleux.dec)

; initialize the output data structure

    distdata = {$
      galaxy:              ' ', $
      ra:                  ' ', $
      dec:                 ' ', $
      cz:                -999.0,$ ; redshift
      litdist:           -999.0,$ ; literature distance
      litdist_err:       -999.0,$ ; error
      litdist_method:       '-',$ ; direct distance method
      litdist_ref:          '-',$ ; reference
      litdist_texref:       '-',$ ; bibtex reference
      modeldist:        -999.0, $ ; model distance
      z_cosmic:         -999.0, $ ; cosmic (Hubble-flow) redshift based on the proper distance
      distance:         -999.0, $ ; final distance
      distance_err:     -999.0, $ ; error
      distance_ref:        ' ', $ ; reference
      distance_texref:     '-', $ ; bibtex reference
      distance_method:     ' '}   ; distance method
    distdata = replicate(distdata,ngalaxy)

    distdata.galaxy = strtrim(veilleux.galaxy,2)
    distdata.ra = veilleux.ra
    distdata.dec = veilleux.dec

    goodz = where(veilleux.z gt -900)
    distdata[goodz].cz = veilleux[goodz].z*light

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
       distdata[modeldist].distance_ref = 'Mould et al. 2000'
       distdata[modeldist].distance_method = 'Infall Model'
    endif

    srt = sort(distdata.ra)
    distdata = distdata[srt]
    
    struct_print, struct_trimtags(distdata,select=['galaxy','ra','dec','cz','distance',$
      'distance_err','distance_ref','distance_method'])

    if keyword_set(write) then begin
       outfile = 'veilleux_distances.fits'
       splog, 'Writing '+veilleux_path+outfile
       mwrfits, distdata, veilleux_path+outfile, /create
       spawn, ['gzip -f '+veilleux_path+outfile], /sh
    endif

return
end    
