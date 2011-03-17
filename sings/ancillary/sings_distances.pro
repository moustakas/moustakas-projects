;+
; NAME:
;       SINGS_DISTANCES
;
; PURPOSE:
;       Compile distances for the SINGS galaxies.
;
; CALLING SEQUENCE:
;       sings_distances, distdata, /write
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
;       J. Moustakas, 2005 Jul 25, U of A
;       jm07jun31nyu - major updates to the distances
;-

pro sings_distances, distdata, write=write

    light = 2.99792458D5 ; speed of light [km/s]

    analysis_path = sings_path(/analysis)

    sings = mrdfits(analysis_path+'sings_ned.fits.gz',1,/silent)
    info = rsex(analysis_path+'sings_ancillary_data.sex')
    ngalaxy = n_elements(sings)
    
    ra = 15.0*im_hms2dec(sings.ra)
    dec = im_hms2dec(sings.dec)

    finaldist = rsex(analysis_path+'sings_distances_v1.2.dat')
    srt = sort(sings[sort(sings.galaxy)].ra) ; alphabetize and then sort by RA
    finaldist = finaldist[srt]
    
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
      murphydist:       -999.0, $ ; distance from E. Murphy
      leedist:          -999.0, $ ; distance from J. Lee (for a restricted set of SINGS galaxies)
      kenn03dist:       -999.0, $ ; distance from Kennicutt et al. 2003
      z_cosmic:         -999.0, $ ; cosmic (Hubble-flow) redshift based on the proper distance
      distance:         -999.0, $ ; final distance
      distance_err:     -999.0, $ ; error
      distance_ref:        ' ', $ ; reference
      distance_texref:     '-', $ ; bibtex reference
      distance_method:     ' '}   ; distance method
    distdata = replicate(distdata,ngalaxy)

    distdata.galaxy = strtrim(sings.galaxy,2)
    distdata.ra     = sings.ra
    distdata.dec    = sings.dec

    goodz = where(sings.z gt -900)
    distdata[goodz].cz = sings[goodz].z*light

; cross-match the distance catalog with the sings galaxies; use a 20"
; positional search radius; search relative to the distance catalog
; rather than the atlas because we want all the multiple matches

    cat = read_distance_catalog()
    raref = 15.0*im_hms2dec(cat.ra)
    decref = im_hms2dec(cat.dec)
    ra = 15.0*im_hms2dec(sings.ra)
    dec = im_hms2dec(sings.dec)

    splog, 'Searching the distance catalog.'
    ntot = im_djs_angle_match(raref,decref,ra,dec,dtheta=20.0/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist,mmax=1)
    good = where(mindx ne -1L,ngood)

;   w = where(strmatch(cat[good].reference,'*shapley*',/fold))
;   good = good[w]
;   
;   srt = sort(mdist[good])
;   niceprint, sings[mindx[good[srt]]].galaxy, finaldist[mindx[good[srt]]].galaxy, $
;     cat[good[srt]].galaxy, cat[good[srt]].distance, finaldist[mindx[good[srt]]].dist, $
;     mdist[good[srt]]*3600.0, cat[good[srt]].reference

; ---------------------------------------------------------------------------
; store the master set of distances
; ---------------------------------------------------------------------------

    ref = finaldist.ref
    ref = repstr(ref,'K03','kara03a')
    ref = repstr(ref,'Ma98','makarova98a')
    ref = repstr(ref,'KLM05','masters05a')
    ref = repstr(ref,'T01','tonry01a')
    ref = repstr(ref,'S96','sharina96a')
    ref = repstr(ref,'F01','freedman01a')
    ref = repstr(ref,'Fe00a','ferrarese00a')
    ref = repstr(ref,'Fe00b','ferrarese00b')
    ref = repstr(ref,'TSB01','tosi01a')
    ref = repstr(ref,'M01','macri01a')
    ref = repstr(ref,'L02','leonard02a')
    ref = repstr(ref,'Jha07','jha07a')
    ref = repstr(ref,'Mei07','mei07a')
    ref = repstr(ref,'S05','seth05a')
    ref = repstr(ref,'DK00','drozdovsky00a')
    ref = repstr(ref,'PG04','pietrzynski04a')
    ref = repstr(ref,'KSH00','kara00a')
;   ref = repstr(ref,'FCJ97','feldmeier97a')
    ref = repstr(ref,'C02','ciardullo02a')
    
    distdata.litdist        = finaldist.dist
    distdata.litdist_err    = finaldist.dist_err
    distdata.litdist_method = finaldist.method
    distdata.litdist_ref    = ref
    distdata.litdist_texref = ref
    
    distdata.murphydist = info.distance_murphy
    distdata.leedist = info.distance_lee
    distdata.kenn03dist = info.distance

; ---------------------------------------------------------------------------
; compute model distances for all objects
; ---------------------------------------------------------------------------

    czgood = where(distdata.cz gt -900.0,nczgood,comp=czbad,ncomp=nczbad)
    mould = mould_distance(distdata[czgood].ra,distdata[czgood].dec,$
      distdata[czgood].cz,object=distdata[czgood].galaxy,/proper)
    
    distdata[czgood].modeldist = mould.distance

; ---------------------------------------------------------------------------
; assign final distances and distance errors
; ---------------------------------------------------------------------------

    litdist = where(distdata.litdist gt -900.0,ndist,comp=modeldist,ncomp=nmodeldist)
    if (ndist ne 0L) then begin
       distdata[litdist].distance        = distdata[litdist].litdist
       distdata[litdist].distance_err    = distdata[litdist].litdist_err
       distdata[litdist].distance_ref    = distdata[litdist].litdist_ref
       distdata[litdist].distance_texref = distdata[litdist].litdist_texref
       distdata[litdist].distance_method = distdata[litdist].litdist_method
    endif
    
    alpha = 0.006 ; [1/Mpc]

    if (nmodeldist ne 0L) then begin
       distdata[modeldist].distance        = distdata[modeldist].modeldist
       distdata[modeldist].distance_err    = distdata[modeldist].modeldist*0.13*exp(-alpha*distdata[modeldist].modeldist)
       distdata[modeldist].distance_ref = 'Mould et al. 2000'
       distdata[modeldist].distance_method = 'Infall Model'
    endif

    struct_print, struct_trimtags(distdata,select=['galaxy',$;'ra','dec','cz',$
      'distance','modeldist','murphydist','kenn03dist','leedist','distance_texref'])

    if keyword_set(write) then begin
       outfile = 'sings_distances.fits'
       splog, 'Writing '+analysis_path+outfile
       mwrfits, distdata, analysis_path+outfile, /create
       spawn, ['gzip -f '+analysis_path+outfile], /sh
    endif

return
end    
