;+
; NAME:
;       ATLAS_DISTANCES
;
; PURPOSE:
;       Derive distances for the spectral atlas.
;
; CALLING SEQUENCE:
;       atlas_distances, /write, /postscript
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       write      - write the results
;       postscript - generate postscript output comparing direct and
;                    model-based distance estimates
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       ATLAS_PATH(), RED, MRDFITS(), IM_HMS2DEC(),
;       READ_DISTANCE_CATALOG(), IM_DJS_ANGLE_MATCH(), READ_96TULLY(),
;       IM_DMOD2D(), PARSE_LEDA(), MOULD_DISTANCE(), STRUCT_PRINT,
;       STRUCT_TRIMTAGS, SPLOG, MWRFITS, IM_OPENCLOSE, IM_STATS(),
;       PAGEMAKER, DJS_PLOT, DJS_OPLOT, LEGEND
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Feb 23, U of A - written
;       jm05jul21uofa - documented and updated
;-

pro atlas_distances, distdata, write=write, postscript=postscript

    if keyword_set(write) then postscript = 1L
    
    light = 2.99792458D5 ; speed of light [km/s]
    
    analysis_path = atlas_path(/analysis)

    red, h100=0.7, omega_0=0.3, omega_lambda=0.7
    H0 = redh100()*100.0
    omega0 = redomega0()
    omega_lambda = redomegal()
    
    q0 = omega0/2.0 - omega_lambda

    atlas = mrdfits(analysis_path+'atlas_ned.fits.gz',1,/silent) ; this must match LEDA!
    leda = atlas_read_leda()
    ngalaxy = n_elements(atlas)

    ra = 15.0*im_hms2dec(atlas.ra)
    dec = im_hms2dec(atlas.dec)

; initialize the output data structure

    distdata = {$
      galaxy:               '', $
;     nedgalaxy:            '', $
      ra:                   '', $
      dec:                  '', $
      cz:                -999.0,$ ; redshift
      litdist:           -999.0,$ ; literature distance
      litdist_err:       -999.0,$ ; error
      litdist_method:        '',$ ; direct distance method
      litdist_ref:           '',$ ; reference
      litdist_texref:        '',$ ; bibtex reference
      modeldist:        -999.0, $ ; model distance
      z_cosmic:         -999.0, $ ; cosmic (Hubble-flow) redshift based on the proper distance
      distance:         -999.0, $ ; final distance
      distance_err:     -999.0, $ ; error
      distance_ref:         '', $ ; reference
      distance_texref:      '', $ ; bibtex reference
      distance_method:      ''}   ; distance method
    distdata = replicate(distdata,ngalaxy)

    distdata.galaxy = strtrim(atlas.galaxy,2)
;   distdata.nedgalaxy = strtrim(atlas.ned_galaxy,2)
    distdata.ra = atlas.ra
    distdata.dec = atlas.dec

    goodz = where(atlas.z gt -900)
;   distdata.cz = atlas.cz
    distdata[goodz].cz = atlas[goodz].z*light

; read the distance catalog

    cat = read_distance_catalog()
    raref = 15.0*im_hms2dec(cat.ra)
    decref = im_hms2dec(cat.dec)

; cross-match the distance catalog with the atlas galaxies; use a 20"
; positional search radius; search relative to the distance catalog
; rather than the atlas because we want all the multiple matches 

    splog, 'Searching the distance catalog.'
    ntot = im_djs_angle_match(raref,decref,ra,dec,dtheta=20.0/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist,mmax=1)
    good = where(mindx ne -1L,ngood)

;   srt = sort(mdist[good])
;   niceprint, atlas[mindx[good[srt]]].galaxy, cat[good[srt]].galaxy, $
;     cat[good[srt]].distance, mdist[good[srt]]*3600.0, cat[good[srt]].reference

    atlasgood = atlas[mindx[good]]
    catgood = cat[good]
;   niceprint, atlas[mindx[good]].galaxy, catgood.galaxy

    galaxy = strtrim(atlasgood.galaxy,2)
    ugalaxy = galaxy[uniq(galaxy,sort(galaxy))]

    atlasgalaxy = strtrim(atlas.galaxy,2)

; Notes on distances: Freedman Cepheid distances are always preferred
; over Ferrarese Cepheid distances    
    
; Notes on multiple distances: 
;   NGC1569 - Shapley gets 1.7 while Karachentsev gets 1.95 Mpc.  The
;             Shapley result was based on Karachentsev et al. (1997),
;             while the new Karachentsev result was based on Makarova
;             & Karachentsev (2003), so take that one.
;
;   NGC2366 - Karachentsev gets 3.19 Mpc from TRGB while Shapley gets
;             3.44 Mpc from the Cepheid distance by Tolstoy et
;             al. (1995)
;


    nunique = n_elements(ugalaxy)
    for i = 0L, nunique-1L do begin

       match = where(galaxy eq ugalaxy[i],nmatch)

; if there is a single match then take that       

       if (nmatch eq 1L) then begin

          bigmatch = where(atlasgalaxy eq galaxy[match[0]],nbigmatch)
          
          distdata[bigmatch].litdist        = catgood[match[0]].distance
          distdata[bigmatch].litdist_err    = catgood[match[0]].distance_err
          distdata[bigmatch].litdist_method = catgood[match[0]].method
          distdata[bigmatch].litdist_ref    = catgood[match[0]].reference
          distdata[bigmatch].litdist_texref = catgood[match[0]].texref

       endif
       
; if there are multiple matches then take the distance according to
; the following priority: Freedman, Ferrarese, Karachentsev, Tonry,
; Whiting 

       if (nmatch gt 1L) then begin

          refs = catgood[match].reference
          methods = catgood[match].method

          indx = where(strmatch(refs,'*Freedman*',/fold) eq 1B,nindx)
;         if (nindx eq 0L) then indx = where(strmatch(refs,'*Ferrarese*',/fold) eq 1B,nindx)
          if (nindx eq 0L) then indx = where(strmatch(refs,'*Karachentsev*',/fold) eq 1B,nindx)
          if (nindx eq 0L) then indx = where(strmatch(refs,'*Tonry*',/fold) eq 1B,nindx)
;         if (nindx eq 0L) then indx = where(strmatch(refs,'*Whiting*',/fold) eq 1B,nindx)
          if (nindx eq 0L) then indx = where(strmatch(refs,'*Shapley*',/fold) eq 1B,nindx)
          if (nindx ne 1L) then message, 'Problem here!'

          bigmatch = where(atlasgalaxy eq galaxy[match[indx[0]]],nbigmatch)

          distdata[bigmatch].litdist        = catgood[match[indx[0]]].distance
          distdata[bigmatch].litdist_err    = catgood[match[indx[0]]].distance_err
          distdata[bigmatch].litdist_method = catgood[match[indx[0]]].method
          distdata[bigmatch].litdist_ref    = catgood[match[indx[0]]].reference
          distdata[bigmatch].litdist_texref = catgood[match[indx[0]]].texref

;         print, string(i,format='(I2)')+'  ', ugalaxy[i], '   '+distdata[bigmatch].galaxy, $
          print, '---------------------------------------------------------------------------'
          print, ugalaxy[i], '   '+distdata[bigmatch].galaxy, $
            distdata[bigmatch].litdist, '   '+distdata[bigmatch].litdist_ref
          niceprint, catgood[match].galaxy, catgood[match].distance, catgood[match].method, catgood[match].reference
;         print, '---------------------------------------------------------------------------'
          print
;         cc = get_kbrd(1)

       endif
       
    endfor

    splog, 'Found '+string(nunique,format='(I3)')+'/'+string(ngalaxy,format='(I3)')+' literature distances.'

; ---------------------------------------------------------------------------
; Ursa Major cluster distances (Tully et al. 1996, AJ, 112, 2471);
; reference for the distance to Ursa Major of 19.8 Mpc is Freedman et
; al. 2001 (compare with 18.6 Mpc from Tully et al. 2000, ApJ, 533,
; 744).  the positions are very good, so use a 10" search radius.  the
; distance error is taken as the "depth" of the cluster, which is
; quoted by Verheijen 2001 as 0.17 mag (=1.6 Mpc)
; ---------------------------------------------------------------------------

    print, format='("Reading the Ursa Major cluster data . . . ",$)'
    umajor = read_96tully()

    ura = 15.0*im_hms2dec(umajor.ra)
    ude = im_hms2dec(umajor.dec)

    ntot = im_djs_angle_match(ra,dec,ura,ude,dtheta=10.0/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    good = where(mindx ne -1L,ngood)
    print, 'found '+string(ngood,format='(I3)')+'/'+string(ngalaxy,format='(I3)')+' galaxies.'

    srt = sort(mdist[good])
    niceprint, distdata[good[srt]].galaxy, umajor[mindx[good[srt]]].galaxy, $
      mdist[good[srt]]*3600.0

; what to do if a literature distance already exists?  here's Janice's
; thinking: "It would depend on the error in the direct distance
; estimate.  If it is smaller than the size of the cluster as inferred
; from its velocity dispersion, then I would take the direct distance
; estimate *if* it is consistent with the cluster membership distance
; and the size of the cluster.  Otherwise I would just take the
; cluster membership distance."

; Having reviewed the literature, I decided that the direct distances
; below are an underestimate of the new Cepheid-calibrated Key Project
; distance to the Ursa Major cluster.  So don't replace any of these.

    yesdist = where(distdata[good].litdist gt -900.0,nyesdist,comp=nodist,ncomp=nnodist)
    if (nyesdist ne 0L) then begin

       niceprint, distdata[good[yesdist]].galaxy, distdata[good[yesdist]].litdist, $
         distdata[good[yesdist]].litdist_err, distdata[good[yesdist]].litdist_ref, $
         distdata[good[yesdist]].litdist_method

;      replace = where((strmatch(distdata[good[yesdist]].litdist_ref,'*shapley*',/fold) eq 1B) or $
;        (strmatch(distdata[good[yesdist]].litdist_ref,'*whiting*',/fold) eq 1B),nreplace)
;      if (nreplace ne 0L) then begin
;         distdata[good[yesdist[replace]]].litdist = 19.8
;         distdata[good[yesdist[replace]]].litdist_err = 1.6
;         distdata[good[yesdist[replace]]].litdist_ref = 'Ursa Major Cluster Membership'
;      endif

    endif

    if (nnodist ne 0L) then begin
       distdata[good[nodist]].litdist = 19.8
       distdata[good[nodist]].litdist_err = 1.6
       distdata[good[nodist]].litdist_ref = 'Tully et al. 1996'
       distdata[good[nodist]].litdist_texref = 'tully96'
       distdata[good[nodist]].litdist_method = 'UM Cluster Membership'
    endif

; ---------------------------------------------------------------------------
; Compile distances from LEDA; these distances are based on kinematic 
; parameters (see LEDA documentation): 
;   
;   logavmm = "The maximum velocity of rotation, logavmm, is expressed
; as the decimal logaritm of the rotational velocity in km/s.  It is a
; weighted average of the measurements from HI observations (logvm)
; with those from optical rotation curves (logvmha)."
;
;   btc = bt - ag - ai - ak.v/10000 = "The corrected apparent
; B-magnitude, where ak is taken from de Vaucouleurs et al (1976, RC2
; p33 rel.25) and v is the heliocentric radial velocity in km/s.
; When unknown, v and t are set to 0.  No calculation is made for btc
; when either bt or ag or ai is unknown.
;
; The final distance modulus is:
;
;     amu = btc + 6.5 logavmm + 6.3
; ---------------------------------------------------------------------------

; assume a fixed distance uncertainty
    
    nodist = where(distdata.litdist lt -900.0,nnodist)
    if (nnodist ne 0L) then begin

       good = where(leda[nodist].mup ne '',ngood)
       if (ngood ne 0L) then begin

          niceprint, leda[nodist[good]].objname, distdata[nodist[good]].galaxy

          dist = im_dmod2d(leda[nodist[good]].mup)
          dist_err = dist*0.10
          
;         distdata[nodist[good]].litdist = dist
;         distdata[nodist[good]].litdist_err = dist_err
;         distdata[nodist[good]].litdist_ref = 'LEDA'
;         distdata[nodist[good]].litdist_method = 'Kinematic Parameters'

       endif

    endif
    
; ---------------------------------------------------------------------------
; Tully (1988) NBG distances; assume a fixed distance uncertainty 
; ---------------------------------------------------------------------------

    nodist = where(distdata.litdist lt -900.0,nnodist)
    flow = flow_distance(distdata.ra,distdata.dec,H0=H0,galaxy=distdata.galaxy)

    good = where(flow[nodist].tully_dist gt -900,ngood)
    if (ngood ne 0L) then begin
    
       srt = sort(distdata[nodist[good]].ra)
       struct_print, flow[nodist[good[srt]]]

       dist = flow[nodist[good]].tully_dist
       dist_err = dist*0.10
       
;      distdata[nodist[good]].litdist = dist
;      distdata[nodist[good]].litdist_err = dist_err
;      distdata[nodist[good]].litdist_ref = 'Tully 1988'
;      distdata[nodist[good]].litdist_method = 'NBG'

    endif
    
; ---------------------------------------------------------------------------
; compute model distances for all objects
; ---------------------------------------------------------------------------

    czgood = where(distdata.cz gt -900.0,nczgood,comp=czbad,ncomp=nczbad)
    mould = mould_distance(distdata[czgood].ra,distdata[czgood].dec,$
      distdata[czgood].cz,object=distdata[czgood].galaxy,/proper,$
      H0=redh100()*100.0,omega0=redomega0(),omega_lambda=redomegal())

    distdata[czgood].modeldist = mould.distance

; ---------------------------------------------------------------------------
; what's left!?!
; ---------------------------------------------------------------------------

    nodist = where(distdata.litdist lt -900.0,comp=yesdist)
;   srt = sort(distdata[nodist].ra)
    srt = reverse(sort(distdata[nodist].cz))
    struct_print, struct_trimtags(distdata[nodist[srt]],select=['galaxy','ra','dec','cz','modeldist'])

    openw, lun, 'distances_needed.dat', /get_lun
    struct_print, struct_trimtags(distdata[nodist[srt]],select=['galaxy','nedgalaxy','ra','dec','cz','modeldist']), lun=lun
    free_lun, lun

; ---------------------------------------------------------------------------
; assign final distances and distance errors; choose an exponential
; model with an e-folding distance of 1/alpha Mpc.  the model
; parameters are fixed according to the numbers in Freedman et
; al. (2001), p. 55:
;
;   x = [3000.0D,30000.0]/70.0 & y = [0.10D,0.01]
;   myfunc = 'P[0]*exp(-P[1]*X)'
;   p = mpfitexpr(myfunc,x,y,yfit=yfit,/quiet)
;   print, p
;         0.12915497    0.0059696651
;
; choose P = [0.13,0.006]
;
; ---------------------------------------------------------------------------

    litdist = where(distdata.litdist gt -900.0,ndist,comp=modeldist)
    distdata[litdist].distance        = distdata[litdist].litdist
    distdata[litdist].distance_err    = distdata[litdist].litdist_err
    distdata[litdist].distance_ref    = distdata[litdist].litdist_ref
    distdata[litdist].distance_texref = distdata[litdist].litdist_texref
    distdata[litdist].distance_method = distdata[litdist].litdist_method
    
    alpha = 0.006 ; [1/Mpc]
    
    distdata[modeldist].distance        = distdata[modeldist].modeldist
    distdata[modeldist].distance_err    = distdata[modeldist].modeldist*0.13*exp(-alpha*distdata[modeldist].modeldist)
    distdata[modeldist].distance_ref = 'Mould et al. 2000'
    distdata[modeldist].distance_texref = 'mould00'
    distdata[modeldist].distance_method = 'Infall Model'
    struct_print, struct_trimtags(distdata,select=['galaxy','ra','dec','cz','distance',$
      'distance_err','distance_ref','distance_method'])

; finally compute the cosmic (Hubble-flow) redshift using the final
; set of proper distances

;   dist = distdata.distance
;   distdata.z_cosmic = 2.0/(1+q0) - sqrt(1.0/(1+q0)^2 - (2*H0*dist)/(light*(1+q0)))

    if keyword_set(write) then begin

       outfile = 'atlas_distances.fits'
       splog, 'Writing '+analysis_path+outfile
       mwrfits, distdata, analysis_path+outfile, /create
       spawn, ['gzip -f '+analysis_path+outfile], /sh

    endif

; ---------------------------------------------------------------------------
; generate QA plots
; ---------------------------------------------------------------------------

    if keyword_set(postscript) then begin
       postthick = 8.0 
    endif else begin
       im_window, 0, /square
       postthick = 2.0
    endelse
       
; --------------------------------------------------
; compare the infall model and literature distances  
; --------------------------------------------------

    psname = 'dist_lit_vs_dist_model.ps'
    im_openclose, analysis_path+psname, postscript=postscript
    
    good = where((distdata.litdist gt -900.0) and (distdata.modeldist gt -900),ngood)

    x = distdata[good].litdist
    y = distdata[good].modeldist

    residuals = 100*(x-y)/x

    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    xtitle = 'Literature Distance [Mpc]'
    ytitle = 'Infall Model Distance [Mpc]'
    ytitle2 = 'Residuals [%]'

    xrange = [min(x)<min(y),max(x)>max(y)]*[0.95,1.05]
    yrange = xrange
    
    xrange2 = xrange
    yrange2 = max(abs(residuals))*[-1.05,1.05]

    plotsym, 0, 1, /fill
    
    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.5,0.3], $
      ymargin=[0.3,1.3], position=pos, height=[6.0,3.4], /normal

    djs_plot, x, y, ps=8, xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, $
      position=pos[*,0], charsize=1.8, xtickname=replicate(' ',10), xsty=3, ysty=3, $
      charthick=postthick, xthick=postthick, ythick=postthick
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    
    djs_plot, x, residuals, ps=8, xtitle=xtitle, ytitle=ytitle2, xrange=xrange2, $
      yrange=yrange2, position=pos[*,1], /noerase, charsize=1.8, xsty=3, ysty=3, $
      charthick=postthick, xthick=postthick, ythick=postthick
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=1.5, $
      charthick=postthick, spacing=0
    
    im_openclose, postscript=postscript, /close    

; --------------------------------------------------
; compare the infall model and Tully distances
; --------------------------------------------------

    psname = 'dist_tully_vs_dist_model.ps'
    im_openclose, analysis_path+psname, postscript=postscript
    
    good = where((flow.tully_dist gt -900.0) and (distdata.modeldist gt -900) and $
      (distdata.litdist lt -900),ngood)

    x = flow[good].tully_dist
    y = distdata[good].modeldist

    residuals = 100*(x-y)/x

    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    xtitle = 'Tully Distance [Mpc]'
    ytitle = 'Infall Model Distance [Mpc]'
    ytitle2 = 'Residuals [%]'

    xrange = [min(x)<min(y),max(x)>max(y)]*[0.95,1.05]
    yrange = xrange
    
    xrange2 = xrange
    yrange2 = max(abs(residuals))*[-1.05,1.05]

    plotsym, 0, 1, /fill
    
    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.5,0.3], $
      ymargin=[0.3,1.3], position=pos, height=[6.0,3.4], /normal

    djs_plot, x, y, ps=8, xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, $
      position=pos[*,0], charsize=1.8, xtickname=replicate(' ',10), xsty=3, ysty=3, $
      charthick=postthick, xthick=postthick, ythick=postthick
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    
    djs_plot, x, residuals, ps=8, xtitle=xtitle, ytitle=ytitle2, xrange=xrange2, $
      yrange=yrange2, position=pos[*,1], /noerase, charsize=1.8, xsty=3, ysty=3, $
      charthick=postthick, xthick=postthick, ythick=postthick
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=1.5, $
      charthick=postthick, spacing=0

    im_openclose, postscript=postscript, /close    
    
stop    

return
end    
