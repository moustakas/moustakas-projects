pro gto_distances, distdata, write=write, postscript=postscript
; jm06dec12nyu - written

    if keyword_set(write) then postscript = 1L
    
    light = 2.99792458D5 ; speed of light [km/s]
    
    analysis_path = gto_path(/ancillary)

    red, h100=0.7, omega_0=0.3, omega_lambda=0.7
    H0 = redh100()*100.0
    omega0 = redomega0()
    omega_lambda = redomegal()
    
    q0 = omega0/2.0 - omega_lambda

    gto = mrdfits(analysis_path+'gto_ned.fits.gz',1,/silent)
    ngalaxy = n_elements(gto)

    ra = 15.0*im_hms2dec(gto.ra)
    dec = im_hms2dec(gto.dec)

; initialize the output data structure

    distdata = {$
      galaxy:               '', $
      ned_galaxy:           '', $
      ra:                   '', $
      dec:                  '', $
      cz:                -999.0,$ ; redshift
      litdist:           -999.0,$ ; literature distance
      litdist_err:       -999.0,$
      litdist_texref:     '...',$ ; reference
      litdist_ref:        '...',$ ; reference
      modeldist:        -999.0, $ ; model distance
;     modeldist_ref:     '...', $ ; model distance
      tullydist:        -999.0, $ ; model distance
;     tullydist_ref:     '...', $ ; model distance
      distance:         -999.0, $ ; final distance
      distance_err:     -999.0, $ ; final distance
      distance_ref:      '...', $ ; reference
      distance_texref:   '...', $ ; latex reference
      distance_method:   '...'}   ; distance method
    distdata = replicate(distdata,ngalaxy)

    distdata.galaxy = strtrim(gto.galaxy,2)
    distdata.ned_galaxy = strtrim(gto.nedgalaxy,2)
    distdata.ra  = gto.ra
    distdata.dec = gto.dec

    goodz = where(gto.z gt -900.0)
    distdata[goodz].cz = gto[goodz].z*light

; read and parse Chad's distance catalog

;   chad1 = rsex(analysis_path+'chad_distances.dat')
    chad1_txt = djs_readlines(analysis_path+'chad_distances.dat')
    chad1_txt = chad1_txt[where(strmatch(chad1_txt,'*#*') eq 0B)]
    nchad1 = n_elements(chad1_txt)

    chad1 = replicate(create_struct('galaxy', '', 'distance', 0.0, 'distance_err', $ 
      0.0, 'texref', '', 'ref', ''),nchad1)
    for ii = 0L, nchad1-1L do begin
       str = strsplit(chad1_txt[ii],' ',/extract)
       chad1[ii].galaxy       = str[0]
       chad1[ii].distance     = str[1]
       chad1[ii].distance_err = str[2]
       chad1[ii].texref       = str[3]
       chad1[ii].ref = strjoin(str[4:n_elements(str)-1],' ')
    endfor

;   niceprint, cmset_op(strtrim(chad1.galaxy,2),'and',/not2,strtrim(gto.galaxy,2))

    match, strtrim(strupcase(gto.galaxy),2), strtrim(strupcase(chad1.galaxy),2), gto_indx, chad_indx

    distdata[gto_indx].litdist        = chad1[chad_indx].distance
    distdata[gto_indx].litdist_err    = chad1[chad_indx].distance_err
    distdata[gto_indx].litdist_texref = chad1[chad_indx].texref
    distdata[gto_indx].litdist_ref    = chad1[chad_indx].ref
    
; ---------------------------------------------------------------------------
; Tully (1988) NBG distances; assume a fixed distance uncertainty 
; ---------------------------------------------------------------------------

    flow = flow_distance(distdata.ra,distdata.dec,H0=H0,galaxy=distdata.galaxy)

    good = where(flow.tully_dist gt -900,ngood)
    if (ngood ne 0L) then begin
    
       srt = sort(distdata[good].ra)
       struct_print, flow[good[srt]]

       distdata[good].tullydist = flow[good].tully_dist
;      distdata[good].tullydist_ref = 'Tully 1988'

    endif
    
; ---------------------------------------------------------------------------
; compute model distances for all objects
; ---------------------------------------------------------------------------

    czgood = where(distdata.cz gt -900.0,nczgood,comp=czbad,ncomp=nczbad)
    mould = mould_distance(distdata[czgood].ra,distdata[czgood].dec,$
      distdata[czgood].cz,object=distdata[czgood].galaxy,/proper,$
      H0=redh100()*100.0,omega0=redomega0(),omega_lambda=redomegal())

    distdata[czgood].modeldist = mould.distance
;   distdata[czgood].modeldist_ref = 'Mould et al. 2000'

; ---------------------------------------------------------------------------
; final distances
; ---------------------------------------------------------------------------
    
    litdist = where(distdata.litdist gt -900.0,ndist,comp=modeldist)
    distdata[litdist].distance        = distdata[litdist].litdist
    distdata[litdist].distance_err    = distdata[litdist].litdist_err
    distdata[litdist].distance_texref = distdata[litdist].litdist_texref
    distdata[litdist].distance_ref    = distdata[litdist].litdist_ref
    
    struct_print, struct_trimtags(distdata,select=['galaxy','ra','dec','cz','distance',$
      'distance_err','distance_ref'])

    if keyword_set(write) then begin

       outfile = 'gto_distances.fits'
       splog, 'Writing '+analysis_path+outfile
       mwrfits, distdata, analysis_path+outfile, /create
       spawn, 'gzip -f '+analysis_path+outfile, /sh

       struct_print, struct_trimtags(distdata,select=['galaxy','ned_galaxy','ra','dec','cz','modeldist']), $
         file=analysis_path+'gto_distances.dat'
       
    endif

; ---------------------------------------------------------------------------
; generate QA plots
; ---------------------------------------------------------------------------

    if keyword_set(postscript) then begin
       postthick = 5.0 
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

    xrange = [(min(x)<min(y))>0.3,max(x)>max(y)]*[0.95,1.1]
    yrange = xrange
    
    xrange2 = xrange
    yrange2 = max(abs(residuals))*[-1.1,1.1]

    plotsym, 0, 1.8, /fill
    
    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.5,0.3], $
      ymargin=[0.3,1.3], position=pos, height=[6.0,3.4], /normal

    djs_plot, x, y, ps=8, xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, $
      position=pos[*,0], charsize=1.8, xtickname=replicate(' ',10), xsty=3, ysty=3, $
      charthick=postthick, xthick=postthick, ythick=postthick, /ylog, /xlog
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick
;   djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    
    djs_plot, x, residuals, ps=8, xtitle=xtitle, ytitle=ytitle2, xrange=xrange2, $
      yrange=yrange2, position=pos[*,1], /noerase, charsize=1.8, xsty=3, ysty=3, $
      charthick=postthick, xthick=postthick, ythick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
;   djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=1.5, $
      charthick=postthick, spacing=0
    
    im_openclose, postscript=postscript, /close    

;; --------------------------------------------------
;; compare the infall model and Tully distances
;; --------------------------------------------------
;
;    psname = 'dist_tully_vs_dist_model.ps'
;    im_openclose, analysis_path+psname, postscript=postscript
;    
;    good = where((flow.tully_dist gt -900.0) and (distdata.modeldist gt -900) and $
;      (distdata.litdist lt -900),ngood)
;
;    x = flow[good].tully_dist
;    y = distdata[good].modeldist
;
;    residuals = 100*(x-y)/x
;
;    stats = im_stats(residuals,/verbose)
;    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
;
;    xtitle = 'Tully/Kraan-Korteweg Distance [Mpc]'
;    ytitle = 'Infall Model Distance [Mpc]'
;    ytitle2 = 'Residuals [%]'
;
;    xrange = [min(x)<min(y),max(x)>max(y)]*[0.95,1.05]
;    yrange = xrange
;    
;    xrange2 = xrange
;    yrange2 = max(abs(residuals))*[-1.05,1.05]
;
;    plotsym, 0, 1.8, /fill
;    
;    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.5,0.3], $
;      ymargin=[0.3,1.3], position=pos, height=[6.0,3.4], /normal
;
;    djs_plot, x, y, ps=8, xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, $
;      position=pos[*,0], charsize=1.8, xtickname=replicate(' ',10), xsty=3, ysty=3, $
;      charthick=postthick, xthick=postthick, ythick=postthick, /xlog, /ylog
;    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick
;;   djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
;    
;    djs_plot, x, residuals, ps=8, xtitle=xtitle, ytitle=ytitle2, xrange=xrange2, $
;      yrange=yrange2, position=pos[*,1], /noerase, charsize=1.8, xsty=3, ysty=3, $
;      charthick=postthick, xthick=postthick, ythick=postthick, /xlog
;    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
;;   djs_oplot, !x.crange, [0,0], line=0, thick=postthick
;
;    legend, textoidl(xstr), /right, /top, box=0, charsize=1.5, $
;      charthick=postthick, spacing=0
;
;    im_openclose, postscript=postscript, /close    
    
return
end    
