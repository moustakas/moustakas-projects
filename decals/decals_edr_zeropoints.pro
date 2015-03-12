pro decals_isedfit_zeropoints, debug=debug, calculate=calculate, $
  makeplots=makeplots

    prefix = 'decals'
    isedfit_dir = getenv('DECALS_DIR')+'/isedfit/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    edr_dir = getenv('HOME')+'/edr/'
    zptdir = isedfit_dir+'zeropoints/'

    dweff = k_lambda_eff(filterlist=decals_filterlist())
    sweff = k_lambda_eff(filterlist=sdss_filterlist())
    
; fitted sample    
    tractor = mrdfits(edr_dir+'edr-specz-dr10.fits',1)
    specz = mrdfits(edr_dir+'edr-specz-dr10.fits',2)
    phot = mrdfits(edr_dir+'edr-specz-dr10.fits',3)
    these = where(specz.z ge 0.05 and specz.z le 0.7 and specz.zwarning eq 0 and $
      strtrim(specz.class,2) eq 'GALAXY' and tractor.brick_primary eq 'T',ngal)
    tractor = tractor[these]
    specz = specz[these]

; isedfit results    
    ised = read_isedfit(isedfit_paramfile,isedfit_dir=isedfit_dir,$
      outprefix='nodecals')

; bin the catalog in sky points
;   hgrid = hogg_histogram(transpose([[tractor.ra],[tractor.dec]]),$
;     transpose([[min(tractor.ra),min(tractor.dec)],$
;     [max(tractor.ra),max(tractor.dec)]]),[25,25])

    brickwid = 0.25
    nbins = [20,20]
    ramin =min(tractor.ra)
    ramax = max(tractor.ra)
    decmin = min(tractor.dec)
    decmax = max(tractor.dec)

    rawid = ramax-ramin
    decwid = decmax-decmin
    nbins = [fix(rawid/(2*brickwid)),fix(decwid/(2*brickwid))]
    splog, 'Number of pixels ', nbins
    splog, 'Width (deg) ', rawid, decwid
    splog, 'Block size (deg) ', [rawid,decwid]/nbins

    if keyword_set(calculate) then begin
       grid = hist_nd(transpose([[tractor.ra],[tractor.dec]]),nbins=nbins,$
         min=[ramin,decmin],max=[ramax,decmax],reverse=rev)
       ntotbins = cmproduct(nbins)

; build several different maps of the zeropoint differences, one for
; each band (grz) 
       maps = {$
         nobj:   intarr(nbins[0],nbins[1]),$
         median: fltarr(nbins[0],nbins[1]),$
         sigma:  fltarr(nbins[0],nbins[1])}
       maps = replicate(maps,3)

;   if keyword_set(debug) then djs_plot, tractor.ra, tractor.dec, psym=3, xsty=3, ysty=3
       for kk = 0, ntotbins-1 do begin
          if (rev[kk+1] gt rev[kk]) and (grid[kk] gt 5) then begin
             indx1 = rev[rev[kk]:rev[kk+1]-1]
             imindx = array_indices(nbins,kk,/dim)
             
             if keyword_set(debug) then begin
                djs_plot, tractor.ra, tractor.dec, psym=3, xsty=3, ysty=3
                djs_oplot, [tractor[indx1].ra], [tractor[indx1].dec], psym=6, color='orange'
             endif
;         print, n_elements(indx1), grid[imindx[0],imindx[1]]
             
             for iband = 0, 2 do begin
                good = where(ised[indx1].maggies[iband] gt 0,ngood)
                zpt = -2.5*alog10(ised[indx1[good]].maggies[iband]/$
                  ised[indx1[good]].bestmaggies[iband])
                maps[iband].nobj[imindx[0],imindx[1]] = ngood
                maps[iband].median[imindx[0],imindx[1]] = djs_median(zpt)
                maps[iband].sigma[imindx[0],imindx[1]] = djsig(zpt)
             endfor
          endif 
       endfor
       im_mwrfits, maps, zptdir+'zpt_maps.fits', /clobber
    endif

; make some cool plots       
    if keyword_set(makeplots) then begin

; zeropoint offset plots
       nband = 8
       zpt = fltarr(nband)
       zpterr = fltarr(nband)
       nobj = lonarr(nband)
       for iband = 0, nband-1 do begin
          good = where(ised.maggies[iband] gt 0,ngood)
          nobj[iband] = ngood
          zpt[iband] = djs_median(-2.5*alog10(ised[good].maggies[iband]/ised[good].bestmaggies[iband]))
          zpterr[iband] = djsig(-2.5*alog10(ised[good].maggies[iband]/ised[good].bestmaggies[iband]))
       endfor

       dnobj = nobj[0:2]
       dzpt = zpt[0:2]
       dzpterr = zpterr[0:2]
       snobj = nobj[3:7]
       szpt = zpt[3:7]
       szpterr = zpterr[3:7]
       
       psfile = zptdir+'qa_specz_zeropoints.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.5, xmargin=[1.3,0.4], width=6.8
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='Filter Wavelength (\AA)', $
         ytitle='\Delta'+'m (AB mag)', $
         xrange=[3200,9500], yrange=0.6*[-1,1]
       djs_oplot, !x.crange, [0,0], line=0, color=cgcolor('grey')
       oploterror, sweff, szpt, szpterr, psym=symcat(16), $
         color=cgcolor('dodger blue'), errcolor=cgcolor('dodger blue'), $
         symsize=2.0
       oploterror, dweff, dzpt, dzpterr, psym=symcat(15), $
         color=cgcolor('firebrick'), errcolor=cgcolor('firebrick'), $
         symsize=2.0
       im_legend, ['SDSS ugriz','DECaLS grz'], psym=[16,15], $
         color=['dodger blue','firebrick'], /right, /top, box=0
       im_plotconfig, /psclose, /pdf, psfile=psfile

; make some residual plots
       iband = 1
       good = where(ised.maggies[iband] gt 0,ngood)
       zpt = -2.5*alog10(ised[good].maggies[iband]/ised[good].bestmaggies[iband])
       djs_plot, ised[good].z, zpt, psym=3
       djs_plot, tractor[good].extinction[1], zpt, psym=3
       djs_plot, -2.5*alog10(ised[good].maggies[5]), zpt, psym=3, xr=[14,23], yr=[-1,1]

stop       
       
; 2D maps       
       ncol = 10
       cgloadct, 3, ncolors=ncol, bottom=1
       band = ['g','r','z']+'-band'
    
; median offset
       range = [0,ceil(max(maps.median)/10.0)*10]*100 ; millimag
       for iband = 0, 2 do begin
          cgimage, maps[iband].median*100, /axes, multimargin=[0,0,0,4], ncolors=ncol, $
            bottom=1, xrange=[ramax,ramin], yrange=[decmin,decmax], oposition=pos, $
            xtitle='R. A. (degrees)', ytitle='Dec (degrees)', /keep_aspect
          cpos = [pos[2]+0.03,pos[1],pos[2]+0.06,pos[3]]
          cgcolorbar, range=range, ncolors=ncol+1, title='Median Offset ('+band[iband]+', mmag)', $
            /vertical, position=cpos, bottom=1, format='(I0)', /right, /discrete
          cc = get_kbrd(1)
       endfor
    
; number of objects    
       range = [0,ceil(max(maps.nobj)/10.0)*10]
       cgimage, maps[1].nobj, /axes, multimargin=[0,0,0,4], ncolors=ncol, $
         bottom=1, xrange=[ramax,ramin], yrange=[decmin,decmax], oposition=pos, $
         xtitle='R. A. (degrees)', ytitle='Dec (degrees)', /keep_aspect
       cpos = [pos[2]+0.03,pos[1],pos[2]+0.06,pos[3]]
       cgcolorbar, range=range, ncolors=ncol+1, title='Number of Galaxies (r-band)', $
         /vertical, position=cpos, bottom=1, format='(I0)', /right, /discrete
    endif
    
stop    

return
end
    
