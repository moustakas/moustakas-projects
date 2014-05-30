;+
; NAME:
;   deep2_zsuccess
; PURPOSE:
;   Build the completeness (redshift success rate) for each field as
;   a function of observed magnitude and color.
; CALLING SEQUENCE:
;   deep2_zsuccess, /clobber
; INPUTS: 
; KEYWORDS:
;   clobber - overwrite an existing output file
; OUTPUTS: 
;   FITS table with everything you need
; COMMENTS:
;   The completeness maps are tied to the 0d version. 
; MODIFICATION HISTORY:
;   J. Moustakas, 2014 May 26 - based on MF_COMPLETENESS
;-

function get_poisson_error, array, nbins=nbins
    poisson_err = im_poisson_limits(array)
    err = reform((poisson_err[1,*]-poisson_err[0,*])/2.0)
    if (n_elements(nbins) ne 0) then err = reform(err,nbins)
return, err
end

function get_zsuccess_and_error, got, attempt, zsuccess_error=zsuccess_error
    zsuccess = attempt*0-1.0
    zsuccess_error = attempt*0-1.0
    gd = where(attempt gt 0.0,ngd)
    if (ngd ne 0L) then begin
       zsuccess[gd] = got[gd]/attempt[gd]
; now compute the error
       fail = attempt[gd]-got[gd] ; hack for cosmos!
;      fail = (attempt[gd]-got[gd])>0 ; hack for cosmos!
       fail_err = get_poisson_error(fail)
       got_err = get_poisson_error(got[gd])
       zsuccess_error[gd] = sqrt((fail*got_err/attempt[gd]^2)^2 + (got[gd]*fail_err/attempt[gd]^2)^2)
    endif
return, zsuccess
end

function smooth_map, zsuccess, zsuccess_err, attempt=attempt, got=got, $
  params=params, verbose=verbose, minsnr=minsnr, smooth_got=got_smooth, $
  smooth_zsuccess_err=zsuccess_smooth_err, smooth_attempt=attempt_smooth, $
  color2=color2

    magaxis = mf_grid_axis(params,0)
    if keyword_set(color2) then begin
       coloraxis = deep2_grid_axis(params,2)
       nb = params.nbins[[0,2]]
    endif else begin
       coloraxis = deep2_grid_axis(params,1)
       nb = params.nbins[[0,1]]
    endelse

    maxcount = 8 ; should be a multiple of 4 (one for each dimension)
    minsnr = 3.0 ; 5.0

    got_smooth = got*0            ; initialize with the original maps
    attempt_smooth = attempt*0
    zsuccess_smooth = zsuccess*0-1.0
    zsuccess_smooth_err = zsuccess_err*0-1.0

; don't smooth pixels with well-measured weights
    msk = (zsuccess/zsuccess_err gt minsnr)*(zsuccess_err gt 0) ; 1=don't smooth
    wmsk = where(msk)
    got_smooth[wmsk] = got[wmsk]
    attempt_smooth[wmsk] = attempt[wmsk]
    zsuccess_smooth[wmsk] = zsuccess[wmsk]
    zsuccess_smooth_err[wmsk] = zsuccess_err[wmsk]

; sort by the most to the least populated pixels
    srt = reverse(sort(attempt))
    indx = array_indices(nb,srt,/dim)
    for ii = 0, cmproduct(nb)-1 do begin
       xx = indx[0,ii]
       yy = indx[1,ii]
       if (msk[xx,yy] eq 0) then begin
          count = 0
          xadd = 0 & xsub = 0 & yadd = 0 & ysub = 0
          doit = 1
          while doit do begin
;         while doit and (count lt maxcount) do begin
             thisattempt = total(attempt[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
               (yy-ysub)>0:(yy+yadd)<(nb[1]-1)]*(msk[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
               (yy-ysub)>0:(yy+yadd)<(nb[1]-1)] eq 0))
             thisgot = total(got[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
               (yy-ysub)>0:(yy+yadd)<(nb[1]-1)]*(msk[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
               (yy-ysub)>0:(yy+yadd)<(nb[1]-1)] eq 0))

             thissuccess = get_zsuccess_and_error(thisgot,thisattempt,zsuccess_error=thissuccess_err)

             if verbose then print, magaxis[xx], coloraxis[yy], count, xx, yy, $
               xadd, xsub, yadd, ysub, thisgot, thisattempt, thissuccess, thissuccess_err, $
               thissuccess/thissuccess_err

             if ((thissuccess ne -1.0) and (thissuccess/thissuccess_err gt minsnr)) or (count eq maxcount-1) then begin
; compute the success rate from the fatter pixel; the +1 and +=
; account for the -1 initialization of the completeness map
                zsuccess_smooth[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)] = zsuccess[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)]*msk[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)]+(msk[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)] eq 0)*thissuccess
                zsuccess_smooth_err[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)] = zsuccess_err[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)]*msk[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)]+(msk[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)] eq 0)*thissuccess_err

                npix = total(msk[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),(yy-ysub)>0:(yy+yadd)<(nb[1]-1)] eq 0)
;               if total(msk[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),(yy-ysub)>0:(yy+yadd)<(nb[1]-1)] gt 0) then stop
                
                attempt_smooth[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)] = attempt[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)]*msk[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)]+(msk[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)] eq 0)*thisattempt/npix
                
                got_smooth[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)] = got[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)]*msk[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)]+(msk[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),$
                  (yy-ysub)>0:(yy+yadd)<(nb[1]-1)] eq 0)*thisgot/npix
;               if count eq maxcount-1 then stop
                msk[(xx-xsub)>0:(xx+xadd)<(nb[0]-1),(yy-ysub)>0:(yy+yadd)<(nb[1]-1)] = 1
                doit = 0        ; well-measured success rate; exit
             endif else begin
                xsub += (count mod 4) eq 0    ; go brighter first
                xadd += (count mod 4) eq 1    ; then go fainter
                yadd += (count mod 4) eq 2    ; then go redder
                ysub += (count mod 4) eq 3    ; finally go bluer
                count++ 
             endelse 
          endwhile 
          if verbose then print
;         print, ii, total(attempt[where(msk)]>0), total(attempt_smooth[where(msk)]>0)
;         if (total(attempt_smooth[where(msk)]>0)-total(attempt[where(msk)]>0) gt 0.01) then stop
       endif 
    endfor 
    
return, zsuccess_smooth
end

function get_fracobj, attempt, nbins=nbins
    cumindex = reverse(sort(attempt))
    cumimage = fltarr(nbins)
    cumimage[cumindex] = total((attempt)[cumindex],/cumu)
    fracobj= cumimage/total(attempt)
return, fracobj
end

function init_outmap, params
; initialize the output completeness map
    arr1 = fltarr(params.nbins)-1.0
    arr2 = fltarr(params.nbins[[0,1]])-1.0
    arr3 = fltarr(params.nbins[[0,2]])-1.0
    arr4 = fltarr(params.nbins[0])-1.0

    out1 = {$
      field:        strtrim(params.field,2),$
;     nall:                    0L,$ ; total number of objects
      nobj:                    0L,$ ; number of primary objects (less AGN and stars)
;     catalog_weight:         0.0,$ ; catalog incompleteness
;     agn_frac:               0.0,$ ; fraction classified as AGN
; e.g., i - unweighted
      attempt_mag:         arr4,$ ; number attempted
      attempt_mag_err:     arr4,$
      got_mag:             arr4,$ ; number got
      got_mag_err:         arr4,$
      zsuccess_mag:        arr4,$ ; redshift success rate 
      zsuccess_mag_err:    arr4,$
; e.g., i - weighted
      attempt_weighted_mag:         arr4,$ ; number attempted
      attempt_weighted_mag_err:     arr4,$
;     got_weighted_mag:             arr4,$ ; number got
;     got_weighted_mag_err:         arr4,$
      zsuccess_weighted_mag:        arr4,$ ; redshift success rate 
      zsuccess_weighted_mag_err:    arr4,$
; e.g., i vs r-i *and* g-i - unweighted
      attempt_color1_color2:         arr1,$ ; number attempted
      attempt_color1_color2_err:     arr1,$
      attempt_color1_color2_smooth:         arr1,$ ; number attempted
      attempt_color1_color2_smooth_err:     arr1,$
      got_color1_color2:             arr1,$ ; number got
      got_color1_color2_err:         arr1,$
      got_color1_color2_smooth:             arr1,$ ; number got
      got_color1_color2_smooth_err:         arr1,$
      zsuccess_color1_color2:        arr1,$ ; redshift success rate 
      zsuccess_color1_color2_err:    arr1,$
      zsuccess_color1_color2_smooth:        arr1,$ ; redshift success rate 
      zsuccess_color1_color2_smooth_err:    arr1,$
      fracobj_color1_color2:         arr1,$ ; cumulative fraction of objects (for mfplot_completeness)
; e.g., i vs r-i *and* g-i - weighted
      attempt_weighted_color1_color2:         arr1,$ ; number attempted
      attempt_weighted_color1_color2_err:     arr1,$
;     got_weighted_color1_color2:             arr1,$ ; number got
;     got_weighted_color1_color2_err:         arr1,$
      zsuccess_weighted_color1_color2:        arr1,$ ; redshift success rate 
      zsuccess_weighted_color1_color2_err:    arr1,$
      fracobj_weighted_color1_color2:         arr1,$ ; cumulative fraction of objects (for mfplot_completeness)
; e.g., i vs r-i - unweighted
      attempt_color1:             arr2,$ ; number attempted
      attempt_color1_err:         arr2,$ 
      attempt_color1_smooth:      arr2,$
      attempt_color1_smooth_err:  arr2,$
      got_color1:                 arr2,$ ; number got
      got_color1_err:             arr2,$
      got_color1_smooth:          arr2,$
      got_color1_smooth_err:      arr2,$
      zsuccess_color1:            arr2,$ ; redshift success rate 
      zsuccess_color1_err:        arr2,$
      zsuccess_color1_smooth:     arr2,$
      zsuccess_color1_smooth_err: arr2,$
      fracobj_color1:             arr2,$ ; cumulative fraction of objects (for mfplot_completeness)
; e.g., i vs r-i - weighted
      attempt_weighted_color1:             arr2,$ ; number attempted
      attempt_weighted_color1_err:         arr2,$ 
      attempt_weighted_color1_smooth:      arr2,$
      attempt_weighted_color1_smooth_err:  arr2,$
;     got_weighted_color1:                 arr2,$ ; number got
;     got_weighted_color1_err:             arr2,$
;     got_weighted_color1_smooth:          arr2,$
;     got_weighted_color1_smooth_err:      arr2,$
      zsuccess_weighted_color1:            arr2,$ ; redshift success rate 
      zsuccess_weighted_color1_err:        arr2,$
      zsuccess_weighted_color1_smooth:     arr2,$
      zsuccess_weighted_color1_smooth_err: arr2,$
      fracobj_weighted_color1:             arr2,$ ; cumulative fraction of objects (for mfplot_completeness)
; e.g., i vs g-i - unweighted
      attempt_color2:             arr3,$ ; number attempted
      attempt_color2_err:         arr3,$
      attempt_color2_smooth:      arr2,$
      attempt_color2_smooth_err:  arr2,$
      got_color2:                 arr3,$ ; number got
      got_color2_err:             arr3,$
      got_color2_smooth:          arr2,$
      got_color2_smooth_err:      arr2,$
      zsuccess_color2:            arr3,$ ; redshift success rate 
      zsuccess_color2_err:        arr3,$
      zsuccess_color2_smooth:     arr3,$
      zsuccess_color2_smooth_err: arr3,$
      fracobj_color2:             arr3,$ ; cumulative fraction of objects (for mfplot_completeness)
; e.g., i vs g-i - weighted
      attempt_weighted_color2:             arr3,$ ; number attempted
      attempt_weighted_color2_err:         arr3,$
      attempt_weighted_color2_smooth:      arr2,$
      attempt_weighted_color2_smooth_err:  arr2,$
;     got_weighted_color2:                 arr3,$ ; number got
;     got_weighted_color2_err:             arr3,$
;     got_weighted_color2_smooth:          arr2,$
;     got_weighted_color2_smooth_err:      arr2,$
      zsuccess_weighted_color2:            arr3,$ ; redshift success rate 
      zsuccess_weighted_color2_err:        arr3,$
      zsuccess_weighted_color2_smooth:     arr3,$
      zsuccess_weighted_color2_smooth_err: arr3,$
      fracobj_weighted_color2:             arr3}  ; cumulative fraction of objects (for mfplot_completeness)
    if (n_elements(nmask) ne 0) then out1 = replicate(out1,nmask)
return, out1
end

pro deep2_zsuccess, clobber=clobber

    catpath = deep2_path(/catalogs)
    outfile = catpath+'deep2_zsuccess.fits'
    if im_file_test(outfile+'.gz',clobber=clobber) then return
    
    field = ['EGS','Fields24']
    nfield = n_elements(field)

    filters = deep2_filterlist()
    
; read the Q34 redshift catalog and the corresponding targeting
; weights for all the fields
    allzcat = read_deep2_zcat(photo=allphoto,weight=allweight,/all)
    egs = where(strmid(strtrim(allzcat.objno,2),0,1) eq '1',$
      negs,comp=fields24,ncomp=nfields24)
    
; for the QAplot
    levels = [0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.99,0.995]
    cannotation = string(levels,format='(F5.3)')
    mincol = 50
    psfile = catpath+'qa_deep2_zsuccess.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], width=6.8, height=6.0
    loadct, 16, /silent
    
; build the completeness maps for the EGS field and for fields 2-4
; separately 
    for ifield = 0, nfield-1 do begin
       splog, 'Working on field '+field[ifield]
       params = get_deep2_completeness_params(field[ifield])
       out = init_outmap(params)

       if ifield eq 0 then index = egs else index = fields24
       nobj = n_elements(index)
       zcat = allzcat[index]
       photo = allphoto[index]
       weight = allweight[index].targ_weight
       deep2_to_maggies, photo, maggies, /unwise
       
; get the selection magnitude and color(s)
       out.nobj = nobj
       colors = get_deep2_completeness_colors(maggies,targ_mag=photo.r,$
         params=params,filterlist=filters)

; distribution of objects targeted in bins of, e.g., r-i, g-i, and i
       nattempt_mag = nobj
       iattempt_mag = lindgen(nattempt_mag)
       iattempt_color1 = where(colors.good1)
       iattempt_color2 = where(colors.good2)
       iattempt_color1_color2 = where(colors.good1 and colors.good2)

; ...unweighted
       out.attempt_mag = hogg_histogram(colors[iattempt_mag].mag,$
         [params.hmin[0],params.hmax[0]],params.nbins[0])
       out.attempt_color1 = hogg_histogram(transpose([[colors[iattempt_color1].mag],$
         [colors[iattempt_color1].color1]]),transpose([[params.hmin[[0,1]]],[params.hmax[[0,1]]]]),$
         params.nbins[[0,1]])
       out.attempt_color2 = hogg_histogram(transpose([[colors[iattempt_color2].mag],$
         [colors[iattempt_color2].color2]]),transpose([[params.hmin[[0,2]]],[params.hmax[[0,2]]]]),$
         params.nbins[[0,2]])
       out.attempt_color1_color2 = hogg_histogram(transpose([[colors[iattempt_color1_color2].mag],$
         [colors[iattempt_color1_color2].color1],[colors[iattempt_color1_color2].color2]]),$
         transpose([[params.hmin],[params.hmax]]),params.nbins)

       out.attempt_mag_err = get_poisson_error(out.attempt_mag,nbins=params.nbins[0])
       out.attempt_color1_err = get_poisson_error(out.attempt_color1,nbins=params.nbins[[0,1]])
       out.attempt_color2_err = get_poisson_error(out.attempt_color2,nbins=params.nbins[[0,2]])
       out.attempt_color1_color2_err = get_poisson_error(out.attempt_color1_color2,nbins=params.nbins[[0,1,2]])
       
; ...weighted
       out.attempt_weighted_mag = hogg_histogram(colors[iattempt_mag].mag,$
         [params.hmin[0],params.hmax[0]],params.nbins[0],weight=weight[iattempt_mag])
       out.attempt_weighted_color1 = hogg_histogram(transpose([[colors[iattempt_color1].mag],$
         [colors[iattempt_color1].color1]]),transpose([[params.hmin[[0,1]]],[params.hmax[[0,1]]]]),$
         params.nbins[[0,1]],weight=weight[iattempt_color1])
       out.attempt_weighted_color2 = hogg_histogram(transpose([[colors[iattempt_color2].mag],$
         [colors[iattempt_color2].color2]]),transpose([[params.hmin[[0,2]]],[params.hmax[[0,2]]]]),$
         params.nbins[[0,2]],weight=weight[iattempt_color2])
       out.attempt_weighted_color1_color2 = hogg_histogram(transpose([[colors[iattempt_color1_color2].mag],$
         [colors[iattempt_color1_color2].color1],[colors[iattempt_color1_color2].color2]]),$
         transpose([[params.hmin],[params.hmax]]),params.nbins,weight=weight[iattempt_color1_color2])

       out.attempt_weighted_mag_err = get_poisson_error(out.attempt_weighted_mag,nbins=params.nbins[0])
       out.attempt_weighted_color1_err = get_poisson_error(out.attempt_weighted_color1,nbins=params.nbins[[0,1]])
       out.attempt_weighted_color2_err = get_poisson_error(out.attempt_weighted_color2,nbins=params.nbins[[0,2]])
       out.attempt_weighted_color1_color2_err = get_poisson_error(out.attempt_weighted_color1_color2,nbins=params.nbins[[0,1,2]])

; corresponding distribution of objects with good redshifts (Q>=3) 
       igot_mag = where(zcat.zquality ge 3,ngot_mag)
       igot_color1 = where(colors.good1 and zcat.zquality ge 3)
       igot_color2 = where(colors.good2 and zcat.zquality ge 3)
       igot_color1_color2 = where(colors.good1 and colors.good2 and $
         zcat.zquality ge 3)

; ...unweighted       
       out.got_mag = hogg_histogram(colors[igot_mag].mag,$
         [params.hmin[0],params.hmax[0]],params.nbins[0])
       out.got_color1 = hogg_histogram(transpose([[colors[igot_color1].mag],$
         [colors[igot_color1].color1]]),transpose([[params.hmin[[0,1]]],[params.hmax[[0,1]]]]),$
         params.nbins[[0,1]])
       out.got_color2 = hogg_histogram(transpose([[colors[igot_color2].mag],$
         [colors[igot_color2].color2]]),transpose([[params.hmin[[0,2]]],[params.hmax[[0,2]]]]),$
         params.nbins[[0,2]])
       out.got_color1_color2 = hogg_histogram(transpose([[colors[igot_color1_color2].mag],$
         [colors[igot_color1_color2].color1],[colors[igot_color1_color2].color2]]),$
         transpose([[params.hmin],[params.hmax]]),params.nbins)

       out.got_mag_err = get_poisson_error(out.got_mag,nbins=params.nbins[0])
       out.got_color1_err = get_poisson_error(out.got_color1,nbins=params.nbins[[0,1]])
       out.got_color2_err = get_poisson_error(out.got_color2,nbins=params.nbins[[0,2]])
       out.got_color1_color2_err = get_poisson_error(out.got_color1_color2,nbins=params.nbins[[0,1,2]])
       
; redshift success rate and error

; ...unweighted       
       out.zsuccess_mag = get_zsuccess_and_error(out.got_mag,out.attempt_mag,zsuccess_error=err)
       out.zsuccess_mag_err = err

       out.zsuccess_color1 = get_zsuccess_and_error(out.got_color1,out.attempt_color1,zsuccess_error=err)
       out.zsuccess_color1_err = err

       out.zsuccess_color2 = get_zsuccess_and_error(out.got_color2,out.attempt_color2,zsuccess_error=err)
       out.zsuccess_color2_err = err
       
       out.zsuccess_color1_color2 = get_zsuccess_and_error(out.got_color1_color2,$
         out.attempt_color1_color2,zsuccess_error=err)
       out.zsuccess_color1_color2_err = err

; ...weighted       

       out.zsuccess_weighted_mag = get_zsuccess_and_error(out.got_mag,out.attempt_weighted_mag,zsuccess_error=err)
       out.zsuccess_weighted_mag_err = err

;      ww = where(out.attempt_weighted_color1-out.got_color1 lt 0)
;      niceprint, out.got_color1[ww], out.attempt_weighted_color1[ww], out.attempt_color1[ww]
       out.zsuccess_weighted_color1 = get_zsuccess_and_error(out.got_color1,out.attempt_weighted_color1,zsuccess_error=err)
       out.zsuccess_weighted_color1_err = err

       out.zsuccess_weighted_color2 = get_zsuccess_and_error(out.got_color2,out.attempt_weighted_color2,zsuccess_error=err)
       out.zsuccess_weighted_color2_err = err
       
       out.zsuccess_weighted_color1_color2 = get_zsuccess_and_error(out.got_color1_color2,$
         out.attempt_weighted_color1_color2,zsuccess_error=err)
       out.zsuccess_weighted_color1_color2_err = err

; get the cumulative fraction of objects in each bin (for the contour
; plot in mfplot_completeness)
       out.fracobj_color1 = get_fracobj(out.attempt_color1,nbins=params.nbins[[0,1]])
       out.fracobj_color2 = get_fracobj(out.attempt_color2,nbins=params.nbins[[0,2]])
       out.fracobj_color1_color2 = get_fracobj(out.attempt_color1_color2,nbins=params.nbins)

       out.fracobj_weighted_color1 = get_fracobj(out.attempt_weighted_color1,nbins=params.nbins[[0,1]])
       out.fracobj_weighted_color2 = get_fracobj(out.attempt_weighted_color2,nbins=params.nbins[[0,2]])
       out.fracobj_weighted_color1_color2 = get_fracobj(out.attempt_weighted_color1_color2,nbins=params.nbins)

; adaptively smooth the completeness map; the idea here is to include
; galaxies from adjacent (high-resolution) pixels until we get a
; minimum S/N ratio on the success rate

; ...unweighted       
       out.zsuccess_color1_smooth = smooth_map(out.zsuccess_color1,out.zsuccess_color1_err,$
         attempt=out.attempt_color1,got=out.got_color1,smooth_zsuccess_err=zsuccess_smooth_err,$
         params=params,verbose=0,minsnr=minsnr,smooth_attempt=attempt_smooth,smooth_got=got_smooth)
       out.zsuccess_color1_smooth_err = zsuccess_smooth_err
       out.got_color1_smooth = got_smooth
       out.attempt_color1_smooth = attempt_smooth

       out.zsuccess_color2_smooth = smooth_map(out.zsuccess_color2,out.zsuccess_color2_err,$
         attempt=out.attempt_color2,got=out.got_color2,smooth_zsuccess_err=zsuccess_smooth_err,$
         params=params,verbose=0,minsnr=minsnr,smooth_attempt=attempt_smooth,smooth_got=got_smooth,$
         /color2)
       out.zsuccess_color2_smooth_err = zsuccess_smooth_err
       out.got_color2_smooth = got_smooth
       out.attempt_color2_smooth = attempt_smooth

;; smooth_map() does not yet support a 3D grid       
;       out.zsuccess_color1_color2_smooth = smooth_map(out.zsuccess_color1_color2,out.zsuccess_color1_color2_err,$
;         attempt=out.attempt_color1_color2,got=out.got_color1_color2,smooth_zsuccess_err=zsuccess_smooth_err,$
;         params=params,verbose=0,minsnr=minsnr,smooth_attempt=attempt_smooth,smooth_got=got_smooth)
;       out.zsuccess_color1_color2_smooth_err = zsuccess_smooth_err
;       out.got_color1_color2_smooth = got_smooth
;       out.attempt_color1_color2_smooth = attempt_smooth

; ...weighted       
       out.zsuccess_weighted_color1_smooth = smooth_map(out.zsuccess_weighted_color1,out.zsuccess_weighted_color1_err,$
         attempt=out.attempt_weighted_color1,got=out.got_color1,smooth_zsuccess_err=zsuccess_weighted_smooth_err,$
         params=params,verbose=0,minsnr=minsnr,smooth_attempt=attempt_weighted_smooth,smooth_got=got_smooth)
       out.zsuccess_weighted_color1_smooth_err = zsuccess_weighted_smooth_err
       out.attempt_weighted_color1_smooth = attempt_weighted_smooth
       
       out.zsuccess_weighted_color2_smooth = smooth_map(out.zsuccess_weighted_color2,out.zsuccess_weighted_color2_err,$
         attempt=out.attempt_weighted_color2,got=out.got_color2,smooth_zsuccess_err=zsuccess_weighted_smooth_err,$
         params=params,verbose=0,minsnr=minsnr,smooth_attempt=attempt_weighted_smooth,smooth_got=got_smooth,$
         /color2)
       out.zsuccess_weighted_color2_smooth_err = zsuccess_weighted_smooth_err
       out.attempt_weighted_color2_smooth = attempt_weighted_smooth

;; test the code -        
;;      wetry = where((colors.mag le 23.2) and (colors.mag ge bright) and $
;;        (colors.color1 ge params.hmin[1]) and (colors.color1 le params.hmax[1]),ntry)
;;      wegot = where(zerod.isgood and (zerod.zprimus_zconf ge 3) and $
;;        (colors.mag le 23.2) and (colors.mag ge bright) and $
;;        (colors.color1 ge params.hmin[1]) and (colors.color1 le params.hmax[1]),ngot)
;
;;       wetry = where((colors.mag le minmag) and (colors.mag ge bright) and colors.good1,ntry)
;;       wegot = where(zerod.isgood and (zerod.zprimus_zconf ge 3) and $
;;         (colors.mag le minmag) and (colors.mag ge bright) and colors.good1,ngot)
;       
;;       wetry = where((colors.mag lt minmag) and (colors.mag gt bright) and colors.good1 and $
;;         (colors.color1 ge params.hmin[1]) and (colors.color1 le params.hmax[1]),ntry)
;;       wegot = where(zerod.isgood and (zerod.zprimus_zconf ge 3) and $
;;         (colors.mag lt minmag) and (colors.mag gt bright) and colors.good1 and $
;;         (colors.color1 gt params.hmin[1]) and (colors.color1 lt params.hmax[1]),ngot)
;
;       minmag = 22.0
;       wetry = where((colors.mag le minmag) and (colors.mag ge bright),ntry)
;       wegot = where(zerod.isgood and (zerod.zprimus_zconf ge 3) and $
;         (colors.mag le minmag) and (colors.mag ge bright),ngot)
;       
;       zw1 = get_mf_zsuccess(field[ifield],mag=colors[wegot].mag,$
;         color1=colors[wegot].color1,color2=colors[wegot].color2,$
;         /noextrap,flag=flg)
;       nall = total(zerod[wetry].targ_weight)
;       print, ngot, ntry, nall, total(1/zw1), total(1/zw1)/nall
;
;       case1 = where(flg eq 1) & case2 = where(flg eq 2) & case3 = where(flg eq 4)
;       print, total(1/zw1[case1]), total(1/zw1[case2]), total(1/zw1[case3]), total(1/zw1), total(1/zw1)/nall
;       
;       miss = where((colors[wegot].color1 lt params.hmin[1]) or $
;         (colors[wegot].color1 gt params.hmax[1]))
;       nmiss = long(total(zerod[wegot[miss]].targ_weight))
;
;;; using just the magnitude lookup table works perfectly!       
;;       wetry = where((colors.mag le minmag) and (colors.mag ge bright),ntry)
;;       wegot = where(zerod.isgood and (zerod.zprimus_zconf ge 3) and $
;;         (colors.mag le minmag) and (colors.mag ge bright),ngot)
;;       test = get_mf_zsuccess(field[ifield],mag=colors[wegot].mag)
;;       print, ngot, ntry, nall, total(1/test), total(1/test)/nall

; make the QAplot
       magrange = [params.hmin[0],params.hmax[0]]
       colorrange = [params.hmin[1],params.hmax[1]]
       magaxis = mf_grid_axis(params,0,/midbin)
       coloraxis = mf_grid_axis(params,1,/midbin)
       mag2daxis = magaxis#(coloraxis*0.0+1.0)
       color2daxis = transpose(coloraxis#(magaxis*0.0+1.0))

; observed - all
       good = where(out.zsuccess_color1 ge 0.0)
       mf_hogg_scatterplot, mag2daxis[good], color2daxis[good], weight=out.zsuccess_color1[good], position=pos, $
         xstyle=1, ystyle=1, xrange=magrange, yrange=colorrange, xnpix=params.nbins[0], ynpix=params.nbins[1], $
         xtitle=textoidl(params.nice_targ_filter+' (AB mag)'), ytitle=textoidl(params.nice_color_filter[0]), $
         color=cgcolor('black',1), exp=2.0, darkest=mincol, /nocontour
       contour, out.fracobj_color1, mag2daxis, color2daxis, levels=levels, $
         /overplot, c_annotation=cannotation, color=cgcolor('black',1)
       im_legend, field[ifield]+' - observed, all', /left, /top, box=0, $
         charsize=1.3, margin=0, textcolor=cgcolor('black',0)
; smoothed - all
       good = where(out.zsuccess_color1_smooth ge 0.0)
       mf_hogg_scatterplot, mag2daxis[good], color2daxis[good], weight=out.zsuccess_color1_smooth[good], position=pos, $
         xstyle=1, ystyle=1, xrange=magrange, yrange=colorrange, xnpix=params.nbins[0], ynpix=params.nbins[1], $
         xtitle=textoidl(params.nice_targ_filter+' (AB mag)'), ytitle=textoidl(params.nice_color_filter[0]), $
         color=cgcolor('black',1), exp=2.0, darkest=mincol, /nocontour
       contour, out.fracobj_color1, mag2daxis, color2daxis, levels=levels, $
         /overplot, c_annotation=cannotation, color=cgcolor('black',1)
       im_legend, field[ifield]+' - smoothed, all', /left, /top, box=0, $
         charsize=1.3, margin=0, textcolor=cgcolor('black',0)
; observed - high S/N
       good = where(out.zsuccess_color1/out.zsuccess_color1_err gt 1.0) ; minsnr
       mf_hogg_scatterplot, mag2daxis[good], color2daxis[good], weight=out.zsuccess_color1[good], position=pos, $
         xstyle=1, ystyle=1, xrange=magrange, yrange=colorrange, xnpix=params.nbins[0], ynpix=params.nbins[1], $
         xtitle=textoidl(params.nice_targ_filter+' (AB mag)'), ytitle=textoidl(params.nice_color_filter[0]), $
         color=cgcolor('black',1), exp=2.0, darkest=mincol, /nocontour
       contour, out.fracobj_color1, mag2daxis, color2daxis, levels=levels, $
         /overplot, c_annotation=cannotation, color=cgcolor('black',1)
       im_legend, field[ifield]+' - observed, S/N>1', /left, /top, box=0, $
         charsize=1.3, margin=0, textcolor=cgcolor('black',0)
; smoothed - high S/N
       good = where(out.zsuccess_color1_smooth/out.zsuccess_color1_smooth_err gt 1.0) ; minsnr
       mf_hogg_scatterplot, mag2daxis[good], color2daxis[good], weight=out.zsuccess_color1_smooth[good], position=pos, $
         xstyle=1, ystyle=1, xrange=magrange, yrange=colorrange, xnpix=params.nbins[0], ynpix=params.nbins[1], $
         xtitle=textoidl(params.nice_targ_filter+' (AB mag)'), ytitle=textoidl(params.nice_color_filter[0]), $
         color=cgcolor('black',1), exp=2.0, darkest=mincol, /nocontour
       contour, out.fracobj_color1, mag2daxis, color2daxis, levels=levels, $
         /overplot, c_annotation=cannotation, color=cgcolor('black',1)
       im_legend, field[ifield]+' - smoothed, S/N>1', /left, /top, box=0, $
         charsize=1.3, margin=0, textcolor=cgcolor('black',0)

;;; fit a 2D model - the bspline code below is fine, but the
;;; completeness maps don't lend themselves well to a simple model 
;;       magaxis = mf_grid_axis(params,0)
;;       color1axis = mf_grid_axis(params,1)
;;       color2axis = mf_grid_axis(params,2)
;;
;;       magcolor1_2daxis = magaxis#(color1axis*0.0+1.0)
;;       color1_2daxis = transpose(color1axis#(magaxis*0.0+1.0))
;;       ivar = 1.0/out.zsuccess_color1_err^2*(out.zsuccess_color1_err gt 0.0)
;;
;;;      gd = where(out.zsuccess_color1*sqrt(ivar) gt 3.0)
;;       gd = where(out.zsuccess_color1 ge 0.0,ngd)
;;       yy = out.zsuccess_color1[gd]
;;       x2 = magcolor1_2daxis[gd]
;;       xx = color1_2daxis[gd]
;;       iv = ivar[gd]
;;       sset = bspline_iterfit(xx,yy,invvar=iv,x2=x2,bkspace=1.0,$
;;         nord=4,npoly=4,xmin=params.hmin[0]+1.5)
;;       yfit = bspline_valu(color1_2daxis,sset,x2=magcolor1_2daxis)
;;       chi2 = total(ivar*(out.zsuccess_color1-yfit)^2)
;;       splog, 'chi^2 ', chi2
;;
;;       for bb = 0, params.nbins[1]-1 do begin
;;          gd = where(out.zsuccess_color1[*,bb] gt -1.0,comp=bad)
;;          ploterror, magaxis[gd], out.zsuccess_color1[gd,bb], $
;;            out.zsuccess_color1_err[gd,bb], psym=6, symsize=2.5, $
;;            xrange=minmax(magaxis), yrange=[-0.1,1.2], xsty=3, ysty=3
;;          djs_oplot, magaxis, yfit[*,bb], color=cgcolor('green'), thick=4, psym=-8
;;          cc = get_kbrd(1)
;;       endfor
;;
;;       for bb = 0, params.nbins[0]-1 do begin
;;          gd = where(out.zsuccess_color1[bb,*] gt -1.0,comp=bad)
;;;         coeff = poly_fit(color1axis[gd],out.zsuccess_color1[bb,gd],2,$
;;;           measure_errors=out.zsuccess_color1_err[bb,gd])
;;;         yfit = poly(color1axis,coeff)
;;          ploterror, color1axis[gd], out.zsuccess_color1[bb,gd], $
;;            out.zsuccess_color1_err[bb,gd], psym=6, symsize=2.5, $
;;            xrange=minmax(color1axis), yrange=[-0.1,1.2], xsty=3, ysty=3
;;          if (bad[0] ne -1) then djs_oplot, magaxis[bad], out.zsuccess_color1[bb,bad], $
;;            psym=6, symsize=2.5, color=cgcolor('blue')
;;          djs_oplot, color1axis, yfit[bb,*], color=cgcolor('red'), thick=4, psym=-8
;;;         djs_oplot, color1axis, yfit, color=cgcolor('red'), thick=4
;;          cc = get_kbrd(1)
;;       endfor

; make a big structure!       
       if (n_elements(allout) eq 0L) then allout = out else $
         allout = [temporary(allout),out]
    endfor ; close field loop

    im_plotconfig, /psclose, psfile=psfile, /pdf
    loadct, 0, /silent
    
; write out
    im_mwrfits, allout, outfile, clobber=clobber

return
end
