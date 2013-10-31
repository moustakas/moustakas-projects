function get_bracket_files, Zmetal, sspinfo=sspinfo
; internal support routine: get the two SSPs that bracket the desired
; metallicity, to allow for interpolation
    nZ = n_elements(Zmetal)
    bracket_files = strarr(2,nZ)

    for ii = 0L, nZ-1 do begin
       indx = findex(sspinfo.Zmetal,Zmetal[ii])
       below = fix(indx)
       above = ceil(indx)
       bracket_files[0,ii] = sspinfo.sspfile[below]
       bracket_files[1,ii] = sspinfo.sspfile[above]
    endfor
    
return, bracket_files
end

function read_and_interpolate, files, Zmetal=Zmetal, age=age, $
  ssppath=ssppath, fits_grid=fits_grid
; internal support routine: read and interpolate the base model file
; onto the desired metallicity grid 

    nZ = n_elements(Zmetal)
    nage = n_elements(age)
    
    files = strtrim(files,2)
    junk1 = file_search(ssppath+files[0],count=c1)
    junk2 = file_search(ssppath+files[1],count=c2)
    if (c1 eq 0) or (c2 eq 0) then message, 'Problem finding files'
    
    if (strtrim(files[0],2) eq strtrim(files[1],2)) then begin
       fits = gz_mrdfits(ssppath+files[0],1,/silent)
       fits = replicate(temporary(fits),nZ) ; identical for all objects
    endif else begin
       fits_grid = [gz_mrdfits(ssppath+files[0],1,/silent),$
         gz_mrdfits(ssppath+files[1],1,/silent)]
       fits = replicate(fits_grid[0],nZ)
       indx = findex(fits_grid.Zmetal,Zmetal)
; this interpolate can be very slow because we're building an
; NPIX,NAGE,nZ elements array
;      t0 = systime(1)
       fits.flux = interpolate(fits_grid.flux,indx)
;      splog, 'Total time (sec) = ', (systime(1)-t0)
       fits.mstar = interpolate(fits_grid.mstar,indx)
       fits.Zmetal = interpolate(fits_grid.Zmetal,indx)
    endelse

; now interpolate in age
    ageindx = findex(fits.age,age*1D9)

    npix = n_elements(fits.wave)
    out = {zmetal: fits.Zmetal, age: age, wave: fits.wave, $
      flux: interpolate(fits.flux,ageindx)}
    
return, out
end

pro plotbcgsfhs_ssps, pdf=pdf, build_models=build_models, clobber=clobber
; jm13oct25siena - construct a grid of SSPs of varying metallicity and
; age; plot the model grids and the data on a plot for the paper

    sample = read_bcgsfhs_sample()
    ncl = n_elements(sample)
    
    bcgpath = bcgsfhs_path()
    if keyword_set(pdf) then begin
       paperpath = bcgsfhs_path()
       suffix = '.ps'
    endif else begin
       paperpath = bcgsfhs_path(/paper)
       suffix = '.eps'
    endelse

    spsmodels = 'fsps_v2.4_miles'
    imf = 'salp'

    outfile = bcgpath+'ssps_'+spsmodels+'_'+imf+'.fits'

; build the SSPs; use the iSEDfit SSP data model and info structure     
    if keyword_set(build_models) then begin
       ssppath = getenv('ISEDFIT_SSP_DIR')+'/'
       info = mrdfits(ssppath+'info_'+spsmodels+'_'+imf+'.fits.gz',1)

; desired output bands, age grid, and metallicities
       filt = ['sdss_u0','sdss_r0','twomass_J','twomass_H']+'.par'
       nfilt = n_elements(filt)

       nage = 8
       nZmetal = 10
       Zmetal = range(0.004,0.03,nZmetal)
       age = range(5.0,12.0,nage)

       out = replicate({Zmetal: 0.0, age: age, ur: age*0.0, $
         rj: age*0.0, jh: age*0.0},nZmetal)
       out.Zmetal = Zmetal

       files = get_bracket_files(Zmetal,sspinfo=info)
       for iZ = 0, nZmetal-1 do begin
          ssp1 = read_and_interpolate(files[*,iZ],Zmetal=Zmetal[iZ],$
            age=age,ssppath=ssppath+spsmodels+'/')
          wave_edges = k_lambda_to_edges(ssp1.wave)

          maggies = fltarr(nfilt,nage)
          for ia = 0, nage-1 do begin
             maggies[*,ia] = reform(k_project_filters(wave_edges,$
               ssp1.flux[*,ia],filterlist=filt))
          endfor
          out[iZ].ur = reform(-2.5*alog10(maggies[0,*]/maggies[1,*]))
          out[iZ].rj = reform(-2.5*alog10(maggies[1,*]/maggies[2,*]))
          out[iZ].jh = reform(-2.5*alog10(maggies[2,*]/maggies[3,*]))
       endfor

       djs_plot, out.jh, out.ur, psym=8, xsty=3, ysty=3
       djs_plot, out.rj, out.ur, psym=8, xsty=3, ysty=3

       im_mwrfits, out, outfile, clobber=clobber
    endif

; read the models
    ssp = gz_mrdfits(outfile,1)
    nZmetal = n_elements(ssp)
    nage = n_elements(ssp[0].age)    
    
; make the plot    
    ageline = 0
    Zmetalline = 0
    
    agecolor = 'firebrick'
    Zmetalcolor = 'forest green'

    xrange = [0.6,1.4]
    yrange = [2.0,3.0]

    psfile = paperpath+'bcg_ssps'+suffix
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.3

    pos = im_getposition(nx=5,ny=3,/landscape,yspace=0.1,xspace=0.1,$
      xmargin=[1.0,0.4])

    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)

       if ic gt 9 then begin
          delvarx, xtickname
       endif else begin
          xtickname = replicate(' ',10)
       endelse

       if (ic mod 5) eq 0 then begin
          delvarx, ytickname
       endif else begin
          ytickname = replicate(' ',10)
       endelse
       
       djs_plot, [0], [0], /nodata, position=pos[*,ic], noerase=ic gt 0, $
         xsty=1, ysty=1, yrange=yrange, xrange=xrange, $
         ytickname=ytickname, xtickname=xtickname, ytickinterval=0.3, $
         xtickinterval=0.3
       im_legend, strupcase(cluster), /left, /top, box=0, margin=0, charsize=1.0

; overplot the model grids       
       for iZ = 0, nZmetal-1 do djs_oplot, ssp[iZ].rj, ssp[iZ].ur, $
         line=Zmetalline, color=cgcolor(Zmetalcolor), thick=2
       for ia = 0, nage-1 do djs_oplot, ssp.rj[ia], ssp.ur[ia], $
         line=ageline, color=cgcolor(agecolor), thick=2
    endfor

    xyouts, min(pos[0,*])-0.05, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
      textoidl('u - r'), orientation=90, align=0.5, charsize=1.4, /norm
    xyouts, (max(pos[2,*])-min(pos[0,*]))/2.0+min(pos[0,*]), min(pos[1,*])-0.09, $
      textoidl('r - J'), align=0.5, charsize=1.4, /norm
    
    im_plotconfig, psfile=psfile, /psclose, pdf=pdf

return
end
    
    
