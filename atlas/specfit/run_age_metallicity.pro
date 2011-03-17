pro run_age_metallicity, object, postscript=postscript
; jm03dec3uofa
; investigate the age-metallicity degeneracy problem for our template
; fitting results

    atlaspath = atlas_path(/atlas1d)
    specfitpath = atlas_path(/specfit)
    atlas = read_integrated(/silent)
    
; to examine just one object crop the list
    
    if n_elements(object) ne 0L then begin

       doit = match_string(object,atlas.galaxy,index=match) ; match by galaxy name
       if match[0] eq -1L then begin
          splog, 'Object '+object+' not found!'
          return
       endif
       atlas = atlas[match] 
       
    endif

    galaxy = strcompress(atlas.galaxy,/remove)
    title = strupcase(galaxy)
;   ngalaxy = 3L & atlas = atlas[0L:ngalaxy-1L]
    ngalaxy = n_elements(galaxy)

    line = read_linepars(linepath=specfitpath,linefile='elinelist.dat')
    linewave = line.wave
    nline = n_elements(linewave)

    ageZ = {$
      Z_TEMPLATE:   0.0,         $ ; best-fitting continuum metallicity
      BEST_INDX:    -1L,         $ ; index of the best metallicity template
      CHI2:         fltarr(3),   $ ; chi2 of the continuum fit [0.05,0.02,0.004]
      OII_3727:     fltarr(2,3), $
      OIII_5007:    fltarr(2,3), $
      H_ALPHA:      fltarr(2,3), $
      OII_3727_EW:  fltarr(2,3), $
      OIII_5007_EW: fltarr(2,3), $
      H_ALPHA_EW:   fltarr(2,3)}
    ageZ = replicate(ageZ,ngalaxy)
    ageZ = struct_addtags(struct_trimtags(atlas,select=['GALAXY','KD*']),ageZ)

    Z = [0.05,0.02,0.004]
    strZ = 'Z = '+['0.05','0.02','0.004']

    psname = 'run_age_metallicity.ps'
    if (not keyword_set(postscript)) then window, 0, xs=550, ys=450
    im_openclose, psname, postscript=postscript

    pagemaker, nx=1, ny=2, position=pos, /normal, yspace=0.0, $
      xspace=0.0, xmargin=[1.1,0.15], ymargin=[0.5,1.0]
    scale = 1E15

    stime0 = systime(1)
    for k = 0L, ngalaxy-1L do begin
    
       fsuper = atlas1d_specfit(atlas[k].driftfile,datapath=atlaspath,specfit=fitsuper,eigenspec=2)
       fsun = atlas1d_specfit(atlas[k].driftfile,datapath=atlaspath,specfit=fitsun,eigenspec=0)
       flmc = atlas1d_specfit(atlas[k].driftfile,datapath=atlaspath,specfit=fitlmc,eigenspec=1)

       specdata = reform([ [fsuper], [fsun], [flmc] ])
       specfit = [ [ [fitsuper] ], [ [fitsun] ], [ [fitlmc] ] ]

       ageZ[k].chi2 = specdata.continuum_chi2
       bestchi2 = min(ageZ[k].chi2,best)
       ageZ[k].best_indx = best

       other = lindgen(3)
       remove, best, other

       ageZ[k].Z_TEMPLATE = Z[best] ; 12+log(O/H) = 10^(ageZ.Z_TEMPLATE/0.02)+8.69
       ageZ[k].OII_3727 = specdata.OII_3727
       ageZ[k].OIII_5007 = specdata.OIII_5007
       ageZ[k].H_ALPHA = specdata.H_ALPHA
       ageZ[k].OII_3727_EW = specdata.OII_3727_EW
       ageZ[k].OIII_5007_EW = specdata.OIII_5007_EW
       ageZ[k].H_ALPHA_EW = specdata.H_ALPHA_EW
       
; generate the plot

       wave = reform(specfit[*,0,best])/(1+atlas[k].z_obj)
       mask = emission_mask(wave,good=good,/telluric)
       flux = reform(specfit[*,1,best])-reform(specfit[*,3,best])

       spec1 = reform(specfit[*,2,best])
       spec2 = reform(specfit[*,2,other[0]])
       spec3 = reform(specfit[*,2,other[1]])

       ytitle = 'Flux Density [10^{-15} '+flam_units()+']'
       xtitle = 'Rest Wavelength ['+angstrom()+']'
       Zlegend = ['Best Fit '+strZ[best],strZ[other]]

       yrange = minmax(scale*flux)*[1,1.3]

       djs_plot, wave, scale*flux, xsty=3, ysty=3, ps=10, xthick=2.0, $
         ythick=2.0, xtickname=replicate(' ',10), yrange=yrange, $
         charsize=1.5, charthick=2.0, position=pos[*,0], thick=3.0
       djs_oplot, wave, scale*spec2, ps=10, color='red', thick=2.0
       djs_oplot, wave, scale*spec3, ps=10, color='green', thick=2.0
       djs_oplot, wave, scale*spec1, ps=10, color='light blue', thick=2.0
       legend, title[k], /left, /top, box=0, charsize=2.0, charthick=2.0
       legend, Zlegend, /right, /top, box=0, charsize=1.3, charthick=2.0, $
         line=0, color=djs_icolor(['light blue','red','green']), thick=4.0
       
       yrange = max(scale*abs(flux[good]-spec1[good]))*[-0.75,0.75]
       
       djs_plot, wave, scale*(flux-spec1), xsty=3, ysty=3, ps=10, xthick=2.0, $
         ythick=2.0, charsize=1.5, charthick=2.0, position=pos[*,1], thick=3.0, $
         xtitle=xtitle, /noerase, color='grey', yrange=yrange
       djs_oplot, wave, scale*(smooth(flux-spec2,8)), ps=10, color='red', thick=4.0
       djs_oplot, wave, scale*(smooth(flux-spec3,8)), ps=10, color='green', thick=4.0
       
       for i = 0L, nline-1L do djs_oplot, linewave[i]*[1,1], !y.crange, line=2, $
         thick=4.0, color='purple'
       
       xyouts, 0.05, 0.52, textoidl(ytitle), /normal, align=0.5, orientation=90, $
         charsize=1.5, charthick=2.0

       if (not keyword_set(postscript)) and (ngalaxy ne 1L) then begin
          splog, 'Press any key to continue.'
          cc = get_kbrd(1)
       endif
       
    endfor
    splog, format='("Total time for RUN_AGE_METALLICITY = ",G0," minutes.")', (systime(1)-stime0)/60.0

    im_openclose, postscript=postscript, /close

    if keyword_set(postscript) then begin
       mwrfits, ageZ, 'age_metallicity.fits', /create
       spawn, ['gzip -f age_metallicity.fits'], /sh
    endif
    
stop

return
end    
