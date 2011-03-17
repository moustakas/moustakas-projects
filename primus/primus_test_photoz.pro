pro primus_test_photoz, extract, slits, oned, oned_photoz, $
  rest, photoz=photoz, nophotoz=nophotoz, ages_photoz=ages_photoz, $
  ages_nophotoz=ages_nophotoz, analyze=analyze, ps=ps, $
  choosedeep=choosedeep
; jm08maynyu - started
; jm08jul04nyu - added support for AGES template fitting
; rerun 0213 is soft-linked to rerun 1300 (the official reduction)
; echo "primus_test_photoz, /photoz" | idl > & primus_photoz.log &
; echo "primus_test_photoz, /nophotoz" | idl > & primus_nophotoz.log &

    path = im_primus_path()
    rerun = '0215' ; true rerun = '1300'
    rootname = 'd02b0061' ; pretty good spectrophtometry, especially for R>22

    if keyword_set(choosedeep) then begin
       primus_read_1dinputs, rootname, rerun, extract=extract, slits=slits
       alldeep = read_deep2(/kcorr)
       match, alldeep.objno, slits.objno, m1, m2
       mb = alldeep[m1].ubvri_absmag[1]
       ub = alldeep[m1].ubvri_absmag[0]-alldeep[m1].ubvri_absmag[1]
       rmag = alldeep[m1].magr
       these = where(ub gt 0.2 and rmag lt 22.0) & thisindex = m2[these]
       plot, mb, ub, ps=4, xr=[-15,-24], yr=[-0.8,0.8], xsty=1, ysty=1
       djs_oplot, mb[these], ub[these], ps=7, color='red', thick=3
;      deep2_plot_kcorrect_sedfits, alldeep[m1[these]]
    endif else thisindex = [363,1658,1679,1146,2033,1830] ; only works for d02b0061

; AGES templates    
    
    if keyword_set(ages_photoz) then begin
       primus_fit_redshift, rootname, rerun, /ages, /simple, /constpsf, $
         /nostar, /noagn, reweight=0, /chooseslit, photoz=1, /test_photoz, $
         /verbose
    endif

stop    
    
    if keyword_set(ages_nophotoz) then begin
       primus_fit_redshift, rootname, rerun, oiitweak=1, /constpsf, $
         /nostar, /noagn, reweight=0, /chooseslit, photoz=0, /test_photoz, $
         /verbose;index=thisindex, 
    endif

; K-correct templates    
    
    if keyword_set(photoz) then begin
       primus_fit_redshift, rootname, rerun, oiitweak=1, /constpsf, $
         /nostar, /noagn, reweight=0, /chooseslit, photoz=1, /test_photoz, $
         /verbose;index=thisindex, 
    endif

    if keyword_set(nophotoz) then begin
       primus_fit_redshift, rootname, rerun, oiitweak=1, /constpsf, $
         /nostar, /noagn, reweight=0, /chooseslit, photoz=0, /test_photoz, $
         /verbose;index=thisindex, 
    endif

; study the plots    
    
    if keyword_set(analyze) then begin

       if (n_elements(extract) eq 0L) or (n_elements(slits) eq 0L) or $
         (n_elements(oned) eq 0L) or (n_elements(oned_photoz) eq 0L) or $
         (n_elements(rest) eq 0L) then begin
          primus_read_1dinputs, rootname, rerun, extract=extract, slits=slits, oned=oned
          primus_read_1dinputs, rootname, rerun, oned=oned_photoz, /photoz
          rest = primus_get_restframe(slits)
       endif
       photoinfo = primus_get_photomaggies(rootname,slits)

; results       

       these = where(total(oned.zmin_gal,1) gt 0.0)
       these_photoz = where(total(oned_photoz.zmin_gal,1) gt 0.0)
       ngal = n_elements(these)

       rmag = slits[these].mag
       snr = total([[extract[these].sn1],[extract[these].sn2]],2)/2.0 ; mean S/N
       ub = rest[these].ubvri_absmag[0]-rest[these].ubvri_absmag[1]
       
       ztrue = slits[these].currz
       zspec = oned[these].zmin_gal[0]
       zphot = oned_photoz[these_photoz].zmin_gal[0]
       zphot2 = oned_photoz[these_photoz].zmin_gal[1] ; second minimum

       dzspec = (zspec-ztrue)/(1.0+ztrue)
       dzphot = (zphot-ztrue)/(1.0+ztrue)

       dzcut = 0.4
       catacut = 0.05
       rmagcut = 23.5
       faintcut = 23.0

       cataspec = where((abs(dzspec) gt catacut) and (rmag lt rmagcut),$
         ncataspec,comp=goodspec)
       cataphot = where((abs(dzphot) gt catacut) and (rmag lt rmagcut),$
         ncataphot,comp=goodphot)
       faint = where(rmag gt faintcut,comp=bright)

       specstats = im_stats(dzspec[goodspec])
       photstats = im_stats(dzphot[goodphot])
       
;      'F(\delta'+'z/(1+z) > 0.05) = '
       labelspec = [$
         '\Delta = '+strtrim(string(specstats.mean,format='(F12.4)'),2)+$
         '\pm'+strtrim(string(specstats.sigma,format='(F12.4)'),2),$
;        '\Delta_{med} = '+strtrim(string(specstats.median,format='(F12.4)'),2),$
         'F(>'+string(catacut,format='(F4.2)')+') = '+$
         strtrim(string(100.0*ncataspec/float(ngal),$
         format='(F12.1)'),2)+'% ('+string(ncataspec,format='(I0)')+'/'+$
         string(ngal,format='(I0)')+')']
       labelphot = [$
         '\Delta = '+strtrim(string(photstats.mean,format='(F12.4)'),2)+$
         '\pm'+strtrim(string(photstats.sigma,format='(F12.4)'),2),$
;        '\Delta_{med} = '+strtrim(string(photstats.median,format='(F12.4)'),2),$
         'F(>'+string(catacut,format='(F4.2)')+') = '+$
         strtrim(string(100.0*ncataphot/float(ngal),$
         format='(F12.1)'),2)+'% ('+string(ncataphot,format='(I0)')+'/'+$
         string(ngal,format='(I0)')+')']
       
       chi2spec = oned[these].chi2min_gal[0]
       chi2phot = oned_photoz[these_photoz].chi2min_gal[0]
       niceprint, ztrue, zspec, zphot, chi2spec, chi2phot
       
; now make some plots 

       if keyword_set(ps) then begin
          postthick1 = 4.0
          postthick2 = 3.0
          postthick3 = 8.0
       endif else begin
          im_window, 0, xratio=0.5, /square
          postthick1 = 2.0
          postthick2 = 2.0
          postthick3 = 2.0
       endelse

       xpage = 8.5 & ypage = 11.0
       psname = path+'primus_test_photoz.ps'
       im_openclose, psname, postscript=keyword_set(ps), $
         xsize=xpage, ysize=ypage

       pagemaker, nx=2, ny=4, xpage=xpage, ypage=ypage, $
         xspace=0.0, yspace=0.55, xmargin=[1.1,0.4], $
         ymargin=[0.1,1.25], width=3.5*[1,1], height=[2.6,1.8,1.8,1.8], $
         /normal, position=pos2
       pagemaker, nx=1, ny=1, xpage=xpage, ypage=ypage, $
         xspace=0.0, yspace=0.0, xmargin=[1.4,0.3], $
         ymargin=[0.7,1.1], width=6.8, height=6.8, $
         /normal, position=pos1
       pagemaker, nx=1, ny=2, xpage=xpage, ypage=ypage, $
         xspace=0.0, yspace=1.0, xmargin=[1.4,0.3], $
         ymargin=[0.7,1.1], width=6.8, height=3.4*[1,1], $
         /normal, position=pos3
       
; compare the spectro+photoz redshifts against DEEP2

       charsize1 = 1.2
       zrange = [-0.1,1.3]
       zaxis = findgen((2.0-(-1.0))/0.01)*0.01+(-1.0)
       residrange = 0.35*[-1,1]
       psize1 = 0.55
       psize2 = 0.15
       
       plot, [0], [0], /nodata, charsize=charsize1, charthick=postthick2, $
         xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
         xrange=zrange, yrange=zrange, xtickinterval=0.2, $
         xtitle='Redshift [DEEP2]', ytitle='Redshift [PRIMUS]', $
         title='No PHOTOZ', $
;        title=strtrim(rootname,2)+' no photoz', $
         position=pos2[*,0]
       oplot, zaxis, zaxis, line=0, thick=2.0
       plotsym, 0, psize1, fill=1, thick=postthick1
       djs_oplot, ztrue, zspec, ps=8, color='dark green'

       plot, [0], [0], /nodata, /noerase, charsize=charsize1, charthick=postthick2, $
         xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
         xrange=zrange, yrange=zrange, xtickinterval=0.2, $
         xtitle='Redshift [DEEP2]', ytitle='', ytickname=replicate(' ',10), $
         title='PHOTOZ', $
;        title=strtrim(rootname,2)+' with photoz', $
         position=pos2[*,1]
       oplot, zaxis, zaxis, line=0, thick=2.0
       oplot, zaxis, zaxis-dzcut, line=2, thick=2.0
;      plotsym, 8, psize1, fill=1
;      djs_oplot, ztrue[bright], zphot[bright], ps=8, color='red'
;      plotsym, 8, psize1, fill=0
;      djs_oplot, ztrue[faint], zphot[faint], ps=8, color='blue'
       plotsym, 8, psize1, fill=1
       djs_oplot, ztrue, zphot, ps=8, color='red'
;      plotsym, 3, psize2, fill=1
;      djs_oplot, ztrue, zphot2, ps=8, color='blue'

;      xyouts, pos2[2,0], pos2[1,0]-0.04, 'Redshift [DEEP2]', /normal, $
;        charsize=charsize1, charthick=postthick2, align=0.5
       
; residuals vs R magnitude

       rmagrange = [19.5,23.8]
       
       plot, [0], [0], /nodata, /noerase, charsize=charsize1, charthick=postthick2, $
         xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
         xrange=rmagrange, yminor=5, yrange=residrange, $
         xtitle='R (mag)', ytitle=textoidl('\delta'+'z/(1+z)'), position=pos2[*,2]
       oplot, !x.crange, [0,0], line=0, thick=2.0
       oplot, !x.crange, catacut*[1,1], line=1, thick=2.0
       oplot, !x.crange, -catacut*[1,1], line=1, thick=2.0
       plotsym, 0, psize1, fill=1
       djs_oplot, rmag, dzspec, ps=8, color='dark green'
;      legend, textoidl(labelspec), /right, /top, box=0, $
;        charsize=1.3, charthick=postthick2, /clear, spacing=1.9, margin=0
       
       plot, [0], [0], /nodata, /noerase, charsize=charsize1, charthick=postthick2, $
         xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
         xrange=rmagrange, yminor=5, yrange=residrange, $
         xtitle='R (mag)', ytitle='', ytickname=replicate(' ',10), position=pos2[*,3]
       oplot, !x.crange, [0,0], line=0, thick=2.0
       oplot, !x.crange, catacut*[1,1], line=1, thick=2.0
       oplot, !x.crange, -catacut*[1,1], line=1, thick=2.0
       plotsym, 8, psize1, fill=1
       djs_oplot, rmag, dzphot, ps=8, color='red'
;      plotsym, 8, psize1, fill=1
;      djs_oplot, rmag[bright], dzphot[bright], ps=8, color='red'
;      plotsym, 8, psize1, fill=0
;      djs_oplot, rmag[faint], dzphot[faint], ps=8, color='blue'
;      legend, textoidl(labelphot), /right, /top, box=0, $
;        charsize=1.3, charthick=postthick2, /clear, spacing=1.9, margin=0
       
; residuals vs S/N

       snrrange = [0,49.0]
       
       plot, [0], [0], /nodata, /noerase, charsize=charsize1, charthick=postthick2, $
         xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
         xrange=snrrange, yminor=5, yrange=residrange, $
         xtitle='<S/N>', ytitle=textoidl('\delta'+'z/(1+z)'), position=pos2[*,4]
       oplot, !x.crange, [0,0], line=0, thick=2.0
       oplot, !x.crange, catacut*[1,1], line=1, thick=2.0
       oplot, !x.crange, -catacut*[1,1], line=1, thick=2.0
       plotsym, 0, psize1, fill=1
       djs_oplot, snr, dzspec, ps=8, color='dark green'
       
       plot, [0], [0], /nodata, /noerase, charsize=charsize1, charthick=postthick2, $
         xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
         xrange=snrrange, yminor=5, yrange=residrange, $
         xtitle='<S/N>', ytitle='', ytickname=replicate(' ',10), position=pos2[*,5]
       oplot, !x.crange, [0,0], line=0, thick=2.0
       oplot, !x.crange, catacut*[1,1], line=1, thick=2.0
       oplot, !x.crange, -catacut*[1,1], line=1, thick=2.0
       plotsym, 8, psize1, fill=1
       djs_oplot, snr, dzphot, ps=8, color='red'
       
; residuals vs U-B color

       ubrange = [-0.5,0.5]
       
       plot, [0], [0], /nodata, /noerase, charsize=charsize1, charthick=postthick2, $
         xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
         xrange=ubrange, yminor=5, yrange=residrange, $
         xtitle='U-B', ytitle=textoidl('\delta'+'z/(1+z)'), position=pos2[*,6]
       oplot, !x.crange, [0,0], line=0, thick=2.0
       oplot, !x.crange, catacut*[1,1], line=1, thick=2.0
       oplot, !x.crange, -catacut*[1,1], line=1, thick=2.0
       plotsym, 0, psize1, fill=1
       djs_oplot, ub, dzspec, ps=8, color='dark green'
       
       plot, [0], [0], /nodata, /noerase, charsize=charsize1, charthick=postthick2, $
         xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
         xrange=ubrange, yminor=5, yrange=residrange, $
         xtitle='U-B', ytitle='', ytickname=replicate(' ',10), position=pos2[*,7]
       oplot, !x.crange, [0,0], line=0, thick=2.0
       oplot, !x.crange, catacut*[1,1], line=1, thick=2.0
       oplot, !x.crange, -catacut*[1,1], line=1, thick=2.0
       plotsym, 8, psize1, fill=1
       djs_oplot, ub, dzphot, ps=8, color='red'
       
       if (not keyword_set(ps)) then cc = get_kbrd(1)

; individual spectra and chi2 curves

;      check = lindgen(ngal) & ncheck = ngal
       check = where(((ztrue-zphot) gt dzcut) and (rmag lt faintcut),ncheck) ; no abs()
;      check = where((strmatch(slits[these].comment,'*qso*',/fold) eq 1B) or $
;       (strmatch(slits[these].comment,'*broad*',/fold) eq 1B),ncheck)
       
;      jj = 5L & thisindex = 127
;      plot, oned_photoz[these_photoz[check[jj]]].zgrid_gal, $
;        oned_photoz[these_photoz[check[jj]]].chi2_gal, ysty=3

       plotscale = 1E17
       light = 2.99792458D18
       weff = photoinfo[0].weff
       filters = strtrim(photoinfo[0].filterlist,2)
       fnu2flam = 10^(-0.4*48.6)*light/weff^2.0
       bandindx = 1L

; pull out the primus spectra and best-fitting model (at the first
; minimum)
       flux = primus2flambda(extract[these[check]],wave=wave,$
         counts2flam=counts2flam)
       modelflux = flux*0.0
       modelflux[*,0,*] = oned_photoz[these_photoz[check]].fmod_gal1[*,0]/$
         (counts2flam[*,0,*]+(counts2flam[*,0,*] eq 0.0))*(counts2flam[*,0,*] ne 0.0)
       modelflux[*,1,*] = oned_photoz[these_photoz[check]].fmod_gal2[*,0]/$
         (counts2flam[*,1,*]+(counts2flam[*,1,*] eq 0.0))*(counts2flam[*,1,*] ne 0.0)
       
       for jj = 0L, ncheck-1L do begin
; calculate some preliminaries and the xyranges
          norm1 = photoinfo[these[check[jj]]].maggies[bandindx]/$
            (k_project_filters(k_lambda_to_edges(wave[*,0,jj]),$
            flux[*,0,jj],filterlist=filters,/silent))[bandindx]
          norm2 = photoinfo[these[check[jj]]].maggies[bandindx]/$
            (k_project_filters(k_lambda_to_edges(wave[*,1,jj]),$
            flux[*,1,jj],filterlist=filters,/silent))[bandindx]
          modelnorm1 = photoinfo[these[check[jj]]].maggies[bandindx]/$
            (k_project_filters(k_lambda_to_edges(wave[*,0,jj]),$
            modelflux[*,0,jj],filterlist=filters,/silent))[bandindx]
          modelnorm2 = photoinfo[these[check[jj]]].maggies[bandindx]/$
            (k_project_filters(k_lambda_to_edges(wave[*,1,jj]),$
            modelflux[*,1,jj],filterlist=filters,/silent))[bandindx]
          yrange = plotscale*[$
            im_min(modelnorm1*modelflux[*,0,jj],sigrej=4.0,/ignorezero)<$
            im_min(modelnorm2*modelflux[*,1,jj],sigrej=4.0,/ignorezero),$
            im_max(modelnorm1*modelflux[*,0,jj],sigrej=4.0)>$
            im_max(modelnorm2*modelflux[*,1,jj],sigrej=4.0)]
          xrange = [3500,9900]
; now make the plot
          djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
            xthick=postthick1, ythick=postthick1, xsty=1, ysty=3, $
            thick=postthick1, xtitle='Observed Wavelength (\AA)', $
            ytitle='Flux (10^{-17} '+flam_units()+')', $
            charsize=1.5, charthick=postthick2, position=pos3[*,0], $
            title='DEEP2/'+string(slits[these[check[jj]]].objno,format='(I0)')+$
            ', z = '+strtrim(string(ztrue[check[jj]],format='(F12.4)'),2)+$
            ', R = '+strtrim(string(rmag[check[jj]],format='(F12.2)'),2)
; PRIMUS spectra
          djs_oplot, wave[*,0,jj], plotscale*norm1*flux[*,0,jj], line=1, $
               color='purple', thick=postthick1
          djs_oplot, wave[*,1,jj], plotscale*norm2*flux[*,1,jj], line=1, $
               color='red', thick=postthick1
; best-fitting model spectra
          djs_oplot, wave[*,0,jj], plotscale*modelnorm1*modelflux[*,0,jj], line=0, $
               thick=postthick1
          djs_oplot, wave[*,1,jj], plotscale*modelnorm2*modelflux[*,1,jj], line=0, $
               thick=postthick1
; BRI photometry          
          plotsym, 8, 2.0, fill=1, thick=postthick3
          djs_oplot, weff, plotscale*photoinfo[these[check[jj]]].maggies*fnu2flam, $
            ps=8, color=djs_icolor('dark green'), thick=postthick1
; P(z) plots
;         pzspec = exp(-0.5*(oned[these[check[jj]]].chi2_gal-chi2spec[check[jj]])
;         pzspec = pzspec/total(pzspec)
;         yrange = [min(oned_photoz[these_photoz[check[jj]]].chi2_gal)<$
;           min(oned[these[check[jj]]].chi2_gal),$
;           max(oned_photoz[these_photoz[check[jj]]].chi2_gal)>$
;           max(oned[these[check[jj]]].chi2_gal)]
;         yrange = yrange*[0.9,1.05]
;         djs_plot, [0], [0], /nodata, xrange=[0.0,1.19], yrange=yrange, $
;           xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
;           thick=postthick1, xtitle='Redshift', ytitle='\chi^{2}', $
;           charsize=1.8, charthick=postthick2, position=pos1[*,0]
;         djs_oplot, zphot[check[jj]]*[1,1], !y.crange, line=0, color='red', thick=postthick3
;         djs_oplot, zspec[check[jj]]*[1,1], !y.crange, line=5, color='dark green', thick=postthick3
;         djs_oplot, ztrue[check[jj]]*[1,1], !y.crange, line=3, color='blue', thick=postthick3
;         djs_oplot, oned_photoz[these_photoz[check[jj]]].zgrid_gal, $
;           oned_photoz[these_photoz[check[jj]]].chi2_gal, $
;           line=0, ps=10, thick=postthick1, color='purple'
;         djs_oplot, oned[these[check[jj]]].zgrid_gal, oned[these[check[jj]]].chi2_gal, $
;           line=2, ps=10, thick=postthick1
; chi2 plots          
          yrange = [min(oned_photoz[these_photoz[check[jj]]].chi2_gal)<$
            min(oned[these[check[jj]]].chi2_gal),$
            max(oned_photoz[these_photoz[check[jj]]].chi2_gal)>$
            max(oned[these[check[jj]]].chi2_gal)]
          yrange = yrange*[0.9,1.05]
          djs_plot, [0], [0], /nodata, /noerase, xrange=[0.0,1.19], yrange=yrange, $
            xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
            thick=postthick1, xtitle='Redshift', ytitle='\chi^{2}', $
            charsize=1.5, charthick=postthick2, position=pos3[*,1]
          djs_oplot, zphot[check[jj]]*[1,1], !y.crange, line=0, color='red', thick=postthick3
          djs_oplot, zspec[check[jj]]*[1,1], !y.crange, line=5, color='dark green', thick=postthick3
          djs_oplot, ztrue[check[jj]]*[1,1], !y.crange, line=3, color='blue', thick=postthick3
          djs_oplot, oned_photoz[these_photoz[check[jj]]].zgrid_gal, $
            oned_photoz[these_photoz[check[jj]]].chi2_gal, $
            line=0, ps=10, thick=postthick1, color='purple'
          djs_oplot, oned[these[check[jj]]].zgrid_gal, oned[these[check[jj]]].chi2_gal, $
            line=2, ps=10, thick=postthick1
          legend, textoidl(['z_{DEEP2}','z_{PRIMUS}','z_{PHOTOZ}']), /left, /top, $
            box=0, charsize=1.5, charthick=postthick2, line=[3,5,0], $
            color=djs_icolor(['blue','dark green','red']), thick=postthick1, $
            /clear
          legend, textoidl(['\chi^{2}_{ PRIMUS}','\chi^{2}_{ PHOTOZ}']), $
            /right, /top, box=0, charsize=1.5, charthick=postthick2, $
            line=[2,0], thick=postthick1, /clear, color=djs_icolor(['default','purple'])
       endfor

;;; distribution of redshift residuals
;;
;;       im_plothist, dzspec, bin=0.01, x1, y1, /noplot
;;       im_plothist, dzphot, bin=0.01, x2, y2, /noplot
;;       yrange = [0.0,max(y1)>max(y2)]
;;
;;       plot, [0], [0], /nodata, charsize=2.0, charthick=postthick2, $
;;         xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
;;         xrange=[-0.5,0.5], yrange=yrange, ytitle='Number of Galaxies', $
;;         xtitle=textoidl('(z_{primus}-z_{deep2})/(1+z_{deep2})'), $
;;         title=strtrim(rootname,2), position=pos1[*,0]
;;       
;;       im_plothist, dzspec, bin=0.01, /overplot, line=0, thick=postthick1, $
;;         color=djs_icolor('dark green')
;;       im_plothist, dzphot, bin=0.01, /overplot, line=2, thick=postthick1, $
;;         color=djs_icolor('red')

;      plot, [0], [0], /nodata, charsize=1.7, charthick=postthick2, $
;        xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
;        xrange=[-0.05,1.1], yrange=[-0.05,1.1], $
;        xtitle='Redshift [DEEP2]', ytitle='Redshift [PRIMUS]', $
;        title=strtrim(rootname,2), position=pos1[*,0]
;      oplot, !x.crange, !y.crange, line=0, thick=2.0
;      im_legend, ['No PHOTOZ','With PHOTOZ'], psym=[108,106], $
;        color=djs_icolor(['dark green','red']), /left, /top, box=0, $
;        charsize=1.3, charthick=postthick2, fill=[0,1], $
;        spacing=2.2, symthick=postthick1, symsize=1.5
;      im_legend, ['No PHOTOZ (first minimum)','With PHOTOZ (first minimum)',$
;        'With PHOTOZ (second minimum)'], psym=[108,106,115], $
;        color=djs_icolor(['dark green','red','blue']), /left, /top, box=0, $
;        charsize=1.3, charthick=postthick2, fill=[0,1,1], $
;        spacing=2.2, symthick=postthick1, symsize=1.5
;      plotsym, 0, 1.5, fill=0, thick=postthick1
;      djs_oplot, ztrue, zspec, ps=8, color='dark green'
;      plotsym, 8, 1.5, fill=1
;      djs_oplot, ztrue, zphot, ps=8, color='red'
;      plotsym, 3, 1.5, fill=1
;      djs_oplot, ztrue, zphot2, ps=8, color='blue'
;      for ii = 0L, ngal-1L do oplot, [ztrue[ii],ztrue[ii]], $
;        [zspec[ii],zphot[ii]], line=5, thick=postthick2
;      for ii = 0L, ngal-1L do oplot, [ztrue,ztrue], $
;        [zphot,zphot2], line=5, thick=postthick2
       
       im_openclose, postscript=keyword_set(ps), /close

stop
       
    endif
       
return
end
