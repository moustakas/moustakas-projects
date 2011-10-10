pro stanford_render_priorplot, sample, noevol=noevol, super=super, xrange=xrange, $
  yrange=yrange, dolegend=dolegend, textcolor=textcolor, ratio=ratio, keynote=keynote
; render the plot that shows the variance in the MF due to different
; prior assumptions where SAMPLE='all', 'quiescent', or 'active'

    mfpath = mf_path(/mfs)
;   superstring = string(super.supergrid,format='(I2.2)')
    nsuper = n_elements(super)

; plotting preferences    
    label = ['BC03','Ma05','FSPS','Pegase']
;   label = [$
;     'BC03-Zrange-CF00-P_{burst}=0.1',$ ; 1
;     'BC03-Zrange-Calz-P_{burst}=0.1',$ ; 2
;     'BC03-Zrange-CF00-P_{burst}=0.0',$ ; 3
;     'BC03-Zrange-CF00-P_{burst}=0.5',$ ; 4
;     'Ma05-Zrange-CF00-P_{burst}=0.1',$ ; 5
;     'FSPS-Zrange-CF00-P_{burst}=0.1',$ ; 6
;     'BC03-Z_{'+sunsymbol()+'}-CF00-P_{burst}=0.1',$ ; 7
;     'BC03-Zrange-No dust-P_{burst}=0.1'] ; 8
      
;   label = 'Prior set '+string(super.supergrid,format='(I0)')
    superfit_color = ['black','red','navy blue','navy blue'];,'navy blue','navy blue','navy blue','navy blue']
    supermf_color = ['powder blue','tomato','dodger blue','orange'];,'forest green','firebrick','orchid','navy']
    superpsym = [6,9,4,1];,5,2,11,1]
    superline = [0,0,0,0];,0,0,0,0]

;; for comparison, read the supergrid-average mass function
;    mfdata = read_mf_vmax('sdss',quiescent=(sample eq 'quiescent'),$
;      active=(sample eq 'active'),/silent)
;    these = where((mfdata.phi gt 0.0) and (mfdata.mass gt xrange[0]) and $
;      (mfdata.mass lt xrange[1]) and (mfdata.phi gt 10^yrange[0]) and $
;      (mfdata.phi lt 10^yrange[1]),nthese)
;    refmass = mfdata.mass[these]
;    refphi = alog10(mfdata.phi[these])
    
    for jj = 0, nsuper-1 do begin
       mfdata = read_mf_vmax('sdss',super=super[jj].supergrid,noevol=noevol,$
         quiescent=(sample eq 'quiescent'),active=(sample eq 'active'),/silent)
       mffit = read_mf_vmax('sdss',super=super[jj].supergrid,noevol=noevol,$
         quiescent=(sample eq 'quiescent'),active=(sample eq 'active'),/bestfit,/silent)
       
       these = where((mfdata.phi gt 0.0) and (mfdata.mass gt xrange[0]) and $
         (mfdata.mass lt xrange[1]) and (mfdata.phi gt 10^yrange[0]) and $
         (mfdata.phi lt 10^yrange[1]) and (mfdata.limit eq 1),nthese)

       mass1 = mfdata.mass[these]
       phi1 = alog10(mfdata.phi[these])
       phierr1 = mfdata.phierr[these]/mfdata.phi[these]/alog(10)

       if (jj eq 0) then begin
          mass = mass1
          phi = phi1
          phierr = phierr1
          bigsuperpsym = replicate(superpsym[jj],nthese)
          bigsupermf_color = replicate(supermf_color[jj],nthese)
       endif else begin
          mass = [mass,mass1]
          phi = [phi,phi1]
          phierr = [phierr,phierr1]
          bigsuperpsym = [bigsuperpsym,replicate(superpsym[jj],nthese)]
          bigsupermf_color = [bigsupermf_color,replicate(supermf_color[jj],nthese)]
       endelse

;      if keyword_set(ratio) then begin ; show the ratio MF
;         good = where((refmass ge min(mass1)) and (refmass le max(mass1)))
;         djs_oplot, refmass[good], refphi-interpol(phi1,mass1,refmass[good]), $
;           psym=symcat(superpsym[jj],thick=3), symsize=0.8, $
;           color=im_color(supermf_color[jj],100)
;      endif
    endfor

; shuffle the indices so the points don't line on top of one another    
    if (keyword_set(ratio) eq 0) then begin
       indx = shuffle_indx(n_elements(mass))
       for pp = 0, n_elements(mass)-1 do begin
          if keyword_set(showerr) then begin
             oploterror, mass[indx[pp]], phi[indx[pp]], phierr[indx[pp]], $
               psym=symcat(bigsuperpsym[indx[pp]],thick=3), symsize=0.8, errthick=!p.thick, $
               color=im_color(bigsupermf_color[indx[pp]],100+jj), $
               errcolor=im_color(bigsupermf_color[indx[pp]],100+jj)
          endif else begin
             plots, mass[indx[pp]], phi[indx[pp]], psym=symcat(bigsuperpsym[indx[pp]],thick=3), $
               symsize=0.8, color=im_color(bigsupermf_color[indx[pp]],100+jj)
          endelse
       endfor
;      djs_oplot, refmass, refphi, psym=symcat(16), color='black', symsize=1.0
    endif 
       
; overplot the best fits
;   djs_oplot, maxis1, alog10(mf_schechter(maxis1,mffit)), thick=8.0, $
;     line=superline[jj], color=im_color(superfit_color[jj],10+jj)

; legend
    if keyword_set(dolegend) then begin
       im_legend, label, /left, /bottom, box=0, charsize=1.4, $
         psym=superpsym, color=supermf_color, symsize=1.2, textcolor=textcolor
    endif
    
return    
end
