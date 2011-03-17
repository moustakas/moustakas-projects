function bootes_stellarlocus_project_kurucz, filterlist
; jm09jul24ucsd - project specified bandpasses onto the Kurucz model
; spectra  
    
    kuruczfile = getenv('IDLSPEC2D_DIR')+'/etc/kurucz_stds_v5.fit'
    flux = mrdfits(kuruczfile,0,hdr,/silent)
    info = mrdfits(kuruczfile,1,/silent)
    keep = where(info.feh eq -1.0 and info.g eq 4.5)
    flux = flux[*,keep]

    nx = n_elements(flux[*,0])
    wave = 10d^(sxpar(hdr, 'CRVAL1') + findgen(nx)*sxpar(hdr, 'CD1_1'))
    nband = n_elements(filterlist)
    nstar = n_elements(flux[0,*])

    lambda = k_lambda_to_edges(wave)
    maggies = fltarr(nband,nstar)
    for i = 0L, nstar-1L do maggies[*,i] = k_project_filters(lambda,$
      flux[*,i],filterlist=filterlist,/silent)
    
return, maggies
end

pro bootes_stellarlocus_scatterplot, xall, yall, xx, yy, xk, yk, $
  xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, $
  type=type, feh=feh, label=label, _extra=extra

; I,II = supergiant; III,IV = giant; V = dwarf    
;   dwarf = where(strmatch(strtrim(type,2),'*V*') or $
;     strmatch(strtrim(type,2),'*IV*'))
;   giant = where(strmatch(strtrim(type,2),'*I*') and $
;     (strmatch(strtrim(type,2),'*IV*') eq 0))

    poordwarf = where(strmatch(strtrim(type,2),'*V*') and $
      (strmatch(strtrim(type,2),'*IV*') eq 0) and (feh lt 0.0))
    richdwarf = where(strmatch(strtrim(type,2),'*V*') and $
      (strmatch(strtrim(type,2),'*IV*') eq 0) and (feh ge 0.0))
    poorgiant = where(strmatch(strtrim(type,2),'*I*') and (feh lt 0.0))
    richgiant = where(strmatch(strtrim(type,2),'*I*') and (feh ge 0.0))

    hogg_scatterplot, xall, yall, xtitle=xtitle, ytitle=ytitle, $
      xrange=xrange, yrange=yrange, /xsty, /ysty, xthick=3.0, ythick=3.0, $
      charsize=1.5, charthick=3.0, /outliers, outcolor=djs_icolor(''), $
      /internal, _extra=extra
    legend, textoidl(label), /left, /top, box=0, charsize=1.2, charthick=3.0

    legend, ['Pickles'], /right, /bottom, box=0, charsize=1.2, $
      charthick=3.0, psym=[6], thick=5.0, color=djs_icolor(['dark green'])
    djs_oplot, xx[richdwarf], yy[richdwarf], psym=6, symsize=1.0, color='dark green', thick=5.0

;   legend, ['Metal-Rich Giant','Metal-Poor Giant','Metal-Rich Dwarf',$
;     'Metal-Poor Dwarf','Kurucz Models'], /right, /bottom, box=0, $
;     charsize=1.0, charthick=3.0, psym=[1,1,6,6,4], thick=5.0, $
;     color=djs_icolor(['red','blue','dark green','cyan blue',''])
;   djs_oplot, xk, yk, psym=4, color='', thick=5.0
;   djs_oplot, xx[richgiant], yy[richgiant], psym=1, symsize=1.3, color='red', thick=5.0
;   djs_oplot, xx[poorgiant], yy[poorgiant], psym=1, symsize=1.3, color='blue', thick=5.0
;   djs_oplot, xx[richdwarf], yy[richdwarf], psym=6, symsize=1.0, color='dark green', thick=5.0
;   djs_oplot, xx[poordwarf], yy[poordwarf], psym=6, symsize=1.0, color='cyan blue', thick=5.0
    
return
end

function bootes_stellarlocus_project_pickles, filterlist, pickles=pickles
; jm09jul24ucsd - project specified bandpasses onto the Pickles+98
;   spectra 
    
    if (n_elements(pickles) eq 0) then pickles = read_98pickles()
    nband = n_elements(filterlist)
    nstar = n_elements(pickles)
    
    lambda = k_lambda_to_edges(pickles[0].wave)
    nstar = n_elements(pickles)
    maggies = fltarr(nband,nstar)
    for i = 0L, nstar-1L do maggies[*,i] = k_project_filters(lambda,$
      pickles[i].flux,filterlist=filterlist,/silent)
    
return, maggies
end

function bootes_stellarlocus_colors, maggies, filterlist=filterlist
; jm09jul28ucsd - support routine
    case filterlist of
       'BRI': begin
          colors = {$
            br: reform(-2.5*alog10(maggies[0,*]/maggies[1,*])), $
            ri: reform(-2.5*alog10(maggies[1,*]/maggies[2,*]))}
       end
    endcase

return, colors
end

pro bootes_stellarlocus, nozpoffset=nozpoffset
; jm09aug25ucsd - 

    common stellarlocus_bootes, iband1, bwband1, rband1, bootes1

    bootespath = getenv('RESEARCHPATH')+'/data/bootes/'
    if keyword_set(nozpoffset) then $
      psfile = bootespath+'bootes_stellarlocus_nozpoffset.ps' else $
      psfile = bootespath+'bootes_stellarlocus_zpoffset.ps'

; read the photometric catalogs
    if (n_elements(iband1) eq 0) then begin
       bootesfile = bootespath+'bootes_I.fits.gz'
       splog, 'Reading '+bootesfile
       iband1 = mrdfits(bootesfile,1)

       bootesfile = bootespath+'bootes_Bw.fits.gz'
       splog, 'Reading '+bootesfile
       bwband1 = mrdfits(bootesfile,1)

       bootesfile = bootespath+'bootes_R.fits.gz'
       splog, 'Reading '+bootesfile
       rband1 = mrdfits(bootesfile,1)
    endif

; select stars
    istar = where($
      (iband1.class_star gt 0.95) and $
      (iband1.flag_duplicate eq 0) and $
      (iband1.flag_subfield eq 1) and $
      (iband1.mag_psf gt 18.0) and (iband1.mag_psf lt 20.5) and $
      (bwband1.mag_psf gt 0.0) and (bwband1.mag_psf lt 90.0) and $
      (rband1.mag_psf gt 0.0) and (rband1.mag_psf lt 90.0),nstar)
    splog, 'N = ', nstar

    tags = tag_names(iband1)
    iband = im_struct_trimtags(iband1[istar],select=tags,newtags='I_'+tags)
    tags = tag_names(bwband1)
    bwband = im_struct_trimtags(bwband1[istar],select=tags,newtags='Bw_'+tags)
    tags = tag_names(rband1)
    rband = im_struct_trimtags(rband1[istar],select=tags,newtags='R_'+tags)
    bootes1 = struct_addtags(struct_addtags(temporary(iband),temporary(rband)),temporary(bwband))

    bootes_to_maggies, bootes1, maggies, ivarmaggies, $
      nozpoffset=nozpoffset, filterlist=filterlist, /psf
    maggies = maggies[0:2,*]
    ivarmaggies = ivarmaggies[0:2,*]
    filterlist = filterlist[0:2]

; cut the sample    
    label = ['NDWFS/BOOTES - BwRI','18<I_{PSF}<20.5','Class Star>0.95']
    if keyword_set(nozpoffset) then $
      label = [label,'No zeropoint offsets'] else $
      label = [label,'\Delta'+'Zpt = [-0.129,-0.066,-0.029]']

; compute all the color combinations
    colors = bootes_stellarlocus_colors(maggies,filterlist='BRI')

; Pickles+98    
    pmaggies = bootes_stellarlocus_project_pickles(filterlist,pickles=pickles)
    pcolors = bootes_stellarlocus_colors(pmaggies,filterlist='BRI')

; Kurucz models
    kmaggies = bootes_stellarlocus_project_kurucz(filterlist)
    kcolors = bootes_stellarlocus_colors(kmaggies,filterlist='BRI')

; make the plots        
    brrange = [-1.0,3.9]
    rirange = [-1.0,2.5]

    im_plotconfig, 0, pos, psfile=psfile
    bootes_stellarlocus_scatterplot, colors.ri, colors.br, pcolors.ri, pcolors.br, $
      kcolors.ri, kcolors.br, xtitle='R-I', ytitle='Bw-R', xrange=rirange, $
      yrange=brrange, type=pickles.type, feh=pickles.feh, label=label, position=pos[*,0]
    im_plotconfig, psfile=psfile, /psclose, /gzip

return    
end
