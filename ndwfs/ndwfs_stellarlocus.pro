function ndwfs_stellarlocus_project_kurucz, filterlist
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

pro ndwfs_stellarlocus_scatterplot, xall, yall, xx, yy, xk, yk, $
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

function ndwfs_stellarlocus_project_pickles, filterlist, pickles=pickles
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

function ndwfs_stellarlocus_colors, maggies, filterlist=filterlist
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

pro ndwfs_stellarlocus, nozpoffset=nozpoffset
; jm09aug25ucsd - 

    common stellarlocus_ndwfs, ndwfs1

    photodir = getenv('RESEARCHPATH')+'/data/ndwfs/'
    if keyword_set(nozpoffset) then $
      psfile = photodir+'ndwfs_stellarlocus_nozpoffset.ps' else $
      psfile = photodir+'ndwfs_stellarlocus_zpoffset.ps'

; read the photometric catalog
    if (n_elements(ndwfs1) eq 0) then begin
       photofile = file_search(photodir+'NDWFS_??_??.fits.gz',count=nfile)
       for ii = 0, nfile-1 do begin
          temp = mrdfits(photofile[ii],1)
          keep = where((temp.bw_flags eq 0) and (temp.r_flags eq 0) and $
            (temp.i_flags eq 0),nkeep)
          splog, nkeep
          temp = temp[keep]
          if (ii eq 0) then ndwfs1 = temporary(temp) else $
            ndwfs1 = [temporary(ndwfs1),temporary(temp)]
       endfor
    endif

; select stars
    istar = where($
      (ndwfs1.bw_mag_auto gt 0.0) and (ndwfs1.bw_mag_auto lt 30.0) and $
      (ndwfs1.r_mag_auto gt 17.0) and (ndwfs1.r_mag_auto lt 21.5) and $
      (ndwfs1.i_mag_auto gt 0.0) and (ndwfs1.i_mag_auto lt 30.0) and $
      (ndwfs1.bw_flags eq 0) and (ndwfs1.r_flags eq 0) and (ndwfs1.i_flags eq 0) and $
      (ndwfs1.r_class_star gt 0.95))
    ndwfs = ndwfs1[istar]

    filterlist = ndwfs_filterlist()
    ndwfs_to_maggies, ndwfs, maggies, maggiesivar, nozpoffset=nozpoffset

; cut the sample    
    label = ['NDWFS - BwRI','17<R<21.5','Class Star>0.95']
    if keyword_set(nozpoffset) then $
      label = [label,'No zeropoint offsets'] else $
      label = [label,'\Delta'+'Zpt = [-0.13,-0.07,-0.04]']

; compute all the color combinations
    colors = ndwfs_stellarlocus_colors(maggies,filterlist='BRI')

; Pickles+98    
    pmaggies = ndwfs_stellarlocus_project_pickles(filterlist,pickles=pickles)
    pcolors = ndwfs_stellarlocus_colors(pmaggies,filterlist='BRI')

; Kurucz models
    kmaggies = ndwfs_stellarlocus_project_kurucz(filterlist)
    kcolors = ndwfs_stellarlocus_colors(kmaggies,filterlist='BRI')

; make the plots        
    brrange = [-1.0,3.9]
    rirange = [-1.0,2.5]

    im_plotconfig, 0, pos, psfile=psfile
    ndwfs_stellarlocus_scatterplot, colors.ri, colors.br, pcolors.ri, pcolors.br, $
      kcolors.ri, kcolors.br, xtitle='R-I', ytitle='Bw-R', xrange=rirange, $
      yrange=brrange, type=pickles.type, feh=pickles.feh, label=label, position=pos[*,0]
    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep

return    
end
