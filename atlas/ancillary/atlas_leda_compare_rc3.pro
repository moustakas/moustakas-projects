pro atlas2d_leda_compare_rc3, postscript=postscript
; jm05may16uofa
; compare data from LEDA to the RC3

    datapath = atlas_path(/ned)
    
    leda = parse_leda('atlas2d_leda.dat',ledapath=datapath)
    ra = 15.0*leda.al2000
    dec = leda.de2000
    
    rc3 = read_rc3()
    rc3ra = 15.0D*im_hms2dec(rc3.ra)
    rc3dec = im_hms2dec(rc3.dec)

    searchrad = 125.0 ; search radius

    splog, 'RC3 search radius = '+string(searchrad,format='(G0.0)')+'".'
    ntot = im_djs_angle_match(ra,dec,rc3ra,rc3dec,dtheta=searchrad/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    match = where(mindx ne -1L,nmatch,comp=nomatch,ncomp=nnomatch)

    niceprint, leda[match].name, leda[match].objname, rc3[mindx[match]].name, $
      rc3[mindx[match]].altname, rc3[mindx[match]].pgc

    subleda = leda[match]
    subrc3 = rc3[mindx[match]]

; generate plots    
    
    if keyword_set(postscript) then begin
       dfpsplot, datapath+'atlas2d_leda_compare_rc3.ps', /square
       postthick = 8.0
    endif else begin
       im_window, 0, xratio=0.4, /square
       postthick = 2.0
    endelse
    
    good = where((subleda.bt ne '') and (subrc3.bt gt -900),ngood)
    plot, subleda[good].bt, subrc3[good].bt, ps=4, xsty=3, ysty=3, $
      xrange=[7.5,16.5], yrange=[7.5,16.5], xthick=postthick, $
      ythick=postthick, charsize=1.8, charthick=postthick, $
      xtitle='B [LEDA]', ytitle='B [RC3]'
    oplot, !x.crange, !y.crange, line=0, thick=postthick
    if (not keyword_set(postscript)) then cc = get_kbrd(1)
    
    good = where((subleda.pa ne '') and (subrc3.pa gt -900),ngood)
    plot, subleda[good].pa, subrc3[good].pa, ps=4, xsty=3, ysty=3, $
      xrange=[0,180], yrange=[0,180], xthick=postthick, $
      ythick=postthick, charsize=1.8, charthick=postthick, $
      xtitle='PA [LEDA]', ytitle='PA [RC3]'
    oplot, !x.crange, !y.crange, line=0, thick=postthick
    if (not keyword_set(postscript)) then cc = get_kbrd(1)
    
    good = where((subleda.incl ne '') and (subrc3.inclination gt -900),ngood)
    plot, subleda[good].pa, subrc3[good].pa, ps=4, xsty=3, ysty=3, $
      xrange=[0,90], yrange=[0,90], xthick=postthick, $
      ythick=postthick, charsize=1.8, charthick=postthick, $
      xtitle='Inclination [LEDA]', ytitle='Inclination [RC3]'
    oplot, !x.crange, !y.crange, line=0, thick=postthick
    if (not keyword_set(postscript)) then cc = get_kbrd(1)
    
    good = where((subleda.t ne '') and (subrc3.t gt -900),ngood)
    plot, subleda[good].t, subrc3[good].t, ps=4, xsty=3, ysty=3, $
      xrange=[-10,15], yrange=[-10,15], xthick=postthick, $
      ythick=postthick, charsize=1.8, charthick=postthick, $
      xtitle='Type [LEDA]', ytitle='Type [RC3]'
    oplot, !x.crange, !y.crange, line=0, thick=postthick
    if (not keyword_set(postscript)) then cc = get_kbrd(1)
    
    if keyword_set(postscript) then dfpsclose

stop    
    
return
end
    
