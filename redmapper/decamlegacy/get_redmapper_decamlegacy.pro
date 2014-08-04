pro get_redmapper_decamlegacy, clobber=clobber
; jm14jul29siena - get the redmapper cluster that are in the 2014 Aug
; observing footprint

    common com_cat, cat

    path = redmapper_path(/decam)    
    
    if n_elements(cat) eq 0L then begin
       catfile = redmapper_path()+'dr8_run_redmapper_v5.10_lgt5_catalog.fits.gz'
       cat = mrdfits(catfile,1)
    endif

    file = path+'14aug_plan.txt'
    readcol, file, ra, dec, pass, order, block, airmass, $
      time, filt, form='D,D,I,I,I,F,F,A', /silent
    isg = where(filt eq 'g')
    isr = where(filt eq 'r')
    isz = where(filt eq 'z')

    djs_plot, ra, dec, psym=8, xsty=3, ysty=3
    djs_oplot, ra[isg], dec[isg], psym=8, color='green'
    djs_oplot, ra[isr], dec[isr], psym=6, color='orange'
    djs_oplot, ra[isz], dec[isz], psym=7, color='blue'
;   cc = get_kbrd(1)
    
    group = spheregroup(ra,dec,10D)
    ugroup = group[uniq(group,sort(group))]

    grp0 = where(group eq 0)
    djs_plot, ra, dec, psym=8, xsty=3, ysty=3
    djs_oplot, ra[grp0], dec[grp0], psym=8, color='green'
;   cc = get_kbrd(1)
    
; only group 0 is expected to have 3-band coverage
    for ig = 0, 0 do begin 
;   for ig = 0, n_elements(ugroup)-1 do begin
       these = where(ugroup[ig] eq group)
       inra = minmax(ra[these])
       indec = minmax(dec[these])

       indx1 = where(cat.ra gt inra[0] and cat.ra lt inra[1] and $
         cat.dec gt indec[0] and cat.dec lt indec[1],ncl)
       if ig eq 0 then indx = indx1 else indx = [indx,indx1]
    endfor

    rich = where(cat[indx].lambda_chisq gt 50)
    this = 0

    psfile = path+'redmapper_decam_14aug.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0
    
    djs_plot, ra[grp0], dec[grp0], psym=6, xrange=[265,225], $
      yrange=[-7,15], xsty=1, ysty=1, xtitle='\alpha_{J2000} (deg)', $
      ytitle='\delta_{J2000} (deg)', position=pos
    djs_oplot, cat[indx].ra, cat[indx].dec, psym=symcat(9), symsize=0.1, $
      color=cgcolor('forest green')
    djs_oplot, cat[indx[rich]].ra, cat[indx[rich]].dec, psym=symcat(15), $
      color='red'
;   djs_oplot, [cat[indx[rich[this]]].ra], [cat[indx[rich[this]]].dec], $
;     psym=symcat(9), color=cgcolor('navy'), symsize=5

    im_legend, ['DECam Pointings (Aug 2014)','redMaPPer clusters','\lambda>50 clusters'], $
      /left, /top, box=0, psym=[6,16,15], color=['black','forest green','red'], $
      charsize=1.4, margin=0
    im_plotconfig, psfile=psfile, /psclose, /pdf

    outfile = path+'redmapper_decam_14aug.fits'
    im_mwrfits, cat[indx[rich]], outfile, clobber=clobber

stop    
    
return
end
    
