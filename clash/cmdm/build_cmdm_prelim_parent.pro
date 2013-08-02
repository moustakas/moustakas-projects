pro build_cmdm_prelim_parent, clobber=clobber
; jm12feb27ucsd - get cutouts of all the cluster members based on the
; PRELIMINARY version of this project!
    
    path = clash_path()+'projects/cmdm/'

    cl = 'macs1206'
    zcl = 0.440
    dz = 0.01

    zcat1 = read_clash_catalog(cl,/redshift)
    cat1 = rsex(path+'data/'+'macsj1206_2009_IC_O1.cat')
    subarucat1 = rsex(path+'data/'+'MACSJ1206_Subaru.cat')

;   djs_plot, subarucat1.ra, subarucat1.dec, psym=3, xsty=3, ysty=3
;   djs_oplot, cat1.xwin_world, cat1.ywin_world, psym=3, color='orange'
;   djs_oplot, zcat1.ra, zcat1.dec, psym=6, symsize=0.5, color='yellow'

; my catalog - 1193 matches, 293 members
; official catalog - 1234 matches (2 without I-band photometry), 298 members
    spherematch, cat1.xwin_world, cat1.ywin_world, zcat1.ra, zcat1.dec, 1D/3600, m1, m2
;   spherematch, cat1.ra, cat1.dec, zcat1.ra, zcat1.dec, 1D/3600, m1, m2

    gal = where(zcat1[m2].z gt 0.4 and zcat1[m2].z lt 0.49 and cat1[m1].flux_radius2 gt 4,ngal)
;   gal = where(zcat1[m2].z gt 0.01 and zcat1[m2].z lt 1 and cat1[m1].flux_radius2 gt 4,ngal)
    mem = where(zcat1[m2[gal]].z lt zcl+dz and zcat1[m2[gal]].z gt zcl-dz,nmem)

    splog, 'Final sample, probable members ', ngal, nmem
    cat = struct_addtags(cat1[m1[gal]],zcat1[m2[gal]])

    
;   subarucat = subarucat1[m1[gal]],zcat1[m2[gal]]

; sort by size    
    srt = reverse(sort(cat.flux_radius2))
    cat = cat[srt]

    im_mwrfits, cat, path+cl+'_cat.fits', clobber=clobber

; some QAplots    
    djs_plot, cat1.mag_auto, cat1.flux_radius2, psym=3, xsty=3, ysty=3, $
      xrange=[16,24], yrange=[0,50]
    djs_oplot, cat.mag_auto, cat.flux_radius2, psym=6, symsize=0.8, color='red'
    djs_oplot, cat1[m1[gal[mem]]].mag_auto, cat1[m1[gal[mem]]].flux_radius2, $
      psym=symcat(16), symsize=0.8, color='yellow'
    djs_oplot, !x.crange, 4*[1,1]
;   djs_plot, cat.mag_auto, cat.class_star, psym=6, xsty=3, ysty=3, yrange=[0,1]

stop    
    
return
end
    
