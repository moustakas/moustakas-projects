pro ages_match_xbootes, xbootes
; jm09apr29nyu - match the AGES and XBOOTES catalogs

; X-ray; apply three separate cuts before matching to the AGES catalog
; (see Brand et al. 2006 and the XBootes README for more details): (1)
; OPTRANK=1 restricts the sample to the 3213/5318 sources with the
; most probable optical counterpart, including the possibility of *no*
; counterpart; (2) FLAG=1 further restricts the sample to those
; 3140/3213 sources *with* an optical counterpart; and (3) BAYPROB>0.5
; imposes a liberal 50% probability that the match is good, leaving a
; total of 3029/3140 objects; finally, we match to the AGES catalog
; (whereas Brand et al. matched to the full NDWFS imaging survey)
; using a 1" search radius, resulting in a sample of 1538
; objects; note that the 1" search radius is pretty generous
; (1508/1538 match to within 0.5"), but we again want to be
; generous in the matching since we would rather throw out a small
; fraction of our sample that *might* have an AGN rather than retain
; them 

    catpath = ages_path(/mycatalogs)

    ages = mrdfits(ages_path(/catalogs)+'catalog.cat.noguidestars.fits.gz',1)
    ages.ra = ages.ra*15.0D
    ngal = n_elements(ages)

    splog, 'Reading '+catpath+'xbootes/xbootes_cat_xray_opt_IR_21jun_v1.0.txt'
    xray = im_read_fmr(catpath+'xbootes/xbootes_cat_xray_opt_IR_21jun_v1.0.txt')

    xkeep = where((xray.optrank eq 1) and (xray.flag eq 1) and $
      (xray.bayprob gt 0.5))
    m1 = im_spherematch(ages,xray[xkeep],ratagname2='oradeg',$
      dectagname2='odecdeg',radius=1.0,match2=m2)
;   spherematch, ages.ra, ages.dec, xray[xkeep].oradeg, $
;     xray[xkeep].odecdeg, searchrad/3600.0, m1, m2, d12
    splog, 'Found '+string(n_elements(m1),format='(I0)')+$
      ' objects with good X-ray photometry'
;   im_plothist, d12*3600.0, ysty=3, bin=0.05, xr=[-0.1,1.5]
    
;   djs_plot, xray[xkeep].oradeg, xray[xkeep].odecdeg, $
;     ps=6, xsty=3, ysty=3, xr=[217,217.5], yr=[34,34.5]
;   djs_oplot, out.ra, out.dec, ps=6, sym=0.4, color='red'
;   djs_oplot, xray[xkeep[m2]].oradeg, xray[xkeep[m2]].odecdeg, ps=6, color='cyan'

    xbootes = im_empty_structure(xray[0],empty_value=-999.0,ncopies=ngal)
    xbootes[m1] = xray[xkeep[m2]]

    moretags = replicate({match: 0, object_position: 0L},ngal)
    xbootes = struct_addtags(moretags,temporary(xbootes))
    xbootes[m1].match = 1
    xbootes[m1].object_position = xkeep[m2]
    
; write out    
    outfile = catpath+'ages_xbootes.fits'
    im_mwrfits, xbootes, outfile, /clobber

return
end
