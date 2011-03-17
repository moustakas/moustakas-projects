pro tkrs_merge_goods, phot
; jm08apr24nyu - merge the TKRS uber-redshift catalog with the GOODS-N
;                BVIz photometric catalog published by Giavalisco et
;                al. (2004), which was retrieved from Vizier; although
;                the TKRS people provide this photometric catalog,
;                they don't give uncertainties, and I want to use
;                Blanton's GOODS_TO_MAGGIES

    path = tkrs_path()

    filterlist = 'goods_'+['acs_'+['f435w','f606w','f775w','f850lp'],$
      ['J','H','Ks']+'_isaac_etc']+'.par'

; read the TKRS+many other redshifts catalog (for details see
; http://tkserver.keck.hawaii.edu/tksurvey/data_products/goods_desc.html)  
    
    tkrs = mrdfits(path+'tkrs_by_ra.fits',1) ; just TKRS
    tkrs_all = mrdfits(path+'goods-n-all.fits',1) 

    spherematch, 15.0D*tkrs.ra, tkrs.dec, $
      tkrs_all.ra, tkrs_all.dec, 5.0/3600.0, t1, t2
    moretags = im_empty_structure(im_struct_trimtags(tkrs[0],$
      select=['ID','ID2','RMAG'],newtags=['ID','GOODS','RMAG']),$
      empty_value=-999,ncopies=n_elements(tkrs_all))
    moretags[t2].id = tkrs[t1].id
    moretags[t2].goods = tkrs[t1].id2
    moretags[t2].rmag = tkrs[t1].rmag
    tkrs_all = struct_addtags(moretags,tkrs_all)
    
    b = mrdfits(path+'goods-n.b.fits',1)
    v = mrdfits(path+'goods-n.v.fits',1)
    i = mrdfits(path+'goods-n.i.fits',1)
    z = mrdfits(path+'goods-n.z.fits',1)

    spherematch, b.raj2000, b.dej2000, $
      tkrs_all.ra, tkrs_all.dec, 3.0/3600.0, m1, m2

    phot = tkrs_all
;   phot = struct_trimtags(tkrs_all,except=['*MAG*'])
    
; build the k-correct-compatible goods catalog; the GOODS magnitudes
; are crap; adopt the TKRS values and a fixed magnitude error
    
    names = ['B','V','I','Z','J','H','K']
    phot = struct_addtags(phot,$
      struct_addtags(mrd_struct(names+'mag_magauto',$
      replicate('-99.0',n_elements(names)),$
      n_elements(tkrs_all)),mrd_struct(names+'magerr_magauto',$
      replicate('-99.0',n_elements(names)),$
      n_elements(tkrs_all))))

    phot[m2].bmag_magauto    = b[m1].magbest
    phot[m2].bmagerr_magauto = b[m1].e_magbest
    phot[m2].vmag_magauto    = v[m1].magbest
    phot[m2].vmagerr_magauto = v[m1].e_magbest
    phot[m2].imag_magauto    = i[m1].magbest
    phot[m2].imagerr_magauto = i[m1].e_magbest
    phot[m2].zmag_magauto    = z[m1].magbest
    phot[m2].zmagerr_magauto = z[m1].e_magbest
 
;   plot, tkrs_all[m2].bmag, b[m1].magbest, ps=4, xsty=3, $
;     ysty=3, xr=[14,40], yr=[14,40]

;   phot.bmag_magauto    = tkrs_all.bmag
;   phot.bmagerr_magauto = 0.02
;   phot.vmag_magauto    = tkrs_all.vmag
;   phot.vmagerr_magauto = 0.02
;   phot.imag_magauto    = tkrs_all.imag
;   phot.imagerr_magauto = 0.02
;   phot.zmag_magauto    = tkrs_all.zmag
;   phot.zmagerr_magauto = 0.02

; add some tags from the TKRS structure - OBSOLETE!

;   spherematch, 15.0D*tkrs.ra, tkrs.dec, $
;     tkrs_all.ra, tkrs_all.dec, 5.0/3600.0, m1, m2
;   phot = struct_addtags(phot,struct_trimtags($
;     im_empty_structure(tkrs[0],ncopies=n_elements(phot),$
;     empty_value=-999.0, empty_string='...'),$
;     except=['RA','DEC','Z']))
;   final = phot[m2]
;   struct_assign, tkrs[m1], final, /nozero
;   phot[m2] = final

; write out just the objects with good redshifts

    outfile = path+'tkrs_zcat.fits'
    splog, 'Writing '+outfile
    phot = phot[where((phot.z gt 0.0) and (finite(phot.z) eq 1B))]
    mwrfits, phot, outfile, /create
    spawn, 'gzip -f '+outfile, /sh
    
return
end
    
