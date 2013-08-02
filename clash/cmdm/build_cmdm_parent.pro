function read_zcat, file
    readcol, file, ra, dec, z, zerr, r, flag, quasar, $
      format='D,D,F,F,F,I,A', /silent, comment='%'
    ngal = n_elements(ra)
    zcat = replicate({zcat_ra: 0D, zcat_dec: 0D, zobj: 0.0, $
      zobjerr: 0.0, zcat_flag: 0},ngal)
    zcat.zcat_ra = ra
    zcat.zcat_dec = dec
    zcat.zobj = z
    zcat.zobjerr = zerr
;   zcat.r = r
    zcat.zcat_flag = flag
return, zcat
end

function read_members, file
; Here is Mauricio description of how these galaxies were picked:
; "First, I identify the red sequence.  Then I select "all" the
; cluster members that are inside of 2 sigma of the red sequence
; (using F814w and F625w) and are inside the cluster redshift +- 0.2
; (it depends, sometimes is 0.1, other 0.3). With these members I
; create a new red sequence, with a new sigma and I repeat the process
; again until to obtain a stable number of members. Then I calculate
; their luminosity function and I obtain the m star of the best fit
; (schechter). Then I select all the members that are brighter than m*
; +2 magnitudes and that are inside of a radius of 120" (the SL
; region). "
    readcol, file, id, ra, dec, a, b, th, mag, lum, $
      format='L,D,D,F,F,F,F,F', /silent, comment='%'
    ngal = n_elements(ra)
    mem = replicate({id: 0L, ra: 0D, dec: 0D, a: 0.0, b: 0.0, $
      th: 0.0, mag: 0.0, lum: 0.0},ngal)
    mem.id = id
    mem.ra = ra
    mem.dec = dec
    mem.a = a
    mem.b = b
    mem.th = th
    mem.mag = mag
    mem.lum = lum
return, mem
end

pro build_cmdm_parent, phot, clobber=clobber
; jm13jul16siena - build the parent sample of objects for the
; CLASH/CMDM project with Doron
    
    path = getenv('IM_PROJECTS_DIR')+'/clash/cmdm/'

; match the Subaru catalog and the spectroscopic redshift catalog,
; keeping just the cluster members 
    allphot = rsex(path+'MACSJ1206_Subaru.cat')
    
    zcat = read_zcat(path+'MACSJ1206_specz_Marc_compilation.cat')
    these = where(zcat.zobj ge 0.4 and zcat.zobj le 0.48,ngal)
    splog, '0.4<z<0.48 cuts: ', ngal, n_elements(zcat)
    zcat = zcat[these]

    spherematch, allphot.ra, allphot.dec, zcat.zcat_ra, $
      zcat.zcat_dec, 3D/3600.0, m1, m2
    phot1 = struct_addtags(allphot[m1],zcat[m2])
    splog, 'Cluster galaxies with photometry ', n_elements(phot1)

; add redshift info to the cluster members identified by
; Mauricio; if there's no spectroscopic redshift then assign
; the cluster redshift, z=0.44
;   hst = rsex(path+'macs1206_ACS_IR.cat')
    mem = read_members(path+'mem_209g_final.pelli')
    spherematch, mem.ra, mem.dec, zcat.zcat_ra, zcat.zcat_dec, 3D/3600.0, m1, m2
    splog, 'Members catalog and matching spectroscopy ', n_elements(mem), n_elements(m1)
    mem_zcat = im_empty_structure(zcat[0],ncopies=n_elements(mem))
    mem_zcat.zobj = 0.44 ; assume cluster redshift
    mem_zcat[m1] = zcat[m2]    

; add the photometry and then merge
    spherematch, allphot.ra, allphot.dec, mem.ra, mem.dec, 3D/3600.0, m1, m2
    splog, 'Members catalog and Subaru photometry ', n_elements(mem), n_elements(m1)
    phot2 = struct_addtags(allphot[m1],mem_zcat[m2])

    phot = [phot1,phot2]
;   splog, 'Total number of objects, including duplicates'
    
; get rid of duplicates
    ing = spheregroup(phot.ra,phot.dec,1D/3600.0)
    phot = phot[uniq(ing,sort(ing))]
    splog, 'Final number of objects ', n_elements(phot)
    
    djs_plot, allphot.ra, allphot.dec, psym=3, xsty=3, ysty=3, $
      xr=[181.45,181.65], yr=[-8.9,-8.7]
    djs_oplot, phot.ra, phot.dec, psym=6, color='blue'
    djs_oplot, zcat.zcat_ra, zcat.zcat_dec, psym=8, color='orange'

; write out    
    im_mwrfits, phot, path+'cmdm_cat.fits', clobber=clobber

return
end
    
