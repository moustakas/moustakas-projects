pro deep2_desi_elg_targeting
; jm14mar01siena - optimize the DESI targeting using DEEP2
; spectroscopy of galaxies in the Groth Strip (Field 1)

    desipath = getenv('IM_PROJECTS_DIR')+'/desi/deep2/'
    isedfit_paramfile = desipath+'desi_deep2_paramfile.par'
    outpath = desipath+'targeting/'

; read the full DEEP2 photometric catalog and the window function     
    deep2phot = mrdfits(deep2_path(/cat)+'pcat_ext.fits.gz',1)
    deep2_to_maggies, deep2phot, mm, ii

stop    
    
    keep = where(where(strmid(strtrim(deep2phot.objno,2),0,1) eq 1) and $
      deep2phot.pgal ge 0.5 and deep2phot.r gt 18 and deep2phot.r lt 23.6,ngal)
    splog, ngal
    phot = deep2phot[keep]
    
    ww = where(phot.zquality lt -90,comp=gg)
    djs_plot, phot.ra_deep, phot.dec_deep, psym=3, ysty=3, xsty=3, $
      xtitle='RA', ytitle='Dec'
    djs_oplot, phot[gg].ra_deep, phot[gg].dec_deep, psym=3, color='orange'

    im_plothist, phot[ww].r, xr=[15,27], bin=0.1
    im_plothist, phot[gg].r, bin=0.1, /over, /fill
    
    
    


stop    
    
    keep = where(deep2phot.r gt 22 and deep2phot.r lt 24)

; read all the relevant files
    cat = read_deep2_zcat(phot=phot)
    ppxf = read_deep2(/ppxf)
    kcorr = read_deep2(/kcorr)
    
    ised = read_isedfit(isedfit_paramfile,isedfit_dir=desipath)
    ikcorr = mrdfits(desipath+'desi_deep2_fsps_v2.4_miles_chab_charlot_sfhgrid01.fits.gz',1)

    select = where(cat.z gt 0.7 and cat.z lt 1.6 and phot.r gt 22 and phot.r lt 24)
    

    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, position=pos, $
      xrange=[-1,2], yr=[-1,2]
    djs_oplot, deep2phot[keep].r-deep2phot[keep].z, $
      deep2phot[keep].g-deep2phot[keep].r, psym=3
    djs_oplot, phot[select].r-phot[select].z, $
      phot[select].g-phot[select].r, psym=3, color='red'
    
stop    
    
return
end
    
