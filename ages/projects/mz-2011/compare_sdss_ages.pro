pro compare_sdss_ages
; jm09apr03nyu - compare the AGES and SDSS selection functions

    ss = read_mz_emline_sample(/mzhii_ancillary,/sdss)
    aa = read_mz_emline_sample(/mzhiiplus_ancillary)

    filt = ['sdss_r0.par','ndwfs_I.par']
    v2ab = k_vega2ab(filterlist=filt,/kurucz,/silent)
    k_load_vmatrix, vm, lam, vname='default.nolines'
    k_projection_table, rm, vm, lam, zvals, filt, zmax=1.0

    k_reconstruct_maggies, ss.coeffs, ss.z, sdssmag, rmatrix=rm, zvals=zvals
    k_reconstruct_maggies, aa.coeffs, aa.z, agesmag, rmatrix=rm, zvals=zvals

    djs_plot, -2.5*alog10(sdssmag[1,*]), -2.5*alog10(sdssmag[0,*]/sdssmag[1,*]), $
      psym=3, xrange=[13,22.0], yrange=[-0.2,1.2], xsty=1, ysty=1, $
      xtitle='I', ytitle='r - I'
    djs_oplot, -2.5*alog10(agesmag[1,*]), -2.5*alog10(agesmag[0,*]/agesmag[1,*]), $
      psym=6, symsize=0.12, color='red'
    

    
return
end
    
