pro ages_sfrm, ised1
; jm09mar09nyu - build the AGES mass functions

; read the data    
    aa1 = read_ages(/ancillary)
    datapath = ages_path(/isedfit)
    isedfile = datapath+'BwRIzK_salp_sfhgrid02_isedfit.fits.gz'
    if (n_elements(ised1) eq 0L) then ised1 = mrdfits(isedfile,1,/silent)

; sample definitions and Vmax parameters
    ifilter = 'ndwfs_I.par'
    vega2ab = k_vega2ab(filterlist=ifilter,/kurucz,/silent)
    ibright = 15.0 
    ifaint = 19.95 
    area = 7.9*!dtor^2 ; [sr]
    sample_zmin = 0.01
    sample_zmax = 0.8
    band_shift = 0.0
    h100 = 0.7
    vname = 'default.nolines'

    imag1 = -2.5*alog10(aa1.abmaggies[2]) ; BwRIzK
    main = where(aa1.main_flag and (aa1.select_absmag ne -999.0) and (imag1 ge ibright) and $
      (imag1 le ifaint) and (aa1.z ge sample_zmin) and (aa1.z le sample_zmax),nmain)
    aa = aa1[main]
    ised = ised1[main]
    imag = imag1[main]

; compute Vmax [h^{-3} Mpc^3]; remember - SELECT_ABSMAG is for h=1 
    lf_calc_vmax, imag, aa.select_absmag, aa.coeffs, ifilter, area, $
      ibright, ifaint, sample_zmin, sample_zmax, vmax=vmax, zmin=zmin, $
      zmax=zmax, band_shift=band_shift, actual_z=aa.z, zvals=zvals, $
      rmatrix=rmatrix, vname=vname

    vvmax = replicate({vmax: 0.0, zmin: 0.0, zmax: 0.0},nmain)
    vvmax.vmax = vmax/h100^3.0 ; h=1 --> h=0.7
    vvmax.zmax = zmax
    vvmax.zmin = zmin
    aa = struct_addtags(aa,vvmax)

; --------------------------------------------------    

    massaxis1 = findgen((12.5-8.0)/0.01+1)*0.01+8.0
    binsize = 0.1
    
    xpage = 8.5 & ypage = 8.8
    pagemaker, nx=2, ny=2, xspace=0, yspace=0.0, width=3.3*[1,1], height=3.3*[1,1], $
      xmargin=[1.5,0.4], ymargin=[1.1,1.1], xpage=xpage, ypage=ypage, $
      position=pos, /normal

    zbin1 = where((aa.z ge 0.01) and (aa.z lt 0.2) and (ised.mass50 gt 9.5))
    mass1 = ised[zbin1].mass50
    weight1 = aa[zbin1].spec_weight*(1.0/aa[zbin1].vmax)
    yy1 = im_hist1d(mass1,weight1,binsize=binsize,obin=xx1,$
      binedge=0,h_err=yy1_err)
    yy1 = yy1/binsize & yy1_err = yy1_err/binsize
    mf_fit_schechter, 10^xx1, yy1, yy1_err, fit1
    
    djs_plot, xx1, yy1, ps=10, /ylog
    djs_oplot, massaxis1, mf_schechter(10^massaxis1,fit1), color='red', thick=3
    

stop    

return
end
    
