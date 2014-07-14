function get_dndm, rmag, weight=weight, faintcut=faintcut, $
  brightcut=brightcut, magaxis=magaxis
; get the number counts    
    area = 0.4342 ; deg^2
    binsize = 0.2
    nbins = ceil((faintcut-brightcut)/binsize)
    magaxis = lindgen(nbins)*binsize+brightcut+binsize/2.0
    dndm = hogg_histogram(rmag,[brightcut,faintcut],$ ; #/0.2 mag/deg^2
      nbins,weight=weight)/area
;   dndm = im_hist1d(rmag,weight,min=brightcut,$ 
;     max=faintcut,binsize=binsize)/area/binsize
return, dndm
end

pro desi_cdr_plots, write_dndz=write_dndz
; jm14jun16siena - build plots for the DESI/CDR section on ELGs

    catpath = deep2_path(/cat)
    cdrpath = getenv('IM_SVNREPOS')+'/desi/cdr/4_Targets/plots/'
    targpath = getenv('IM_PROJECTS_DIR')+'/desi/targeting/'
    nzpath = getenv('DESIMODEL')+'/data/targets/'
    mbrownpath = getenv('IM_DATA_DIR')+'/mbrownatlas/'

    brightcut = 18.5
    faintcut = 24.0

;   magcut2 = 23.0
    oiisnrcut = 1.0
    oiicut1 = 8D-17 ; [erg/s/cm2]
    area = 0.4342 ; deg^2

    phot = mrdfits(targpath+'deep2egs-phot.fits.gz',1)
    zcat = mrdfits(targpath+'deep2egs-oii.fits.gz',1)
    stars = mrdfits(targpath+'deep2egs-stars.fits.gz',1)

; --------------------------------------------------
; check how our grz color box compares with M. Brown's
; templates at z=0.6-2 
    gal = [$
      'NGC_0337',$
      'NGC_0695',$
      'NGC_3079',$
      'Mrk_33',$
      'UGCA_219',$
      'NGC_3521',$ ; red
      'NGC_3690',$
      'NGC_4125',$ ; red
      'NGC_4138',$
      'NGC_4552',$ ; red
      'NGC_4725',$ ; red
      'NGC_5256',$
      'CGCG_049-057',$ ; red
      'NGC_5953',$
      'IC_4553',$ ; red
      'NGC_6090',$
      'NGC_6240',$
      'II_Zw_096']+'_spec.dat'
    nuvmu = [0.96,1.17,1.31,0.59,0.01,2.02,0.72,3.44,$
      1.83,3.01,2.52,1.35,3.17,1.89,2.35,0.82,1.69,0.49]
    keep = where(nuvmu lt 1.5,ngal)
    gal = gal[keep]

;   gal = ['II_Zw_096','NGC_6090','NGC_3690','UGCA_219',$
;     'Mrk_33','NGC_0337']+'_spec.dat'
;   ngal = n_elements(gal)
;   zvals = range(0.1,2,10)
    zvals = range(0.1,2,20)
    this = where(zvals eq 0.6)
    loz = where(zvals lt 0.6,nloz)
    hiz = where(zvals ge 0.6 and zvals le 1.6,nhiz)
    vhiz = where(zvals gt 1.6,nvhiz)
;   hiz = where(zvals ge 0.6,nhiz,comp=loz,ncomp=nloz)

; ##########
; stars    
    psfile = cdrpath+'deep2-grz-models-stars.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
      xmargin=[1.5,0.4], charsize=1.7
    
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.2], yrange=[-0.3,2]

    stars1 = stars[where(stars.cfhtls_r lt 23,nstar)]
    djs_oplot, stars1.cfhtls_r-stars1.cfhtls_z, symsize=0.1, $
      stars1.cfhtls_g-stars1.cfhtls_r, psym=6, $
      color=cgcolor('grey')
    
    cgloadct, 27, ncolors=ngal+1, bottom=0, /brewer
    for ii = 0, ngal-1 do begin
       if file_test(mbrownpath+gal[ii]) eq 0 then message, 'Help!'
       readcol, mbrownpath+gal[ii], wave, flux, $
         format='F,F,X,X', comment='#', /silent
       npix = n_elements(wave)
       k_projection_table, rmatrix, reform(flux,npix,1), wave, $
         zvals, cfhtls_filterlist()
       gr = -2.5*alog10(rmatrix[*,0,1]/rmatrix[*,0,2])
       rz = -2.5*alog10(rmatrix[*,0,2]/rmatrix[*,0,4])
;      cgoplot, rz, gr, psym=-symcat(16), color=ii+1
       cgoplot, rz[loz], gr[loz], psym=-symcat(16), color=ii+1
       cgoplot, rz[hiz], gr[hiz], psym=-symcat(15), color=ii+1
       cgoplot, rz[vhiz], gr[vhiz], psym=-symcat(17), color=ii+1
       cgoplot, [rz[loz[nloz-1]],rz[hiz[0]]], [gr[loz[nloz-1]],gr[hiz[0]]], $
         line=0, color=ii+1
       cgoplot, [rz[hiz[nhiz-1]],rz[vhiz[0]]], [gr[hiz[nhiz-1]],gr[vhiz[0]]], $
         line=0, color=ii+1
       plots, rz[0], gr[0], psym=symcat(9), symsize=2, color=ii+1
;      plots, rz[this], gr[this], psym=symcat(6), symsize=2, color=ii+1
    endfor
    cgloadct, 0

; overplot the proposed color box
    rzaxis1 = range(0.2,0.9,100)
    rzaxis2 = range(0.9,1.4,100)
    int1 = -0.1 & slope1 = 1.0
    int2 = 1.7 & slope2 = -1.0
    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[int1,slope1])], thick=8
    djs_oplot, rzaxis1, poly(rzaxis1,[int1,slope1]), thick=8
    djs_oplot, rzaxis2, poly(rzaxis2,[int2,slope2]), thick=8
    djs_oplot, 1.4*[1,1], [!y.crange[0],poly(1.4,[int2,slope2])], thick=8
    djs_oplot, [0.0,0.0], [!y.crange[0],0.1], line=2, thick=8
    djs_oplot, [0.0,0.2], [0.1,0.1], line=2, thick=8

    im_legend, ['z=0.1','z=0.1-0.5','z=0.6-1.6','z=1.7-2'], /left, /top, box=0, $
      psym=[16,16,15,17], margin=0, spacing=2.5, symsize=[1.0,1.0,1.0,1.3]*1.6
    plots, 0.228, 0.8723, psym=symcat(9,thick=6), symsize=3, /norm;, $
;     color=cgcolor('grey'), 
    im_legend, 'Stars (r<23)', /right, /top, box=0, psym=6, color=cgcolor('grey')
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

; ##########
; z<22 galaxies
    psfile = cdrpath+'deep2-grz-models-z22gals.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
      xmargin=[1.5,0.4], charsize=1.7
    
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.2], yrange=[-0.3,2]

    z22gals = phot[where(phot.cfhtls_z gt 20 and phot.cfhtls_z lt 22,nz22gals)]
    djs_oplot, z22gals.cfhtls_r-z22gals.cfhtls_z, symsize=0.1, $
      z22gals.cfhtls_g-z22gals.cfhtls_r, psym=6, $
      color=cgcolor('grey')
    
    cgloadct, 27, ncolors=ngal+1, bottom=0, /brewer
    for ii = 0, ngal-1 do begin
       if file_test(mbrownpath+gal[ii]) eq 0 then message, 'Help!'
       readcol, mbrownpath+gal[ii], wave, flux, $
         format='F,F,X,X', comment='#', /silent
       npix = n_elements(wave)
       k_projection_table, rmatrix, reform(flux,npix,1), wave, $
         zvals, cfhtls_filterlist()
       gr = -2.5*alog10(rmatrix[*,0,1]/rmatrix[*,0,2])
       rz = -2.5*alog10(rmatrix[*,0,2]/rmatrix[*,0,4])
;      cgoplot, rz, gr, psym=-symcat(16), color=ii+1
       cgoplot, rz[loz], gr[loz], psym=-symcat(16), color=ii+1
       cgoplot, rz[hiz], gr[hiz], psym=-symcat(15), color=ii+1
       cgoplot, rz[vhiz], gr[vhiz], psym=-symcat(17), color=ii+1
       cgoplot, [rz[loz[nloz-1]],rz[hiz[0]]], [gr[loz[nloz-1]],gr[hiz[0]]], $
         line=0, color=ii+1
       cgoplot, [rz[hiz[nhiz-1]],rz[vhiz[0]]], [gr[hiz[nhiz-1]],gr[vhiz[0]]], $
         line=0, color=ii+1
       plots, rz[0], gr[0], psym=symcat(9), symsize=2, color=ii+1
;      plots, rz[this], gr[this], psym=symcat(6), symsize=2, color=ii+1
    endfor
    cgloadct, 0

; overplot the proposed color box
    rzaxis1 = range(0.2,0.9,100)
    rzaxis2 = range(0.9,1.4,100)
    int1 = -0.1 & slope1 = 1.0
    int2 = 1.7 & slope2 = -1.0
    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[int1,slope1])], thick=8
    djs_oplot, rzaxis1, poly(rzaxis1,[int1,slope1]), thick=8
    djs_oplot, rzaxis2, poly(rzaxis2,[int2,slope2]), thick=8
    djs_oplot, 1.4*[1,1], [!y.crange[0],poly(1.4,[int2,slope2])], thick=8
    djs_oplot, [0.0,0.0], [!y.crange[0],0.1], line=2, thick=8
    djs_oplot, [0.0,0.2], [0.1,0.1], line=2, thick=8

    im_legend, ['z=0.1','z=0.1-0.5','z=0.6-1.6','z=1.7-2'], /left, /top, box=0, $
      psym=[16,16,15,17], margin=0, spacing=2.5, symsize=[1.0,1.0,1.0,1.3]*1.6
    plots, 0.228, 0.8723, psym=symcat(9,thick=6), symsize=3, /norm
    im_legend, 'Galaxies (20<z<22)', /right, /top, box=0, psym=15, color=cgcolor('grey')
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

; ##########
; r<23 galaxies
    psfile = cdrpath+'deep2-grz-models-r23gals.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
      xmargin=[1.5,0.4], charsize=1.7
    
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.2], yrange=[-0.3,2]

    r23gals = phot[where(phot.cfhtls_r gt 21 and phot.cfhtls_r lt 23,nr23gals)]
    djs_oplot, r23gals.cfhtls_r-r23gals.cfhtls_z, symsize=0.1, $
      r23gals.cfhtls_g-r23gals.cfhtls_r, psym=6, $
      color=cgcolor('grey')
    
    cgloadct, 27, ncolors=ngal+1, bottom=0, /brewer
    for ii = 0, ngal-1 do begin
       if file_test(mbrownpath+gal[ii]) eq 0 then message, 'Help!'
       readcol, mbrownpath+gal[ii], wave, flux, $
         format='F,F,X,X', comment='#', /silent
       npix = n_elements(wave)
       k_projection_table, rmatrix, reform(flux,npix,1), wave, $
         zvals, cfhtls_filterlist()
       gr = -2.5*alog10(rmatrix[*,0,1]/rmatrix[*,0,2])
       rz = -2.5*alog10(rmatrix[*,0,2]/rmatrix[*,0,4])
;      cgoplot, rz, gr, psym=-symcat(16), color=ii+1
       cgoplot, rz[loz], gr[loz], psym=-symcat(16), color=ii+1
       cgoplot, rz[hiz], gr[hiz], psym=-symcat(15), color=ii+1
       cgoplot, rz[vhiz], gr[vhiz], psym=-symcat(17), color=ii+1
       cgoplot, [rz[loz[nloz-1]],rz[hiz[0]]], [gr[loz[nloz-1]],gr[hiz[0]]], $
         line=0, color=ii+1
       cgoplot, [rz[hiz[nhiz-1]],rz[vhiz[0]]], [gr[hiz[nhiz-1]],gr[vhiz[0]]], $
         line=0, color=ii+1
       plots, rz[0], gr[0], psym=symcat(9), symsize=2, color=ii+1
;      plots, rz[this], gr[this], psym=symcat(6), symsize=2, color=ii+1
    endfor
    cgloadct, 0

; overplot the proposed color box
    rzaxis1 = range(0.2,0.9,100)
    rzaxis2 = range(0.9,1.4,100)
    int1 = -0.1 & slope1 = 1.0
    int2 = 1.7 & slope2 = -1.0
    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[int1,slope1])], thick=8
    djs_oplot, rzaxis1, poly(rzaxis1,[int1,slope1]), thick=8
    djs_oplot, rzaxis2, poly(rzaxis2,[int2,slope2]), thick=8
    djs_oplot, 1.4*[1,1], [!y.crange[0],poly(1.4,[int2,slope2])], thick=8
    djs_oplot, [0.0,0.0], [!y.crange[0],0.1], line=2, thick=8
    djs_oplot, [0.0,0.2], [0.1,0.1], line=2, thick=8

    im_legend, ['z=0.1','z=0.1-0.5','z=0.6-1.6','z=1.7-2'], /left, /top, box=0, $
      psym=[16,16,15,17], margin=0, spacing=2.5, symsize=[1.0,1.0,1.0,1.3]*1.6
    plots, 0.228, 0.8723, psym=symcat(9,thick=6), symsize=3, /norm
    im_legend, 'Galaxies (21<r<23)', /right, /top, box=0, psym=15, color=cgcolor('grey')
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

; --------------------------------------------------
; run elgtune
    infile = targpath+'deep2egs-oii.fits.gz'
    cfilter = ['cfhtls_g','cfhtls_r','cfhtls_z']
    mlimit = [99,99,99]
    o2tag = 'oii_3727'
    o2limit = 8D-17
    ndensity = 2400
    zpow = 3.0
    elgtune, infile=infile, cfilter=cfilter, mlimit=mlimit, $
      o2tag=o2tag, o2limit=o2limit, ndensity=ndensity, zpow=zpow

; --------------------------------------------------
; compare r- and z-band limited dN/dz distributions
    rlimit = 22.9
    zlimit1 = 22.0
    zlimit2 = 23.0
    zlimit3 = 23.5
    
    hiz_r = desi_get_hizelg(zcat,magcut=rlimit,sigma_kms=zcat.sigma_kms)
    hiz_z1 = desi_get_hizelg(zcat,magcut=zlimit1,sigma_kms=zcat.sigma_kms,/zbandlimit)
    hiz_z2 = desi_get_hizelg(zcat,magcut=zlimit2,sigma_kms=zcat.sigma_kms,/zbandlimit)
    hiz_z3 = desi_get_hizelg(zcat,magcut=zlimit3,sigma_kms=zcat.sigma_kms,/zbandlimit)

    oiibright_r = where(zcat[hiz_r].oii_3727_err ne -2.0 and $
      zcat[hiz_r].oii_3727 gt oiicut1)
    oiibright_z1 = where(zcat[hiz_z1].oii_3727_err ne -2.0 and $
      zcat[hiz_z1].oii_3727 gt oiicut1)
    oiibright_z2 = where(zcat[hiz_z2].oii_3727_err ne -2.0 and $
      zcat[hiz_z2].oii_3727 gt oiicut1)
    oiibright_z3 = where(zcat[hiz_z3].oii_3727_err ne -2.0 and $
      zcat[hiz_z3].oii_3727 gt oiicut1)

    zmin = 0.0
    zmax = 2.0
    zbin = 0.1
    nzbin = ceil((zmax-zmin)/zbin)
    zhist = lindgen(nzbin)*zbin+zmin+zbin/2.0
    nz_hiz_r = hogg_histogram(zcat[hiz_r].zbest,[zmin,zmax],$
      nzbin,weight=zcat[hiz_r].final_weight)/area
    nz_hiz_z1 = hogg_histogram(zcat[hiz_z1].zbest,[zmin,zmax],$
      nzbin,weight=zcat[hiz_z1].final_weight)/area
    nz_hiz_z2 = hogg_histogram(zcat[hiz_z2].zbest,[zmin,zmax],$
      nzbin,weight=zcat[hiz_z2].final_weight)/area
    nz_hiz_z3 = hogg_histogram(zcat[hiz_z3].zbest,[zmin,zmax],$
      nzbin,weight=zcat[hiz_z3].final_weight)/area

    nz_oiibright_r = hogg_histogram(zcat[hiz_r[oiibright_r]].zbest,[zmin,zmax],$
      nzbin,weight=zcat[hiz_r[oiibright_r]].final_weight)/area
    nz_oiibright_z1 = hogg_histogram(zcat[hiz_z1[oiibright_z1]].zbest,[zmin,zmax],$
      nzbin,weight=zcat[hiz_z1[oiibright_z1]].final_weight)/area
    nz_oiibright_z2 = hogg_histogram(zcat[hiz_z2[oiibright_z2]].zbest,[zmin,zmax],$
      nzbin,weight=zcat[hiz_z2[oiibright_z2]].final_weight)/area
    nz_oiibright_z3 = hogg_histogram(zcat[hiz_z3[oiibright_z3]].zbest,[zmin,zmax],$
      nzbin,weight=zcat[hiz_z3[oiibright_z3]].final_weight)/area

; normalize to unity
    nz_oiibright_r *= 1/total(nz_oiibright_r)
    nz_oiibright_z1 *= 1/total(nz_oiibright_z1)
    nz_oiibright_z2 *= 1/total(nz_oiibright_z2)
    nz_oiibright_z3 *= 1/total(nz_oiibright_z3)

; make the plot
    psfile = cdrpath+'deep2-elg-dndz-zlimit.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
      xmargin=[1.5,0.4], charsize=1.7
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='Redshift', ytitle='Fraction', $
      xrange=[zmin,zmax], yrange=[0.0,0.4]
;   djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
;     [0,nz_hiz,0], psym=10, thick=8, line=0, color='black'
    djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
      [0,nz_oiibright_r,0], psym=10, thick=8, line=0, color=cgcolor('dodger blue')
    djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
      [0,nz_oiibright_z1,0], psym=10, thick=8, line=5, color=cgcolor('firebrick')
    djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
      [0,nz_oiibright_z2,0], psym=10, thick=8, line=3, color=cgcolor('forest green')
    djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
      [0,nz_oiibright_z3,0], psym=10, thick=8, line=4, color=cgcolor('orange')

    im_legend, ['r<22.9','z<22.0','z<23.0','z<23.5'], spacing=2.0, /right, /top, box=0, $
      line=[0,5,3,4], color=['dodger blue','firebrick','forest green','orange'], $
      thick=10

    im_plotconfig, psfile=psfile, /psclose, /pdf
    
; --------------------------------------------------
; dN/dz of z>0.6 ELGs with bright [OII]; also writes out the DEEP2/EGS
; dN/dz (see /write_dndz keyword)  
    zbandlimit = 0
    if keyword_set(zbandlimit) then magcut1 = 23.0 else magcut1 = 22.0

    hiz = desi_get_hizelg(zcat,magcut=magcut1,sigma_kms=zcat.sigma_kms,zbandlimit=zbandlimit)
    nhiz = n_elements(hiz)
    oiibright = where(zcat[hiz].oii_3727_err ne -2.0 and zcat[hiz].oii_3727 gt oiicut1,noiibright)
    oiifaint = where(zcat[hiz].oii_3727_err ne -2.0 and zcat[hiz].oii_3727 lt oiicut1,noiifaint)

;    oiibright = where(zcat[hiz].oii_3727[1] ne -2.0 and $ ; oiisnr gt oiisnrcut and $
;      zcat[hiz].oii_3727_2_ew[0]/zcat[hiz].oii_3727_2_ew[1] gt 1.0 and $
;      zcat[hiz].oii_3727[0] gt oiicut1,noiibright)
;    oiifaint = where(zcat[hiz].oii_3727[1] ne -2.0 and $ ; oiisnr gt oiisnrcut and $
;      zcat[hiz].oii_3727_2_ew[0]/zcat[hiz].oii_3727_2_ew[1] gt 1.0 and $
;      zcat[hiz].oii_3727[0] lt oiicut1,noiifaint)
;    oiinone = where(zcat[hiz].oii_3727[1] eq -2 or $
;      zcat[hiz].oii_3727_2_ew[0]/zcat[hiz].oii_3727_2_ew[1] le 1.0,noiinone)
    
    zmin = 0.0
    zmax = 2.0
    zbin = 0.1
    nzbin = ceil((zmax-zmin)/zbin)
    zhist = lindgen(nzbin)*zbin+zmin+zbin/2.0
    nz = hogg_histogram(zcat.zbest,[zmin,zmax],nzbin,$
      weight=zcat.final_weight)/area
    nz_hiz = hogg_histogram(zcat[hiz].zbest,[zmin,zmax],$
      nzbin,weight=zcat[hiz].final_weight)/area
    nz_oiibright = hogg_histogram(zcat[hiz[oiibright]].zbest,[zmin,zmax],$
      nzbin,weight=zcat[hiz[oiibright]].final_weight)/area
    nz_oiifaint = hogg_histogram(zcat[hiz[oiifaint]].zbest,[zmin,zmax],$
      nzbin,weight=zcat[hiz[oiifaint]].final_weight)/area

; extrapolate linearly
    anchor = where(zhist gt 1.1 and zhist lt 1.4)
    extrap = where(zhist gt 1.4 and zhist lt 1.8)

;   plot, zhist, alog10(nz_hiz), psym=8, xr=[1.0,1.6], symsize=3
    nz_hiz_extrap = nz_hiz
;   nz_hiz_extrap[extrap] = interpol(nz_hiz[anchor],$
;     zhist[anchor],zhist[extrap])>nz_hiz[extrap]
    nz_hiz_extrap[extrap] = 10^poly(zhist[extrap],$
      linfit(zhist[anchor],alog10(nz_hiz[anchor])))
;   niceprint, zhist, nz_hiz_extrap, nz_hiz & print

;   plot, zhist, alog10(nz_oiibright), psym=8, xr=[1.0,1.6], symsize=3
    nz_oiibright_extrap = nz_oiibright
;   nz_oiibright_extrap[extrap] = interpol(nz_oiibright[anchor],$
;     zhist[anchor],zhist[extrap])>nz_oiibright[extrap]
    nz_oiibright_extrap[extrap] = 10^poly(zhist[extrap],$
      linfit(zhist[anchor],alog10(nz_oiibright[anchor])))
;   niceprint, zhist, nz_oiibright_extrap, nz_oiibright & print
    
    nz_oiifaint_extrap = nz_oiifaint
    nz_oiifaint_extrap[extrap] = interpol(nz_oiifaint[anchor],$
      zhist[anchor],zhist[extrap])>nz_oiifaint[extrap]
;   niceprint, zhist, nz_oiifaint_extrap, nz_oiifaint

; renormalize to the number of targets from the CFHTLS-DEEP survey
    if keyword_set(zbandlimit) then cfhtlsfile = 'zdndm-cfhtls.txt' else $
      cfhtlsfile = 'rdndm-cfhtls.txt'
    readcol, targpath+cfhtlsfile, cfhtlsmag, dndm_true, dndmerr_true, $
      dndm_degraded, dndmerr_degraded, fraction, /silent
;   niceprint, cfhtlsmag, dndm_true, dndm_degraded
    dndm_degraded_cen = interpol(dndm_degraded,cfhtlsmag,cfhtlsmag-0.05/2) ; center bin
    cfhtlscounts = interpol(dndm_degraded_cen,cfhtlsmag-0.05/2,magcut1)

; normalize to the desired target density
    cfhtlscounts = 2400.0 ; 3000.
  
    norm = cfhtlscounts/total(nz_hiz_extrap)
    nz_hiz_extrap *= norm
    nz_oiibright_extrap *= norm
    nz_oiifaint_extrap *= norm

;   splog, total(nz_hiz_extrap), total(nz_oiibright_extrap), $
;     total(nz_oiibright_extrap)/total(nz_hiz_extrap)

; get the fraction of targets in the z=0.85-1.35 redshift range that
; have a sufficiently high [OII] flux
    hizref = where(zcat[hiz].zbest gt 0.85 and zcat[hiz].zbest lt 1.35)
    hizoii = where(zcat[hiz[oiibright]].zbest gt 0.85 and zcat[hiz[oiibright]].zbest lt 1.35)
    numer = total(zcat[hiz[oiibright[hizoii]]].final_weight)/area
    denom = total(zcat[hiz[hizref]].final_weight)/area
    fracgood = numer/denom
    splog, 'Fraction of targets with [OII]>8x10^-17 ', numer, denom, fracgood
    
; optionally write out dN/dz
    if keyword_set(write_dndz) then begin
       outfile = nzpath+'nz_elg_deep2.dat'
       openw, lun, outfile, /get_lun
       printf, lun, '# Distribution of ELG targets from DEEP2/EGS'
       printf, lun, '# Selection by John Moustakas on 2014 July 11:'
       printf, lun, '#   r < 23.0'
       printf, lun, '#   0.2 < (r-z) < 1.4'
       printf, lun, '#   (g-r) < (r-z) - 0.1'
       printf, lun, '#   (g-r) < 1.7 - (r-z)'
       printf, lun, '#'
       printf, lun, '# The normalization is the one determined from'
       printf, lun, '# the DEEP2 dataset itself and should *not* be used'
       printf, lun, '# for projections (see nz_elg.dat).'
       printf, lun, '#'
       printf, lun, '# Total number per sq deg per dz=0.1 redshift bin'
       printf, lun, '#'
       printf, lun, '# z_lo z_hi  N_elg'
       for jj = 0, n_elements(zhist)-1 do printf, lun, zhist[jj]-zbin/2, $
         zhist[jj]+zbin/2, nz_oiibright_extrap[jj]
       free_lun, lun
    endif

; make the plot
    psfile = cdrpath+'deep2-elg-dndz.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
      xmargin=[1.5,0.4], charsize=1.7
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='Redshift', ytitle='dN/dz (gal / 0.1dz / deg^{2})', $
      xrange=[zmin,zmax], yrange=[0.0,1100]
;     xrange=[zmin,zmax], yrange=[0.0,6500]
;   djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
;     [0,nz,0], psym=10, thick=6, line=0
    djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
      [0,nz_hiz_extrap,0], psym=10, thick=8, line=0, color='black'
    djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
      [0,nz_oiibright_extrap,0], psym=10, thick=6, line=5, color=cgcolor('firebrick')
;   djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
;     [0,nz_oiibright_cor,0], psym=10, thick=8, line=5, color='blue'
    djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
      [0,nz_oiifaint_extrap,0], psym=10, thick=6, line=4, color=cgcolor('forest green')

    magstr = string(magcut1,format='(F4.1)')
    im_legend, ['18.5<r<'+magstr], spacing=2.0, /left, /top, box=0, $
      position=[0.17,0.92], /norm, charsize=1.7

    xyouts, 0.63, 0.88, 'In grz color box:', /norm, align=0.0, charsize=1.3
    im_legend, ['All galaxies','[OII]>8\times10^{-17}','[OII]<8\times10^{-17}'], $
      color=['black','firebrick','forest green'], /left, /top, box=0, thick=6, $
      charsize=1.3, spacing=2.0, margin=0, line=[0,5,4], pspacing=1.5, $
      position=[0.65,0.85], /norm

;   im_legend, ['zcatParent','grz Sample','grz & F([OII])>8\times10^{-17}'], $
;     line=[0,3,5], pspacing=1.7, color=['','orange','blue'], /right, /top, box=0, $
;     charsize=1.4, spacing=2.0, margin=0, thick=6
    im_plotconfig, psfile=psfile, /psclose, /pdf

; --------------------------------------------------
; write out the final dN/dz for ELGs based on the median of the DEEP2
; and COSMOS distributions, normalized to the same target density (see
; the targets.dat file)
    ntarget_elg = 2400.0          ; number of potential targets
    nobs_elg = 0.8*ntarget_elg    ; number that will get a fiber (fiber completeness)
    success_elg = 0.7             ; number with [OII]>8x10^-17 (target completeness)
    
    zbin = 0.1
    readcol, nzpath+'nz_elg_deep2.dat', zlo_deep2, zhi_deep2, dndz_deep2, /silent
    readcol, nzpath+'nz_elg_vvds.dat', zlo_vvds, zhi_vvds, dndz_vvds, /silent
    readcol, nzpath+'nz_elg_cosmos.dat', zlo_cosmos, zhi_cosmos, dndz_cosmos, /silent

; shift the redshift bin, crop, and renormalize
    keep = where(zhi_cosmos le 2.0)
    zlo_cosmos = zlo_cosmos[keep]
    zhi_cosmos = zhi_cosmos[keep]
    dndz_cosmos = dndz_cosmos[keep]

    dndz_deep2 *= 1.0/total(dndz_deep2)
    dndz_cosmos *= 1.0/total(dndz_cosmos)
    dndz_vvds *= 1.0/total(dndz_vvds)

; get the median dN/dz and write out the final n(z) after normalizing 
    zlo = zlo_deep2
    zhi = zhi_deep2
    dndz = djs_median([[dndz_cosmos],[dndz_deep2]],2)
    dndz *= success_elg*nobs_elg/total(dndz)

    if keyword_set(write_dndz) then begin
       outfile = nzpath+'nz_elg.dat'
       openw, lun, outfile, /get_lun
       printf, lun, '# Distribution of ELG targets based on the median of '
       printf, lun, '# the DEEP2 and COSMOS n(z) determinations.'
       printf, lun, '# File written by John Moustakas on 2014 July 11:'
       printf, lun, '#   r < 23.0'
       printf, lun, '#   0.2 < (r-z) < 1.4'
       printf, lun, '#   (g-r) < (r-z) - 0.1'
       printf, lun, '#   (g-r) < 1.7 - (r-z)'
       printf, lun, '#'
       printf, lun, '# Total number per sq deg per dz=0.1 redshift bin'
       printf, lun, '#'
       printf, lun, '# The target densities are normalized to the nobs_elg*success_elg'
       printf, lun, '# values specified in the targets.dat file'
       printf, lun, '#'
       printf, lun, '# z_lo z_hi  N_elg'
       for jj = 0, n_elements(zhist)-1 do printf, lun, zlo[jj], zhi[jj], dndz[jj]
       free_lun, lun
    endif
    
; make the plot of the final distribution for the CDR
    zbin_deep2 = (zhi_deep2-zlo_deep2)/2+zlo_deep2
    zbin_cosmos = (zhi_cosmos-zlo_cosmos)/2+zlo_cosmos
    zbin_vvds = (zhi_vvds-zlo_vvds)/2+zlo_vvds
    psfile = cdrpath+'elg-dndz.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
      xmargin=[1.5,0.4], charsize=1.7
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='Redshift', ytitle='dN/dz (gal / 0.1dz / deg^{2})', $
      xrange=[zmin,zmax], yrange=[0.0,800]
    djs_oplot, [(zbin_deep2[0]-zbin/2.0),zbin_deep2,(zbin_deep2[nzbin-1]+zbin/2.0)], $
      [0,dndz_deep2,0], psym=10, thick=6, line=3, color='dodger blue'
    djs_oplot, [(zbin_cosmos[0]-zbin/2.0),zbin_cosmos,(zbin_cosmos[nzbin-1]+zbin/2.0)], $
      [0,dndz_cosmos,0], psym=10, thick=6, line=5, color=cgcolor('firebrick')
;   djs_oplot, [(zbin_vvds[0]-zbin/2.0),zbin_vvds,(zbin_vvds[nzbin-1]+zbin/2.0)], $
;     [0,dndz_vvds,0], psym=10, thick=8, line=3, color=cgcolor('forest green')
    djs_oplot, [(zbin_deep2[0]-zbin/2.0),zbin_deep2,(zbin_deep2[nzbin-1]+zbin/2.0)], $
      [0,dndz,0], psym=10, thick=8, line=0, color='black'

    magstr = string(magcut1,format='(F4.1)')
    im_legend, ['18.5<r<'+magstr], spacing=2.0, /left, /top, box=0, $
      margin=0, /norm, charsize=1.7       ; position=[0.17,0.92], 

    im_legend, ['DEEP2/EGS','COSMOS','Median'], /right, /top, box=0, $
      line=[3,5,0], color=['dodger blue','firebrick','black'], $
      pspacing=1.9, thick=8
;   im_legend, ['DEEP2/EGS','COSMOS','VVDS'], /right, /top, box=0, $
;     line=[0,5,3], color=['black','orange','forest green'], $
;     pspacing=1.9, thick=8
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

    
stop
    
; --------------------------------------------------
; gr vs rz coded by [OII] strength
    zmin = 0.6 ; 0.8
    zmax = 1.6
    magcut1 = 23.0
    
    all = where(zcat.cfhtls_r lt magcut1,nall)
    loz = where(zcat.zbest lt zmin and zcat.cfhtls_r lt magcut1,nloz)
    oiibright = where(zcat.zbest ge zmin and $ ; zcat.zbest lt zmax and $
      zcat.cfhtls_r lt magcut1 and zcat.oii_3727_err ne -2.0 and $
      zcat.oii_3727 gt oiicut1,noiibright)
    oiifaint = where(zcat.zbest ge zmin and $ ; zcat.zbest lt zmax and $
      zcat.cfhtls_r lt magcut1 and zcat.oii_3727_err ne -2.0 and $
      zcat.oii_3727 lt oiicut1,noiifaint)
    oiinone = where(zcat.zbest ge zmin and $ ; zcat.zbest lt zmax and $ ; all at z>1.4
      zcat.cfhtls_r lt magcut1 and zcat.oii_3727_err eq -2,noiinone)

; for testing
    oiibright_vhiz = where(zcat.zbest ge 1.2 and $ ; zcat.zbest lt zmax and $
      zcat.cfhtls_r lt magcut1 and zcat.oii_3727_err ne -2.0 and $
      zcat.oii_3727[0] gt oiicut1)
    splog, nall, nloz, noiibright, noiifaint, noiinone, $
      nloz+noiibright+noiifaint+noiinone

    psfile = cdrpath+'deep2-elg-grz-oii.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
      xmargin=[1.5,0.4], charsize=1.7
    
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.2], yrange=[-0.3,2]
    djs_oplot, zcat[loz].cfhtls_r-zcat[loz].cfhtls_z, $
      zcat[loz].cfhtls_g-zcat[loz].cfhtls_r, psym=symcat(16), symsize=0.3
    djs_oplot, zcat[oiifaint].cfhtls_r-zcat[oiifaint].cfhtls_z, $
      zcat[oiifaint].cfhtls_g-zcat[oiifaint].cfhtls_r, psym=symcat(6), $
      color=cgcolor('forest green'), symsize=0.5
    djs_oplot, zcat[oiibright].cfhtls_r-zcat[oiibright].cfhtls_z, $
      zcat[oiibright].cfhtls_g-zcat[oiibright].cfhtls_r, psym=symcat(5), $
      color=cgcolor('firebrick'), symsize=0.3
    djs_oplot, zcat[oiibright_vhiz].cfhtls_r-zcat[oiibright_vhiz].cfhtls_z, $
      zcat[oiibright_vhiz].cfhtls_g-zcat[oiibright_vhiz].cfhtls_r, psym=symcat(7), $
      color=cgcolor('blue'), symsize=0.3
;   djs_oplot, zcat[oiinone].cfhtls_r-zcat[oiinone].cfhtls_z, $
;     zcat[oiinone].cfhtls_g-zcat[oiinone].cfhtls_r, psym=symcat(7), $
;     color=cgcolor('blue'), symsize=0.3

    hiz_all = desi_get_hizelg(zcat,magcut=magcut1)
    hiz_oiibright = desi_get_hizelg(zcat[oiibright]);,magcut=magcut1)
    hiz_oiifaint = desi_get_hizelg(zcat[oiifaint]);,magcut=magcut1)
;   hiz_oiinone = desi_get_hizelg(zcat[oiinone]);,magcut=magcut1)
    hiz_loz = desi_get_hizelg(zcat[loz]);,magcut=magcut1)
    splog, 'All ', n_elements(hiz_all), $
      ' Bright ', n_elements(hiz_oiibright), 1.0*n_elements(hiz_oiibright)/n_elements(hiz_all), $
      ' Faint ', n_elements(hiz_oiifaint), 1.0*n_elements(hiz_oiifaint)/n_elements(hiz_all), $
;     ' None ', 1.0*n_elements(hiz_oiinone)/n_elements(hiz_all), $
      ' Low-z ', n_elements(hiz_loz), 1.0*n_elements(hiz_loz)/n_elements(hiz_all)
;   help, hiz_all, hiz_oiibright, hiz_oiifaint, hiz_oiinone, loz
;   djs_oplot, zcat[hiz_all].cfhtls_r-zcat[hiz_all].cfhtls_z, $
;     zcat[hiz_all].cfhtls_g-zcat[hiz_all].cfhtls_r, psym=7, symsize=0.2, $
;     color='orange'
;    djs_oplot, zcat[loz[hiz_loz]].cfhtls_r-zcat[loz[hiz_loz]].cfhtls_z, $
;      zcat[loz[hiz_loz]].cfhtls_g-zcat[loz[hiz_loz]].cfhtls_r, psym=7, symsize=1.1, $
;      color='blue'
    
;; Mostek       
;    rzaxis = range(0.2,1.2,500)
;    djs_oplot, rzaxis, poly(rzaxis,[-0.08,0.68]), line=5, thick=6
;    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[-0.08,0.68])], line=5, thick=6
;    djs_oplot, 1.3*[1,1], [!y.crange[0],poly(1.3,[-0.08,0.68])], line=5, thick=6
    
;; proposed
;    rzaxis = range(0.1,1.1,500)
;    int = 0.0 & slope = 0.9
;    djs_oplot, [!x.crange[0],0.3], 0.1*[1,1], line=0, thick=6
;    djs_oplot, [0.3,1.3], poly([0.3,1.3],[int,slope]), line=0, thick=6
;    djs_oplot, 1.3*[1,1], [!y.crange[0],poly(1.3,[int,slope])], line=0, thick=6

; proposed
    rzaxis1 = range(0.2,0.9,100)
;   rzaxis1 = range(0.2,1.4,100)
    rzaxis2 = range(0.9,1.4,100)
    int1 = -0.1 & slope1 = 1.0
    int2 = 1.7 & slope2 = -1.0
    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[int1,slope1])], thick=6
;   djs_oplot, [!x.crange[0],0.2], poly(0.2,[int1,slope1])*[1,1], line=0, thick=6
;   djs_oplot, [-0.2,0.3], poly(0.3,[int1,slope1])*[1,1], line=0, thick=6
;   djs_oplot, 0.1*[1,1], [!y.crange[0],poly(0.1,[int1,slope1])], line=0, thick=6
    djs_oplot, rzaxis1, poly(rzaxis1,[int1,slope1]), line=0, thick=6
;   djs_oplot, rzaxis1, poly(rzaxis1,[int1,slope1])<poly(0.9,[int1,slope1]), line=0, thick=6
    djs_oplot, rzaxis2, poly(rzaxis2,[int2,slope2]), line=0, thick=6
;   djs_oplot, [0.9,!x.crange[1]], poly(0.9,[int1,slope1])*[1,1], line=0, thick=6
;   djs_oplot, [0.9,1.4], poly(0.9,[int1,slope1])*[1,1], line=0, thick=6
    djs_oplot, 1.4*[1,1], [!y.crange[0],poly(1.4,[int2,slope2])], line=0, thick=6

;   djs_oplot, zcat[oiibright[hiz]].cfhtls_r-zcat[oiibright[hiz]].cfhtls_z, $
;     zcat[oiibright[hiz]].cfhtls_g-zcat[oiibright[hiz]].cfhtls_r, psym=symcat(16), $
;     symsize=2, color='green'
    
;; Johan's proposed cuts    
;    rzaxis = range(0.2,1.25,500)
;    int = 0.1 & slope = 0.55/1.25
;    djs_oplot, rzaxis, slope*rzaxis+int, line=5, thick=6
;    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[int,slope])], line=5, thick=6
;    djs_oplot, 1.25*[1,1], [!y.crange[0],poly(1.25,[int,slope])], line=5, thick=6
    
;   rmag-zmag>0.2
;   rmag-zmag<1.25
;   gmag-rmag<0.55*(rmag-zmag)/1.25+0.1
;   gmag-rmag>0

    magstr = string(magcut1,format='(G0)')
    im_legend, ['18.5<r<'+magstr], spacing=2.0, /left, /top, box=0, $
      position=[0.17,0.92], /norm, charsize=1.7

;   im_legend, ['0.75<z<1.45'], spacing=2.0, /left, /top, box=0, $
;     position=[0.18,0.855], /norm, charsize=1.4
;   im_legend, ['[OII]>8\times10^{-17}',$
;     '[OII]<8\times10^{-17}','[OII] unmeasured'], /left, /top, $
;     box=0, psym=[5,6,7], position=[0.21,0.8], /norm, $
;     color=['firebrick','forest green','blue'], charsize=1.2
;   im_oplot_box, 0.9, 0.5, 0.0, xoff=0.0, yoff=1.48

    zstr = string(zmin,format='(F3.1)')
    im_legend, ['z<'+zstr+'','z>'+zstr+'; [OII]>8\times10^{-17}',$
      'z>'+zstr+'; [OII]<8\times10^{-17}'], /left, /top, $
      box=0, psym=[16,5,6], position=[0.2,0.86], /norm, $
      color=['','firebrick','forest green'], charsize=1.5
;   im_legend, ['z<0.8','z>0.8; [OII]>8\times10^{-17}',$
;     'z>0.8; [OII]<8\times10^{-17}',$
;     'z>0.8; [OII] unmeasured'], $
;     /left, /top, $
;     box=0, psym=[16,5,6,7], position=[0.2,0.86], /norm, $
;     color=['','firebrick','forest green','blue'], charsize=1.3
;   im_legend, ['z<0.75','0.75<z<1.45; [OII]>8\times10^{-17}',$
;     '0.75<z<1.45; [OII]<8\times10^{-17}',$
;     '0.75<z<1.45; [OII] unmeasured'], /left, /top, $
;     box=0, psym=[16,5,6,7], position=[0.2,0.86], /norm, $
;     color=['','firebrick','forest green','blue'], charsize=1.4
    
    xyouts, 1.4, 0.35, 'In color box:', align=0.0, charsize=1.5, /data
    im_legend, [$
      'N='+strtrim(n_elements(hiz_loz),2)+' ('+$
      string(round(100.0*n_elements(hiz_loz)/n_elements(hiz_all)),format='(I0)')+'%)',$
      'N='+strtrim(n_elements(hiz_oiibright),2)+' ('+$
      string(round(100.0*n_elements(hiz_oiibright)/n_elements(hiz_all)),format='(I0)')+'%)',$
      'N='+strtrim(n_elements(hiz_oiifaint),2)+' ('+$
      string(round(100.0*n_elements(hiz_oiifaint)/n_elements(hiz_all)),format='(I0)')+'%)'], $
      /left, /bottom, box=0, psym=[16,5,6], position=[0.72,0.25], /norm, $
      color=['','firebrick','forest green'], charsize=1.4

    im_plotconfig, psfile=psfile, /psclose, /pdf

; --------------------------------------------------
; cumulative number counts

; assume that the objects with formal flux limits above our [OII] cut
; are detections (this is a small number of objects)    
    hiz = desi_get_hizelg(zcat.ugriz)
    oiibright = where(zcat[hiz].oii_3727[1] ne -2.0 and $ ; oiisnr gt oiisnrcut and $
      zcat[hiz].oii_3727_2_ew[0]/zcat[hiz].oii_3727_2_ew[1] gt 1.0 and $
      zcat[hiz].oii_3727[0] gt oiicut1,noiibright)
    oiifaint = where(zcat[hiz].oii_3727[1] ne -2.0 and $ ; oiisnr gt oiisnrcut and $
      zcat[hiz].oii_3727_2_ew[0]/zcat[hiz].oii_3727_2_ew[1] gt 1.0 and $
      zcat[hiz].oii_3727[0] lt oiicut1,noiifaint)
    oiinone = where(zcat[hiz].oii_3727[1] eq -2 or $
      zcat[hiz].oii_3727_2_ew[0]/zcat[hiz].oii_3727_2_ew[1] le 1.0,noiinone)

    dndm_phot = get_dndm(phot.cfhtls_r,faintcut=faintcut,$
      brightcut=brightcut,magaxis=magaxis)
    dndm = get_dndm(zcat.cfhtls_r,weight=zcat.final_weight,$
      faintcut=faintcut,brightcut=brightcut,magaxis=magaxis)
    dndm_hiz = get_dndm(zcat[hiz].cfhtls_r,weight=zcat[hiz].final_weight,$
      faintcut=faintcut,brightcut=brightcut)
    dndm_oiibright = get_dndm(zcat[hiz[oiibright]].cfhtls_r,$
      weight=zcat[hiz[oiibright]].final_weight,$
      faintcut=faintcut,brightcut=brightcut)
    
; get the *ratio* of the number of bright-to-faint [OII] sources so
; that we can correct for the objects with missing [OII] measurements
; (see plot below)     
    dndm_oiifaint = get_dndm(zcat[hiz[oiifaint]].cfhtls_r,$
      weight=zcat[hiz[oiifaint]].final_weight,$
      faintcut=faintcut,brightcut=brightcut)
    dndm_oiinone = get_dndm(zcat[hiz[oiinone]].cfhtls_r,$
      weight=zcat[hiz[oiinone]].final_weight,$
      faintcut=faintcut,brightcut=brightcut)

;   denom = dndm_oiibright+dndm_oiifaint
;   dndm_oiicor = dndm_oiinone*dndm_oiibright/(denom+(denom eq 0))*(denom ne 0)
;   dndm_oiibright_cor = dndm_oiibright + dndm_oiicor    
    
    psfile = cdrpath+'deep2-elg-dndm.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
      xmargin=[1.5,0.4], charsize=1.7
    
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='r (AB mag)', ytitle='log N(<r) (gal / deg^{2})', $
      xrange=[brightcut-0.5,faintcut+0.5], yrange=alog10([5,2E5]) ;, /ylog
    djs_oplot, magaxis, alog10(total(dndm_phot,/cumu)), line=0, thick=6, color='grey'
;     line=2, thick=6, color=cgcolor('dark grey')
;   djs_oplot, magaxis, alog10(total(dndm,/cumu)), line=0, thick=6
    
    ww = where(total(dndm_hiz,/cumu) gt 0)
    djs_oplot, magaxis[ww], alog10((total(dndm_hiz,/cumu))[ww]), line=0, thick=8
    
    ww = where(total(dndm_oiibright,/cumu) gt 0)
    djs_oplot, magaxis[ww], alog10((total(dndm_oiibright,/cumu))[ww]), color=cgcolor('firebrick'), $
      line=5, thick=8           ; psym=symcat(6,thick=4), symsize=1.3

    ww = where(total(dndm_oiifaint,/cumu) gt 0)
    djs_oplot, magaxis[ww], alog10((total(dndm_oiifaint,/cumu))[ww]), color=cgcolor('forest green'), $
      line=4, thick=8           ; psym=symcat(6,thick=4), symsize=1.3

    readcol, '/Users/ioannis/research/projects/desi/targeting/counts.txt', skip=1, $
      mag, num, numerr, /silent
    oploterror, mag, alog10(num), numerr/num/alog(10), psym=8, color='blue', symsize=0.8
    
;   ww = where(total(dndm_oiinone,/cumu) gt 0)
;   djs_oplot, magaxis[ww], alog10((total(dndm_oiinone,/cumu))[ww]), color=cgcolor('dark grey'), $
;     line=1, thick=6           ; psym=symcat(6,thick=4), symsize=1.3

;    ww = where(total(dndm_oiibright_cor,/cumu) gt 0)
;    djs_oplot, magaxis[ww], alog10((total(dndm_oiibright_cor,/cumu))[ww]), $
;      color=cgcolor('red'), line=5, thick=8 ; psym=symcat(6,thick=4), symsize=1.3
    
    numcut = 3000
    rcut = interpol(magaxis,total(dndm_hiz,/cumu),numcut)
;   numcut = 2400
;   rcut = interpol(magaxis,total(dndm_oiibright_cor,/cumu),numcut)
;   rcut = 22.7
    oiinumcut = long(interpol(total(dndm_oiibright,/cumu),magaxis,rcut))
    splog, rcut, numcut, oiinumcut, 1.0*oiinumcut/numcut
;   oiinumcut = long(interpol(total(dndm_oiibright_cor,/cumu),magaxis,rcut))
    djs_oplot, [!x.crange[0],rcut], alog10(numcut)*[1,1], line=1, thick=6
    djs_oplot, [!x.crange[0],rcut], alog10(oiinumcut)*[1,1], line=1, thick=6
    djs_oplot, rcut*[1,1], [!y.crange[0],alog10(numcut)], line=1, thick=6
    
    xyouts, 18.3, alog10(numcut*1.2), textoidl(string(numcut,$
      format='(I0)')+' gal/deg^{2}'), align=0.0, /data, charsize=1.5
    xyouts, 18.3, alog10(oiinumcut*0.6), textoidl(string(oiinumcut,$
      format='(I0)')+' gal/deg^{2}'), align=0.0, /data, charsize=1.5
    xyouts, rcut+0.1, 2.3, textoidl('r='+string(rcut,$
      format='(F4.1)')), align=0.0, orientation=270, /data, charsize=1.5
    
    im_legend, ['All Galaxies'], /right, /top, box=0, thick=6, $
      charsize=1.3, spacing=2.0, margin=0, line=0, pspacing=1.5, color='grey'
    xyouts, 0.21, 0.88, 'In grz color box:', /norm, align=0.0, charsize=1.3
    im_legend, ['All galaxies','[OII]>8\times10^{-17}','[OII]<8\times10^{-17}'], $
      color=['black','firebrick','forest green'], /left, /top, box=0, thick=6, $
      charsize=1.3, spacing=2.0, margin=0, line=[0,5,4], pspacing=1.5, $
      position=[0.22,0.85], /norm
;   im_legend, ['PhotParent','zcatParent','grz Color Box',$
;     'grz & F([OII])>8\times10^{-17}'], $
;     color=['dark grey','','orange','blue'], /left, /top, box=0, thick=6, $
;     charsize=1.3, spacing=2.0, margin=0, line=[2,0,3,5], pspacing=1.5
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

; --------------------------------------------------
; gr vs rz - stellar contamination
    allphot = mrdfits(catpath+'deep2.pcat_ext.fits.gz',1)
    keep = where(allphot.zquality ge 3 and allphot.g gt 0 and $
      allphot.r gt 0 and allphot.z gt 0 and allphot.gerr lt 1 and $
      allphot.rerr lt 1 and allphot.zerr lt 1 and $
      allphot.badflag eq 0 and allphot.pgal ge 1)
    allphot = deep2_get_ugriz(allphot[keep])

    stars = mrdfits(targpath+'deep2egs-photstars.fits.gz',1)
    stars = stars[where(stars.cfhtls_r lt magcut1,nstar)]
    
    loz = where(allphot.zhelio lt 0.6 and allphot.cfhtls_r lt magcut1,nloz)
    hiz = where(allphot.zhelio gt 0.6 and allphot.zhelio lt 1.2 and $
      allphot.cfhtls_r lt magcut1,nhiz)
    vhiz = where(allphot.zhelio gt 1.2 and allphot.zhelio lt 1.5 and $
      allphot.cfhtls_r lt magcut1,nvhiz)
    
    w1 = where(allphot.zhelio gt 0.6 and allphot.cfhtls_r lt magcut1,nw1)
    w2 = where(allphot.zhelio gt 1.2 and allphot.cfhtls_r lt magcut1,nw2)
    splog, 'Fraction of z>0.6 in grz box: ', n_elements(desi_get_hizelg($
      allphot[w1].ugriz,magcut=magcut1))/float(nw1)
    splog, 'Fraction of z>1.2 in grz box: ', n_elements(desi_get_hizelg($
      allphot[w2].ugriz,magcut=magcut1))/float(nw2)
    
    psfile = cdrpath+'deep2-elg-stars-grz.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
      xmargin=[1.5,0.4], charsize=1.7
    
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.2], yrange=[-0.3,2]
    djs_oplot, allphot[loz].cfhtls_r-allphot[loz].cfhtls_z, $
      allphot[loz].cfhtls_g-allphot[loz].cfhtls_r, psym=symcat(16), symsize=0.3
    djs_oplot, allphot[hiz].cfhtls_r-allphot[hiz].cfhtls_z, $
      allphot[hiz].cfhtls_g-allphot[hiz].cfhtls_r, psym=symcat(9,thick=2), $
      color=cgcolor('firebrick'), symsize=0.3
    djs_oplot, allphot[vhiz].cfhtls_r-allphot[vhiz].cfhtls_z, $
      allphot[vhiz].cfhtls_g-allphot[vhiz].cfhtls_r, psym=symcat(6), $
      color=cgcolor('forest green'), symsize=0.5
    djs_oplot, stars.cfhtls_r-stars.cfhtls_z, symsize=0.15, $
      stars.cfhtls_g-stars.cfhtls_r, psym=symcat(7), color='blue'
       
;; Mostek       
;    rzaxis = range(0.2,1.2,500)
;    djs_oplot, rzaxis, poly(rzaxis,[-0.08,0.68]), line=5, thick=6
;    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[-0.08,0.68])], line=5, thick=6
;    djs_oplot, 1.3*[1,1], [!y.crange[0],poly(1.3,[-0.08,0.68])], line=5, thick=6
;    
;; proposed
;    rzaxis = range(0.1,1.1,500)
;    int = -0.08 & slope = 1.0
;    djs_oplot, rzaxis, poly(rzaxis,[int,slope]), line=0, thick=6
;    djs_oplot, 0.1*[1,1], [!y.crange[0],poly(0.1,[int,slope])], line=0, thick=6
;    djs_oplot, [1.1,1.8], poly(1.1,[int,slope])*[1,1], line=0, thick=6
;    djs_oplot, 1.8*[1,1], [!y.crange[0],poly(1.1,[int,slope])], line=0, thick=6
    
    im_legend, ['Stars','z<0.6','0.6<z<1.2','1.2<z<1.5'], $
      /left, /top, box=0, psym=[7,16,9,6], position=[0.22,0.86], /norm, $
      color=['blue','','firebrick','forest green'], charsize=1.5
    im_legend, ['18.5<r<23'], spacing=2.0, /left, /top, box=0, margin=0
    
;   im_legend, ['Mostek','Proposed'], /right, /bottom, box=0, $
;     line=[0,5], pspacing=1.7, thick=6
    im_plotconfig, psfile=psfile, /psclose, /pdf


    
stop       

; --------------------------------------------------
; Figure 3.4 - example ELG spectrum
    version = desi_deep2_template_version()
    templatepath = getenv('IM_PROJECTS_DIR')+'/desi/templates/'+version+'/'

    templatefile = templatepath+'desi_deep2elg_templates_'+version+'.fits.gz'
    templates = mrdfits(templatefile,0,hdr)
    wave = 10D^make_wave(hdr)
    info = mrdfits(templatefile,1)

;   ww = where(info.z gt 0.99 and info.z lt 1.01 and info.oii_3727 gt 1D-16 and $
;     info.sigma_kms gt 50 and info.sigma_kms lt 85,nww)
;   for ii = 0, nww-1 do begin
;      djs_plot, wave, templates[*,ww[ii]], xsty=3, ysty=3, $ ; /ylog, $
;        psym=-8, xrange=[3707,3747]      ; , xrange=[2000,10000]
;      splog, info[ww[ii]].z, info[ww[ii]].oii_3727, info[ww[ii]].oii_3727_ew, info[ww[ii]].logmstar, $
;        info[ww[ii]].objno, info[ww[ii]].sigma_kms
;      cc = get_kbrd(1)
;   endfor

    this = where(info.objno eq 13057544L)
    flux = templates[*,this]
    flux = flux/interpol(flux,wave,5500)
;   desi_quicksim, wave, flux, model='elg', simdat=simdat

;   im_lineid_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;     xrange=[1500,1E4], yrange=[0.3,200], xtitle='Rest-Frame Wavelength (\AA)', $
;     ytitle='Relative Flux', /ylog
    psfile = cdrpath+'ELG-deep2-example.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=4.5, width=6.7, $
      xmargin=[1.3,0.5], thick=4
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=[1500,1E4], yrange=[0.35,200], xtitle='Rest-Frame Wavelength (\AA)', $
      ytitle='Relative Flux', /ylog
    djs_oplot, wave, flux, psym=10, color=cgcolor('firebrick'), thick=3

;   zrange = where(wave*(1+info[this].z) gt 6550 and wave*(1+info[this].z) lt 9800)
;   djs_oplot, wave[zrange], flux[zrange], psym=10, color=cgcolor('dodger blue'), thick=3
    
    xyouts, 3710, 115, textoidl('[OII]'), /data, align=0.5, charsize=1.3
    xyouts, 4890, 45, textoidl('H\beta+[OIII]'), /data, align=0.5, charsize=1.3
    xyouts, 6560, 90, textoidl('H\alpha+[NII]'), /data, align=0.5, charsize=1.3
    
    xr = 3727+10*[-1,1] & get_element, wave, xr, ww
    djs_plot, [0], [0], /nodata, /noerase, position=[0.21,0.67,0.35,0.9], $ ; /noerase, 
      xsty=9, ysty=9, xrange=xr, yrange=[0,1.1], /norm, xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), xticks=0, yticks=0
    djs_oplot, wave, flux/max(flux[ww[0]:ww[1]]), psym=10, color=cgcolor('dodger blue')
    
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
stop

; --------------------------------------------------
; Figure 3.6 - grz color-color plot coded by redshift

    gr = allphot.cfhtls_g-allphot.cfhtls_r
    rz = allphot.cfhtls_r-allphot.cfhtls_z

    grrange = [-0.5,1.8]
    rzrange = [-0.5,2.2]
    
    grzhist = hist_nd(transpose([[rz],[gr]]),0.1,min=[rzrange[0],grrange[0]],$
      max=[rzrange[1],grrange[1]],reverse_ind=ri)
    sgrzhist = smooth(grzhist,3)
    zmap = grzhist*0.0;-1
    for ii = 0, n_elements(grzhist)-1 do begin
       	if ri[ii] ne ri[ii+1] then begin
           if sgrzhist[ii] gt 3 then $
             zmap[ii] = djs_median(allphot[ri[ri[ii]:ri[ii+1]-1]].zhelio)
        endif
    endfor

    psfile = cdrpath+'grz-zweight.ps'
    cgPS_Open, psfile
     
    cgLoadCT, 16
    TVLCT, cgColor('white', /Triple), 0
    TVLCT, r, g, b, /Get
    palette = [ [r], [g], [b] ]

    zmax = 1.5
    zcrap = where(zmap lt 0,comp=zgood)
    zimage = bytscl(zmap,min=0.0,max=zmax)
;   zimage[zcrap] = cgcolor('white')
    
    cgImage, zimage, XRange=rzrange, YRange=grrange, $
      /Axes, Palette=palette, XTitle='r - z', YTitle='g -r', $
      Position=[0.125, 0.125, 0.9, 0.8]
;   cgContour, grzhist, LEVELS=max(sgrzhist)*[0.25,0.5,0.75,0.9,0.95], /OnImage, $
;     C_Colors=['Tan','Tan', 'Brown'], $
;     C_Annotation=['0.25', '0.5', '0.75'], $
;     C_Thick=thick, C_CharThick=thick
    cgColorbar, Position=[0.125, 0.875, 0.9, 0.925], Title='Mean Redshift', $
      Range=[0,zmax], NColors=254, Bottom=1, $ ; OOB_Low='white', $
      TLocation='Top'
    cgloadct, 0

; Mostek       
    rzaxis = range(0.2,1.3,500)
    djs_oplot, rzaxis, poly(rzaxis,[-0.08,0.68]), line=0, thick=6
    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[-0.08,0.68])], line=0, thick=6
    djs_oplot, 1.3*[1,1], [!y.crange[0],poly(1.3,[-0.08,0.68])], line=0, thick=6

    cgPS_Close
    
    
stop    
    
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.2], yrange=[-0.3,1.8]
    djs_oplot, zcat[hiz].cfhtls_r-zcat[hiz].cfhtls_z, $
      zcat[hiz].cfhtls_g-zcat[hiz].cfhtls_r, psym=symcat(5), $
      color=cgcolor('firebrick'), symsize=0.2


    
; --------------------------------------------------
; Figure 3.10 - [OII] flux vs redshift from DEEP2


return
end
    

    


    
;;; --------------------------------------------------
;;; gr vs rz - stellar contamination
;;    stars = mrdfits(targpath+'deep2egs-photstars.fits.gz',1)
;;    stars = stars[where(stars.cfhtls_r lt magcut1,nstar)]
;;    
;;    loz = where(zcat.zbest lt 0.6 and zcat.cfhtls_r lt magcut1,nloz)
;;    hiz = where(zcat.zbest gt 0.6 and zcat.zbest lt 1.2 and $
;;      zcat.cfhtls_r lt magcut1,nhiz)
;;    vhiz = where(zcat.zbest gt 1.2 and zcat.zbest lt 1.6 and $
;;      zcat.cfhtls_r lt magcut1,nvhiz)
;;    
;;    w1 = where(zcat.zbest gt 0.6 and zcat.cfhtls_r lt magcut1,nw1)
;;    w2 = where(zcat.zbest gt 1.2 and zcat.cfhtls_r lt magcut1,nw2)
;;    splog, 'Fraction of z>0.6 in grz box: ', n_elements(desi_get_hizelg($
;;      zcat[w1].ugriz,magcut=magcut1))/float(nw1)
;;    splog, 'Fraction of z>1.2 in grz box: ', n_elements(desi_get_hizelg($
;;      zcat[w2].ugriz,magcut=magcut1))/float(nw2)
;;    
;;    psfile = cdrpath+'deep2-elg-stars-grz.ps'
;;    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
;;      xmargin=[1.5,0.4], charsize=1.7
;;    
;;    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;;      xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.2], yrange=[-0.3,1.8]
;;    djs_oplot, zcat[loz].cfhtls_r-zcat[loz].cfhtls_z, $
;;      zcat[loz].cfhtls_g-zcat[loz].cfhtls_r, psym=symcat(16), symsize=0.3
;;    djs_oplot, zcat[hiz].cfhtls_r-zcat[hiz].cfhtls_z, $
;;      zcat[hiz].cfhtls_g-zcat[hiz].cfhtls_r, psym=symcat(5), $
;;      color=cgcolor('firebrick'), symsize=0.2
;;    djs_oplot, zcat[vhiz].cfhtls_r-zcat[vhiz].cfhtls_z, $
;;      zcat[vhiz].cfhtls_g-zcat[vhiz].cfhtls_r, psym=symcat(6), $
;;      color=cgcolor('forest green'), symsize=0.5
;;    djs_oplot, stars.cfhtls_r-stars.cfhtls_z, symsize=0.1, $
;;      stars.cfhtls_g-stars.cfhtls_r, psym=symcat(7), color='blue'
;;       
;;;; Mostek       
;;;    rzaxis = range(0.2,1.2,500)
;;;    djs_oplot, rzaxis, poly(rzaxis,[-0.08,0.68]), line=5, thick=6
;;;    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[-0.08,0.68])], line=5, thick=6
;;;    djs_oplot, 1.3*[1,1], [!y.crange[0],poly(1.3,[-0.08,0.68])], line=5, thick=6
;;;    
;;;; proposed
;;;    rzaxis = range(0.1,1.1,500)
;;;    int = -0.08 & slope = 1.0
;;;    djs_oplot, rzaxis, poly(rzaxis,[int,slope]), line=0, thick=6
;;;    djs_oplot, 0.1*[1,1], [!y.crange[0],poly(0.1,[int,slope])], line=0, thick=6
;;;    djs_oplot, [1.1,1.8], poly(1.1,[int,slope])*[1,1], line=0, thick=6
;;;    djs_oplot, 1.8*[1,1], [!y.crange[0],poly(1.1,[int,slope])], line=0, thick=6
;;    
;;    im_legend, ['Stars','z<0.6','0.6<z<1.2','1.2<z<1.6'], $
;;      /left, /top, box=0, psym=[7,16,5,6], position=[0.22,0.86], /norm, $
;;      color=['blue','','firebrick','forest green'], charsize=1.5
;;    im_legend, ['18.5<r<23'], spacing=2.0, /left, /top, box=0, margin=0
;;    
;;;   im_legend, ['Mostek','Proposed'], /right, /bottom, box=0, $
;;;     line=[0,5], pspacing=1.7, thick=6
;;    im_plotconfig, psfile=psfile, /psclose, /pdf

;;    rmin = 19.5
;;    rmax = 24.0
;;    rbin = 0.2
;;    nrbins = long((rmax-rmin)/rbin+1)
;;    raxis = lindgen(nrbins)*rbin+rmin+rbin/2.0
;;    
;;    rzmin = -0.5
;;    rzmax = 2.2
;;    rzbin = 0.15
;;    nrzbins = long((rzmax-rzmin)/rzbin+1)
;;    rzaxis = lindgen(nrzbins)*rzbin+rzmin+rzbin/2.0
;;
;;    grmin = -0.5
;;    grmax = 2.2
;;    grbin = 0.15
;;    ngrbins = long((grmax-grmin)/grbin+1)
;;    graxis = lindgen(ngrbins)*grbin+grmin+grbin/2.0
;;
;;    zmin = 0.6                  ; 0.8
;;    zmax = 1.6
;;    zbin = 0.05
;;    nzbins = long((zmax-zmin)/zbin+1)
;;    zaxis = lindgen(nzbins)*zbin+zmin+zbin/2.0
;;    
;;    ewoiisnr = zcat.oii_3727_2_ew[0]/zcat.oii_3727_2_ew[1]
;;    refindx = where(zcat.zbest gt zmin and zcat.zbest lt zmax and $
;;      zcat.oii_3727[1] ne -2 and ewoiisnr ge 1.0,nrefindx)
;;    rref = zcat[refindx].cfhtls_r
;;    rzref = zcat[refindx].cfhtls_r-zcat[refindx].cfhtls_z
;;    grref = zcat[refindx].cfhtls_g-zcat[refindx].cfhtls_r
;;    zref = zcat[refindx].zbest
;;    oiiref = zcat[refindx].oii_3727
;;
;;    noneindx = where(zcat.zbest gt zmin and zcat.zbest lt zmax and $
;;      (zcat.oii_3727[1] eq -2 or ewoiisnr lt 1.0),nnoneindx)
;;    rnone = zcat[noneindx].cfhtls_r
;;    rznone = zcat[noneindx].cfhtls_r-zcat[noneindx].cfhtls_z
;;    grnone = zcat[noneindx].cfhtls_g-zcat[noneindx].cfhtls_r
;;    znone = zcat[noneindx].zbest
;;    
;;    histref = hist_nd(transpose([$
;;      [rref],$ ; r-band first
;;      [(rzref>rzmin)<rzmax],$ ; r-z
;;      [(grref>grmin)<grmax],$ ; g-r
;;      [zref]]),$              ; redshift
;;      [rbin,rzbin,grbin,zbin],min=[rmin,rzmin,grmin,zmin],$
;;      max=[rmax,rzmax,grmax,zmax],reverse_ind=riref)
;;    histnone = hist_nd(transpose([$
;;      [rnone],$ ; r-band first
;;      [(rznone>rzmin)<rzmax],$ ; r-z
;;      [(grnone>grmin)<grmax],$ ; g-r
;;      [znone]]),$              ; redshift
;;      [rbin,rzbin,grbin,zbin],min=[rmin,rzmin,grmin,zmin],$
;;      max=[rmax,rzmax,grmax,zmax],reverse_ind=rinone)
;;
;;;    for ii = 0, n_elements(histref)-1 do begin
;;;       if riref[ii] ne riref[ii+1] then begin
;;;print, ii & stop
;;;       endif
;;;    endfor
;;
;;;    for ii = 0, n_elements(histnone)-1 do begin
;;;       indx = array_indices(histnone,ii)
;;;       print, raxis[indx[0]], rzaxis[indx[1]], $
;;;         graxis[indx[2]], zaxis[indx[3]]
;;;    endfor
;;
;;    refgood = where(histref gt 0)
;;    for ii = 0, n_elements(histnone)-1 do begin
;;       if rinone[ii] ne rinone[ii+1] then begin
;;;         print, rinone[rinone[ii]:rinone[ii+1]-1]
;;;         print, rnone[rinone[rinone[ii]:rinone[ii+1]-1]], $
;;;           znone[rinone[rinone[ii]:rinone[ii+1]-1]], $
;;;           grnone[rinone[rinone[ii]:rinone[ii+1]-1]], $
;;;           rznone[rinone[rinone[ii]:rinone[ii+1]-1]]
;;;         keepgoing = 1
;;;         jj = ii
;;          if riref[ii] ne riref[ii+1] then begin
;;;            print, rref[riref[riref[ii]:riref[ii+1]-1]], $
;;;              zref[riref[riref[ii]:riref[ii+1]-1]], $
;;;              grref[riref[riref[ii]:riref[ii+1]-1]]
;;
;;; randomly pick one of the "reference" galaxies
;;             theseref = riref[riref[ii]:riref[ii+1]-1]
;;             for kk = 0, histnone[ii]-1 do begin
;;                thisnone = (rinone[rinone[ii]:rinone[ii+1]-1])[kk]
;;                thisref = theseref[fix(randomu(seed,1)*n_elements(theseref))]
;;                zcat[noneindx[thisnone]].oii_3727 = zcat[refindx[thisref]].oii_3727           ; flux
;;                zcat[noneindx[thisnone]].oii_3727_2_ew = zcat[refindx[thisref]].oii_3727_2_ew ; EW
;;;               if thisnone eq 16 then stop
;;             endfor
;;          endif else begin
;;             indx = array_indices(histnone,ii)
;;;            print, raxis[indx[0]], rzaxis[indx[1]], $
;;;              graxis[indx[2]], zaxis[indx[3]]
;;             
;;; find the "nearest" galaxy with a good [OII] measurement in the ND
;;; space                 
;;             for kk = 0, histnone[ii]-1 do begin
;;                thisnone = (rinone[rinone[ii]:rinone[ii+1]-1])[kk]
;;                dist = sqrt((rref-rnone[thisnone])^2+(rzref-rznone[thisnone])^2+$
;;                  (grref-grnone[thisnone])^2+(zref-znone[thisnone])^2)
;;
;;                mindist = min(dist,thisref)
;;
;;;                print, rref[thisref], raxis[indx[0]], rnone[thisnone]
;;;                print, rzref[thisref], rzaxis[indx[1]], rznone[thisnone]
;;;                print, grref[thisref], graxis[indx[2]], grnone[thisnone]
;;;                print, zref[thisref], zaxis[indx[3]], znone[thisnone]
;;;                print
;;
;;                zcat[noneindx[thisnone]].oii_3727 = zcat[refindx[thisref]].oii_3727           ; flux
;;                zcat[noneindx[thisnone]].oii_3727_2_ew = zcat[refindx[thisref]].oii_3727_2_ew ; EW
;;;               print, thisnone
;;;               if thisnone eq 16 then stop
;;;               if (rinone[rinone[ii]:rinone[ii+1]-1])[kk] eq 16 then stop
;;             endfor
;;
;;;             while keepgoing do begin
;;;                jj++
;;;                if riref[jj] ne riref[jj+1] then begin
;;;                   keepgoing = 0
;;;                endif
;;;             endwhile
;;          endelse
;;;         print, rnone[rinone[rinone[ii]:rinone[ii+1]-1]], rref[riref[riref[jj]:riref[jj+1]-1]], $
;;;           znone[rinone[rinone[ii]:rinone[ii+1]-1]], zref[riref[riref[jj]:riref[jj+1]-1]], $
;;;           grnone[rinone[rinone[ii]:rinone[ii+1]-1]], grref[riref[riref[jj]:riref[jj+1]-1]], $
;;;           rznone[rinone[rinone[ii]:rinone[ii+1]-1]], rzref[riref[riref[jj]:riref[jj+1]-1]]
;;;         print & print, jj, oiiref[*,riref[riref[jj]:riref[jj+1]-1]] & print
;;
;;;; randomly pick one of the "reference" galaxies
;;;          theseref = riref[riref[jj]:riref[jj+1]-1]
;;;          for kk = 0, histnone[ii]-1 do begin
;;;             thisnone = (rinone[rinone[ii]:rinone[ii+1]-1])[kk]
;;;             thisref = theseref[fix(randomu(seed,1)*n_elements(theseref))]
;;;             zcat[noneindx[thisnone]].oii_3727 = zcat[refindx[thisref]].oii_3727 ; flux
;;;             zcat[noneindx[thisnone]].oii_3727_2_ew = zcat[refindx[thisref]].oii_3727_2_ew ; EW
;;;;            if zcat[refindx[thisref]].oii_3727[1] eq -2 then stop
;;;
;;;; print this as a diagnostic             
;;;             print, zcat[noneindx[thisnone]].zbest, zcat[refindx[thisref]].zbest, $
;;;               zcat[noneindx[thisnone]].cfhtls_r, zcat[refindx[thisref]].cfhtls_r
;;;
;;;;            if thisnone eq 100 then stop
;;;;            if (rinone[rinone[ii]:rinone[ii+1]-1])[kk] eq 100 then stop
;;;          endfor
;;;         if histref[jj] gt 1 then stop
;;;         if total(rinone[rinone[ii]:rinone[ii+1]-1] eq 100) gt 0 then stop
;;       endif
;;    endfor
