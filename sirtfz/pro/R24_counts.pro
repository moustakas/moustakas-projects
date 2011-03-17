pro R24_counts, flimit=flimit, sigma=sigma, noevol=noevol
; jm01aug6uofa
; calculate the surface density of 24mu source counts using the
; evolutionary models of Dole et al and at R using Huang et al 2001
; (CADIS deep survey)

    if not keyword_set(flimit) then flimit = 3D-5 ; limiting flux at 24mu
    if not keyword_set(sigma) then sigma = 3.0    ; sigma-detection
    
    path = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='lib')
    if keyword_set(noevol) then cfile = 'dole_counts_0024.noevol.dat' else $
      cfile = 'dole_counts_0024.dat'

; flux is in log flux [Jy] and counts is cumulative source counts > flux [gal/sr]
    
    readcol, path+cfile, logflux, counts, format='F,D,X', /silent
    flux = 10.0^logflux ; [Jy]
    
    sdensity = counts*(!dpi/180D)^2 ; [gal/deg^2]
    sdensity_flimit = interpol(sdensity,flux,sigma*flimit) ; interpolate the integral counts

    print, format='("Surface density [gal/deg^2] at 24mu: ",G0.0,".")', sdensity_flimit

; R-band counts

    rmin = 16.0 & rmax = 23.0 & dr = 0.5
    rmag = (findgen((rmax-rmin+dr)/dr))*dr+rmin

    rcounts = [0.69,1.39,1.17,1.39,1.87,2.10,2.44,2.68,$ ; #/0.5mag/deg^2
               2.86,3.11,3.27,3.47,3.60,3.78,4.05]
    fit = poly_fit(rmag,rcounts,2)

;   rmagnew = congrid(rmag,10*n_elements(rmag),/interp)
    rmagnew = (findgen((26-rmin+0.05)/0.05))*0.05+rmin
    newcounts = poly(rmagnew,fit)

    rdensity = cumulate((10^newcounts)*0.05) ; cumulative LF [gal/deg^2]
    rdensity_26 = total((10^newcounts)*0.05) ; [gal/deg^2] down to R=26
    
    print, format='("Surface density [gal/deg^2] in R: ",G0.0,".")', rdensity

; ----------------------------------------------------------------------
    
    beam = 2.0*206265*24E-6/0.85/3600.0 ; 24mu beamsize (deg)
    p24 = 1.0 - exp(-!pi*sdensity*beam^2)
    pr = 1.0 - exp(-!pi*rdensity*beam^2)

; ----------------------------------------------------------------------

    rdensity_26_beam = rdensity_26*beam^2.0
    sdensity_flimit_beam = sdensity_flimit*beam^2.0

    ps_open, '24mu_counts', /ps_fonts
    device, /inches, /times
    
    djs_plot, logflux, p24, xsty=3, ysty=3, xtitle='log Flux [Jy]', $
      ytitle='Probability', title='24 \mu'+'m Counts (Dole)', charsize=2.0, $
      charthick=2.0, thick=2.0, xthick=2.0, ythick=2.0,       yrange=[0,1]
    oplot, alog10(3.0*flimit)*[1,1], !y.crange
    legend, ['Limiting flux 90 '+textoidl('\mu')+'Jy','Beam '+strn(2.0*206265*24E-6/0.85,length=4)+'"'], $
      /right, /top, box=0, charsize=2.0, thick=2.0

    ps_close


    ps_open, 'r_counts', /ps_fonts
    device, /inches, /times
    
    djs_plot, rmagnew, pr, xsty=3, ysty=3, xtitle='R Magnitude', $
      ytitle='Probability', title='R-Band Counts (Huang)', charsize=2.0, $
      charthick=2.0, thick=2.0, xthick=2.0, ythick=2.0, xr=reverse(minmax(rmagnew)), $
      yrange=[0,1]

    ps_close


stop
    
return
end
