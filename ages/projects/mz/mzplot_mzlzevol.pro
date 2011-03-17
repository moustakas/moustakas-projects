pro mzplot_mzlzevol, ps=ps
; jm09mar27nyu - plot the evolution of the LZ and MZ relations
; jm10oct10ucsd - major update    

    mzpath = ages_path(/projects)+'mz/'
    pspath = ages_path(/papers)+'mz/FIG_MZ/'
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

; read the data    
    zbins = mz_zbins(nzbins)
    agesancillary = read_mz_sample(/mzhii_ancillary)
    agesmass = read_mz_sample(/mzhii_mass)
    agesohdust = read_mz_sample(/mzhii_log12oh)

    mzevol = mrdfits(mzpath+'mzevol.fits.gz',1)
    lzevol = mrdfits(mzpath+'lzevol_B.fits.gz',1)

    mzlocal = mrdfits(mzpath+'mzlocal_sdss_brokenpl.fits.gz',1)
    lzlocal = mrdfits(mzpath+'lzlocal_sdss.fits.gz',1)
    ncalib = n_elements(lzlocal)

; --------------------------------------------------
; AGES/MZ and LZ evolution
    massrange1 = [8.1,12.0]
    magrange1 = [-16.3,-23.7]
    ohrange1 = [8.35,9.37] ; [8.3,9.39]

    localline = 0
    localcolor = 'black'
    evolline = 5
    evolcolor = 'firebrick'
    
    for ii = 0, ncalib-1 do begin
       t04 = 0 & m91 = 0 & kk04 = 0
       case ii of
          0: t04 = 1
          1: m91 = 1
          2: kk04 = 1
       endcase
       if keyword_set(t04) then calib = 't04'
       if keyword_set(m91) then calib = 'm91'
       if keyword_set(kk04) then calib = 'kk04'

       ainfo = mzlz_grab_info(agesohdust,agesancillary,agesmass,$
         t04=t04,m91=m91,kk04=kk04)
       ohtitle = mzplot_ohtitle(t04=t04,m91=m91,kk04=kk04)
       
       psfile = pspath+'mzevol_'+calib+suffix
       mzplot_sixpanel, ainfo.z, ainfo.mass, ainfo.oh, ainfo.weight, $
         psfile=psfile, xtitle=mzplot_masstitle(), ytitle=ohtitle, /ages, $
         xrange=massrange1, yrange=ohrange1, npix=10, mzlocal=mzlocal[ii], $
         mzevol=mzevol[ii], localline=localline, localcolor=localcolor, $
         evolline=evolline, evolcolor=evolcolor
       
       psfile = pspath+'lzevol_'+calib+suffix
       mzplot_sixpanel, ainfo.z, ainfo.mb_ab, ainfo.oh, ainfo.weight, $
         psfile=psfile, xtitle=mzplot_mbtitle(), ytitle=ohtitle, /ages, $
         xrange=magrange1, yrange=ohrange1, npix=10, lzlocal=lzlocal[ii], $
         lzevol=lzevol[ii], localline=localline, localcolor=localcolor, $
         evolline=evolline, evolcolor=evolcolor
    endfor
       
stop
    
return
end
