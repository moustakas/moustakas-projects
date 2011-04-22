pro mzplot_mzlzevol, ps=ps
; jm09mar27nyu - plot the evolution of the LZ and MZ relations
; jm10oct10ucsd - major update    

    mzpath = mz_path()
    qapath = mzpath+'qaplots/'
    if keyword_set(ps) then begin
       pspath = qapath
       suffix = '.ps'
    endif else begin
       pspath = mz_path(/paper)
       suffix = '.eps'
    endelse

; read the data    
    zbins = mz_zbins(nzbins)
    agesancillary = read_mz_sample(/mzhii_ancillary)
    agesmass = read_mz_sample(/mzhii_mass)
    agesohdust = read_mz_sample(/mzhii_log12oh)
    agesohnodust = read_mz_sample(/mzhii_log12oh,/nodust)

    mzevol = mrdfits(mzpath+'mzevol.fits.gz',1)
    lzevol = mrdfits(mzpath+'lzevol_B.fits.gz',1)

; --------------------------------------------------
; AGES/MZ and LZ evolution
    massrange1 = [8.1,12.0]
    magrange1 = [-16.3,-23.7]
    ohrange1 = [8.35,9.37] ; [8.3,9.39]

    localline = 0
    localcolor = 'black'
    evolline = 5
    evolcolor = 'firebrick'
    
    for ii = 2, 2 do begin
       case ii of
          0: begin
             t04 = 1 & m91 = 0 & kk04 = 0
             calib = 't04'
             ohrange1 = [8.3,9.35]
          end
          1: begin
             t04 = 0 & m91 = 1 & kk04 = 0
             calib = 'm91'
             ohrange1 = [8.3,9.2]
          end
          2: begin
             t04 = 0 & m91 = 0 & kk04 = 1
             calib = 'kk04'
             ohrange1 = [8.55,9.35]
          end
       endcase
       ohtitle = mzplot_ohtitle(t04=t04,m91=m91,kk04=kk04)

       mztest = mrdfits(mzpath+'mzlocal_sdss_ews_'+calib+'.fits.gz',1)
       mzlocal = mrdfits(mzpath+'mzlocal_sdss_fluxcor_'+calib+'.fits.gz',1)
       lzlocal = mrdfits(mzpath+'lzlocal_sdss_fluxcor_'+calib+'.fits.gz',1)

;      ainfo = mzlz_grab_info(agesohnodust,agesancillary,agesmass,$
;        t04=t04,m91=m91,kk04=kk04,/nolimit,/flux,zmin=0.05,zmax=0.15)
       ainfo = mzlz_grab_info(agesohdust,agesancillary,agesmass,$
         t04=t04,m91=m91,kk04=kk04,/nolimit)
       
       psfile = pspath+'mzevol_'+calib+suffix
       mzplot_sixpanel, ainfo.z, ainfo.mass, ainfo.oh, ainfo.weight, $
         oh_err=ainfo.oh_err, psfile=psfile, xtitle=mzplot_masstitle(), $
         ytitle=ohtitle, /ages, xrange=massrange1, yrange=ohrange1, npix=10, $
         mzlocal=mzlocal, mzevol=mzevol[ii], localline=localline, $
         localcolor=localcolor, evolline=evolline, evolcolor=evolcolor, $
         postscript=keyword_set(ps), mztest=mztest
       
       psfile = pspath+'lzevol_'+calib+suffix
       mzplot_sixpanel, ainfo.z, ainfo.mb_ab, ainfo.oh, ainfo.weight, $
         oh_err=ainfo.oh_err, psfile=psfile, xtitle=mzplot_mbtitle(), $
         ytitle=ohtitle, /ages, xrange=magrange1, yrange=ohrange1, npix=10, $
         lzlocal=lzlocal, lzevol=lzevol[ii], localline=localline, $
         localcolor=localcolor, evolline=evolline, evolcolor=evolcolor, $
         postscript=keyword_set(ps)
    endfor
       
stop
    
return
end
