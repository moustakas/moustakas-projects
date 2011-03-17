pro eigen_candidates, write=write, doplot=doplot

    outpath = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='jacoby_atlas/eigen')
    
    light = 2.99792458D5 ; speed of light [km/s]

    eigenres = 4.5 ; Jacoby atlas resolution [Angstrom]
    
; read the jacoby atlas and sort by spectral type

    jhc = read_jacoby_atlas()

    srt = sort(strcompress(jhc.sp_type_lit,/remove))
    jhc = jhc[srt]
    
; only keep normal stars that have a radial velocity measurement 

    keep = where(strmatch(jhc.comment,'*normal star*') eq 1B or $
      strmatch(jhc.comment,'*cluster star*') eq 1B and $
      jhc.rv gt -900.0)
    jhc = jhc[keep]
    
; remove all O stars (lifetimes are too short), and GKM main-sequence
; stars (too faint)

    littype = strcompress(jhc.sp_type_lit,/remove)
    jhctype = strcompress(jhc.sp_type,/remove)

    keep = where((strmatch(littype,'O*') eq 0B) and (strmatch(littype,'G*V*') eq 0B) and $
      (strmatch(littype,'K*V*') eq 0B) and (strmatch(littype,'M*V*') eq 0B))
    jhc = jhc[keep]

; remove the following stars

    remstars = [$
      'HD12027',  $ ; pixel spike
      'HD23733',  $ ; stitching problem(s)
      'HD167451', $ ; pixel spike
      'HD14357',  $ ; small pixel spike
      'HD240344', $ ; strange near H-beta
      'HD13505',  $ ; stitching problem(s)
      'HD94028',  $ ; unwanted spectral type (F4V)
      'SAO102986',$ ; unwanted spectral type (F6V)
      'BD+051668',$ ; unwanted spectral type (M5V)
      'HD37767',  $ ; stitching problem(s)
      'HD39136',  $ ; stitching problem(s)
      'HD240296', $ ; stitching problem(s)
      'HD15866',  $ ; stitching problem(s)
      'SAO20899', $ ; stitching problem(s)
      'HD192832', $ ; stitching problem(s)
      'SAO87716', $ ; stitching problem(s)
      'HD842',    $ ; stitching problem(s)
      'HD186293', $ ; stitching problem(s)
      'HD1069',   $ ; stitching problem(s)
      'HD17647',  $ ; unwanted spectral type (G1V)
      'HD22193',  $ ; unwanted spectral type (G6V)
      'HD107214', $ ; unwanted spectral type (G2)
      'SAO57199', $ ; unwanted spectral type?? (F6V)
      'HD6111',   $ ; unwanted spectral type?? (F8V)
      'HD31084',  $ ; unwanted spectral type?? (F9V)
      'HD83140',  $ ; unwanted spectral type?? (F3IV)
      'HD260655', $ ; unwanted spectral type?? (M0V/K7)
      'BD+630137' $ ; unwanted spectral type?? (M1V/K5)
      ]

    doit = match_string(remstars,jhc.star,/exact,index=indx)
    keep = lindgen(n_elements(jhc))
    remove, indx, keep

    jhc = jhc[keep]

    srt = sort(jhc.jhc_id)
    jhc = jhc[srt]

    if keyword_set(write) then begin

       openw, lun1, outpath+'README', /get_lun
       printf, lun1, '# Written by EIGEN_CANDIDATES [J. Moustakas '+im_today()+']'
       printf, lun1, '# The Jacoby, Hunter, & Christian (1984) stellar spectral atlas contains 162 stars.  Stars were '
       printf, lun1, '# removed for the following reasons: [1] no radial velocity measurement; [2] O stars (emission '
       printf, lun1, '# lines and/or lifetimes too short; [3] GKM main sequence stars (too faint); [4] an additional '
       printf, lun1, '# 27 stars with either a "stitching" problem, funny pixel spikes, or unwanted spectral types '
       printf, lun1, '# (e.g., MV).  The five columns in the table below have the following definition, respectively: '
       printf, lun1, '# '
       printf, lun1, '#     JHC_ID      - identification number in the original atlas; '
       printf, lun1, '#     STAR        - star name in the original atlas (SIMBAD compliant); '
       printf, lun1, '#     SP_TYPE     - spectral type in the original atlas; '
       printf, lun1, '#     SP_TYPE_LIT - SIMBAD (literature) spectral type; '
       printf, lun1, '#     RV          - heliocentric radial velocity (SIMBAD) [km/s]'
       printf, lun1, '# '
       struct_print, struct_trimtags(jhc,select=['jhc_id','star','sp_type','sp_type_lit','rv']), $
         lun=lun1, /no_head
       free_lun, lun1
       
    endif else begin

       print_struct, jhc, ['jhc_id','star','sp_type','sp_type_lit','rv']

    endelse

    starcount = n_elements(jhc)
    
    crval1 = float(jhc[0].wave[0])
    cd1_1 = float(jhc[0].wave[1]-jhc[0].wave[0])
    refwave = float(jhc[0].wave)
    npix = n_elements(refwave)

    for i = 0L, starcount-1L do begin

; correct for radial velocity
       
       flux = jhc[i].flux
       wave = jhc[i].wave
       rv = jhc[i].rv
       
       restwave = wave/sqrt((1+rv/light)/(1-rv/light)) ; rest wavelength
       starflux = interpol(flux,restwave,refwave)
       
       if keyword_set(doplot) then begin
          djs_plot, refwave, starflux, xsty=3, ysty=3, charsize=2.0, charthick=2.0, $
            xtitle='Wavelength ('+angstrom()+')', ytitle='Relative f_{\lambda}', $
            xthick=2.0, ythick=2.0, yrange=[-0.1,1.1]
          legend, strcompress([string(jhc[i].jhc_id,format='(I0)'),jhc[i].star,$
            jhc[i].sp_type],/remove), /right, /top, box=0, $
            charsize=2.0, charthick=2.0
          cc = get_kbrd(1)
       endif

       mkhdr, header, float(starflux)
       sxaddpar, header, 'OBJECT', strcompress(jhc[i].star,/remove), ' JHC84 star name'
       sxaddpar, header, 'JHC_ID', strcompress(jhc[i].jhc_id,/remove), ' JHC84 ID number'
       sxaddpar, header, 'RV', float(jhc[i].rv), ' radial velocity [km/s]'
       sxaddpar, header, 'RV_ERR', float(jhc[i].rv_err), ' radial velocity error [km/s]'
       sxaddpar, header, 'SPTYP', strcompress(jhc[i].sp_type,/remove), ' JHC84 spectral type'
       sxaddpar, header, 'SPTYPLIT', strcompress(jhc[i].sp_type_lit,/remove), ' SIMBAD spectral type'
       sxaddpar, header, 'CRVAL1', crval1, ' central wavelength of first pixel'
       sxaddpar, header, 'CD1_1', cd1_1, ' dispersion [Angstrom/pixel]'
       sxaddpar, header, 'CRPIX1', 0, ' starting pixel (0-indexed)'
       sxaddpar, header, 'CTYPE1', 'LINEAR'
       sxaddpar, header, 'DC-FLAG', 0, ' log-linear flag'
       sxaddpar, header, 'EIGENRES', eigenres, 'resolution [Angstrom]'

; write out

       outname = strlowcase(strcompress(jhc[i].star,/remove))+'.fits'
       if keyword_set(write) then begin
          print, 'Writing '+outpath+outname+'.'
          writefits, outpath+outname, float(starflux), header
       endif

    endfor 

return
end
