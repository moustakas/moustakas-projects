pro choose_ancillary_decamdata
; jm15mar25siena - prioritize the set of DECam observations
; that aren't DECaLS data

    path = getenv('IM_PROJECTS_DIR')+'/decals/ancillary/'
    cat = mrdfits(path+'noaoarchive_decalsplus_allbands.fits',1)
;   cat = mrdfits(path+'noaoarchive_decalsplus_grz.fits',1)

    tilepath = '/Users/ioannis/desi/decam/code/dtiling/trunk/'
    alltile = mrdfits(tilepath+'decam-tiles.fits',1)

; --------------------------------------------------
; get the previously observed DECam tiles that overlap the DEEP2/Field
; 3 field
    dd = read_deep2_zcat()

    ww = where(strmid(strtrim(dd.objno,2),0,1) eq 3)
    spherematch, dd[ww].ra, dd[ww].dec, cat.ra, $
      cat.dec, 1D, m1, m2, maxmatch=0
    field3 = m2[uniq(m2,sort(m2))]

    djs_plot, dd[ww].ra, dd[ww].dec, psym=6, xsty=3, ysty=3, $
      xrange=[355,350], yrange=[-1,1], $
      xtitle='RA', ytitle='Dec', title='DEEP2/Field-3'
    djs_oplot, cat.ra, cat.dec, psym=8, symsize=3, color='red'
    djs_oplot, cat[field3].ra, cat[field3].dec, $
      psym=7, symsize=3, color='yellow'
    struct_print, cat[field3]

    ww = where(strmid(strtrim(dd.objno,2),0,1) eq 4)
    spherematch, dd[ww].ra, dd[ww].dec, cat.ra, $
      cat.dec, 1D, m1, m2, maxmatch=0
    field4 = m2[uniq(m2,sort(m2))]

    djs_plot, dd[ww].ra, dd[ww].dec, psym=6, xsty=3, ysty=3, $
      xrange=[39,35], yrange=[-1,2], $
      xtitle='RA', ytitle='Dec', title='DEEP2/Field-4'
    djs_oplot, cat.ra, cat.dec, psym=8, symsize=3, color='red'
    djs_oplot, cat[field4].ra, cat[field4].dec, $
      psym=7, symsize=3, color='yellow'
    struct_print, cat[field4]

    out = cat[[field3,field4]]
    im_mwrfits, out, path+'noaoarchive_decalsplus_highpriority.fits'


; --------------------------------------------------

    mwrfits, cat[field3], path+'decam-deep2-field2.fits', /create


    field4 = where(cat.ra gt 36.5 and cat.ra lt 38.1 and $
      cat.dec gt 0.3 and cat.dec lt 0.9) ; DEEP2/Field 4
    field3 = where(cat.ra gt 351.0 and cat.ra lt 354.0 and $
      cat.dec gt -0.11 and cat.dec lt 0.45) ; DEEP2/Field 3
    cosmos = where(cat.ra gt 149.0 and cat.ra lt 151.0 and $
      cat.dec gt 1.45 and cat.dec lt 2.95) ; COSMOS

    these = [field4,field3,cosmos]
    out = cat[these]

    djs_plot, cat.ra, cat.dec, psym=6, xsty=3, ysty=3
    djs_oplot, out.ra, out.dec, psym=6, color='red'

    gg = where(strtrim(strmid(cat.filter,0,2),2) eq 'g')
    djs_plot, cat.ra, cat.dec, psym=6, xsty=3, ysty=3
    djs_oplot, out.ra, out.dec, psym=6, color='red'
      

; make sure the cosmos tiles are in DR1
    dr1 = mrdfits(path+'dr1/decals-bricks-in-dr1.fits',1)
    djs_plot, dr1.ra, dr1.dec, psym=6, xsty=3, ysty=3
    djs_oplot, cat[cosmos].ra, cat[cosmos].dec, psym=6, color='red'


; --------------------------------------------------
; get the DECaLS tile numbers corresponding to the extragalactic deep
; fields we need to be targeting this spring

; find out which Bootes field tiles to observe
    dd = mrdfits(path+'ages_photometry_v4.0.fits.gz',1)
    print, djs_median(dd.ra), djs_median(dd.dec)

    spherematch, dd.ra, dd.dec, alltile.ra, alltile.dec, 1D, m1, m2, maxmatch=0
    these = m2[uniq(m2,sort(m2))]
    struct_print, alltile[these]

; overlay the clusters from Brodwin+
    cra = ['14:29:15.16','14:32:29.18','14:26:09.51','14:25:03.44',$
      '14:26:30.42','14:34:30.44','14:29:18.51','14:32:38.38',$
      '14:25:19.33','14:33:51.14','14:34:46.33','14:32:18.31','14:25:18.50',$
      '14:38:08.71','14:31:08.06','14:32:24.16']
    cdec = ['+33:57:08.5','+33:32:36.0','+34:03:41.1','+35:20:10.4','+33:39:33.2',$
      '+34:27:12.3','+34:37:25.8','+34:36:49.0','+34:28:38.2','+33:25:51.1',$
      '+35:19:33.5','+32:53:07.8','+32:50:40.5','+34:14:19.2','+34:59:43.3',$
      '+32:50:03.0']

    djs_plot, dd.ra, dd.dec, psym=6, xsty=3, ysty=3, $
      xrange=[222,214], yrange=[31.5,36.5], $
      xtitle='RA', ytitle='Dec', title='Bootes/AGES'
    djs_oplot, alltile.ra, alltile.dec, psym=8, symsize=3, color='red'
    djs_oplot, alltile[these].ra, alltile[these].dec, $
      psym=7, symsize=3, color='yellow'
;   tvcircle, 1.0, alltile[these].ra, alltile[these].dec, $
;     /data, color=djs_icolor('red')
    djs_oplot, 15D*hms2dec(cra), hms2dec(cdec), psym=symcat(15), $
      color=cgcolor('forest green'), symsize=2


    mwrfits, alltile[these], path+'decam-ages.fits', /create
    
; find out which DEEP2/Field 2 tiles to observe
    dd = read_deep2_zcat()
    ww = where(strmid(strtrim(dd.objno,2),0,1) eq 2)

    spherematch, dd[ww].ra, dd[ww].dec, alltile.ra, $
      alltile.dec, 1D, m1, m2, maxmatch=0
    these = m2[uniq(m2,sort(m2))]

    djs_plot, dd[ww].ra, dd[ww].dec, psym=6, xsty=3, ysty=3, $
      xrange=[250,256], yrange=[34,36], $
      xtitle='RA', ytitle='Dec', title='DEEP2/Field-2'
    djs_oplot, alltile.ra, alltile.dec, psym=8, symsize=3, color='red'
    djs_oplot, alltile[these].ra, alltile[these].dec, $
      psym=7, symsize=3, color='yellow'
    tvcircle, 1.0, alltile[these].ra, alltile[these].dec, $
      /data, color=djs_icolor('red')
    struct_print, alltile[these]

    mwrfits, alltile[these], path+'decam-deep2-field2.fits', /create

; find out which DEEP2/Field 2 tiles to observe
    dd = read_deep2_zcat()
    ww = where(strmid(strtrim(dd.objno,2),0,1) eq 2)

    spherematch, dd[ww].ra, dd[ww].dec, alltile.ra, $
      alltile.dec, 1D, m1, m2, maxmatch=0
    these = m2[uniq(m2,sort(m2))]

    djs_plot, dd[ww].ra, dd[ww].dec, psym=6, xsty=3, ysty=3, $
      xrange=[250,256], yrange=[34,36], xtitle='RA', ytitle='Dec', $
      title='DEEP2/Field-2'
    djs_oplot, alltile.ra, alltile.dec, psym=8, symsize=3, color='red'
    djs_oplot, alltile[these].ra, alltile[these].dec, $
      psym=7, symsize=3, color='yellow'
    tvcircle, 1.0, alltile[these].ra, alltile[these].dec, $
      /data, color=djs_icolor('red')
    struct_print, alltile[these]

    mwrfits, alltile[these], path+'decam-deep2-field2.fits', /create

; find out which cosmos tiles to observe
    dd = mrdfits(path+'cosmos_phot_20060103.fits.gz',1)
    racen = djs_median(dd.ra)
    deccen = djs_median(dd.dec)

    spherematch, dd.ra, dd.dec, alltile.ra, $
      alltile.dec, 0.8D, m1, m2, maxmatch=0
    these = m2[uniq(m2,sort(m2))]
    these2 = [0,4,7]
;   these2 = [0,3,4,6,7]

    custom = alltile[these[these2]]
    custom.ra = custom.ra-custom[1].ra+racen
    custom.dec = custom.dec-custom[1].dec+deccen

    djs_plot, dd.ra, dd.dec, psym=6, xsty=3, ysty=3, $
      xrange=[152,148], yrange=[0,4], $
      xtitle='RA', ytitle='Dec', title='COSMOS'
    djs_oplot, alltile.ra, alltile.dec, psym=8, symsize=3, color='red'
    djs_oplot, alltile[these].ra, alltile[these].dec, $
      psym=7, symsize=3, color='yellow'
;   djs_oplot, alltile[these[these2]].ra, alltile[these[these2]].dec, $
;     psym=6, symsize=3, color=cgcolor('dodger blue')
;   tvcircle, 1.0, alltile[these].ra, alltile[these].dec, $
;     /data, color=djs_icolor('red')
;   tvcircle, 1.0, alltile[these[these2]].ra, alltile[these[these2]].dec, $
;     /data, color=cgcolor('dodger blue')
    djs_oplot, custom.ra, custom.dec, $
      psym=6, symsize=3, color=cgcolor('forest green')
    tvcircle, 1.0, custom.ra, custom.dec, $
      /data, color=cgcolor('forest green')
    struct_print, alltile[these]
    print
    struct_print, custom

    mwrfits, custom, path+'decam-cosmos-custom.fits', /create

; --------------------------------------------------
; get DEEP2 field centers and widths
    dd = read_deep2_zcat()
    field = ['1','2','3','4']
    for ii = 0, 3 do begin
       ww = where(strmid(strtrim(dd.objno,2),0,1) eq field[ii])
       ra = djs_median(dd[ww].ra)
       dec = djs_median(dd[ww].dec)
       dra = max(dd[ww].ra)-min(dd[ww].ra)
       ddec = max(dd[ww].dec)-min(dd[ww].dec)
       print, 'DEEP2-Field-'+field[ii], ra, dec, dra, ddec
    endfor

stop

return
end
