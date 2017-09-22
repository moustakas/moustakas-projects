pro noao_report2015, clean=clean, trilogy=trilogy
; jm15jun24siena - make a pretty picture of the interacting galaxy UGC12589 for
; the cover of the monthly NOAO Currents issue


; It would be nice to show more than just the galaxy aspects of our survey. That
; might induce a broader audience to use our data … How about a four-panel
; figure that has the following objects: Interacting system:
    
; http://www.legacysurvey.org/viewer/?ra=351.26&dec=0

; Bow shock:
; http://legacysurvey.org/viewer/?ra=325.6956&dec=1.0094&zoom=14&layer=decals-dr1j
; 
; Star cluster:
; http://legacysurvey.org/viewer/?ra=114.5974&dec=21.5719&zoom=12&layer=decals-dr1j
; 
; Nice galaxy:
; http://legacysurvey.org/viewer/?ra=324.5394&dec=8.9602&zoom=14&layer=decals-dr1j
; or
; http://legacysurvey.org/viewer/?ra=324.9098&dec=8.9448&zoom=14&layer=decals-dr1j

; starcluster: runbrick --threads 24 --no-write --stage image_coadds --radec 114.5974 21.5719 --width 2500 --height 2500
; nicegal1:    runbrick --threads 24 --no-write --stage image_coadds --radec 324.5394 8.9602 --width 2500 --height 2500 --pixscale 0.157
; nicegal2:    runbrick --threads 24 --no-write --stage image_coadds --radec 324.9098 8.9448 --width 2500 --height 2500 --pixscale 0.104
; bowshock:    runbrick --threads 24 --no-write --stage image_coadds --radec 325.6784 0.9916 --width 2500 --height 2500 --pixscale 0.115 --gpsf

    band = ['g','r','z']
    obj = ['starcluster','nicegal1','nicegal2','bowshock']
    brick = ['custom-114597p21571','custom-324539p08960','custom-324909p08944','custom-325678p00991']
    noiselum = ['0.2','0.25','0.15','0.3']
    
    obj = ['nicegal2']
    brick = ['custom-324909p08944']
    noiselum = ['0.15']
    
; then interpolate over the saturated star
    if keyword_set(clean) then begin
       for io = 0, n_elements(obj)-1 do begin
          for ii = 0, n_elements(band)-1 do begin
             im = mrdfits('coadd/cus/'+brick[io]+'/decals-'+brick[io]+'-image-'+band[ii]+'.fits',0,hdr)
;            fim = djs_maskinterp(im,im eq 0,iaxis=1)
;            sim = gauss_smooth(fim,2.0,/edge_truncate)
             if strmatch(obj[io],'*bowshock*') eq 0 then begin
                sim = djs_maskinterp(im,im eq 0,iaxis=1)
                if strmatch(obj[io],'*star*') then begin
                   sim[1000:1300,2200:*] = im[1000:1300,2200:*]
                   sim[1570:*,400:730] = im[1570:*,400:730]
                endif
             endif else begin
                sim = im
             endelse
;            if strmatch(obj[io],'*bowshock*') eq 0 then begin
;               sim = djs_maskinterp(im,im eq 0,iaxis=1)
;            endif else begin
;               sim = im
;            endelse
             mwrfits, sim, obj[io]+'-'+band[ii]+'.fits', hdr, /create
          endfor
       endfor
    endif

; call trilogy
    if keyword_set(trilogy) then begin
      for io = 0, n_elements(obj)-1 do begin
          openw, lun, 'trilogy.in', /get_lun
          printf, lun, 'B'
          printf, lun, obj[io]+'-'+band[0]
          printf, lun, ''
          printf, lun, 'G'
          printf, lun, obj[io]+'-'+band[1]
          printf, lun, ''
          printf, lun, 'R'
          printf, lun, obj[io]+'-'+band[2]
          printf, lun, ''
          printf, lun, 'indir .'
          printf, lun, 'outname '+obj[io]+'-grz'
          printf, lun, 'noiselum '+noiselum[io]
          printf, lun, 'satpercent  0.001'
;         printf, lun, 'stampsize 2000'
;         printf, lun, 'samplesize 2500'
          printf, lun, 'legend 0'
          printf, lun, 'deletetests 1'
          printf, lun, 'show 0'
          printf, lun, 'testfirst 0'
          free_lun, lun
          spawn, 'trilogy trilogy.in'
       endfor
    endif

return
end
