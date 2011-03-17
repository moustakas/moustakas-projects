pro write_99veilleux, vout
; jm03sep15uofa
; jm05may03uofa - updated

; the veilleux sample consists of two papers, one from 1995 and the
; other from 1999.      
    
; ---------------------------------------------------------------------------
; 1995 sample:  the galaxy names and classes are in
; VEILLEUX95_CLASS.DAT (typed by hand).  the galaxy list sent to NED
; was written by WRITE_V95_NED_INPUT.PRO and lives in
; VEILLEUX95_NED_INPUT.TXT.  the NED basic data are contained in
; VEILLEUX95.NED which we will write to VEILLEUX95_NED.FITS.GZ below. 

; parse and read the NED results; ned screws up on several objects so
; we need to hand-type those coordinates and redshifts

    parse_ned_byname, 'veilleux95.ned', v95_ned, $
      inputnedfile='veilleux95_ned_input.txt', outfile='veilleux95_ned.fits'

    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS14341+3017NW*'))].ra  = '14:36:16.8'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS14341+3017NW*'))].dec = '+30:04:14'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS14341+3017NW*'))].z   = 0.034390

    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS14341+3017SE*'))].ra  = '14:36:16.8'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS14341+3017SE*'))].dec = '+30:04:14'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS14341+3017SE*'))].z   =  0.034591
    
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15391+3214SW*'))].ra  = '15:41:05.4'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15391+3214SW*'))].dec = '+32:04:51'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15391+3214SW*'))].z   =  0.053303 

    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15391+3214NE*'))].ra  = '15:41:05.4'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15391+3214NE*'))].dec = '+32:04:51'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15391+3214NE*'))].z   =  0.052603

    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15543+4158NW*'))].ra  = '15:56:04.8'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15543+4158NW*'))].dec = '+41:49:31'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15543+4158NW*'))].z   =  0.134393

    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15543+4158SE*'))].ra  = '15:56:05.4'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15543+4158SE*'))].dec = '+41:49:23'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15543+4158SE*'))].z   =  0.134893

    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15569+2807W*'))].ra  = '15:58:58.9'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15569+2807W*'))].dec = '+27:59:13'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15569+2807W*'))].z   =  0.052336

    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15569+2807E*'))].ra  = '15:58:58.9'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15569+2807E*'))].dec = '+27:59:13'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS15569+2807E*'))].z   =  0.051736 

    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS21549-1206NW*'))].ra  = '21:57:35.4'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS21549-1206NW*'))].dec = '-11:51:47'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS21549-1206NW*'))].z   =  0.050368 

    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS21549-1206SE*'))].ra  = '21:57:35.4'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS21549-1206SE*'))].dec = '-11:51:47'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS21549-1206SE*'))].z   =  0.050935

    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS22204-0214NW*'))].ra  = '22:22:59.2'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS22204-0214NW*'))].dec = '-01:59:17'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS22204-0214NW*'))].z   =  0.139463

    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS22204-0214SE*'))].ra  = '22:22:59.2'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS22204-0214SE*'))].dec = '-01:59:17'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS22204-0214SE*'))].z   =  0.139563 

    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS22279-1112NW*'))].ra  = '22:30:34.4'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS22279-1112NW*'))].dec = '-10:57:29'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS22279-1112NW*'))].z   =  0.086893

    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS22279-1112SE*'))].ra  = '22:30:34.4'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS22279-1112SE*'))].dec = '-10:57:29'
    v95_ned[where(strmatch(v95_ned.galaxy,'*IRAS22279-1112SE*'))].z   =  0.087660
    
; read the classifications

    readcol, 'veilleux95_class.dat', galaxy, class, format='A,A', $
      comment='#', /silent
    ngalaxy = n_elements(galaxy)

;   niceprint, galaxy, v95_ned.galaxy

    class = strcompress(class,/remove)
    
    v95_class = {class:  ''}
    v95_class = replicate(v95_class,ngalaxy)

    for i = 0L, ngalaxy-1L do begin

       if strmatch(class[i],'S2*') or strmatch(class[i],'S1*') then $
         v95_class[i].class = strmid(class[i],0,2) else $
         v95_class[i].class = strmid(class[i],0,1)

    endfor

    v95 = struct_addtags(struct_trimtags(v95_ned,select=['GALAXY',$
      'NEDGALAXY','RA','DEC','Z']),v95_class)

;   struct_print, v95
    
; ---------------------------------------------------------------------------
; 1999 sample
; ---------------------------------------------------------------------------

; parse and read the NED results and the Vizier classification table 2

    parse_ned_literature, 'veilleux99.ned', outfile='veilleux99_ned.fits'
    v99_ned = mrdfits('veilleux99_ned.fits.gz',1,/silent)

    v99_class = mrdfits('veilleux99_class.fits',1,/silent)
    ngalaxy = n_elements(v99_class)
    
;   niceprint, v99_class.x_raj2000, v99_ned.ra, v99_class.x_dej2000, $
;     v99_ned.dec, v99_class.name, v99_ned.galaxy

    v99 = {galaxy: '', nedgalaxy: '', ra: '', dec: '', z: -999.0, class: ''}
    v99 = replicate(v99,ngalaxy)

    v99.galaxy = 'IRAS'+v99_class.name
    v99.nedgalaxy = v99_ned.galaxy
    v99.ra = v99_ned.ra
    v99.dec = v99_ned.dec
    v99.z = v99_ned.z
    v99.class = v99_class.class

    unknown = where(strcompress(v99.class,/remove) eq '',nunknown)
    if nunknown ne 0L then v99[unknown].class = 'U'
    
; ---------------------------------------------------------------------------
; concatenate everything and write out    
; ---------------------------------------------------------------------------
    
    vout = struct_append(v95,v99)
    srt = sort(im_hms2dec(vout.ra))
    vout = vout[srt]
    
    mwrfits, vout, '99veilleux.fits', /create
    spawn, ['gzip -f 99veilleux.fits'], /sh

return
end    
