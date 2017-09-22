pro parse_dr3_zeropoints
; jm16mar21siena 

; In preparation for the tractor runs for DR3, I have computed the zero points
; for the following data:
;   (1) All DECaLS data (i.e., obtained by us under propid 2014B-0404) - 8,224
;     frames
;   (2) All NonDECaLS data (i.e., obtained by others) with existing CP version >
;     3 reductions - 10,403 frames  

    dr2 = mrdfits(getenv('LEGACYPIPE_DIR')+'/decals-ccds-annotated.fits',1)
    pt1 = mrdfits(getenv('IM_PROJECTS_DIR')+'/decals/zeropoints/decals-zpt-dr3pt1.fits',1)

    upropid = dr2[uniq(dr2.propid,sort(dr2.propid))].propid
    for ii = 0, n_elements(upropid)-1 do print, upropid[ii], long(total(dr2.propid eq upropid[ii]))
; 2012B-0001       44225
; 2012B-0003        3111
; 2013A-0351        3172
; 2013A-0360       21655
; 2013A-0386          61
; 2013A-0529         915
; 2013A-0611        1037
; 2013A-0613        2440
; 2013A-0614        5856
; 2013A-0704         183
; 2013A-0716        1342
; 2013A-0717        2501
; 2013A-0719        5185
; 2013A-0723          61
; 2013A-0737        2989
; 2013A-0741       80703
; 2013A-9999        2940
; 2013B-0502        3111
; 2013B-0613         122
; 2013B-0616         549
; 2014A-0191         540
; 2014A-0239        1380
; 2014A-0255         300
; 2014A-0339       14880
; 2014A-0429        3300
; 2014B-0404      311820

    keep = where(strtrim(dr2.propid,2) eq '2014B-0404')
    dr2 = dr2[keep]

; this should just be 2014B-0404, but it's not    
    upropid = pt1[uniq(pt1.propid,sort(pt1.propid))].propid
    for ii = 0, n_elements(upropid)-1 do print, upropid[ii], long(total(pt1.propid eq upropid[ii]))
; 2013B-0440         780
; 2014B-0404      491160
; 2015B-0187        1500

    keep = where(strtrim(pt1.propid,2) eq '2014B-0404')
    pt1 = pt1[keep]

; figure out which files are in DR3 but *not* DR2
    dr2file = repstr(strtrim(file_basename(dr2.image_filename),2),'.fz','')
    udr2 = uniq(dr2file,sort(dr2file))
    dr2file = dr2file[udr2]
    
    pt1file = strtrim(pt1.filename,2)
    upt1 = uniq(pt1file,sort(pt1file))
    pt1file = pt1file[upt1]

    djs_plot, dr2[udr2].ra, dr2[udr2].dec, psym=4, symsize=0.2, xsty=3, ysty=3, $
      xrange=[-10,360], xtitle='RA', ytitle='Dec'
    djs_oplot, pt1[upt1].ra, pt1[upt1].dec, psym=4, symsize=0.2, color='orange'

    help, dr2file, pt1file, cmset_op(dr2file,'and',pt1file), $
      cmset_op(dr2file,'and',/not2,pt1file), cmset_op(dr2file,'and',/not1,pt1file)
    
    


return
end
