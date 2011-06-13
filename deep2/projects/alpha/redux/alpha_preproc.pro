pro alpha_preproc, basepath
; jm09jan19nyu - pre-process the data; reject cosmic rays and use a
; running median to get rid of some electronic noise evident in the
; raw data that is not dealt with properly in the pipeline reductions;
; called by ALPHA_REDUCE_ALL

    if (n_elements(basepath) eq 0L) then begin
       splog, 'Please specify BASEPATH'
       return
    endif
    
    if (file_test(basepath+'preproc',/dir) eq 0L) or $
      (file_test(basepath+'rawdata',/dir) eq 0L) then $
        message, 'Go forth and make the requisite directories'

; assumes [r,b] prefix for RAWLIST and adopts a [p] prefix for the
; output PROCLIST

    rawlist = file_search(basepath+'rawdata/r????.fits.gz',count=nraw)    
    proclist = basepath+'preproc/p'+strmid(file_basename(rawlist),1)
    if (nraw eq 0L) then begin
       splog, 'No files found!'
       stop
    endif

; forage the headers to identify the science exposures; don't
; process standard stars (which we toss out by cutting on exposure
; time); at least one object was mislabeled as a 'COMP' so don't
; cut on EXPTYPE
    info = iforage(rawlist)
;   struct_print, struct_trimtags(info,sel=['file','exptype','object','exptime'])
    obj = where((info.exptime gt 500.0) and strtrim(info.exptype,2) eq 'Object',nobj)
    struct_print, struct_trimtags(info[obj],sel=['file','exptype','object',$
      'exptime','enoise','egain'])

; CR-rejection parameters; they are not optimized, but based on
; Dennis' experiments; MEDWIDTH is the width of the
; median-smoothing filter 
    sigclip = 10.0 
    objlim = 0.01
    sigfrac = 0.01
    medwidth = 301L

;   for ii = 18L, nobj-1L do begin
    for ii = 0L, nobj-1L do begin
; read the image and build a basic inverse variance map
       splog, 'Reading '+rawlist[obj[ii]]
       img = float(readfits(rawlist[obj[ii]],hdr))
       sz = size(img,/dim)

; median-filter the image to get rid of the electronic noice
       splog, 'Median-filtering'
       fimg = img*0.0
       for jj = 0L, sz[1]-1L do for kk = 0L, sz[0]-1L do fimg[kk,jj] = img[kk,jj] - $
         median(img[(kk-medwidth/2L)>0L:(kk+medwidth/2L)<(sz[0]-1L),jj])
;      kk = 1000L & jj = 2000L 
;      djs_plot, img[*,jj], ysty=3 & djs_oplot, fimg[*,jj], color='red'
;      mwrfits, fimg, basepath+'rawdata/junk1.fits', /create
       
; now over-scan subtract (in order to get the statistics right in
; LA_COSMIC) and then reject cosmic rays; note that we pass LA_COSMIC
; the trimmed image (removing the overscan regions) and then chop the
; image into 1024 blocks to speed up the code 
       splog, 'Overscan-subtracting'
       mike_suboscan, fimg, hdr, ovimg, 1, 1, imtype, nobiasrow=nobiasrow, $
         debug=0, redarc=redarc, svbad=svbad, silent=silent
       invvar = vmap_init(ovimg,rdnoise=info[obj[ii]].enoise,gain=info[obj[ii]].egain)

       splog, 'Rejecting cosmic-rays'
       ovimg_trim = ovimg[0L:sz[0]-128L-1L,0L:sz[1]-128L-1L]
       ila_cosmic, ovimg_trim, outlist=outimg_trim, masklist=outmask, sigclip=sigclip, $
         objlim=objlim, sigfrac=sigfrac, /zeroindexed, readn=info[obj[ii]].enoise, $
         gain=info[obj[ii]].egain, /isbig, blocksize=1024
       outimg = ovimg
       outimg[0L:sz[0]-128L-1L,0L:sz[1]-128L-1L] = outimg_trim

; simple sky subtraction
;      sky, ovimg, skymode, skysig, readnoise=info[obj[ii]].enoise, /silent
;      imnosky = ovimg - skymode
;      reject_cr, imnosky, invvar, [0.496,0.246], rejects, $
;        nrejects=nrejects, c2fudge=c2fudge, niter=10L
;      splog, 'Identified '+string(nrejects,format='(I0)')+' cosmic rays'
;      invvar[rejects] = 0.0
;      outimg = float(img)
;      outimg = djs_maskinterp(img,(invvar le 0.0),iaxis=0,/const)
;      xatv, outimg, /bl

       splog, 'Writing '+proclist[obj[ii]]
       mwrfits, outimg, proclist[obj[ii]], hdr, /create

    endfor
    
return
end
    
