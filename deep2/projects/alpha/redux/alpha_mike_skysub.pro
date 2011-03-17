;+ 
; NAME:
; mike_skysub   
;     Version 1.2
;
; PURPOSE:
;    Sky Subtract image and add sky subtracted frame to the final
;    image (e.g. Final/f_mike0020.fits).  The program also outputs
;    information related to the sky fit in the Sky directory.  The
;    user has control over which order to subtract (ORDR=), the
;    wavelength image to use (AIMG=).  The program preforms the
;    following steps:
;
;    1.  Offset the Order structure by the appropriate shift
;    2.  Subtract scattered light (mike_fitgap)
;     -- LOOP ON ORDERS --
;    3.  After masking the object, the code finds the sky pixels and
;    normalizes them by the slit profile.
;    4.  It bspline fit to the sky
;    5.  It assesses the fit and adds in breakpoints in a somewhat
;    smart way.
;    6.  It performs a second bspline fit
;    7.  The coefficients of the fit are saved and the image is
;    subtracted 
;    8.  The code compares the extracted sky spectrum against the UVES
;    sky spectrum to guess at a shift.  This is saved to the Object
;    structure. 
;     -- END OF LOOP
;    9.  The final sky subtracted image is written to the fits file
;    (extension #2)
;
; CALLING SEQUENCE:
;   
;  mike_skysub, mike, setup, obj_id, side, [exp], /CHK, /STD, ORDR=,
;  AIMG=, /USEOLD, /FCHK
;
; INPUTS:
;   mike    -  MIKE structure
;   setup   -  Setup ID
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   side    -  Blue (1) or Red (2) side
;   [exp]   -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;   A sky subtracted image added as extension #2 to the fits file.
;   A file in the directory Sky/ describing the fit
;
;
; OPTIONAL KEYWORDS:
;  /STD     - Sky subtract standard star
;  /CHK     - Show steps along the way
;  AIMG=    - Use alternate Wavelength map for sky subtraction
;             (string filename)
;  ORDR=    - Sky subtract a single order only (e.g.  [5L,5L])
;  /CLOBBER - Overwrite any previos sky image
;  /USEOLD  - Overwrite only the new orders into the old sky sub image
;  /FCHK    - Show the final sky subtracted image
;  /DOALL   - ??
;
; OPTIONAL OUTPUTS:
;  SCATTERED_LIGHT= -- Scattered light image.  It may be useful to
;                      save and use this if you are subtracting individual
;                      orders
;
; COMMENTS:
;
; EXAMPLES:
;   mike_skysub, mike, 1L, 1L, 1L, /FCHK
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_fitgap
;  mike_ordermask
;  bspline_iterfit
;  mike_getfil
;  mike_skyshift
;  mike_skysub_write
;
; REVISION HISTORY:
;   16-May-2003 Written by JXP
;   ??-2003     Modified extensively by SMB
;   26-Jun-2004 Update arc shifting
;-
;------------------------------------------------------------------------------

pro mike_skysub_write, mm, bset, skyfil, flg_ordr
  ;; Write Sky info (bset or spectrum)
  if flg_ordr EQ 0 then begin  ;; ENTIRE SPECTRUM
      if mm EQ 0 then mwrfits, bset, skyfil, /create, /silent $
      else mwrfits, bset, skyfil, /silent 
  endif else begin  ;; SELECT ORDERS
      a = findfile(skyfil+'*',count=na)
      if na EQ 0 then stop
      spawn, 'cp '+skyfil+' tmp_sky.fits'
      case mm of 
          0: begin
              mwrfits, bset, skyfil, /create,/silent 
              for i=1L,nordr-1 do begin
                  tmp_sky = xmrdfits('tmp_sky.fits', i+1, /silent)
                  mwrfits, tmp_sky, skyfil, /silent 
              endfor
          end
          else: begin
              tmp_sky = xmrdfits('tmp_sky.fits', 1, /silent)
              mwrfits, tmp_sky, skyfil, /silent, /create
              ;; MIDDLE
              for i=1L,mm-1 do begin
                  tmp_sky = xmrdfits('tmp_sky.fits', i+1, /silent)
                  mwrfits, tmp_sky, skyfil, /silent 
              endfor
              ;; NEW
              mwrfits, bset, skyfil, /silent 
              ;; REST
              for i=mm+1L,nordr-1 do begin
                  tmp_sky = xmrdfits('tmp_sky.fits', i+1, /silent)
                  mwrfits, tmp_sky, skyfil, /silent 
              endfor
          endelse
      endcase
      spawn, '\rm -f tmp_sky.fits'
  endelse

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro alpha_mike_skysub, mike, setup, obj_id, side, exp, CHK=chk, STD=std, $
                 ORDR=ordr, AIMG=aimg, USEOLD=useold, FCHK=fchk, $
                 DEBUG=debug, DOALL=doall, SCATTERED_LIGHT=iscattered_light, $
  use_objmask=use_objmask, _extra=extra ; jm09jan01nyu

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'mike_skysub, mike, setup, obj_id, side, [exp], /CHK, /STD, '
      print, '      ORDR=, AIMG=, /USEOLD, /FCHK, /DOALL [v1.2]'
      return
  endif 
  
  common xcommon_color, r_vector, g_vector, b_vector
  xicmm_colors

;  Optional Keywords

;  Setup and setup
  case side of 
      1: ostep = 1
      2: ostep = -1
      else: stop
  endcase

;  Find all relevant obj
  if not keyword_set( STD ) then begin
      indx = where(mike.flg_anly NE 0 AND mike.side EQ side[0] AND $
                   mike.obj_id EQ obj_id AND mike.setup EQ setup AND $
                   (strtrim(mike.type,2) EQ 'OBJ' OR $
                    strtrim(mike.type,2) EQ 'STD'), nindx)
      if nindx EQ 0 then begin
          print, 'mike_skysub: No images to sky subtract!', obj_id
          return
      endif
  endif else begin
      indx = obj_id[0]
      nindx = 1L
  endelse
      
;  Exposures
  if n_elements(exp) EQ 0 then exp = lindgen(nindx)

  skylinefile = getenv('MIKE_DIR')+'/pro/Skysub/UVES_SKY/uves.linelist'
  if (findfile(skylinefile))[0] NE '' then begin
     print, 'MIKE_SKYSUB: Reading ', +strcompress(skylinefile)
     readcol, skylinefile[0], uvesindx, uveswave, uvesint, uvesfwhm, uvesflux,$
           format='(i,d,f,f,f)'
     uveslog = alog10(uveswave)
     nuves= n_elements(uveslog)
  endif

;  Read in order structure
  ordr_str = mike_getfil('ordr_str', setup, side=side)
  sv_ostr = ordr_str

  if (where(tag_names(ordr_str) EQ 'PROFILE0'))[0] NE -1 then begin
     if total(ordr_str.profile0) EQ 0 then $
       message, 'Need to run mike_slitflat to obtain slit profile'
  endif
  nordr = n_elements(ordr_str)

; Orders
  if not keyword_set( ORDR ) then begin  ;; Orders to sky subtract
      ordr = [ordr_str[0].order,ordr_str[nordr-1].order]
      flg_ordr = 0
  endif else begin
      if keyword_set( USEOLD ) then flg_ordr = 1 else flg_ordr = 0
      if n_elements(ordr) NE 2 then begin
          print, 'mike_skysub: ORDR must be a 2 element array'
          return
      endif
  endelse

  clr = getcolor(/load)


;  Loop

  for q=0L,n_elements(exp)-1 do begin
    
      thisone = indx[exp[q]]
      ;;  Read Arc
      if not keyword_set(AIMG) then arc_img = mike[thisone].arc_img $
      else arc_img = aimg
      ;; READ
      if x_chkfil(arc_img+'*') EQ 0 then begin
          print, 'mike_skysub: Arc file doesnt exist!', arc_img
          stop
          return
      endif
      print, 'mike_skysub: Reading Wave image file ', arc_img
      img_arc = xmrdfits(arc_img, /silent) 
      sz_arc = size(img_arc, /dimensions)

      ;; xoff
      xyoff = mike[thisone].arc_xyoff

      ;; Offset ordr_str
      ordr_str = sv_ostr
      shft = mike_shifti(xyoff, OSTR=ordr_str)

      ;; Open Obj file
;      objfil = mike[thisone].obj_fil
      objfil = mike_getfil('obj_fil', setup, SUBFIL=mike[thisone].img_root, /name)
      objstr = xmrdfits(objfil, 1, STRUCTYP='mikeobjstrct', /silent)
      if keyword_set( DOALL ) then begin
          a = where((objstr.flg_anly MOD 2) NE 1, na)
          if na NE 0 then objstr[a].flg_anly = objstr[a].flg_anly + 1
      endif

      ;; IMG+VAR Fil 
;      imgfil = strtrim(mike[thisone].img_final,2)
      imgfil = mike_getfil('fin_fil', setup, SUBFIL=mike[thisone].img_root, /name)
      skyfil = mike_getfil('sky_fil', setup, SUBFIL=mike[thisone].img_root, /name)
      scattfil = mike_getfil('scatt_fil', setup, $
                             SUBFIL=mike[thisone].img_root, /name)
      print, 'mike_skysub: Reading Image files... ', imgfil
      img = xmrdfits(imgfil, 0, head, /silent)
      ivar = xmrdfits(imgfil, 1, /silent)
      sz_img = size(img, /dimensions)

;     if (n_elements(linemask) ne 0L) then img = img - linemask ; jm09jan25nyu
      
      ;; Mask
      msk = lonarr(sz_img[0],sz_img[1])
      model = fltarr(sz_img[0],sz_img[1])
      maskimage = mike_ordermask(sz_img[0], sz_img[1], ordr_str, trim=0.)

      ;; Inter-order light fit
      nxbkpt = side EQ 1 ? 4 : 12
      nybkpt = side EQ 1 ? 4 : 8

      qafil = mike_getfil('qa_skysub', setup, SUBFIL=mike[thisone].img_root) ; jm08apr13nyu

      alpha_x_echskysub, img, ivar, img_arc, objstr, ordr_str, img_new, $
        CHK=chk, $
        ORDR=ordr, USEOLD=useold, FCHK=fchk, $
        DEBUG=debug, DOALL=doall, SCATTERED_LIGHT=scattered_light, $
        NXBKPT=nxbkpt, NYBKPT=nybkpt, OSTEP=ostep, SKYFIL=skyfil, $
        COLBIN=mike[thisone].colbin, QAFIL=qafil, MEDSCATT=medscatt, $
        SCATTTRIM=scatttrim, OBJTRIM=objtrim, $
        use_objmask=use_objmask, _extra=extra ; jm09jan01nyu

      ;; Ouptut Obj Struct (with sky arrays)
      mwrfits, objstr, objfil, /create
      spawn, 'gzip -f '+objfil

      ;; CHK
      print, 'mike_skysub: Writing output to: ', imgfil
      if keyword_set( CHK ) or keyword_set( FCHK ) then $
        xatv, img_new, WVIMG=10^img_arc, /block, min=-50., max=50.

      ;; Ouptut New Image
      mwrfits, img, imgfil, head, /create, /silent
      mwrfits, ivar, imgfil, /silent
      mwrfits, img_new, imgfil, /silent

      ;; Scattered light image
      print, 'mike_skysub: Writing scattered image: ', scattfil
      mwrfits, scattered_light, scattfil, /create, /silent

      ;; COMPRESS
      print, 'mike_skysub: Compressing...'
      spawn, 'gzip -f '+imgfil
      spawn, 'gzip -f '+scattfil
  endfor
  
;  DONE
  print, 'mike_skysub: All done! '

  return
end


;      if not keyword_set( ISCATTERED_LIGHT ) then begin
;          print, 'Fitting inter-order light with mike_fitgap', $
;            nxbkpt, ' x ', nybkpt, ' bkpts..', format='(a,i3,a,i2,a)'
;          scattered_light = mike_fitgap(img, ivar, maskimage, $
;                                        nxbkpt=nxbkpt, nybkpt=nybkpt)
;          print, 'mike_skysub: Min, max scattered: ', $
;            min(scattered_light), max(scattered_light), $
;            format='(a,f8.2, f8.2)'
;          if min(scattered_light) LT -50. then begin
;              print, 'mike_skysub: I suspect you have a problem here!, '
;              print, ' If you continue, scattered light will have 0 set as minimum'
;              xatv, scattered_light, /block
;              stop  ;; Problem likely
;              scattered_light = scattered_light > 0
;          endif
;      endif else scattered_light = iscattered_light
;
;      img_sub = img - scattered_light
;;
;      obj_temp = ordr_str
;      slit_length = obj_temp.rhedg - obj_temp.lhedg
;      obj_temp.lhedg = (objstr.trace[0:sz_img[1]-1] - $
;                      (objstr.aper[0] ## replicate(1,sz_img[1])) * $
;                         slit_length/2.)  > ordr_str.lhedg
;      obj_temp.rhedg = (objstr.trace[0:sz_img[1]-1] + $
;                      (objstr.aper[1]  ## replicate(1,sz_img[1])) * $
;                         slit_length/2.)  < ordr_str.rhedg
;
;      mask_temp_obj = mike_ordermask(sz_img[0], sz_img[1], obj_temp, $
;                                      trim=0., /noextra)
;;
;;      Object pixels are order number + 10000
;;
;      maskimage = maskimage + (mask_temp_obj GT 0)*10000L
;
;      ;; Final image
;      ;; if flg_ordr EQ 1 then img_new = xmrdfits(imgfil, 2, /silent) $ 
;      ;; else img_new = fltarr(sz_img[0],sz_img[1])
;;
;      img_new = img_sub
;      usemask_img = maskimage*0
;      fwhm = fltarr(nordr)
;      fwhm_sig = fltarr(nordr)
;
;      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      ;; Loop on Orders
;      for qq=ordr[0],ordr[1],ostep do begin
;          ;;
;          mm = where(ordr_str.order EQ qq, nmm)
;          if nmm EQ 0 then stop else mm = mm[0]
;          if objstr[mm].order NE ordr_str[mm].order then stop
;          if (objstr[mm].flg_anly MOD 2) NE 1 then begin
;              print, 'mike_skysub: Skipping order ', qq
;              continue
;          endif
;          ;; 
;          print, 'mike_skysub: ', qq, FORMAT='(a, i3, $)'
;
;          ;; Order
;          nord = 4L
;
;                  
;          inorder = where(maskimage mod 10000L EQ qq AND img_arc GT 3. $
;                         AND ivar GT 0.)
;          wave_sort = inorder[sort(img_arc[inorder])]
;          xstart = 1.0d*(wave_sort mod sz_img[0])
;          ystart = wave_sort /   sz_img[0]
;          slit_length = ordr_str[mm].rhedg[ystart] - ordr_str[mm].lhedg[ystart]
;;
;;           Ordrcen is just center of flat-field order
;;
;          ordrcen  =  (ordr_str[mm].lhedg[ystart] + $
;                       ordr_str[mm].rhedg[ystart])/2.0d ; Offset
;
;          ywave = x_qckwav(xstart - ordrcen, ystart, ordr_str[mm].arc_m)
;          slit_frac  = frac_order(ordr_str[mm], xstart, ywave)
;                  
;          profile = mike_profile_return(slit_frac, ystart, ordr_str[mm])
;              
;          img_norm = (img_sub)[wave_sort] / $
;                (profile + (profile EQ 0)) * (profile GT 0)
;              
;          img_norm_ivar = ivar[wave_sort] * profile^2 * (profile GT 0)
;          wave_norm =  img_arc[wave_sort]
;              
;          sky_pix = (maskimage[wave_sort] LT 200 AND profile GT 0.3)
;          sky = where(sky_pix)
;          obj_pix = (maskimage[wave_sort] GT 10000L)
;              
;          nsky_pix = total(s;ky_pix)
;          nobj_pix = total(obj_pix)
;             
;
;          evn = (1.2 * nsky_pix / (max(ystart) - min(ystart) + 1)) + 1 
;             
;          if nsky_pix LT 0.2*nobj_pix then $
;            print, 'MIKE_SKYSUB: sky pixel to object pixel ratio is lower than 0.2'
;          
;          if nobj_pix LT 2*sz_img[1] then $
;            print, 'MIKE_SKYSUB: Warning, maybe a short order '+string(nobj_pix)
;          
;          mean_object_flux = djs_avsigclip(img_sub[where(obj_pix)], sigrej=10)
;          bkpt_s = 0    
;          
;;
;;               I don't like how rejection is taking place...
;;               it seems like bspline_iterfit is using once rejected
;;                reject the pixel for all iterations.  Not sure if this
;;                is right.
;;
;          bset_prof = bspline_iterfit(wave_norm[sky], img_norm[sky],everyn=evn,$
;                                      nord=nord, upper=5., lower=5., $
;                                      INVVAR= img_norm_ivar[sky], /groupbadpix,$
;                                      maxrej= 1, bkpt=bkpt_s, maxiter= 15L, $
;                                      OUTMASK=pixmsk, yfit=yfit, /silent)
;              
;          usemask = pixmsk
;          ;; calculate chi^2 in sky pixels
;              
;          img_fit = bspline_valu(wave_norm, bset_prof)
;          img_fit = img_fit+ (yfit[0] - img_fit) * $
;                       (wave_norm LT wave_norm[sky[0]])
;          img_fit = img_fit+ (yfit[nsky_pix-1] - img_fit) * $
;                       (wave_norm GT wave_norm[sky[nsky_pix-1]])
;          img_model = img_fit*profile + scattered_light[wave_sort]
;          chi = (img[wave_sort]-img_model)*  sqrt(ivar[wave_sort])
;          badsky = where(pixmsk EQ 0)
;             
;;
;;   Calculate first moment of chi^2 in slit_frac
;;
;          numer = slit_frac[sky] * ((chi[sky] > (-5.0)) < 5)
;          guess_pix_shift = total(numer)/total(abs(chi) < 5) $
;                * (slit_length[sz_img[1]/2]/2.0)
;              
;          use_shift = (guess_pix_shift > (-2.0)) < 2.0
;          print, ', shift guess: ', use_shift, format='(a,f7.3,$)'
;              
;          bkpt_n = bkpt_s
;              
;          ;; 
;          ;; need a better check for more bkpts
;          ;;
;          min_threshold =  10.0*((median(yfit)>0)+5) > (3.0*mean_object_flux)
;          print, ', Min thresh', min_threshold, format='(a,f8.2)'
;
;          need_more_bkpts = where(chi[sky]^2 GT 16 AND  yfit GT min_threshold )
;              
;          if need_more_bkpts[0] NE -1 then begin 
;              spot = need_more_bkpts * 0
;              
;              for ib = 0, n_elements(need_more_bkpts) - 1 do $
;                spot[ib] = (where(wave_norm[sky[need_more_bkpts[ib]]] GT bkpt_s $
;                                  AND wave_norm[sky[need_more_bkpts[ib]]] LE bkpt_s[1:*]))[0]
;              
;              fill_in_spots = spot[uniq(spot)]
;              last_check = where(fill_in_spots NE -1, nspots)
;              if nspots GT 0 then begin
;                  
;                  fill_in_spots = fill_in_spots[last_check] 
;                  
;                  ;;
;                  ;; Add two bkpts between every set where large chi^2 is found
;                  ;; 
;                  new_bkpts = [bkpt_s[fill_in_spots+1]/2. + $
;                               bkpt_s[fill_in_spots]/2.]
;                                ;new_bkpts = [2*bkpt_s[fill_in_spots+1]/3. + $
;                                ;                  bkpt_s[fill_in_spots]/3., $
;                                ;                  bkpt_s[fill_in_spots+1]/3. + $
;                                ;                  2*bkpt_s[fill_in_spots]/3.]
;                  
;                  bkpt_n = [bkpt_s, new_bkpts]
;                  bkpt_n = bkpt_n[sort(bkpt_n)]
;                  print, 'Added ', n_elements(fill_in_spots),' more bkpoints'
;                 if  n_elements(fill_in_spots) GT 100 then stop
;              endif
;              
;              bset_prof = bspline_iterfit(wave_norm[sky],img_norm[sky],$
;                                          nord=nord, upper=5., lower=5., $
;                                          INVVAR= img_norm_ivar[sky], $
;                                          /groupbadpix, maxrej=1, bkpt=bkpt_n, $
;                                          maxiter= 25L, OUTMASK=pixmsk2, $
;                                          yfit=yfit2, /silent)
;          
;              usemask = pixmsk2
;              ;; chi^2 test in sky pixels
;              
;              
;              img_fit = bspline_valu(wave_norm, bset_prof)
;              img_fit = img_fit+ (yfit2[0] - img_fit) * $
;                       (wave_norm LT wave_norm[sky[0]])
;              img_fit = img_fit+ (yfit2[nsky_pix-1] - img_fit) * $
;                       (wave_norm GT wave_norm[sky[nsky_pix-1]])
;              img_model = img_fit*profile+scattered_light[wave_sort]
;              chi = (img[wave_sort]-img_model)* sqrt(ivar[wave_sort])
;              badsky = where(pixmsk2 EQ 0)
;          endif
;          
;          usemask_img[wave_sort[sky]] = usemask
;          model[wave_sort] = img_model
;          img_new[wave_sort] = img[wave_sort] - img_model
;          
;          ;; Write the fit to skyfil
;          if keyword_set(bset_prof) then begin
;              skystr = { $
;                         fullbkpt: bset_prof.fullbkpt, $
;                         bkmask: bset_prof.bkmask, $
;                         nord: bset_prof.nord, $
;                         coeff: bset_prof.coeff, $
;                         icoeff: bset_prof.icoeff, $
;                         badsky: badsky, order: qq $
;                       }
;              mike_skysub_write, mm, skystr, skyfil, flg_ordr
;          endif
;          ;; Create/Write sky spectrum 
;          ywave = mike_qw(xstart-ordrcen, ystart, ordr_str[mm].arc_m, $
;                          arc_slope=arc_slope, slit_dist=slit_dist)
;;
;;
;;;         Attempt FWHM calculation...
;;         this could be slow
;;
;;        hmm = mike_rectify(img_new, obj_temp[mm].lhedg-1, obj_temp[mm].rhedg+1)
;;;        width = (size(hmm,/dimen))[0]
;;        hmmt = total(hmm,1) 
;;        x = (findgen(width) - width/2) # replicate(1,sz_img[1])
;;        goodt = where(hmmt GT 0, ngood)
;;        if ngood GT 10 then begin
;;          djs_iterstat, (total(hmm * x,1))[goodt]/hmmt[goodt], mean=mn
;;          djs_iterstat, (total(hmm * (x-mn)^2,1))[goodt]/hmmt[goodt], $
;;                 mean=var, sigma=sig
;;          fwhm[mm] = sqrt(var)*2.358
;;          fwhm_sig[mm] = sig/2./2.358/ sqrt(ngood)
;;        endif
;;        print, 'I think the spatial FWHM is ',fwhm[mm] ,$
;;                ' pixels, but do not believe this'
;
;          raw_wave_set = bspline_iterfit(ywave, 1.0d*img_arc[inorder], $
;                                         invvar=1.0d*(img_arc[inorder] GT 3), $
;                                         nbkpt=10, /silent)
;          
;          nrow = (size(img))[2]
;          objstr[mm].sky_wv[0:nrow-1] = bspline_valu(dindgen(nrow),raw_wave_set)
;          objstr[mm].sky[0:nrow-1] = $
;            bspline_valu(objstr[mm].sky_wv[0:nrow-1],bset_prof)
;          
;          if keyword_set(uveslog) then $
;            sky_shift = mike_skyshift(bset_prof, uveslog, uvesflux, $
;                                      objstr[mm].sky_wv[0:nrow-1] ) 
;          
;          
;          objstr[mm].sky_wv[0:nrow-1] = 10^objstr[mm].sky_wv[0:nrow-1] 
;          objstr[mm].skyshift = sky_shift
;          
;
;          
;          if keyword_set( CHK ) then begin
;              xrange = [min(objstr[mm].sky_wv[0:nrow-1])-10, $
;                        max(objstr[mm].sky_wv[0:nrow-1])+10]
;              yrange = [-3., max(objstr[mm].sky[50:nrow-50])]
;              
;              ;;  Top Panel sky fit for this order
;              
;              pos = [0.15,0.5,0.95,0.95]
;              plot, objstr[mm].sky_wv[0:nrow-1], objstr[mm].sky[0:nrow-1], $
;                yrange=yrange, ys=1, xrange=xrange, /xs, pos=pos, $
;                /nodata, ymar=[0,2], xchars=0.0001, ytitle='SkySub fit'
;              oplot, 10^wave_norm[sky], img_norm[sky]*usemask,ps=3
;              oplot, objstr[mm].sky_wv[0:nrow-1], $
;                objstr[mm].sky[0:nrow-1], color=clr.green
;              if need_more_bkpts[0] NE -1 AND keyword_set( new_bkpts ) then $
;                oplot, 10^new_bkpts, new_bkpts*0 + yrange[1]*0.3, ps=1, $
;                    color=clr.red
;              
;              ;;  Bottom Panel sky fit for this order
;              
;              pos = [0.15,0.15,0.95,0.5]
;              djs_plot, objstr[mm].sky_wv[0:nrow-1], objstr[mm].sky[0:nrow-1], $
;                /nodata, yr=[-10, 10], /ys, ymar=[4,0],  ychars=1.5, $
;                xrange=xrange, /xs, pos=pos, /noerase, $
;                ytitle='\chi (\lambda)', xtitle='Wavelength(\AA)'
;              oplot, 10^wave_norm[sky], chi[sky], ps=3
;              
;              if badsky[0] NE -1 then $
;                oplot, 10^wave_norm[sky[badsky]], chi[sky[badsky]], $
;                ps=1, color=clr.green, syms=0.3
;          endif
;
;          ;; DEBUG
;          if keyword_set( DEBUG ) then stop
;
;      endfor
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;      ;; QA
;;      save, objstr, sz_img, ivar, img_new, img_arc, maskimage, filename='tmp.idl'
;;      stop
;;      restore, 'tmp.idl'
;
;      qafil = mike_getfil('qa_skysub', setup, SUBFIL=mike[thisone].img_root)
;      !p.multi=[0,1,1]
;      print, 'mike_skysub: QA file in ', qafil+'.gz'
;      x_psopen, qafil, /maxs
;      clr = getcolor(/load)
;
;      ;; Skyshift
;      plot, objstr.order, objstr.skyshift, xtitle='Order', $
;        ytitle='Skyshift (pix)', charsize=1.5, ps=1, color=clr.black, $
;        background=clr.white, xmargin=[11,2], ymargin=[5,2]
;      oplot, [-99, 1000.], [0.,0.], linestyle=1, color=clr.green
;
;      ;; Sky lines
;      nlin = 0L
;      skylin = dblarr(100,2)
;      skylin[nlin,*] = [3600.0d, 3601.0] & nlin = nlin + 1
;      skylin[nlin,*] = [3800.0d, 3801.0] & nlin = nlin + 1
;      skylin[nlin,*] = [4400.0d, 4401.0] & nlin = nlin + 1
;      skylin[nlin,*] = [5576.9d, 5577.7] & nlin = nlin + 1
;      skylin[nlin,*] = [5895.d, 5897]    & nlin = nlin + 1
;      skylin[nlin,*] = [6299.9, 6300.7]  & nlin = nlin + 1
;      skylin[nlin,*] = [6306.0, 6308.0]  & nlin = nlin + 1
;      skylin[nlin,*] = [6363.1, 6364.5]  & nlin = nlin + 1
;      skylin[nlin,*] = [6947.8, 6950.1]  & nlin = nlin + 1
;      skylin[nlin,*] = [7244.0, 7246.0]  & nlin = nlin + 1
;      skylin[nlin,*] = [7275.6, 7277.3]  & nlin = nlin + 1
;      skylin[nlin,*] = [7711.1, 7714.1]  & nlin = nlin + 1
;
;
;      if side EQ 2 then begin
;          img_mn = -100.
;          img_mx = 100.
;      endif else begin
;          img_mn = -50.
;          img_mx = 50.
;      endelse
;      ntv = 50L
;      tv_image = bytarr(ntv,ntv)
;      gd = where(maskimage GT 0L and img_arc GT 0.)
;      mxwv = max(10^img_arc[gd], min=mnwv)
;
;      isz = round(18L*2./mike[thisone].colbin)
;
;      for qq=0L,nlin-1 do begin
;          if skylin[qq,0] GT mxwv or skylin[qq,1] LT mnwv then continue
;          qqm = qq MOD 3
;          ;; new page
;          if qqm EQ 0L then $
;            plot, [0], [0], /nodata, xrange=[5,10.], yrange=[0,1], $
;            xstyle=5, ystyle=5
;          ;; colors
;          state = { $
;                    ncolors: 0L, $
;                    brightness: 0.5, $
;                    contrast: 0.5 $
;                  }
;          loadct, 0, /silent
;          ncolors = !d.table_size - 9
;          state.ncolors=ncolors
;          r_vector = bytarr(ncolors)
;          g_vector = bytarr(ncolors)
;          b_vector = bytarr(ncolors)
;          ximgd_getct, state, 0, /CLR
;          ;; Invert
;          r_vector = reverse(r_vector)
;          g_vector = reverse(g_vector)
;          b_vector = reverse(b_vector)
;          ximgd_stretchct, state
;
;          ;; 2D image
;          wcen = mean(skylin[qq,*])
;          ;; Find pixel in image nearest the wavelength and in object aperture
;          diff = abs(alog10(wcen)-img_arc) + (maskimage LT 10000L)
;          srt = sort(diff)
;          ;; And neareset the center of the image
;          diff2 = abs( (srt[0:99]/sz_img[0]) - sz_img[1]/2.)
;          srt2 = sort(diff2)
;          yp = round(median(srt[srt2[0:10]]/sz_img[0]))
;          pordr = median( maskimage[srt[srt2[0:10]]] mod 10000L)
;          oordr = where(objstr.order EQ pordr)
;
;          if oordr[0] NE -1 then begin
;            ;; Center on obj
;            xp = objstr[oordr].trace[yp]
;
;            ;; Make image
;            x1 = (xp - isz) > 0L
;            x2 = (xp + isz) < (sz_img[0]-1)
;            y1 = (yp - isz) > 0L
;            y2 = (yp + isz) < (sz_img[1]-1)
;  
;            display_image = bytscl(img_new[x1:x2,y1:y2], $
;                                   min=img_mn, max=img_mx, /nan, $
;                                   top=ncolors-1) + 8B
;            tv_image = congrid(display_image, ntv, ntv)
;            ;; TV
;            tv, tv_image, 0.5+qqm*3.3, 4.25, xsize=3.0, ysize=3.0, /inches
;
;          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;          ;; Sky fit
;          ;; Read in fit
;          bset_prof = xmrdfits(skyfil, oordr[0]+1, /silent)
;
;          ;; Grab sky pixels
;          clr = getcolor(/load)
;          subm = maskimage[x1:x2,y1:y2]
;          suba = img_arc[x1:x2,y1:y2]
;          subi = img[x1:x2,y1:y2]
;          subn = img_new[x1:x2,y1:y2]
;          subv = ivar[x1:x2,y1:y2]
;          ;; Focus on sky line
;          inorder = where(subm mod 10000L EQ pordr AND $
;                          suba GT alog10(skylin[qq,0]) AND $
;                          suba LT alog10(skylin[qq,1]) )
;
;          ;; Make profile
;          sz2 = size(subm, /dimensions)
;          xstart = inorder mod sz2[0]  + x1
;          ystart = inorder /   sz2[0]  + y1
;
;          ordrcen  =  (ordr_str[oordr].lhedg[ystart] + $
;                       ordr_str[oordr].rhedg[ystart])/2.0 
;          slit_length = ordr_str[oordr].rhedg[ystart] - $
;                        ordr_str[oordr].lhedg[ystart]
;
;          slit_frac  = 2.0* (xstart - ordrcen) / slit_length
;          profile = mike_profile_return(slit_frac, ystart, ordr_str[oordr])
;
;          ;; Sky pixels
;          sky_pix = (subm[inorder] LT 200 AND profile GT 0.3)
;          sky = inorder[where(sky_pix, nsky)]
;          
;          ;; Sort
;          srt = sort(suba[sky])
;          swv = suba[sky[srt]]
;          sky_fit = bspline_valu(swv, bset_prof)
;          skyres = subn[sky[srt]]
;          srtres = sort(skyres)
;
;          ;; Max min
;          djs_iterstat, skyres, sigma=sig
;          mxs = max(sqrt(sky_fit))
;          mn = -4.*sig < (-3*mxs)
;          mx = 4.*sig > 3*mxs
;
;          ;; Plot 
;          plot, 10^swv, skyres, color=clr.black, /noerase,$
;            background=clr.white, pos=[0.1+qqm*0.33, 0.1, 0.33*(qqm+1), 0.45], $
;            yrange=[mn,mx], ystyle=1, xticks=2, psym=1, xtitle='Wavelength', $
;            ytitle='Residuals' 
;
;;          mxfit = max(sky_fit)
;          oplot, 10^swv, sqrt(sky_fit), color=clr.blue
;
;          ;; Break points
;          oplot, 10^bset_prof.fullbkpt, fltarr(n_elements(bset_prof.fullbkpt)), $
;            psym=2, color=clr.green, symsize=0.5
;          ;; Rej points
;          rej = where(subv[sky[srt]] LE 0., nrej)
;          if nrej NE 0 then $
;            oplot, 10^swv[rej], skyres[rej], color=clr.red, psym=2
;       endif
;      endfor
;
;      ;; Close
;      x_psclose
;      !p.multi=[0,1,1]
;      spawn, 'gzip -f '+qafil

