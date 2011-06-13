pro alpha_reduce_night, datapath, setup=setup, side=side, clobber=clobber, $
  flat=flat, arc=arc, slitflat=slitflat, proc=proc, emlines=emlines, $
  dotrace=dotrace, skysub=skysub, extract=extract, calibrate=calibrate, $
  coadd=coadd, dostandards=dostandards, makesens=makesens, linlist=linlist, $
  nycoeff=nycoeff, nocoeff=nocoeff, sigrej_2darc=sigrej_2darc
; jm09jan06nyu - wrapper to reduce one night of alpha data; very
; modular 
    
    sz_img = [2048L,4096L] ; stupidly hard-coded

    if (n_elements(datapath) eq 0L) then datapath = './'
    if (n_elements(setup) eq 0L) then setup = 1 ; default setup
    if (n_elements(side) eq 0L) then side = 2   ; red
    
    if (file_test(datapath+'mike.fits',/regular) eq 0L) then begin
       splog, 'MIKE structure not found'
       return
    endif
    mike = mike_ar(datapath+'mike.fits')

    pushd, datapath
    
; ##################################################
; grab all the objects, and cross-match MIKE against the parent
; catalog of observed objects

    objindx = where(mike.setup EQ setup AND mike.type EQ 'OBJ',nobjindx)
    objindx = objindx[uniq(mike[objindx].obj_id,sort(mike[objindx].obj_id))]
;   objindx = objindx[5]
    obj = mike[objindx].obj_id
    nobj = n_elements(obj)

    info = alpha_match_observed_catalog(mike[objindx],$
      setup=setup,side=side,/verbose)

; identify the standards and build an INFO structure    
    stdindx = where((mike.setup EQ setup AND mike.type EQ 'STD'),nstd)
    if (nstd gt 0L) then begin
       std = mike[stdindx].obj_id
       stdinfo = struct_trimtags(mike[stdindx],select=['obj_id',$
         'obj','ra','dec','exp','am','img_root'])
       morestdinfo = alpha_find_standard(stdinfo)
       stdinfo = struct_addtags(temporary(stdinfo),morestdinfo)
       struct_print, stdinfo
    endif else splog, 'No standard stars observed'

;   niceprint, info.img_root, mike[objindx].img_root
;   niceprint, mike[objindx].frame, mike[objindx].obj, mike[objindx].obj_id

; ##################################################
; process the flats and optionally verify

    if keyword_set(flat) then begin
; (1) generate a stacked, normalized milky flat; (2) combine the trace
; flats for order and slit tracing; (3) generate a smooth model of the
; order curvature using the trace flat previously created
       if keyword_set(stepbystep) then begin
          mike_mkmflat, mike, setup, side, smooth=0, clobber=clobber
          mike_mktflat, mike, setup, side, clobber=clobber
          mike_edgeflat, mike, setup, side, inter=inter, chk=chk, clobber=clobber
       endif else begin ; do it all in one step
          mike_allflat, mike, setup, side, clobber=clobber
       endelse
;      spawn, 'gv QA/Flats01/qa_trcflt_01R.ps.gz &'
;      xatv, 'Flats/Flat_R_01_M.fits' 
;      xatv, 'Flats/Flat_R_01_T.fits' 
;      xatv, 'Flats/Flat_R_01_T.fits' 
;      mike_chktrcflat, mike, setup, side, /nostop, /fit
    endif    
    
; ##################################################
; process the arcs

    if keyword_set(arc) then begin
       if keyword_set(stepbystep) then begin
; pre-process the arcs and do the tracing (do not fit and do not
; generate the 2D image)
;         guess_arc = mike_getfil('guess_arc',side=side,/name,chkfil=chkfil,sz=[2048,4096])
;         mike_allarc, mike, setup, side, fits='mike.fits', $
;           clobber=clobber, chk=chk, /nowav, /noimg
       endif else begin
;         arcs = where(mike.type EQ 'ARC' AND mike.flg_anly NE 0 $
;           AND mike.setup EQ setup AND mike.side EQ side)
;         exp = where(strmatch(mike[arcs].img_root,'*0053*'))
;         sigrej_2darc = 2.5
          mike_allarc, mike, setup, side, fits=datapath+'mike.fits', $
            clobber=clobber, chk=chk, linlist=linlist, nycoeff=nycoeff, $
            nocoeff=nocoeff, sig2drej=sigrej_2darc;, exp=exp
       endelse
    endif

; ##################################################
; build the slit profile
    if keyword_set(slitflat) then begin
       mike_slitflat, mike, setup, side, clobber=clobber, chk=chk
    endif

; ##################################################
; reduce the standards
    if keyword_set(dostandards) then begin
       if keyword_set(proc) then begin
          for kk = 0L, nstd-1L do mike_proc, mike, setup=setup, $
            obj=std[kk], side=side, clobber=clobber, /std
       endif
       if keyword_set(dotrace) then begin
          for kk = 0L, nstd-1L do mike_fntobj, mike, setup, $
            std[kk], side, chk=chk, /std
       endif
       if keyword_set(skysub) then begin
          for kk = 0L, nstd-1L do mike_skysub, mike, setup, $
            stdindx[kk], side, chk=chk, /std ; note STDINDX!
       endif
       if keyword_set(extract) then begin
          for kk = 0L, nstd-1L do mike_box, mike, setup, $
            std[kk], side, chk=chk, ochk=ochk, reschk=reschk, /std, $
            /skipskysub, nohelio=0, novac=0
       endif
    endif

; ##################################################
; build the sensitivity function
    if keyword_set(makesens) then begin
       for kk = 0L, nstd-1L do $
         mike_calibstd, mike, stdindx[kk], $ ; note STDINDX!
         esofil=strtrim(stdinfo[kk].esofil,2) 
    endif 

; ##################################################
; initialize the inverse variance map, overscan-subtract, and divide
; by the flat-field; the output is written in /Final
    if (keyword_set(dostandards) eq 0) and keyword_set(proc) then begin
       for jj = 0L, nobj-1L do mike_proc, mike, setup=setup, $
         obj=obj[jj], side=side, clobber=clobber
    endif
    
; ##################################################
; identify the orders of interest (i.e., those containing the emission
; lines)
    if (keyword_set(dostandards) eq 0) and keyword_set(emlines) then begin
       qafile = datapath+'QA/qa_trace_emlines.ps'
       alpha_trace_emlines, info, mike[objindx], qafile=qafile
    endif

; ##################################################
; find and trace the object in the slit; OBJAPER is the aperture to
; mask for sky subtraction (26 unbinned pixels when /STD, and 20
; pixels otherwise); FWIDTH is the fraction of the slit width to use
; when tracing (default 0.25 = 1/4 slit width)
    if (keyword_set(dostandards) eq 0) and keyword_set(dotrace) then begin
       for jj = 0L, nobj-1L do begin
          mike_fntobj, mike, setup, obj[jj], side, $
            objaper=objaper, fwidth=fwidth, chk=chk
          objstrfil = mike_getfil('obj_fil',setup,$
            SUBFIL='Extract/Obj_'+mike[objindx[jj]].img_root,/name)
          objstr = mike_getfil('obj_fil',setup,$
            SUBFIL='Extract/Obj_'+mike[objindx[jj]].img_root)
; read "my" object structure file written out by /EMLINES and
; overwrite the default trace structure; this modified trace does not
; get used by alpha_mike_skysub because we build and pass an object
; mask; note that I have to force the trace to be in the center of the
; order for orders that we don't care about because of a bug in
; x_fntobj, which would otherwise crash the sky-subtraction routine 
          ordr_str = mike_getfil('ordr_str',setup,side=side)
          objstr.trace[0L:sz_img[1]-1L] = (ordr_str.rhedg-ordr_str.lhedg)/2.0+$
            ordr_str.lhedg
; test (ut080415/r0036.fits)
;         mm=14 & djs_plot, objstr[mm].trace[0:4095], ysty=3
;         djs_oplot, ordr_str[mm].lhedg, color='red'
;         djs_oplot, ordr_str[mm].rhedg, color='blue'
; now modify the orders we care about          
          myobjstrfil = repstr(objstrfil,'Obj_','myObj_')+'.gz'
          myobjstr = mrdfits(myobjstrfil,1,/silent)
          for kk = 0L, n_elements(myobjstr)-1L do begin
             indx = where(myobjstr[kk].ordr eq objstr.order)
             objstr[indx].trace[0L:sz_img[1]-1L] = myobjstr[kk].trace
             objstr[indx].aper = myobjstr[kk].aper
          endfor
          mwrfits, objstr, objstrfil, /create ; overwrite!
       endfor 
    endif

; ##################################################
; sky-subtract
    if (not keyword_set(dostandards)) and keyword_set(skysub) then begin
       for jj = 0L, nobj-1L do begin
; read and build the emission-line mask
          imgfil = mike_getfil('fin_fil',setup,/name,$
            SUBFIL=mike[objindx[jj]].img_root)
          maskfile = repstr(imgfil,'.fits','_emlines.fits.gz')
          splog, 'Reading '+maskfile
          linemask = mrdfits(maskfile,0,/silent)
          linefit = mrdfits(maskfile,1,/silent)
          objfil = mike_getfil('obj_fil',setup,$
            SUBFIL=mike[objindx[jj]].img_root,/name)
          myobjfil = repstr(objfil,'Obj_','myObj_')+'.gz'
          splog, 'Reading '+objfil
          objstr = xmrdfits(objfil,1,STRUCTYP='mikeobjstrct',/silent)
          splog, 'Reading '+myobjfil
          myobjstr = mrdfits(myobjfil,1,/silent)
          these_ordrs = myobjstr.ordr
          ordr_str = mike_getfil('ordr_str',setup,side=side)
; in order for x_echskysub to work properly, we also have to mask out
; object pixels in every order, not just the emission lines; since we
; typically don't detect continuum, mask out some fraction of the slit
; in all the orders, but, centered on the object; this code is all in
; x_echskysub
          obj_temp = ordr_str
          frac = 0.5 ; 0.2
          aper = cmreplicate([frac,frac],n_elements(ordr_str)) ; +/-20% of the slit width
          slit_length = obj_temp.rhedg - obj_temp.lhedg
          obj_temp.lhedg = (objstr.trace[0:sz_img[1]-1] - $
            (aper[0,*] ## replicate(1,sz_img[1]))*slit_length/2.0) > ordr_str.lhedg
          obj_temp.rhedg = (objstr.trace[0:sz_img[1]-1] + $
            (aper[1,*] ## replicate(1,sz_img[1]))*slit_length/2.0) < ordr_str.rhedg
          mask_temp_obj = x_ordermask(sz_img[0],sz_img[1],obj_temp,trim=0.0)
; now build the emission-line mask; two possibilities: (1) mask out
; pixels that are 0.1% above the background or (2) mask out the full
; width of the slit at the position of each line to be conservative in
; the sky-subtraction 
          objmask1 = (linemask gt 0.001) ; simple method
;;        ordermask = x_ordermask(sz_img[0],sz_img[1],ordr_str,trim=0.0)
;;        objmask1 = fix(linemask*0.0)
;;        for gg = 0L, n_elements(linefit)-1L do begin
;;           if linefit[gg].goodfit then begin
;;              indx = where(ordr_str.order eq linefit[gg].ordr)
;;              y1 = floor(linefit[gg].ypos-linefit[gg].height/2.0)
;;              y2 = ceil(linefit[gg].ypos+linefit[gg].height/2.0)
;;              yy = transpose(indgen(sz_img[1]) # (intarr(sz_img[0])+1))
;;              these = where((yy ge y1) and (yy le y2) and (ordermask eq ordr_str[indx].order))
;;              objmask1[these] = objmask1[these] or 1
;;              atv, (objmask1 gt 0) or (mask_temp_obj gt 0), /bl
;;           endif
;;        endfor
; mask the Gaussian and the central pixels
;         objmask = (objmask1 gt 0) or (mask_temp_obj gt 0)
; just mask the central continuum pixels
          objmask = (mask_temp_obj gt 0)
;         atv, objmask, /bl
; now sky-subtract!
          alpha_mike_skysub, mike, setup, obj[jj], side, chk=chk, $
            use_objmask=objmask;, $
;           sigrej=3.0, nord=3L, evn=1000L, ordr=[69,62], debug=0
;           ordr=[these_ordrs[0],these_ordrs[n_elements(these_ordrs)-1]]
       endfor
    endif

; ##################################################
; extract 1D spectra
    if (keyword_set(dostandards) eq 0) and keyword_set(extract) then begin
       for jj = 0L, nobj-1L do begin
          splog, 'Performing boxcar extraction'
          objfil = mike_getfil('obj_fil',setup,$
            SUBFIL=mike[objindx[jj]].img_root,/name)
          myobjfil = repstr(objfil,'Obj_','myObj_')+'.gz'
          splog, 'Reading '+myobjfil
          myobjstr = mrdfits(myobjfil,1,/silent)
          these_ordrs = myobjstr.ordr
          nthese_ordrs = n_elements(these_ordrs)
          base_aper = myobjstr[0].aper ; aperture width figured out in /EMLINES
; do not mask cosmic rays as the algorithm in x_extechopt
          mike_box, mike, setup, obj[jj], side, chk=chk, ochk=ochk, reschk=reschk, $
            /boxonly, ordrs=these_ordrs, nohelio=0, novac=0, base_aper=base_aper, $
            /nocrmask ; do not mask cosmic rays (too simplistic!)
; this code will coadd multiple exposures of the same object, but we
; probably want to do this ourselves in post-processing          
       endfor
    endif

; ##################################################
; flux-calibrate the objects and standard stars
    if keyword_set(calibrate) then begin
       if (nstd gt 0L) then begin
          fluxfil = mike_getfil('sens_fil',setup,side=side,$
            subfil=mike[stdindx[0]].img_root,/name)
          for kk = 0L, nstd-1L do mike_flux, mike, setup, $
            stdindx[kk], side, fluxfil=fluxfil, /std
          for jj = 0L, nobj-1L do mike_flux, mike, setup, $
            obj[jj], side, fluxfil=fluxfil, /boxcar
       endif else begin
          splog, 'Unable to flux-calibrate! No standards observed'
          return
       endelse
    endif

; ##################################################
; combine multiple exposures of the same object, which must be run on
; even a single object
;   if keyword_set(combine) then begin
;      for jj = 0L, nobj-1L do begin
;         mike_combspec, mike, setup, obj[jj], side, /useboxflux
;      endfor
;   endif
    
; ##################################################
; read the extracted 1D spectra and stitch the orders together;  also
; coadd the 1D spectra of common objects  
    if keyword_set(coadd) then begin
; unfluxed spectra
       if (nobj gt 0L) then begin
          qafile = datapath+'spec1d/qa_spec1d.ps'
          alpha_coadd_spec1d, mike[objindx], info, side=side, $
            datapath=datapath, qafile=qafile, fluxed=0
;; fluxed spectra
;          if (nstd gt 0L) then begin
;             qafile = datapath+'spec1d/qa_spec1d_fluxed.ps'
;             alpha_coadd_spec1d, mike[objindx], info, side=side, $
;               datapath=datapath, qafile=qafile, fluxed=1
;          endif
       endif
;; standards - do just the fluxed spectra
;       if (nstd gt 0L) then begin
;          qafile = datapath+'spec1d/qa_spec1d_std.ps'
;          alpha_coadd_spec1d, mike[stdindx], stdinfo, side=side, $
;            datapath=datapath, qafile=qafile, fluxed=1, /std
;       endif
    endif 
    popd

return
end
