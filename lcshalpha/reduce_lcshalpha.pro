;+
; NAME:
;   REDUCE_LCSHALPHA
;
; PURPOSE:
;   Reduce all the KPNO/0.9-meter data from AYCamp/2010.
;;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;   inspect - 
;   mkbias - 
;   mkdflat - 
;   assign_calib - 
;   proc - 
;   clobber - 
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 May 31, UCSD
;-

function get_zpt, night, filter
; pull out the right zeropoint
    zptfile = getenv('LCSHALPHA_DATA')+'2010/kp09m/'+$
      night+'/Std/landolt_zpt_'+filter+'.fits.gz'
    if file_test(zptfile) then begin
       table = mrdfits(zptfile,1,/silent)
       return, table.zpt
    endif else return, -1
end

function zeropoint_func, mtrue, param, minst=minst, $
  airmass=airmass, color=color
; m0 = m + Z - k1*AM - k2*C - k3*C*AM
;
;   m0 = Standard magnitude
;   m  = instrumental magnitude
;   Z  = zero point
;   k1 = extinction coefficient [mag/airmass]
;   AM = airmass
;   k2 = color term coefficient
;   C  = Standard system color of the object
;   k3 = second order extinction coefficient

    return, mtrue - param[0] + param[1]*airmass - $
      param[2]*color - param[3]*color*airmass
    
end

pro lcshalpha_unpack, struct, indir=indir, outdir=outdir, $
  night=night, noweight=noweight
    struct_print, struct_trimtags(struct,$
      select=['img_root','obj','filter'])
    nobj = n_elements(struct)
    for ii = 0, nobj-1 do begin
       thisfile = indir+'f_'+strtrim(struct[ii].img_root,2);+'.gz'
       outfile = outdir+strtrim(struct[ii].target,2)+'_'+$
         strmid(struct[ii].img_root,0,5)+'_'+$
         string(struct[ii].frame,format='(I4.4)')+'_'+$
         strtrim(struct[ii].filter,2)+'.fits';.gz'
       if (file_test(thisfile) eq 0) then begin
          splog, 'File '+thisfile+' not found!'
       endif else begin
          spawn, 'rsync -auv '+thisfile+' '+outfile, /sh
; update the header with the zeropoint, if it exists, and do some
; other necessary evil to the headers
          hdr = headfits(outfile)
          sxdelpar, hdr, 'FILTER1'
          sxdelpar, hdr, 'FILTER2'
          sxdelpar, hdr, 'FILTERS'
;         zpt = get_zpt(night,strtrim(struct[ii].filter,2))
;         if (zpt ne -1) then begin
;            sxaddpar, hdr, 'MAGZERO', zpt, ' magnitude zeropoint', $
;              before='HISTORY'
;         endif
          modfits, outfile, 0, hdr, exten_no=0
          if (keyword_set(noweight) eq 0) then begin
             spawn, 'rsync -auv '+repstr(thisfile,'.fits',$
               '.weight.fits')+' '+repstr(outfile,$
               '.fits','.weight.fits'), /sh
          endif
       endelse
    endfor
return
end    

pro reduce_lcshalpha, night=night, inspect=inspect, mkbias=mkbias, $
  mkdflat=mkdflat, assign_calib=assign_calib, proc=proc, astrom=astrom, $
  zeropoints=zeropoints, makelinks=makelinks, sextractor=sextractor, $
  solve=solve, unpack_projects=unpack_projects, clobber=clobber

    datapath = getenv('LCSHALPHA_DATA')+'/'
    splog, 'Data path '+datapath
    tel = 'kp09m'
    ccd = 'S2KB'

    sexpath = getenv('LCSHALPHA_DIR')+'/data/'
    sexconfig = sexpath+'default.sex'
    sexparam = sexpath+'default.sex.param'
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'

    http = ''
;   http = 'http://sdss.physics.nyu.edu/ioannis/'
    xslsex = http+'sex.xsl'

    pixscale = 0.6 ; pixel scale [arcsec/pixel]
    photaper_arcsec = findgen(15)+1.0         ; [diameter,arcsec]
    photaper_pixel = photaper_arcsec/pixscale ; [diameter,pixel]
    photaper = strjoin(strtrim(string(photaper_pixel,format='(F10.2)'),2),',')

     if (n_elements(night) eq 0) then night = $
      ['20110404']
    nnight = n_elements(night)
    structfile = night+'_imgstrct.fits'

;   targetinfo = rsex(datapath+'jun10_targetlist.txt')

; #########################
; preliminary reductions

; inspect all the files, removing problematic or irrelevant exposures 
    if keyword_set(inspect) then begin 
       for inight = 0, nnight-1 do begin
          pushd, datapath+night[inight]
          mosaic_img_strct, allstruct, ccd=ccd, tel=tel, $
            outfil='all_'+structfile[inight], clobber=clobber
; map the filter wheel number to filter name and update the structure
          for iobj = 0, n_elements(allstruct)-1 do begin
             head = headfits(strtrim(allstruct[iobj].rootpth,2)+$
               strtrim(allstruct[iobj].img_root,2));+'.gz')
             filt1 = strtrim(sxpar(head,'filter1'),2)
             filt2 = strtrim(sxpar(head,'filter2'),2)
             case filt1 of
                '7': begin
                   case filt2 of
                      '0': allstruct[iobj].filter = 'Ha6580'
                      '1': allstruct[iobj].filter = 'Ha6620'
                      '2': allstruct[iobj].filter = 'Ha6660'
                      '3': allstruct[iobj].filter = 'Ha6709'
                      '4': allstruct[iobj].filter = 'Ha6740'
                      '5': allstruct[iobj].filter = 'Hb'
                      '6': allstruct[iobj].filter = 'OIII'
                      '7': allstruct[iobj].filter = 'none'
                      else: message, 'Unknown filter!'
                   endcase 
                end
                '0': allstruct[iobj].filter = 'U'
                '1': allstruct[iobj].filter = 'B'
                '2': allstruct[iobj].filter = 'V'
                '3': allstruct[iobj].filter = 'R'
                '4': allstruct[iobj].filter = 'I'
                '5': allstruct[iobj].filter = 'OII'
                '6': allstruct[iobj].filter = 'Na'
                else: message, 'Unknown filter!'
             endcase 
          endfor 
; throw out crap objects
          case night[inight] of
             '20110404': begin
                keep = lindgen(n_elements(allstruct))
; fix some headers
;               allstruct[where(strmatch(allstruct.img_root,'*16jun2010.0074.fits*'))].obj = $
;                 '2059baboquivari,Bfilter'
             end
             else: message, 'Update me!'
          endcase
          struct = allstruct[keep]
          struct.obj = strtrim(struct.obj,2)
          nobj = n_elements(struct)
; match with the object list to add the formal object name
;         struct = struct_addtags(struct,replicate({target: ''},nobj))
;         for jj = 0, n_elements(targetinfo)-1 do begin
;            match = where(strmatch(struct.img_root,'*'+targetinfo[jj].file+'*'))
;            if (match[0] ne -1) then struct[match].target = targetinfo[jj].target
;         endfor
          struct_print, struct_trimtags(struct,select=['img_root',$
            'obj','target','type','exp','am','filter','naxis1','naxis2']), $
            file=repstr(structfile[inight],'.fits','.txt')
          mwrfits, struct, structfile[inight], /create
          popd
       endfor
    endif

; build the master bias frames, including the 512x512 mini-biases 
    if keyword_set(mkbias) then begin
       for inight = 0, nnight-1 do begin
          pushd, datapath+night[inight]
          struct = mrdfits(structfile[inight],1)
          aycamp_img_mkbias, struct, clobber=clobber
          allbias = file_search('Bias/Bias?x?.fits',count=nbias)
          for ibias = 0, nbias-1 do begin
             info = file_info(allbias[ibias])
;            if (info.symlink eq 0) then begin
;               splog, 'Reading '+allbias[ibias]
;               bias2048 = mrdfits(allbias[ibias],0,hdr,/silent)
;               bias512 = bias2048[768:1279,768:1279] ; see the observing log
;               outfile = repstr(allbias[ibias],'.fits','_512.fits')
;               splog, 'Writing '+outfile
;               mwrfits, bias512, outfile, hdr, /create
;            endif 
          endfor 
          popd
       endfor
    endif 

; build the master dome flat-fields, including the 512x512
; mini-flat-fields 
    if keyword_set(mkdflat) then begin
       for inight = 0, nnight-1 do begin
          pushd, datapath+night[inight]
          struct = mrdfits(structfile[inight],1)
          aycamp_img_mkdflat, struct, clobber=clobber
          alldflat = file_search('Flats/DFlat_*?x?.fits',count=ndflat)
          for idflat = 0, ndflat-1 do begin
             info = file_info(alldflat[idflat])
;            if (info.symlink eq 0) then begin
;               splog, 'Reading '+alldflat[idflat]
;               dflat2048 = mrdfits(alldflat[idflat],0,hdr,/silent)
;               dflat512 = dflat2048[768:1279,768:1279] ; see the observing log
;               outfile = repstr(alldflat[idflat],'.fits','_512.fits')
;               splog, 'Writing '+outfile
;               mwrfits, dflat512, outfile, hdr, /create
;            endif
          endfor 
          popd
       endfor 
    endif

; #########################
; assign each science image a bias frame and dome flat-field
    if keyword_set(assign_calib) then begin
       pushd, datapath
       for inight = 0, nnight-1 do begin
          biasfiles = file_search(night[inight]+'/Bias/Bias*.fits',count=nbias)
          if (nbias gt 0) then begin
             biasinfo1 = aycamp_forage(biasfiles)
             if (n_elements(biasinfo) eq 0) then $
               biasinfo = biasinfo1 else $
                 biasinfo = [biasinfo,biasinfo1]
          endif
          dflatfiles = file_search(night[inight]+'/Flats/DFlat*.fits',count=ndflat)
          if (ndflat gt 0) then begin
             dflatinfo1 = aycamp_forage(dflatfiles)
             if (n_elements(dflatinfo) eq 0) then $
               dflatinfo = dflatinfo1 else $
                 dflatinfo = [dflatinfo,dflatinfo1]
          endif 
       endfor 

       for inight = 0, nnight-1 do begin
          thisstructfile = night[inight]+'/'+night[inight]+'_imgstrct.fits'
          struct = mrdfits(thisstructfile,1,/silent)
          obj = where(((struct.type eq 'OBJ') or (struct.type eq 'STD')) and $
            (struct.flg_anly ne 0),nobj)
          for iobj = 0, nobj-1 do begin
             bin = strtrim(struct[obj[iobj]].cbin,2)+'x'+$
               strtrim(struct[obj[iobj]].rbin,2)
             filt = strtrim(struct[obj[iobj]].filter,2)
             naxis2 = strtrim(struct[obj[iobj]].naxis2,2)
; do the necessary bias frames and dome flats exist?  if not, then
; create a soft link pointing to the appropriate file from another
; night (if it exists)
             if (naxis2 eq 512) then suffix = '_512' else suffix = ''
             biasfile = 'Bias/Bias'+bin+suffix+'.fits'
             dflatfile = 'Flats/DFlat_'+filt+bin+suffix+'.fits'
             if (file_test(night[inight]+'/'+biasfile) eq 0) then begin
                if (n_elements(biasinfo) ne 0) then begin
                   these = where((biasinfo.cbin eq struct[obj[iobj]].cbin) and $
                     (biasinfo.rbin eq struct[obj[iobj]].rbin) and $
                     (biasinfo.naxis2 eq struct[obj[iobj]].naxis2),nthese)
                   if (nthese eq 0) then begin
                      splog, 'No master bias frame found for binning '+$
                        bin+' and axis size '+strtrim(naxis2,2)+'x'+strtrim(naxis2,2)+'!'
                   endif else begin
                      mindiff = min(abs(biasinfo[these].jdmean-struct[obj[iobj]].date),mindx)
                      pushd, night[inight]+'/Bias'
                      spawn, 'ln -sfv ../../'+strtrim(biasinfo[these[mindx]].file,2)+$
                        ' '+file_basename(biasfile), /sh
                      popd
                      struct[obj[iobj]].bias_fil = biasfile ; assigned bias frame
                   endelse
                endif
             endif else struct[obj[iobj]].bias_fil = biasfile ; assigned bias frame

; find flats with the right binning and filter name
             if (file_test(night[inight]+'/'+dflatfile) eq 0) then begin
                if (n_elements(dflatinfo) ne 0) then begin
                   these = where((dflatinfo.cbin eq struct[obj[iobj]].cbin) and $
                     (dflatinfo.rbin eq struct[obj[iobj]].rbin) and $
                     (strtrim(dflatinfo.filter,2) eq filt) and $
                     (dflatinfo.naxis2 eq struct[obj[iobj]].naxis2),nthese)
                   if (nthese eq 0) then begin
                      splog, 'No master dome-flat found for binning '+$
                        bin+', axis size '+strtrim(naxis2,2)+'x'+$
                        strtrim(naxis2,2)+' and filter '+filt+'!'
                   endif else begin
                      mindiff = min(abs(dflatinfo[these].jdmean-struct[obj[iobj]].date),mindx)
                      pushd, night[inight]+'/Flats'
                      spawn, 'ln -sfv ../../'+strtrim(dflatinfo[these[mindx]].file,2)+$
                        ' '+file_basename(dflatfile), /sh
                      popd
                      struct[obj[iobj]].flat_fil = dflatfile ; assigned dome-flat field
                   endelse
                endif
             endif else struct[obj[iobj]].flat_fil = dflatfile ; assigned dome-flat field
          endfor ; close object loop
          splog, '### Night '+night[inight]
          splog, 'Updating '+thisstructfile
          mwrfits, struct, thisstructfile, /create
          struct_print, struct_trimtags(struct,select=$
            ['img_root','obj','type','exp',$
            'am','filter','naxis1','naxis2']), $
            file=repstr(thisstructfile,'.fits','.txt')
; test info for the screen
          if (nobj ne 0) then begin
             struct_print, struct_trimtags(struct[obj],select=$
               ['obj','filter','naxis2','cbin','rbin','bias_fil','flat_fil'])
             print
          endif
       endfor    ; close night loop
       popd
    endif    

; #########################
; overscan-subtract, bias-subtract, and flat-field all the science
; images
    if keyword_set(proc) then begin
       for inight = 0, nnight-1 do begin
          pushd, datapath+night[inight]
          struct = mrdfits(structfile[inight],1)
          obj = where(((struct.type eq 'OBJ') or (struct.type eq 'STD')) and $
            (struct.flg_anly ne 0),nobj)
          if (nobj ne 0) then begin
; check for missing bias frames and dome flats
             miss = where(strtrim(struct[obj].bias_fil eq '') or $
               strtrim(struct[obj].flat_fil eq ''),nmiss)
             if (nmiss ne 0) then begin
                splog, 'The following objects are missing bias and/or dome flat-fields!'
                stop
             endif
             aycamp_img_proc, struct[obj], outroot=outroot, $
               /nozip, clobber=clobber
          endif
          popd
       endfor 
    endif 

;; #########################
;    if keyword_set(astrom) then begin
;       for inight = 0, nnight-1 do begin
;          pushd, datapath+night[inight]
;          allfiles = file
;          python ~/local/bin/anet.py
;          
;stop          
;          popd
;       endfor
;    endif       
       
; #########################
; derive photometric zeropoints
    if keyword_set(zeropoints) then begin
; read the Stetson (2000) photometric catalog    
       stetson = read_00stetson() 
       for inight = 0, nnight-1 do begin
          pushd, datapath+night[inight]
; ---------------
; copy the standard-star files and derive rough astrometric headers
          if keyword_set(astrom) then begin
             struct = mrdfits(night[inight]+'_imgstrct.fits',1,/silent)
             std = where(strmatch(struct.obj,'*landolt*',/fold) or $
               strmatch(struct.obj,'*standard*',/fold),nstd)
             if (nstd eq 0) then begin
                splog, 'No standards observed for night '+night[inight]+'!'
             endif else begin
                struct = struct[std]
                pushd, 'Std'
                alloutfile = strarr(nstd)
                for istd = 0, nstd-1 do begin
                   infile = '../Final/f_'+strtrim(struct[istd].img_root,2);+'.gz'
                   if file_test(infile) then begin
                      outfile = strtrim(struct[istd].target,2)+'_'+$
                        string(struct[istd].frame,format='(I4.4)')+'_'+$
                        strtrim(struct[istd].filter,2)+'.fits';.gz'
                      spawn, 'rsync -auv '+infile+' '+outfile, /sh
                      alloutfile[istd] = outfile
;                     spawn, 'ln -svf '+infile+' '+outfile, /sh
                   endif else message, 'Fix this'
                endfor
; call astrometry.net if no header exists
                for istd = 0, nstd-1 do begin
                   hdr = headfits(alloutfile[istd])
                   anjob = sxpar(hdr,'AN_JOBID',count=job)
                   if (job eq 0) or keyword_set(clobber) then $
                     spawn, 'python ${HOME}/local/bin/anet.py '+alloutfile[istd], /sh
                endfor 
                popd
             endelse
          endif 
;; ---------------
;; build soft links
;          if keyword_set(makelinks) then begin
;             struct = mrdfits(night[inight]+'_imgstrct.fits',1,/silent)
;             std = where(strmatch(struct.obj,'*landolt*',/fold) or $
;               strmatch(struct.obj,'*standard*',/fold),nstd)
;             if (nstd eq 0) then begin
;                splog, 'No standards observed for night '+night[inight]+'!'
;             endif else begin
;                struct = struct[std]
;                pushd, 'Std'
;                for istd = 0, nstd-1 do begin
;                   infile = '../Final/f_'+strtrim(struct[istd].img_root,2);+'.gz'
;                   if file_test(infile) then begin
;                      outfile = strtrim(struct[istd].target,2)+'_'+$
;                        string(struct[istd].frame,format='(I4.4)')+'_'+$
;                        strtrim(struct[istd].filter,2)+'.fits';.gz'
;                      spawn, 'ln -svf '+infile+' '+outfile, /sh
;                   endif else stop
;                endfor
;                popd
;             endelse
;          endif
; ---------------
; build SE catalogs
          if keyword_set(sextractor) then begin
             pushd, 'Std'
             imagelist = file_search('*.fits',count=nimage)
;            aycamp_do_sextractor, imagelist

             catlist = file_basename(repstr(imagelist,'.fits','.cat'))
       
; initialize the SE configuration parameters
             config = init_sex_config(nimage)
             configfile = 'sex.config'
             
             config.catalog_name = catlist
             config.catalog_type = 'FITS_LDAC'
             config.parameters_name = sexparam
             config.filter_name = sexconv
             config.starnnw_name = sexnnw
          
             config.detect_thresh = 20.0
             config.analysis_thresh = 20.0
             config.detect_minarea = 10.0
             config.deblend_mincont = 0.1
             config.interp_type = 'NONE'
             
             config.seeing_fwhm = 1.0
;            config.mag_zeropoint = 0.0

             config.checkimage_type = 'NONE' ; SEGMENTATION
             config.phot_apertures = photaper
             im_mwrfits, config, configfile+'.fits', /clobber

             im_sex, imagelist, config, silent=silent ; do not pass CONFIGFILE
             popd
          endif 

; ---------------
; solve for the BVRI zeropoints
          if keyword_set(solve) then begin
             pushd, 'Std'
             info = aycamp_forage('*.fits')
             allfilt = strtrim(info.filter,2)
             filt = allfilt[uniq(allfilt,sort(allfilt))]
             nfilt = n_elements(filt)
             for ifilt = 0, nfilt-1 do begin
                zpt = {filter: filt[ifilt], night: night[inight], $
                  nallstar: 0, nstar: 0, zpt: 0.0, zpt_err: 0.0, $
                  extcoeff: 0.0, extcoeff_err: 0.0, colorcoeff: 0.0, $
                  colorcoeff_err: 0.0}

                these = where(filt[ifilt] eq allfilt,nimage)
                imagelist = strtrim(info[these].file,2)
                catlist = repstr(imagelist,'.fits','.cat')

                delvarx, info
                for ii = 0, nimage-1 do begin
                   cat = mrdfits(catlist[ii],2)
                   spherematch, cat.xwin_world, cat.ywin_world, $
                     stetson.ra, stetson.dec, 5.0/3600.0, m1, m2
                   if (m1[0] ne -1L) then begin
                      stet = stetson[m2]
                      info1 = struct_trimtags(cat[m1],select=[$
                        'XWIN_IMAGE','YWIN_IMAGE','XWIN_WORLD','YWIN_WORLD',$
                        'MAG_AUTO','MAGERR_AUTO','MAG_APER','MAGERR_APER',$
                        'FWHM_WORLD','FLAGS'])
                      moretags = replicate({file: file_basename(imagelist[ii]), $
                        date: '', exptime: 0.0, airmass: 0.0},n_elements(info1))
                      info1 = struct_addtags(struct_addtags(moretags,info1),stet)
                      info1.date = sxpar(headfits(imagelist[ii]),'DATE-OBS')
                      info1.exptime = sxpar(headfits(imagelist[ii]),'EXPTIME')
                      info1.airmass = sxpar(headfits(imagelist[ii]),'AIRMASS')
                      if (n_elements(info) eq 0L) then info = info1 else info = [info,info1]
                   endif 
                endfor ; close image loop
                outfile = 'landolt_'+filt[ifilt]+'.fits'
                im_mwrfits, info, outfile, /clobber

; compute the zeropoint
                case filt[ifilt] of
                   'B': begin
                      these = where((info.mag_b lt 90.0) and (info.mag_v lt 90.0),nthese)
                      mtrue = info[these].mag_b & color = info[these].mag_b-info[these].mag_v
                      zpt.extcoeff = 0.25
                   end
                   'V': begin
                      these = where((info.mag_r lt 90.0) and (info.mag_v lt 90.0),nthese)
                      mtrue = info[these].mag_r & color = info[these].mag_v-info[these].mag_r
                      zpt.extcoeff = 0.15
                   end
                   'R': begin
                      these = where((info.mag_r lt 90.0) and (info.mag_v lt 90.0),nthese)
                      mtrue = info[these].mag_r & color = info[these].mag_v-info[these].mag_r
                      zpt.extcoeff = 0.10
                   end
                   'I': begin
                      these = where((info.mag_i lt 90.0) and (info.mag_r lt 90.0),nthese)
                      mtrue = info[these].mag_i & color = info[these].mag_r-info[these].mag_i
                      zpt.extcoeff = 0.08
                   end 
                   else: message, 'Fix this'
                endcase
                
; iteratively reject outliers
                minst = info[these].mag_auto + 2.5*alog10(info[these].exptime)

                djs_iterstat, mtrue-minst, median=zpt_guess, $
                  mask=mask, sigrej=2.0, maxiter=100
                goodstars = where(mask,nstar)
                if (nstar eq 0) then message, 'This is very bad'

                zpt.nallstar = nthese
                zpt.nstar = nstar

                minst = minst[goodstars]
                mtrue = mtrue[goodstars]
                color = color[goodstars]
                airmass = info[these[goodstars]].airmass
                minst_err = info[these[goodstars]].magerr_auto

                parinfo = replicate({value: 0.0D, fixed: 0},4) ; 4-parameter fit
                parinfo[0].value = zpt_guess & parinfo[0].fixed = 0 ; zeropoint
                
                parinfo[1].value = zpt.extcoeff & parinfo[1].fixed = 1
;               parinfo[1].value = 0.1D & parinfo[1].fixed = 0 ; linear extinction coeff
;               parinfo[2].value = colorcoeff[iq,ib] & parinfo[2].fixed = 1 ; color term

                parinfo[2].value = 0.0D & parinfo[2].fixed = 1 ; color term
;               parinfo[2].value = 0.1D & parinfo[2].fixed = 0 ; color term
                parinfo[3].value = 0.0D & parinfo[3].fixed = 1 ; quadratic extinction coeff

                functargs = {minst: minst, airmass: airmass, color: color}
                param = mpfitfun('zeropoint_func',mtrue,minst,minst_err,$
                  functargs=functargs,perror=perror,quiet=1,parinfo=parinfo,$
                  dof=dof,bestnorm=bestnorm,yfit=minst_best)
;               niceprint, minst_best, minst, mtrue-param[0]

                zpt.zpt = param[0]
                zpt.zpt_err = djsig(minst-minst_best,sigrej=2.0)
                zpt.extcoeff = param[1]
;               zpt.extcoeff_err = extcoeff_err[ib]
                zpt.colorcoeff = param[2] ; note!
                zpt.colorcoeff_err = 0.0

;; generate a QA plot
;               djs_plot, [0], [0], /nodata, xthick=2.0, ythick=2.0, $
;                 charsize=2.0, charthick=2.0, xsty=1, ysty=1, $
;                 xrange=[-0.3,2.0], xtitle='R-I', ytitle='I-i', $
;                 yrange=minmax(mtrue-minst)+[-0.1,0.1]
;               plotsym, 8, 1.5, fill=1
;               djs_oplot, color, mtrue-minst, ps=8
;               cc = get_kbrd(1)

; write out the FITS table for this bandpass
                outfile = 'landolt_zpt_'+filt[ifilt]+'.fits'
                im_mwrfits, zpt, outfile, /clobber
                splog, 'Filter '+filt[ifilt]
                struct_print, zpt
             endfor    ; close filter loop
          endif 
       endfor ; close night                  
    endif 

; #########################
; unpack the various projects into dedicated directories
    if keyword_set(unpack_projects) then begin

       for inight = 0, nnight-1 do begin
          struct = mrdfits(datapath+night[inight]+'/'+night[inight]+$
            '_imgstrct.fits',1,/silent)
          obj = where(((struct.type eq 'OBJ') or (struct.type eq 'STD')) and $
            (struct.flg_anly ne 0),nobj)
          if (nobj ne 0) then begin
             struct = struct[obj]
             struct_print, struct_trimtags(struct,$
               select=['img_root','obj','filter'])
          endif
; -------------------------
; crescent
          project = 'crescent'
          outdir = datapath+'projects/'+project+'/'
          if (file_test(outdir,/dir) eq 0) then $
            spawn, 'mkdir -p '+outdir, /sh
          these = where($
            strmatch(struct.obj,'*crescent*',/fold),nthese)
          if (nthese ne 0) then begin
             lcshalpha_unpack, struct[these], outdir=outdir, $
               indir=datapath+night[inight]+'/Final/', $
               night=night[inight]
          endif
; -------------------------
; PNe
          project = 'pne'
          outdir = datapath+'projects/'+project+'/'
          if (file_test(outdir,/dir) eq 0) then $
            spawn, 'mkdir -p '+outdir, /sh
          these = where($
            strmatch(struct.obj,'*abell43*',/fold) or $
            strmatch(struct.obj,'*ic4593*',/fold) or $
            strmatch(struct.obj,'*abell39*',/fold) or $
            strmatch(struct.obj,'*abell72*',/fold),nthese)
          if (nthese ne 0) then begin
             lcshalpha_unpack, struct[these], outdir=outdir, $
               indir=datapath+night[inight]+'/Final/', $
               night=night[inight]
          endif
; -------------------------
; star-forming regions
          project = 'sfregions'
          outdir = datapath+'projects/'+project+'/'
          if (file_test(outdir,/dir) eq 0) then $
            spawn, 'mkdir -p '+outdir, /sh
          these = where($
            strmatch(struct.obj,'*m20*',/fold) or $
            strmatch(struct.obj,'*lagoon*',/fold) or $
            strmatch(struct.obj,'*omega*',/fold) or $
            strmatch(struct.obj,'*trifid*',/fold) or $
            strmatch(struct.obj,'*eagle*',/fold),nthese)
          if (nthese ne 0) then begin
             lcshalpha_unpack, struct[these], outdir=outdir, $
               indir=datapath+night[inight]+'/Final/', $
               night=night[inight]
          endif
; -------------------------
; bulgedisk
          project = 'bulgedisk'
          outdir = datapath+'projects/'+project+'/'
          if (file_test(outdir,/dir) eq 0) then $
            spawn, 'mkdir -p '+outdir, /sh
          these = where($
            strmatch(struct.obj,'*m85*',/fold) or $
            strmatch(struct.obj,'*m51*',/fold) or $
            strmatch(struct.obj,'*m101*',/fold) or $
            strmatch(struct.obj,'*ngc4374*',/fold),nthese)
          if (nthese ne 0) then begin
             lcshalpha_unpack, struct[these], outdir=outdir, $
               indir=datapath+night[inight]+'/Final/', $
               night=night[inight]
          endif
; -------------------------
; interacting galaxies
          project = 'intergals'
          outdir = datapath+'projects/'+project+'/'
          if (file_test(outdir,/dir) eq 0) then $
            spawn, 'mkdir -p '+outdir, /sh
          these = where($
            strmatch(struct.obj,'*NGC4038*',/fold) or $
            strmatch(struct.obj,'*NGC*4676*',/fold) or $
            strmatch(struct.obj,'*NGC4520*',/fold) or $
            strmatch(struct.obj,'*NGC4546*',/fold),nthese)
          if (nthese ne 0) then begin
             lcshalpha_unpack, struct[these], outdir=outdir, $
               indir=datapath+night[inight]+'/Final/', $
               night=night[inight]
          endif
; -------------------------
; open/globular clusters
          project = 'clusters'
          outdir = datapath+'projects/'+project+'/'
          if (file_test(outdir,/dir) eq 0) then $
            spawn, 'mkdir -p '+outdir, /sh
          these = where($
            strmatch(struct.obj,'*M71*',/fold) or $
            strmatch(struct.obj,'*M39*',/fold),nthese)
          if (nthese ne 0) then begin
             lcshalpha_unpack, struct[these], outdir=outdir, $
               indir=datapath+night[inight]+'/Final/', $
               night=night[inight]
          endif
; -------------------------
; baboquivari
          project = 'baboquivari'
          outdir = datapath+'projects/'+project+'/'
          if (file_test(outdir,/dir) eq 0) then $
            spawn, 'mkdir -p '+outdir, /sh
          these = where(strmatch(struct.obj,'*babo*',/fold),nthese)
          if (nthese ne 0) then begin
             lcshalpha_unpack, struct[these], outdir=outdir, $
               indir=datapath+night[inight]+'/Final/', $
               night=night[inight]
          endif
; -------------------------
; standard-star fields (don't need the weight maps)
          project = 'landolt'
          outdir = datapath+'projects/'+project+'/'
          if (file_test(outdir,/dir) eq 0) then $
            spawn, 'mkdir -p '+outdir, /sh
          these = where(strmatch(struct.obj,'*landolt*',/fold) or $
            strmatch(struct.obj,'*standard*',/fold),nthese)
          if (nthese ne 0) then begin
             lcshalpha_unpack, struct[these], outdir=outdir, $
               indir=datapath+night[inight]+'/Final/', $
               night=night[inight], /noweight
          endif
       endfor
    endif

return
end
