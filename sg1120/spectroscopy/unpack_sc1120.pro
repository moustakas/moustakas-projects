;+
; NAME:
;       UNPACK_SC1120
;
; PURPOSE:
;       Unpack the raw SC1120 standard-star observations. 
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;       Trim the relevant portions of each quadrant where the star was
;       illuminated.  Unfortunately, the arc lamps are saturated at
;       the position of the standards and the dome flats, so trim a
;       different region for those objects.  
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Nov 05, U of A
;-

pro unpack_sc1120, debug=debug, wfits=wfits

    rootpath = '/home/ioannis/research/projects/sc1120/'

    cwdpath = rootpath+'redux/'
    rawpath = rootpath+'raw/'
    datapath = rootpath+'data/'

    pushd, rawpath
    flist = file_search('*.fits*',count=fcount)
    popd

    if keyword_set(wfits) then begin
       splog, 'Delete all FITS files from '+datapath+'Q1/: [Y/N]?' & cc = get_kbrd(1)
       if strupcase(cc) eq 'Y' then spawn, ['/bin/rm -f '+datapath+'Q1/*.fits'], /sh
       splog, 'Delete all FITS files from '+datapath+'Q2/: [Y/N]?' & cc = get_kbrd(1)
       if strupcase(cc) eq 'Y' then spawn, ['/bin/rm -f '+datapath+'Q2/*.fits'], /sh
       splog, 'Delete all FITS files from '+datapath+'Q3/: [Y/N]?' & cc = get_kbrd(1)
       if strupcase(cc) eq 'Y' then spawn, ['/bin/rm -f '+datapath+'Q3/*.fits'], /sh
       splog, 'Delete all FITS files from '+datapath+'Q4/: [Y/N]?' & cc = get_kbrd(1)
       if strupcase(cc) eq 'Y' then spawn, ['/bin/rm -f '+datapath+'Q4/*.fits'], /sh
    endif

    for i = 0L, fcount-1L do begin

       if keyword_set(wfits) then begin
          image = readfits(rawpath+flist[i],h,/silent) 
       endif else begin
;         image = readfits(rawpath+flist[i],h,/silent) 
          splog, 'Reading '+flist[i]+'.'
          h = headfits(rawpath+flist[i])
       endelse

; parse the FITS header into a data structure and then trim to a list
; of useful tag names

       hinfo = vlt_header_forage(h)
       hinfo = create_struct({rawfile: flist[i]},hinfo)

       subhinfo = struct_trimtags(hinfo,select=['RAWFILE','DATE_OBS','EXPTIME','RA','DEC',$
         'EQUINOX','UTC','OBS_NAME','OBS_TARG_NAME','TPL_ID','TPL_NAME','DPR_CATG','DPR_TYPE',$
         'SEQ_SPEC_TARG','INS_MODE','INS_SLIT_NAME','INS_SLIT_RA','INS_SLIT_DEC','INS_SLIT_WID',$
         'OCS_CON_QUAD','INS_GRIS1_DISP','INS_GRIS1_WLEN','ORIGFILE'])
       if (i eq 0L) then bighinfo = subhinfo else bighinfo = struct_append(bighinfo,subhinfo)

       quad = strtrim(subhinfo.ocs_con_quad)
       mask = strmid(flist[i],8,1) ; NOT GENERAL!

; MOS: standard stars

       if (strtrim(subhinfo.dpr_type,2) eq 'STD') then begin

          newh = update_vlt_header(h,hinfo,/object)
          
          stdname = strtrim(strupcase(subhinfo.obs_targ_name),2)
          outfile = 'a.Q'+quad+'-mask'+mask+'-'+stdname+'.fits'
          
       endif

; MOS: dome flats

       if (strtrim(subhinfo.dpr_type,2) eq 'FLAT,LAMP') then begin
          
          newh = update_vlt_header(h,hinfo,/domeflat)
          outfile = 'a.Q'+quad+'-mask'+mask+'-domeflat.fits'
          
       endif

; MOS: arc lamps

       if (strtrim(subhinfo.dpr_type,2) eq 'WAVE,LAMP') then begin

          newh = update_vlt_header(h,hinfo,/arclamp)
          outfile = 'a.Q'+quad+'-mask'+mask+'-arc.fits'
          
       endif
       
; MOS: bias

       if (strtrim(subhinfo.dpr_type,2) eq 'BIAS') then begin

          newh = update_vlt_header(h,hinfo,/bias)
          outfile = 'a.Q'+quad+'-mask'+mask+'-bias.fits'
          
       endif

; define the trim regions

       case quad of
          '1': if (mask eq 1) or (mask eq 2) then $
            if strmatch(flist[i],'*arc*') then xslit = 1268L else xslit = 1118L else $
            if strmatch(flist[i],'*arc*') then xslit = 1272L else xslit = 1122L
          '2': if strmatch(flist[i],'*arc*') then xslit = 1137L else xslit = 987L
          '3': if strmatch(flist[i],'*arc*') then xslit = 1202L else xslit = 1051L
          '4': if strmatch(flist[i],'*arc*') then xslit = 1261L else xslit = 1111L
       endcase

       imtrim = [0,4095,xslit+25*[-1,1]]

; ---------------------------------------------------------------------------       
; trim testing
; ---------------------------------------------------------------------------       

       if keyword_set(debug) and (strmatch(flist[i],'*arc*') or strmatch(flist[i],'*flat*')) then begin
          
          ncols = 2148
          colaxis = findgen(ncols)
          
          pixscale = 0.205      ; [arcsec/pixel]
          slitwidth = 10.0      ; [arcsec]

          npix = round(slitwidth/pixscale)+1L
          
          refprofile = fltarr(ncols)
          refprofile[ncols/2L-npix/2L:ncols/2L+npix/2L-1L] = 1.0
          
          profile = total(image[*,2000:2100],2)
          profile = profile / max(profile)
          profile[600:ncols-1L] = profile[600:ncols-1L]*20 ; boost the inner slits

          lagmax = ncols-1.0
          lagmin = 0.0
          dlag = 1.0
          lags = findgen((lagmax-lagmin)/dlag+1)*dlag-lagmax/2.0

          nfind = 8L
          
          corr = 1.0 - djs_correlate(refprofile,profile,lags)
          xslit = find_nminima(corr,colaxis,nfind=nfind,ypeak=yslit)

;         plot, colaxis, 1-corr, ps=-4
          plot, colaxis, profile, xsty=3, ysty=3, ps=10, title=flist[i]
          for ifind = 0L, nfind-1L do djs_oplot, xslit[ifind]*[1,1], !y.crange, color='red'

          print, 'Q'+quad+'-'+mask, string(xslit[sort(xslit)],format='(I4)')
;         cc = get_kbrd(1)

       endif
          
; ---------------------------------------------------------------------------
       
;      mm2pix = 0.015 ; [mm/pixel]
;
;      objinfo = struct_trimtags(hinfo,select=['*SLIT?_OBJ_RA'])
;      xpos = fltarr(8)
;      for i = 0L, 7L do xpos[i] = hms2dec(objinfo.(i)) * 15.0 / 3600.0 / 0.205
;
;      slitinfo = struct_trimtags(hinfo,select='*SLIT?_X')
;      slitx = fltarr(8)
;      for i = 0L, 7L do slitx[i] = slitinfo.(i) / mm2pix * 0.205
;      
;      speclength = float(hinfo.ins_adf_grism_spectlen) ; spectrum length [pixels]
;
;      refpix = float(hinfo.crpix1) ; [pixel]
;      slitx = float(hinfo.ins_slit5_x) / mm2pix
;      slity = float(hinfo.ins_slit5_y) / mm2pix
;      slitdimx = float(hinfo.ins_slit5_dimx) / mm2pix
;      slitdimy = float(hinfo.ins_slit5_dimy) / mm2pix
;      print, slitx, slity, slitdimx, slitdimy
;      
;      plot, total(image,2)
;      djs_oplot, slity*[1,1], !y.crange, color='red'

       if keyword_set(wfits) then begin

          timage = transpose(image)
          itrim, timage, header=h, trim=imtrim ; note transposition!
          
          splog, 'Writing '+datapath+'Q'+quad+'/'+outfile+'.'
          mwrfits, timage, datapath+'Q'+quad+'/'+outfile, newh, /create
;         spawn, ['gzip -f '+datapath+outfile], /sh
          
       endif

    endfor 
    print

stop
    
return
end
