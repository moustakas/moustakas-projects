pro unpack_ediscs_spec1d, outinfo, debug=debug
; jm09jun09nyu - unpack the EDisCS 1D spectra
; jm10apr28ucsd - do not exclude the left- and right most 50 pixels;
;   also store the pixel size in the output data structure

    datapath = ediscs_path(/bamford_spec1d)
    outpath = ediscs_path(/spec1d)
    catpath = ediscs_path(/catalogs)
    mycatpath = ediscs_path(/mycatalogs)

    infofile = mycatpath+'ediscs_spec1d_info.fits'
    
    if (file_test(outpath,/directory) ne 0L) then begin
       splog, 'Remove all FITS files from '+outpath+' [Y/N]?'
       cc = get_kbrd(1)
       if (strupcase(cc) eq 'Y') then begin
          flist = file_search(outpath+'*.fits*')
          rmfile, flist
;         spawn, ['/bin/rm '+outpath+'*.fits*'], /sh
       endif
    endif

    thisfile = catpath+'combined_spec_table_selected_2007-08-22.fits.gz'
    splog, 'Reading '+thisfile
    info = mrdfits(thisfile,1)
    ngal = n_elements(info)

    noname = where(strcompress(info.spec_ediscsid,/remove) eq '',nnoname)
    if (nnoname ne 0) then info[noname].spec_ediscsid = $
      'NONAME'+string(lindgen(nnoname)+1,format='(I2.2)')

    fluxfile = strtrim(info.spec_name,2)+'.fits'
    ferrfile = repstr(fluxfile,'_sci_','_sig_')
    outfile = repstr(repstr(repstr(repstr(repstr($
      fluxfile,'___','_'),'_sci_','_'),'_1D',''),'.FT',''),'_T','')

    cluster = strarr(ngal)
    for igal = 0, ngal-1 do cluster[igal] = $
      repstr(strmid(outfile[igal],0,strpos(outfile[igal],'_','')),'.fits')

; construct the output "info" structure; PHOTSPLIT indicates two
; spectra, but a single photometric object; these objects have an "a"
; or a "b" in the old ID numbers
    photsplit = strcompress(info.spec_old_ediscsid_u,/remove) ne ''

    outinfo = {ediscs_id: -1, galaxy: '', specfile: '', $
      ra: -999.0D, dec: -999.0D, ebv_mw: 0.0, z: 0.0, zquality: 0.0, $
      minwave: 0.0D, maxwave: 0.0D, dwave: 0.0D, original_memberflag: '', memberflag: '', $
      targetflag: 0, photsplit: 0, type: '', comment: '', $
      cluster: '', cluster_subname: '', cluster_fullname: '', cluster_z: -999.0, $
      cluster_sigma: -999.0, cluster_sigma_p: -999.0, cluster_sigma_m: -999.0, $
      cluster_nmember: 0L, cluster_bcg: '', cluster_fieldname: ''}
    outinfo = replicate(outinfo,ngal)

    outinfo.ediscs_id = lindgen(ngal)
    outinfo.galaxy = strtrim(info.spec_ediscsid,2)
    outinfo.specfile = outfile
    outinfo.cluster = cluster
    outinfo.z = info.spec_z
    outinfo.original_memberflag = strtrim(info.spec_mem,2)
;   outinfo.memberflag = strtrim(info.spec_mem,2)
    outinfo.targetflag = fix(info.spec_tar)
    outinfo.photsplit = fix(photsplit)
    outinfo.type = strtrim(info.spec_type,2)
    outinfo.comment = strtrim(info.spec_comment,2)
    
    good = where(strmatch(info.spec_ra,'*99:99:99.9*') eq 0B,ngoodcoords)
    outinfo[good].ra = im_hms2dec(info[good].spec_ra)*15.0D
    outinfo[good].dec = im_hms2dec(info[good].spec_dec)
;   outinfo.ra = strtrim(info.spec_ra,2)
;   outinfo.dec = strtrim(info.spec_dec,2)

; additional cluster properties
    cl = rsex(catpath+'ediscs_clusters.sex')
    groups = rsex(catpath+'ediscs_group_catalog.sex')
    
; recover the average MW reddening values
    readcol, catpath+'E_BV_corrections.v2.2.tab', $
      ebv_cluster, ebv_la, ebv_la, ebv_mw, $
      format='A,F,F,F', comment='#', /silent
    ebv_cluster = strcompress(repstr(ebv_cluster,'-',''),/remove)
    for ii = 0, n_elements(ebv_cluster)-1 do begin
       if strmatch(ebv_cluster[ii],'*1054*') then mmax = 8 else mmax = 6
       match = where(strmid(ebv_cluster[ii],0,mmax) eq strtrim(outinfo.cluster,2),nmatch)
       if (nmatch ne 0) then outinfo[match].ebv_mw = ebv_mw[ii]
    endfor

; loop on each object    
    for igal = 0, ngal-1 do begin
       print, format='("EDisCS galaxy ",I0,"/",I0,A0,$)', $
         igal+1, ngal, string(13b)

; assign group/cluster membership; first check if the object is
; in Bianca's group catalog; if not, then use the official
; membership flag 
       match = where(strtrim(outinfo[igal].galaxy,2) eq strtrim(groups.galaxy,2),nmatch)
       if (nmatch ne 0) then begin
          subname = strtrim(strmid(groups[match].cluster,13),2)
          outinfo[igal].memberflag = '1'+subname
       endif else outinfo[igal].memberflag = outinfo[igal].original_memberflag
       
; add cluster/group properties *just* for members
       if strmatch(outinfo[igal].memberflag,'*1*') then begin
          subname_ref = strtrim(repstr(cl.cluster_subname,'-',''),2)
          subname = strtrim(strmid(outinfo[igal].memberflag,1),2)
          match = where((subname eq subname_ref) and $
            strtrim(outinfo[igal].cluster,2) eq strtrim(cl.cluster,2),nmatch)
          if (nmatch ne 1) then message, 'Fix me'
          outinfo[igal].cluster_fullname  = strtrim(cl[match].cluster_fullname,2)
          outinfo[igal].cluster_subname   = repstr(cl[match].cluster_subname,'-','')
          outinfo[igal].cluster_z         = cl[match].z
          outinfo[igal].cluster_sigma     = cl[match].sigma
          outinfo[igal].cluster_sigma_p   = cl[match].sigma_p
          outinfo[igal].cluster_sigma_m   = cl[match].sigma_m
          outinfo[igal].cluster_nmember   = cl[match].nmember
          outinfo[igal].cluster_bcg       = cl[match].bcg
          outinfo[igal].cluster_fieldname = cl[match].fieldname
       endif
          
; parse the spectra       
       if (file_test(datapath+fluxfile[igal]) eq 0L) then $
         message, 'Spectrum '+datapath+fluxfile[igal]+' not found'
       if (file_test(datapath+ferrfile[igal]) eq 0L) then $
         splog, 'Error spectrum '+datapath+ferrfile[igal]+' not found'
       
       flux = mrdfits(datapath+fluxfile[igal],0,h,/silent)
       ferr = mrdfits(datapath+ferrfile[igal],0,/silent)
       wave = make_wave(h,cd1_1=dwave,crval1=minwave)

; a negative error means a bad pixel that has been interpolated; set
; the errors in these pixels to a large number
       neg = where(ferr le 0.0,nneg)
       if (nneg ne 0L) then ferr[neg] = 1D16
;      ferr = sqrt(ferr^2)

       flux = flux[10:n_elements(flux)-10]
       ferr = ferr[10:n_elements(flux)-10]
       wave = wave[10:n_elements(flux)-10]
       npix = n_elements(flux)

; ---------------------------------------------------------------------------
; clean up cosmic rays
; ---------------------------------------------------------------------------
;       emask = emission_mask(wave,z=info[igal].spec_z,good=good,bad=bad,$
;         width=maskwidth,qsowidth=qsowidth,/telluric,/bluemask,/qsomask,$
;         spectrum=flux,/cosmic,sigrej_cosmic=4.0,crmask=crmask,skymask=skymask)
;;      mask = (crmask + skymask) gt 1B
;       mask = crmask
;       good = where(mask eq 1B,comp=bad)
;       flux = djs_maskinterp(flux,(mask eq 0B),xval=wave)
; ---------------------------------------------------------------------------
          
; funky edges
       case outfile[igal] of
          'cl1018_m05_s07_1.fits': begin
;            debug = 1
             good = where((wave gt 5650.0) and (wave lt 8700.0),ngood)
          end
          'cl1018_m07_s07_1.fits': begin
;            debug = 1
             good = where((wave gt 4250.0) and (wave lt 7390.0),ngood)
          end
          else: good = lindgen(npix)
       endcase
       flux = flux[good]
       ferr = ferr[good]
       wave = wave[good]

       zero = where(ferr eq 0.0,nzero)
       if (nzero ne 0L) then begin
          splog, 'Zero uncertainty in spectrum '+fluxfile[igal]
       endif

       kl = k_lambda(wave,/odonnell,r_v=3.1)
       flux = flux*10^(0.4*outinfo[igal].ebv_mw*kl) ; correct for foreground Galactic extinction
       ferr = ferr*10^(0.4*outinfo[igal].ebv_mw*kl)

       outinfo[igal].minwave = min(wave)
       outinfo[igal].maxwave = max(wave)
       outinfo[igal].dwave = dwave

       mkhdr, header, float(flux), /extend
       sxdelpar, header, 'COMMENT'
       sxaddpar, header, 'OBJECT', strtrim(outinfo[igal].galaxy,2)
       sxaddpar, header, 'GALAXY', strtrim(outinfo[igal].galaxy,2)
       sxaddpar, header, 'CRVAL1', float(min(wave)), ' wavelength at CRPIX1'
       sxaddpar, header, 'CRPIX1', float(1.0), ' reference pixel number'
       sxaddpar, header, 'CD1_1',  float(dwave), ' dispersion [Angstrom/pixel]'
       sxaddpar, header, 'CDELT1', float(dwave), ' dispersion [Angstrom/pixel]'
       sxaddpar, header, 'CTYPE1', 'LINEAR', ' projection type'
;      sxaddpar, header, 'Z', float(info[igal].spec_z), ' spectroscopic redshift'
;      sxaddpar, header, 'VDISP', 100.0

       if keyword_set(debug) then begin
;         fstats = im_stats(1D17*flux,sigrej=5.0)
;         yrange = [fstats.minrej,fstats.maxrej]
          yrange = minmax(flux*1D17)
;         yramge = minmax(flux[good]*1D17)
;         yrange = fstats.median+fstats.sigma_rej*[-3.0,5.0]
;         ploterror, wave, smooth(1D17*flux,3), smooth(1D17*ferr,3), $
;           ps=10, xsty=3, ysty=3, thick=1.0, charthick=2.0, charsize=1.5, $
          plot, wave, 1D17*flux, ps=10, xsty=3, ysty=3, thick=1.0, charthick=2.0, charsize=1.5, $
            xthick=2.0, ythick=2.0, title=outinfo[igal].galaxy+' ['+repstr(outfile[igal],'.fits','')+'] z='+$
            strtrim(string(info[igal].spec_z,format='(F12.3)'),2), $
            xtitle='Observed Wavelength', ytitle='Relative Flux', yrange=yrange
;            if (bad[0] ne -1L) then djs_oplot, wave[bad], 1D17*flux[bad], ps=4, color='green'
          cc = get_kbrd(1)
       endif

       splog, 'Writing '+outpath+outfile[igal]
       mwrfits, float(flux), outpath+outfile[igal], header, /create
       mwrfits, float(ferr), outpath+outfile[igal], /silent
    endfor 

;   ww = where(strmatch(outinfo.memberflag,'*1a') and strtrim(outinfo.cluster,2) eq 'cl1103')
    ww = where(strtrim(outinfo.memberflag,2) ne strtrim(outinfo.original_memberflag,2)) 
;   ww = where(strmatch(outinfo.memberflag,'*1*'))
    struct_print, struct_trimtags(outinfo[ww],select=['galaxy','cluster',$
      'cluster_fullname','cluster_subname','original_memberflag','memberflag',$
      'cluster_z','z'])

    im_mwrfits, outinfo, infofile, /clobber

return
end
    
    
