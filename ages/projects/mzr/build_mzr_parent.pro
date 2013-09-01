pro build_mzr_parent, sdss=sdss, clobber=clobber
; jm13aug28siena - build the AGES and SDSS MZR parent samples; based
;   on MZ-2011 code

    mzrpath = mzr_path()
    isedpath = mz_path(/isedfit)
    vagcpath = getenv('VAGC_REDUX')+'/'

    if keyword_set(sdss) then begin
       splog, '#########################'
       splog, 'Building the SDSS comparison sample'

; POSTSTR/35: -17>Mr>-24; 0.033<z<0.25; note that if POSTSTR changes
; then SDSS_ISEDFIT needs to be rerun
; 
; area=2.1187566 sr or 6955.47 deg^2 (see
; $VAGC_LSS/bsafe/35/lss/README.dr72bsafe35) 
;
; total volume = jhnvol(0.033,0.25)*6955.47*3600.0 = 7.23E+08 Mpc^3  
       vagc = mz_get_vagc(sample=sample,letter=letter,$
         poststr=poststr,zminmax=zminmax)
       post = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/postlss)
       massoh = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/mpamassoh)

       maggies = mz_get_maggies(post,/sdss,ivarmaggies=ivarmaggies, $
         filterlist=filterlist)

; apply the window
       windowfile = mzrpath+'dr72bsafe35.ply'
       keep = where(im_is_in_window(windowfile,ra=post.ra,$
         dec=post.dec,polyid=allpolyid) and (post.z ge zminmax[0]) and $
         (post.z le zminmax[1]),ngal)
       polyid = allpolyid[keep]

       post = post[keep]
       massoh = massoh[keep]
       
       maggies = maggies[*,keep]
       ivarmaggies = ivarmaggies[*,keep]

; read and write out the stellar masses; just use the fiducial grid (01) 
       ised = mrdfits(isedpath+'sdss_bc03_chab_charlot_sfhgrid01.fits.gz',1,row=keep)
       im_mwrfits, ised, mzrpath+'sdss_parent_mass.fits', /clobber

; add fgot as the statistical weight
       post = struct_addtags(temporary(post),replicate({vagc_object_position: 0L, $
         final_weight: 1.0},ngal))
       read_mangle_polygons, windowfile, win
       post.vagc_object_position = keep
       post.final_weight = 1.0/win[polyid].weight
       
; build the final catalog by keeping just the tags we want
       sample = struct_trimtags(post,select=['object_position','ra','dec',$
         'z','vagc_object_position','final_weight'])
       sample = struct_addtags(temporary(sample),struct_trimtags(massoh,$
         except=['ra','dec','z','oh_entropy']))
       sample = struct_addtags(temporary(sample),struct_trimtags(sdssphot[keep],$
         select=['modelflux*','petroflux*','fiberflux*','extinction','petror*']))

;      sample = struct_addtags(temporary(sample),$
;        im_struct_trimtags(vmax_noevol,$
;        select=['zmin_local','zmax_local','vmax'],$
;        newtags=['zmin','zmax','vmax']+'_noevol'))
;      sample = struct_addtags(temporary(sample),$
;        im_struct_trimtags(vmax_evol,$
;        select=['zmin_local','zmax_local','vmax'],$
;        newtags=['zmin','zmax','vmax']+'_evol'))
;      sample.vmax_noevol = sample.vmax_noevol/h100^3.0 ; h=1-->h=0.7
;      sample.vmax_evol = sample.vmax_evol/h100^3.0 ; h=1-->h=0.7
       
; compute K-corrections
       splog, 'Computing K-corrections'
       kcorr = mz_kcorrect(sample.z,maggies,ivarmaggies,filterlist=filterlist)
       sample = struct_addtags(temporary(sample),kcorr)

       im_mwrfits, sample, mzrpath+'sdss_parent.fits', /clobber

    endif else begin
; area = 8.013 deg^2 (i.e., print, ages_survey_area());
; total volume = jhnvol(0.05,0.75)*8.0135*3600.0 = 1.51329E+07 Mpc^3
; at z=0.05-0.15 the total volume is = jhnvol(0.05,0.15)*8.0135*3600.0
; = 1.86578E+05 Mpc^3

       splog, '#########################'
       splog, 'Building the AGES sample'

       phot = read_ages(/photo)
       ppxf = read_ages(/ppxf)

; identify main sample galaxies and apply additional sample cuts:
; 0.05<z<0.75; 15<Itot<19.95; reject unfluxed plates and require that
; the spectrum has been fitted by PPXF
       ifaint = mz_ifaint(ibright=ibright,select_filter=select_filter,$ 
         select_vega2ab=iv2ab,/vega)
       ages_to_maggies, phot, allmaggies, allivarmaggies, filterlist=filterlist

       sample_zmin = 0.05
       sample_zmax = 0.75
       nminphot = 3
       rejplates = [$
         104,$ ; average fluxing vector
         106,110,209,310,311,$ ; unfluxed
         604]  ; very crummy 
       
       ippxf = intarr(n_elements(phot))
       ippxf[ppxf.ages_id] = 1 ; fitted

       iparent = (phot.imain eq 1) and (phot.i_tot le ifaint) and (phot.i_tot ge ibright) and $
         (ippxf eq 1) and (phot.z ge sample_zmin) and (phot.z le sample_zmax) and $
         (total(allmaggies gt 0,1) ge nminphot) and $
         (phot.pass ne rejplates[0]) and (phot.pass ne rejplates[1]) and $
         (phot.pass ne rejplates[2]) and (phot.pass ne rejplates[3]) and $
         (phot.pass ne rejplates[4]) and (phot.pass ne rejplates[5]) and $
         (phot.pass ne rejplates[6])
       parent = where(iparent,ngal)
       splog, 'Ngal = ', ngal

       sample = phot[parent]
       maggies = allmaggies[*,parent]
       ivarmaggies = allivarmaggies[*,parent]

; in addition to the corrections for spectroscopic incompleteness
; (~2.1%), fiber incompleteness (~4.3%), and sparse-sampling, we also
; need to account for the catalog incompleteness (catalog_weight) and
; for the plates that were rejected above (field_weight); from
; Eisenstein, the catalog incompleteness is ~3%-6%, so assume an
; average value of 4%; finally, compute the final galaxy weight as the
; product of all the various selection terms        
       moretags = replicate({catalog_weight: 1.04, field_weight: 1.0, $
         final_weight: 1.0},ngal)
       sample = struct_addtags(temporary(sample),moretags)
       
       field_weight = ages_upweight_field(rejplates,field=field)
       for ii = 0, n_elements(field)-1 do begin
          these = where(field[ii] eq sample.field,nthese)
          if (nthese ne 0) then sample[these].field_weight = field_weight[ii]
       endfor
       sample.final_weight = sample.spec_weight*sample.target_weight*$
         sample.fiber_weight*sample.catalog_weight*1.0/sample.field_weight ; note 1/FIELD_WEIGHT!
pick up here       
       
; compute K-corrections
       splog, 'Computing K-corrections'
       kcorr = mz_kcorrect(sample.z,maggies,$
         ivarmaggies,filterlist=filterlist)
       sample = struct_addtags(temporary(sample),kcorr)

; write out       
       im_mwrfits, sample, mzrpath+'ages_parent.fits', /clobber
       im_mwrfits, outised, mzrpath+'ages_parent_mass.fits', /clobber
    endelse 

return
end
    
