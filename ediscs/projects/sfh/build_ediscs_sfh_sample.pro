pro build_ediscs_sfh_sample, clobber=clobber
; jm10may03ucsd - build the EDisCS/SFH samples; use the smoothed
; indices measured by EDISCS_SFH_INDICES

; output file names    
    sfhpath = ediscs_path(/projects)+'sfh/'
    fieldfile = sfhpath+'ediscs_sfh_field.fits'
    clusterfile = sfhpath+'ediscs_sfh_cluster.fits'

    allfieldfile = repstr(fieldfile,'.fits','_all.fits')
    allclusterfile = repstr(clusterfile,'.fits','_all.fits')

; read the data
    ppxf = read_ediscs(/ppxf)
    kcorr = read_ediscs(/kcorr)
    info = read_ediscs(/spec1d)
    indices = mrdfits(sfhpath+'ediscs_sfh_smooth_indices.fits.gz',1)

; define the parent sample; exclude "false" clusters and clusters
; whose spectroscopy was likely biased (see Sec 5.1 in Bo's paper) 
    parent = where($
      (strtrim(info.cluster,2)+strtrim(info.cluster_subname,2) ne 'cl1119') and $
      (strtrim(info.cluster,2)+strtrim(info.cluster_subname,2) ne 'cl1238') and $
      (strtrim(info.cluster,2)+strtrim(info.cluster_subname,2) ne 'cl1037a') and $ ; z=0.43
      (strtrim(info.cluster,2)+strtrim(info.cluster_subname,2) ne 'cl1103') and $ ; z=0.96
      (ppxf.continuum_snr gt 2.0) and $
      (info.minwave/(1.0+info.z) le 3700.0) and $ ; blue side of [OII]
      (info.maxwave/(1.0+info.z) ge 4370.0) and $ ; red side of Hg_A
      (strmatch(info.galaxy,'*name*',/fold) eq 0),ngal)
    splog, 'Full parent sample = ', ngal

    sample = ppxf[parent]
    sample = struct_addtags(struct_trimtags(sample,$
      except=['indices','d4000*','lick*','dbalmer*']),$
      struct_trimtags(indices[parent],except=['ediscs_id']))
    sample = struct_addtags(temporary(sample),$
      struct_trimtags(info[parent],except=['ediscs_id',$
      'galaxy','specfile','ra','dec','z','cluster']))
    sample = struct_addtags(temporary(sample),$
      struct_trimtags(kcorr[parent],$
      except=['ediscs_id','galaxy','z']))
    
; split into cluster and field samples
    cluster = where(strtrim(sample.cluster_fullname,2) ne '' and $
      strmatch(sample.memberflag,'*1*'),ncluster)
    splog, 'N(cluster) = ', ncluster

    allcl = strtrim(sample[cluster].cluster_fullname,2)
    uindx = uniq(allcl,sort(allcl))
    cl = allcl[uindx]
    target = strtrim(sample[cluster[uindx]].cluster,2) ; targeted cluster
    zclust = sample[cluster[uindx]].cluster_z
    niceprint, target, zclust, cl

    delvarx, field
    for ii = 0, n_elements(zclust)-1 do begin
       field1 = where((strtrim(sample.cluster,2) eq target[ii]) and $
         (strmatch(sample.memberflag,'*1*') eq 0) and $
         (abs(sample.z-zclust[ii]) lt 0.2),nfield1)
;      print, sample[field1].memberflag & stop
       if (nfield1 eq 0) then splog, 'No field galaxies for '+$
         'cluster '+cl[ii]+' (z='+strtrim(zclust[ii],2)+')' else begin
          if (n_elements(field) eq 0) then field = field1 else $
            field = [field,field1]
       endelse
    endfor
    field = field[uniq(field,sort(field))] ; remove duplicates
    nfield = n_elements(field)
    splog, 'N(field) = ', nfield

; write out the full samples
    allfield_sample = sample[field]
    allcluster_sample = sample[cluster]
    im_mwrfits, allfield_sample, allfieldfile, clobber=clobber
    im_mwrfits, allcluster_sample, allclusterfile, clobber=clobber

; now also write out "cleaned" samples where crummy spectra have been
; rejected
;   toss = 
    field_sample = allfield_sample
    cluster_sample = allcluster_sample
    im_mwrfits, field_sample, fieldfile, clobber=clobber
    im_mwrfits, cluster_sample, clusterfile, clobber=clobber
    
return
end
