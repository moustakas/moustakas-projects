pro trgb_keck_redux, objname, ni, nv, seq_num0, night, datafile


IF (night eq 1) THEN BEGIN
   RESTORE, filename = '/deep3/marc/trgb/data/22dec97/drk-flt-msk.dat'
   LRIS_PROC, ni, seq_num0, drk_oa2, siflt_om1, imgodf, header, /FITS, RVAL = 1.08
   LRIS_PROC, nv, seq_num0+ni, drk_oa2, tvflt_om1, imgodf, header, /FITS, RVAL = 1.08
ENDIF

IF (night eq 2) THEN BEGIN
    RESTORE,filename = '/deep3/marc/trgb/data/23dec97/drk-flt-msk.dat'
   LRIS_PROC, ni, seq_num0, dark, siflt_om2, imgodf, header, /FITS, RVAL = 1.075
   LRIS_PROC, nv, seq_num0+ni, dark, tvflt_om2, imgodf, header, /FITS, RVAL = 1.075
ENDIF


iname = objname+'_I'
iname = STRCOMPRESS(iname,/REMOVE_ALL)
READCUBE, ni, iname, seq_num0, icube, iheadcube

vname = objname+'_V'
vname = STRCOMPRESS(vname,/REMOVE_ALL)
READCUBE, nv, vname, seq_num0+ni, vcube, vheadcube

csize = SIZE(icube)
cube = FLTARR(csize[1], csize[2], ni+nv)
cube[*, *, 0:ni-1] = icube
cube[*, *, ni+nv-1] = vcube

hsize = SIZE(iheadcube)
header =  STRARR(hsize[1], ni+nv)
header[*, 0:ni-1] = iheadcube
header[*, ni+nv-1] =  vheadcube


IMAGE_ALIGN,cube,datafile,mask,out,outmask,/FITS,seq_num0,header


  return
end
