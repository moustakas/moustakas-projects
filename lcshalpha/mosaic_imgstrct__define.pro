pro mosaic_imgstrct__define
; define the imaging inspection structure
    tmp = {mosaic_imgstrct, $
      frame: 0,   $             ; FRAME Number
      flg_anly: 0,$             ; Analysis flag 0=Don't Analyse, 2=bias sub, 4=scatt light
      obj: '', $               ; Object Name
      naxis1: 0, $              ; number of columns
      naxis2: 0, $              ; number of rows
      type: '',   $            ; ObjTyp: OBJ,STD,DRK,ZRO,DFLT
      exp: 0.0,   $             ; Exposure time
      filter: '', $            ; Image Filter
      cbin: 0, $                ; Column bin
      rbin: 0, $                ; Row bin
      am:   0.,   $             ; Airmass
      ccd: '',    $            ; CCD
      tel: '',    $            ; Telescope
      gain: 0.0,   $            ; Gain
      readno: 0.0, $            ; Read Noise
      date: 0.0d,  $            ; JD Date of Obs
      UT: '',     $            ; UT
      RA: 0D,     $             ; RA
      DEC: 0D,    $             ; DEC
      pixscale: 0D, $         ; pixel scale [arcsec/pix]
      equinox: 0.,$            ; EQUINOX
      rootpth: '',$            ; Path of the Root
      img_root: '',$           ; Root name (usually in Raw directory)
      bias_fil: '', $          ; Name of the master Bias image file (fits)
      flat_fil: '', $          ; Name of the master Dome Flat image file (fits)
      twiflat_fil: ''$         ; Name of the master Twilight Flat image file (fits)
      }
return
end
  
         
