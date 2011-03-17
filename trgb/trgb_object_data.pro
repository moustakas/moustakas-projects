pro writetxt, s, fname=fname

	openw, lun1, fname, /get_lun
        printf, lun1, '# '+systime()
        printf, lun1, '# -----------------------------------------------------'$
          +'-----------------------------------------------------'
        printf, lun1, '# Comments: '
        printf, lun1, '#   v  : heliocentric velocity (km/s, NED)'
        printf, lun1, '#   z  : heliocentric redshift (km/s, NED)'
        printf, lun1, '#   cz : Local Group redshift (km/s, IRAS predicted)'
        printf, lun1, '#   v_p: Local Group peculiar velocity (km/s, IRAS predicted)'
        printf, lun1, '# -----------------------------------------------------'$
          +'-----------------------------------------------------' 
        printf, lun1, '#        Galaxy      Type     Longitude    Latitude     A_V     A_I    E(V-I)   v'$
          +'       z        cz     v_p  '
        printf, lun1, '# -----------------------------------------------------'$
          +'-----------------------------------------------------'
        printf, lun1, '# '
        for k = 0L, n_elements(s)-1L do $
          printf, lun1, s[k].truename, s[k].type, s[k].longitude, s[k].latitude, $
          s[k].a_v, s[k].a_i, s[k].e_v_i, s[k].v_helio, s[k].z_helio, s[k].cz_lg, $
          s[k].pec_vel_pred, format='(A16,2x,A8,1x,2F12.5,1x,3F8.3,1x,F6.1,1x,F9.6,1x,2F7.1)'
        free_lun, lun1

return
end

pro trgb_object_data
;+
; NAME:
;	TRGB_OBJECT_DATA
;
; PURPOSE:
;	Create and write out an array of data structures and a text
;	file of information for each Keck and HST object. 
;
; OUTPUTS:
;	Writes two idl save sets to /deep1/ioannis/trgb/results:
;	hst_object.dat and  keck_object.dat.  The save sets are
;	restored by TRGB_READATA.  The routine also creates text files
;	of the information in the structure in the same directory as
;	above, but with the suffix .txt.
;
; PROCEDURES USED:
;	SWRITE()
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 20, UCB
;	jm00sep8uofa, added WRITETXT subroutine
;-

; ----------------------------------------------------------------------------------------------------
; master template structure
; ----------------------------------------------------------------------------------------------------

	template = {name:	'Object Data Structure', $
                    object:		strarr(1), $	; object name
                    truename:		strarr(1), $	; object name for plots
                    color_info:		strarr(1), $	; is there color information? (Yes/No)
; ----------------------------------------------------------------------------------------------------
; NED database information
; ----------------------------------------------------------------------------------------------------
		    type:		strarr(1), $	; morphological type
                    longitude:		fltarr(1), $	; Galactic longitude (2000)
                    latitude:		fltarr(1), $	; Galactic latitude (2000)
                    a_v:		fltarr(1), $	; V-band extinction
                    a_i:		fltarr(1), $	; I-band extinction
                    e_v_i:		fltarr(1), $ 	; E(V-I) color excess
                    v_helio:		fltarr(1), $	; heliocentric radial velocity (km/s)
                    z_helio:		dblarr(1), $	; heliocentric redshift
; ----------------------------------------------------------------------------------------------------
; IRAS map information
; ----------------------------------------------------------------------------------------------------
                    cz_lg:		fltarr(1), $	; local group redshift * 3x10^5 (km/s)
                    pec_vel_pred:	fltarr(1)}	; IRAS predicted local group peculiar velocity (km/s)

        light = 2.99E5	; speed of light (km/s)

        nhst = 18L      ; total number of HST objects
        nkeck = 12L	; total number of KECK objects

	hst_object = replicate(template,nhst)
	keck_object = replicate(template,nkeck)

; HST data:

        hst_object[0].object = 'ugc07577'
        hst_object[0].truename = 'UGC 07577'
        hst_object[0].color_info = 'NO'
        hst_object[0].type = 'Im'
        hst_object[0].longitude = 137.7509
        hst_object[0].latitude = 72.9451
        hst_object[0].a_v = 0.068
        hst_object[0].a_i = 0.040
        hst_object[0].e_v_i = 0.028
        hst_object[0].v_helio = 196.
        hst_object[0].z_helio = 0.000654
        hst_object[0].cz_lg = 290.
        hst_object[0].pec_vel_pred = 36.
        
        hst_object[1].object = 'ugc03476'
        hst_object[1].truename = 'UGC 03476'
        hst_object[1].color_info = 'YES'
        hst_object[1].type = 'Im'
        hst_object[1].longitude = 180.7657
        hst_object[1].latitude = 10.5128
        hst_object[1].a_v = 0.784
        hst_object[1].a_i = 0.459
        hst_object[1].e_v_i = 0.325
        hst_object[1].v_helio = 469.
        hst_object[1].z_helio = 0.001564
        hst_object[1].cz_lg = 577.
        hst_object[1].pec_vel_pred = -246.

        hst_object[2].object = 'ugc03698'
        hst_object[2].truename = 'UGC 03698'
        hst_object[2].color_info = 'NO'
        hst_object[2].type = 'Im'
        hst_object[2].longitude = 172.9687
        hst_object[2].latitude = 21.6196
        hst_object[2].a_v = 0.325
        hst_object[2].a_i = 0.190
        hst_object[2].e_v_i = 0.135
        hst_object[2].v_helio = 420.
        hst_object[2].z_helio = 0.001401
        hst_object[2].cz_lg = 592.
        hst_object[2].pec_vel_pred = -222.

        hst_object[3].object = 'ugc03755'
        hst_object[3].truename = 'UGC 03755'
        hst_object[3].color_info = 'YES'
        hst_object[3].type = 'Im'
        hst_object[3].longitude = 206.0171
        hst_object[3].latitude = 9.7140
        hst_object[3].a_v = 0.295
        hst_object[3].a_i = 0.172
        hst_object[3].e_v_i = 0.123
        hst_object[3].v_helio = 314.
        hst_object[3].z_helio = 0.001047
        hst_object[3].cz_lg = 163.
        hst_object[3].pec_vel_pred = -89.

        hst_object[4].object = 'ugc03860'
        hst_object[4].truename = 'UGC 03860'
        hst_object[4].color_info = 'NO'
        hst_object[4].type = 'Im'
        hst_object[4].longitude = 177.8063
        hst_object[4].latitude = 18.6608
        hst_object[4].a_v = 0.194
        hst_object[4].a_i = 0.114
        hst_object[4].e_v_i = 0.080
        hst_object[4].v_helio = 354.
        hst_object[4].z_helio = 0.001181
        hst_object[4].cz_lg = 0.
        hst_object[4].pec_vel_pred = 0.

        hst_object[5].object = 'ugc03974'
        hst_object[5].truename = 'UGC 03974'
        hst_object[5].color_info = 'NO'
        hst_object[5].type = 'IB(s)m'
        hst_object[5].longitude = 203.1003
        hst_object[5].latitude = 18.5408
        hst_object[5].a_v = 0.111
        hst_object[5].a_i = 0.065
        hst_object[5].e_v_i = 0.046
        hst_object[5].v_helio = 272.
        hst_object[5].z_helio = 0.000907
        hst_object[5].cz_lg = 138.
        hst_object[5].pec_vel_pred = -56.

        hst_object[6].object = 'ugc04115'
        hst_object[6].truename = 'UGC 04115'
        hst_object[6].color_info = 'NO'
        hst_object[6].type = 'IAm'
        hst_object[6].longitude = 207.0078
        hst_object[6].latitude = 20.8981
        hst_object[6].a_v = 0.143
        hst_object[6].a_i = 0.084
        hst_object[6].e_v_i = 0.059
        hst_object[6].v_helio = 338.
        hst_object[6].z_helio = 0.001127
        hst_object[6].cz_lg = 0.
        hst_object[6].pec_vel_pred = 0.

        hst_object[7].object = 'ngc1313'
        hst_object[7].truename = 'NGC 1313'
        hst_object[7].color_info = 'NO'
        hst_object[7].type = 'SB(s)d'
        hst_object[7].longitude = 283.3588
        hst_object[7].latitude = -44.6442
        hst_object[7].a_v = 0.362
        hst_object[7].a_i = 0.212
        hst_object[7].e_v_i = 0.140
        hst_object[7].v_helio = 475
        hst_object[7].z_helio = 0.001584
        hst_object[7].cz_lg = 71.
        hst_object[7].pec_vel_pred = -11.

        hst_object[8].object = 'ngc2683'
        hst_object[8].truename = 'NGC 2683'
        hst_object[8].color_info = 'NO'
        hst_object[8].type = 'SA(rs)b'
        hst_object[8].longitude = 190.4565
        hst_object[8].latitude = 38.7609
        hst_object[8].a_v = 0.109
        hst_object[8].a_i = 0.064
        hst_object[8].e_v_i = 0.045
        hst_object[8].v_helio = 411.
        hst_object[8].z_helio = 0.001371
        hst_object[8].cz_lg = 389.
        hst_object[8].pec_vel_pred = -139.

        hst_object[9].object = 'ngc2903'
        hst_object[9].truename = 'NGC 2903'
        hst_object[9].color_info = 'NO'
        hst_object[9].type = 'SB(s)d'
        hst_object[9].longitude = 208.7118
        hst_object[9].latitude = 44.5398
        hst_object[9].a_v = 0.103
        hst_object[9].a_i = 0.060
        hst_object[9].e_v_i = 0.043
        hst_object[9].v_helio = 556.
        hst_object[9].z_helio = 0.001855
        hst_object[9].cz_lg = 391.
        hst_object[9].pec_vel_pred = -111.

        hst_object[10].object = 'ngc6503'
        hst_object[10].truename = 'NGC 6503'
        hst_object[10].color_info = 'NO'
        hst_object[10].type = 'SA(s)cd'
        hst_object[10].longitude = 100.5735
        hst_object[10].latitude = 30.6398
        hst_object[10].a_v = 0.106
        hst_object[10].a_i = 0.062
        hst_object[10].e_v_i = 0.044
        hst_object[10].v_helio = 62.
        hst_object[10].z_helio = 0.000207
        hst_object[10].cz_lg = 543.
        hst_object[10].pec_vel_pred = -113.

        hst_object[11].object = 'ic342'
        hst_object[11].truename = 'IC 342'
        hst_object[11].color_info = 'YES'
        hst_object[11].type = 'SAB(rs)cd'
        hst_object[11].longitude = 138.1730
        hst_object[11].latitude = 10.5803
        hst_object[11].a_v = 1.849
        hst_object[11].a_i = 1.082
        hst_object[11].e_v_i = 0.767
        hst_object[11].v_helio = 31.
        hst_object[11].z_helio = 0.00010
        hst_object[11].cz_lg = 512.
        hst_object[11].pec_vel_pred = -56.

        hst_object[12].object = 'iras1248'
        hst_object[12].truename = 'IRAS 12483-1311'
        hst_object[12].color_info = 'NO'
        hst_object[12].type = 'SB?'
        hst_object[12].longitude = 302.7468
        hst_object[12].latitude = 49.4149
        hst_object[12].a_v = 0.180
        hst_object[12].a_i = 0.106
        hst_object[12].e_v_i = 0.074
        hst_object[12].v_helio = 346.
        hst_object[12].z_helio = 0.00115
        hst_object[12].cz_lg = 0.
        hst_object[12].pec_vel_pred = 0.

        hst_object[13].object = 'ngc5457'
        hst_object[13].truename = 'NGC 5457'
        hst_object[13].color_info = 'YES'
        hst_object[13].type = 'SAB(rs)cd'
        hst_object[13].longitude = 102.0373
        hst_object[13].latitude = 59.7716
        hst_object[13].a_v = 0.028
        hst_object[13].a_i = 0.017
        hst_object[13].e_v_i = 0.011
        hst_object[13].v_helio = 241.
        hst_object[13].z_helio = 0.00080
        hst_object[13].cz_lg = 0.
        hst_object[13].pec_vel_pred = 0.

        hst_object[14].object = 'ngc5128'
        hst_object[14].truename = 'NGC 5128'
        hst_object[14].color_info = 'YES'
        hst_object[14].type = 'S0 pec'
        hst_object[14].longitude = 309.5158
        hst_object[14].latitude = 19.4173
        hst_object[14].a_v = 0.381
        hst_object[14].a_i = 0.223
        hst_object[14].e_v_i = 0.158
        hst_object[14].v_helio = 547.
        hst_object[14].z_helio = 0.00183
        hst_object[14].cz_lg = 0.
        hst_object[14].pec_vel_pred = 0.

        hst_object[15].object = 'ngc3621'
        hst_object[15].truename = 'NGC 3621'
        hst_object[15].color_info = 'YES'
        hst_object[15].type = 'SA(s)d'
        hst_object[15].longitude = 281.2115
        hst_object[15].latitude = 26.1003
        hst_object[15].a_v = 0.266
        hst_object[15].a_i = 0.156
        hst_object[15].e_v_i = 0.110
        hst_object[15].v_helio = 727.
        hst_object[15].z_helio = 0.00243
        hst_object[15].cz_lg = 0.
        hst_object[15].pec_vel_pred = 0.

        hst_object[16].object = 'ngc1705'
        hst_object[16].truename = 'NGC 1705'
        hst_object[16].color_info = 'YES'
        hst_object[16].type = 'SA0-pec'
        hst_object[16].longitude = 261.0788
        hst_object[16].latitude = -38.7427
        hst_object[16].a_v = 0.027
        hst_object[16].a_i = 0.016
        hst_object[16].e_v_i = 0.011
        hst_object[16].v_helio = 628.
        hst_object[16].z_helio = 0.00210
        hst_object[16].cz_lg = 0.
        hst_object[16].pec_vel_pred = 0.

        hst_object[17].object = 'ugc06456'
        hst_object[17].truename = 'UGC 06456'
        hst_object[17].color_info = 'YES'
        hst_object[17].type = 'Pec'
        hst_object[17].longitude = 127.8360
        hst_object[17].latitude = 11.3999
        hst_object[17].a_v = 0.119
        hst_object[17].a_i = 0.070
        hst_object[17].e_v_i = 0.049
        hst_object[17].v_helio = -100.
        hst_object[17].z_helio = -0.00033
        hst_object[17].cz_lg = 0.
        hst_object[17].pec_vel_pred = 0.

; save it:

        swrite, hst_object, '/deep1/ioannis/trgb/results/hst_object.dat'
        writetxt, hst_object, fname = '/deep1/ioannis/trgb/results/hst_object.txt'

; KECK data

        keck_object[0].object = 'SextansB'
        keck_object[0].truename = 'Sextans B'
        keck_object[0].color_info = 'YES'
        keck_object[0].type = 'ImIV-V'
        keck_object[0].longitude = 233.2001
        keck_object[0].latitude = 43.7838
        keck_object[0].a_v = 0.105
        keck_object[0].a_i = 0.062
        keck_object[0].e_v_i = 0.043
        keck_object[0].v_helio = 301.
        keck_object[0].z_helio = 0.001
        keck_object[0].cz_lg = -30.
        keck_object[0].pec_vel_pred = 7.
        
        keck_object[1].object = 'leoi'
        keck_object[1].truename = 'Leo I'
        keck_object[1].color_info = 'YES'
        keck_object[1].type = 'E;dSph'
        keck_object[1].longitude = 225.9823
        keck_object[1].latitude = 49.1104
        keck_object[1].a_v = 0.120
        keck_object[1].a_i = 0.070
        keck_object[1].e_v_i = 0.050
        keck_object[1].v_helio = 168.
        keck_object[1].z_helio = 0.00056
        keck_object[1].cz_lg = -99.
        keck_object[1].pec_vel_pred = 19.

        keck_object[2].object = 'NGC2903'
        keck_object[2].truename = 'NGC 2903'
        keck_object[2].color_info = 'YES'
        keck_object[2].type = 'SB(s)d'
        keck_object[2].longitude = 208.7110
        keck_object[2].latitude = 44.5403
        keck_object[2].a_v = 0.103
        keck_object[2].a_i = 0.060
        keck_object[2].e_v_i = 0.043
        keck_object[2].v_helio = 556.
        keck_object[2].z_helio = 0.00186
        keck_object[2].cz_lg = 391.
        keck_object[2].pec_vel_pred = -111.

        keck_object[3].object = 'ic2574'
        keck_object[3].truename = 'IC 2574'
        keck_object[3].color_info = 'YES'
        keck_object[3].type = 'SAB(s)m'
        keck_object[3].longitude = 140.2102
        keck_object[3].latitude = 43.6049
        keck_object[3].a_v = 0.120
        keck_object[3].a_i = 0.070
        keck_object[3].e_v_i = 0.050
        keck_object[3].v_helio = 57.
        keck_object[3].z_helio = 0.00019
        keck_object[3].cz_lg = 357.
        keck_object[3].pec_vel_pred = -49.

        keck_object[4].object = 'NGC3109'
        keck_object[4].truename = 'NGC 3109'
        keck_object[4].color_info = 'YES'
        keck_object[4].type = 'SB(s)m'
        keck_object[4].longitude = 262.1006
        keck_object[4].latitude = 23.0701
        keck_object[4].a_v = 0.221
        keck_object[4].a_i = 0.129
        keck_object[4].e_v_i = 0.092
        keck_object[4].v_helio = 403.
        keck_object[4].z_helio = 0.00134
        keck_object[4].cz_lg = -150.
        keck_object[4].pec_vel_pred = 31.

        keck_object[5].object = 'I342'
        keck_object[5].truename = 'IC 342'
        keck_object[5].color_info = 'YES'
        keck_object[5].type = 'SAB(rs)cd'
        keck_object[5].longitude = 138.1730
        keck_object[5].latitude = 10.5803
        keck_object[5].a_v = 1.849
        keck_object[5].a_i = 1.082
        keck_object[5].e_v_i = 0.767
        keck_object[5].v_helio = 31.
        keck_object[5].z_helio = 0.00010
        keck_object[5].cz_lg = 512.
        keck_object[5].pec_vel_pred = -56.

        keck_object[6].object = 'NGC1560'
        keck_object[6].truename = 'NGC 1560'
        keck_object[6].color_info = 'YES'
        keck_object[6].type = 'SA(s)d'
        keck_object[6].longitude = 138.3669
        keck_object[6].latitude = 16.0221
        keck_object[6].a_v = 0.624
        keck_object[6].a_i = 0.365
        keck_object[6].e_v_i = 0.259
        keck_object[6].v_helio = -36.
        keck_object[6].z_helio = -0.00012
        keck_object[6].cz_lg = 420.
        keck_object[6].pec_vel_pred = -67.

        keck_object[7].object = 'NGC2366'
        keck_object[7].truename = 'NGC 2366'
        keck_object[7].color_info = 'YES'
        keck_object[7].type = 'IB(s)m'
        keck_object[7].longitude = 146.4192
        keck_object[7].latitude = 28.5350
        keck_object[7].a_v = 0.120
        keck_object[7].a_i = 0.070
        keck_object[7].e_v_i = 0.050
        keck_object[7].v_helio = 100.
        keck_object[7].z_helio = 0.00033
        keck_object[7].cz_lg = 377.
        keck_object[7].pec_vel_pred = -80.

        keck_object[8].object = 'HOLMBERGI'
        keck_object[8].truename = 'Holmberg I'
        keck_object[8].color_info = 'YES'
        keck_object[8].type = 'IAB(s)m'
        keck_object[8].longitude = 140.7294
        keck_object[8].latitude = 38.6620
        keck_object[8].a_v = 0.159
        keck_object[8].a_i = 0.093
        keck_object[8].e_v_i = 0.066
        keck_object[8].v_helio = 143.
        keck_object[8].z_helio = 0.00048
        keck_object[8].cz_lg = 480.
        keck_object[8].pec_vel_pred = -83.

        keck_object[9].object = 'holmbergii'
        keck_object[9].truename = 'Holmberg II'
        keck_object[9].color_info = 'YES'
        keck_object[9].type = 'Im'
        keck_object[9].longitude = 144.2833
        keck_object[9].latitude = 32.6910
        keck_object[9].a_v = 0.107
        keck_object[9].a_i = 0.063
        keck_object[9].e_v_i = 0.044
        keck_object[9].v_helio = 157.
        keck_object[9].z_helio = 0.00052
        keck_object[9].cz_lg = 511.
        keck_object[9].pec_vel_pred = -97.

        keck_object[10].object = 'HolmbergIX'
        keck_object[10].truename = 'Holmberg IX'
        keck_object[10].color_info = 'YES'
        keck_object[10].type = 'Im'
        keck_object[10].longitude = 141.9769
        keck_object[10].latitude = 41.0521
        keck_object[10].a_v = 0.263
        keck_object[10].a_i = 0.154
        keck_object[10].e_v_i = 0.109
        keck_object[10].v_helio = 46.
        keck_object[10].z_helio = 0.00015
        keck_object[10].cz_lg = 364.
        keck_object[10].pec_vel_pred = -56.
        
        keck_object[11].object = 'NGC2976'
        keck_object[11].truename = 'NGC 2976'
        keck_object[11].color_info = 'YES'
        keck_object[11].type = 'SAc pec'
        keck_object[11].longitude = 143.9143
        keck_object[11].latitude = 40.9022
        keck_object[11].a_v = 0.230
        keck_object[11].a_i = 0.135
        keck_object[11].e_v_i = 0.095
        keck_object[11].v_helio = 3.
        keck_object[11].z_helio = 0.00001
        keck_object[11].cz_lg = 319.
        keck_object[11].pec_vel_pred = -47.

        swrite, keck_object, '/deep1/ioannis/trgb/results/keck_object.dat'
        writetxt, keck_object, fname = '/deep1/ioannis/trgb/results/keck_object.txt'

return
end

