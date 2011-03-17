
PRO aitoff_position, POSITION=p, XMARGIN=XMargin, YMARGIN=YMargin

IF n_elements( XMargin ) EQ 0 THEN XMargin = !X.margin
IF n_elements( YMargin ) EQ 0 THEN YMargin = !Y.margin
IF n_elements( XMargin ) EQ 1 THEN XMargin = [XMargin, XMargin]
IF n_elements( YMargin ) EQ 1 THEN YMargin = [YMargin, YMargin]

; --Total number of plots on page
ncols = !P.multi[1]
nrows = !P.multi[2]
IF ncols LE 0 THEN ncols = 1
IF nrows LE 0 THEN nrows = 1
nplots = ncols * nrows

; --Number of plots remaining
nrem = !P.multi[0]
IF nrem LE 0 THEN nrem = nplots

; --(row,col) index of next plot, bottom left is (1,1)
row = (nrem-1) / ncols + 1
col = nrem MOD ncols
IF col EQ 0 THEN col = ncols
col = ncols - col + 1

dx = !D.x_size / float(ncols)
dy = !D.y_size / float(nrows)

XMarginp = XMargin * !D.x_ch_size / float(ncols)
YMarginp = YMargin * !D.y_ch_size / float(nrows)

XMid = (col - 0.5)*dx + (XMarginp[0] - XMarginp[1]) * 0.5
YMid = (row - 0.5)*dy + (YMarginp[0] - YMarginp[1]) * 0.5

dxp = dx - total(XMarginp)
dyp = dy - total(YMarginp)

IF dyp/dxp LE 0.5 THEN BEGIN
    dxp = 2*dyp
ENDIF ELSE BEGIN
    dyp = dxp/2.
ENDELSE

x0 = (XMid - dxp/2.) / float(!d.x_size)
x1 = (XMid + dxp/2.) / float(!d.x_size) 
y0 = (YMid - dyp/2.) / float(!d.y_size)
y1 = (YMid + dyp/2.) / float(!d.y_size)

p = [x0,y0,x1,y1]

END
