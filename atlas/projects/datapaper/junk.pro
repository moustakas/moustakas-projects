pro junk

    arm_plotconfig, psfile='junk.eps', /landscape, xpage=8.5, ypage=7.55, coord=pos
    plot, findgen(10), position=pos[*,0]
    dfpsclose

return
end
