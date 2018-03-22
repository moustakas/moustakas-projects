pro astr330_s18_grades, alldata, test=test, sendit=sendit, final=final
; jm18feb05siena - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/330-S18/grades/'
    
    date = '18mar05' ; update this
    semester = 'Spring 2018'
    class = 'ASTRO330 - Astrophysics Seminar I'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'astr330_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Participation','Written Summary','Journal Club']
    weight = [0.20,0.40,0.40]

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final, droplowest=droplowest
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
