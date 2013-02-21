pro phys392_s13_grades, alldata, test=test, sendit=sendit
; jm13feb06siena - parse the grades for this class 

    path = '~/Dropbox/Teaching/392-S13/grades/'
    
    date = '13feb16' ; update this
    semester = 'Spring 2013'
    class = 'Physics 392: Introductory Astrophysics II'
    
; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'phys392_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)
    nstudent = n_elements(data)

;; as a special case for this class alone, assign everybody 100% on
;; their participation grade
;    assign = [assign,'Participation']
;    part = replicate({participation: 100.0, participation_details: 'APOD Discussion', $
;      participation_date: 'various dates', participation_points: 100.0},nstudent)
;    data = struct_addtags(data,part)
    
; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Computational','Quizzes','Presentations','Final Exam']
    weight = [0.3,0.25,0.25,0.2]

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata
    struct_print, alldata

return
end
