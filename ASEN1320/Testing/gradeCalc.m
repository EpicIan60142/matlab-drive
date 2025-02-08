% I used sample grades here for privacy. You should have more inputs than this,
% but it shouldn't matter how many values you have in the vector.
% So just put in as many as you have.

quizGrades = [15 15 15 15 13.5 12 14 14 14 14 14]; % exclude quiz 0
recitationGrades = [10 10 10 10 10 10 10 10 10 10 10 10]; % use as many as are graded
hwGrades = [20 20 20 20 20 20 20 20 20 20]; % exclude hw 0
interviewGrades = [100 100 100]; % use as many as are graded

quizAvg = sum(quizGrades)/(length(quizGrades)*15);
recitationAvg = sum(recitationGrades)/(length(recitationGrades)*10);
hwAvg = sum(hwGrades)/(length(hwGrades)*20);
interviewAvg = sum(interviewGrades)/(length(interviewGrades)*100);

weightQuiz = .15*quizAvg;
weightRecitation = .1*recitationAvg;
weightHW = .3*hwAvg;
weightInterview = .25*interviewAvg;

currentGrade = 100*(weightQuiz + weightRecitation + weightHW + weightInterview)/.8;
fprintf('Your current grade (excluding final project) is %3.2f%%.\n', currentGrade)