function [Q U]=polardecomp(F)

C=transpose(F)*F;
[evec eval]=eig(C); eval=sort([eval(1,1) eval(2,2) eval(3,3)],'descend');
%C=eval(1,1)*dyad(evec(:,1),evec(:,1))+eval(2,2)*dyad(evec(:,2),evec(:,2))+eval(3,3)*dyad(evec(:,3),evec(:,3))
%U=sqrt(C);
U=[sqrt(eval(1)) 0 0;0 sqrt(eval(2)) 0; 0 0 sqrt(eval(3))];
Q=F/U;