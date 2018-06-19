function R = tetrOrientations(aC, aT, cT, Ui, Uj, Rtwin, Rhabit)

%% variant i
m1 = Rhabit * Rtwin * Ui * aC*[1;0;0];
m2 = Rhabit * Rtwin * Ui * aC*[0;1;0];
m3 = Rhabit * Rtwin * Ui * aC*[0;0;1];

M = [m1 m2 m3];

% Fix order of lattice vectors so that order is m1, m2, m3 and a = b < c
CHECKS = [norm(m1) norm(m2) norm(m3);
    ang(m1,m2)*180/pi ang(m2,m3)*180/pi ang(m3,m1)*180/pi];
[temp ind] = sort(CHECKS(1,:),2);
m1 = M(:,ind(1));
m2 = M(:,ind(2));
m3 = M(:,ind(3));

% Calculate R via polar decomposition
M_unrotated = [aT 0 0;
    0 aT 0;
    0 0 cT];
M_rotated = [m1 m2 m3];
F = M_rotated / M_unrotated;
[R(:,:,1), Utemp] = polardecomp(F);


%% variant j
m1 = Rhabit * Uj * aC*[1;0;0];
m2 = Rhabit * Uj * aC*[0;1;0];
m3 = Rhabit * Uj * aC*[0;0;1];

M = [m1 m2 m3];

% Fix order of lattice vectors so that order is m1, m2, m3 and a = b < c
CHECKS = [norm(m1) norm(m2) norm(m3);
    ang(m1,m2)*180/pi ang(m2,m3)*180/pi ang(m3,m1)*180/pi];
[temp ind] = sort(CHECKS(1,:),2);
m1 = M(:,ind(1));
m2 = M(:,ind(2));
m3 = M(:,ind(3));

% Calculate R via polar decomposition
M_unrotated = [aT 0 0;
    0 aT 0;
    0 0 cT];
M_rotated = [m1 m2 m3];
F = M_rotated / M_unrotated;
[R(:,:,2), Utemp] = polardecomp(F);
