% function [nI, aI, shearI, nII, aII, shearII, bplusI, mplusI, bminusI, mminusI, shapestrainI, bplusII, mplusII, bminusII, mminusII, shapestrainII, Ui, Uj, lamI, lamII, Q] = cubic_to_tetr(aC, aT, cT, i, j)
function [nI, nII, RtwinI, RtwinII, mplusI, mminusI, mplusII, mminusII, RhabitplusI, RhabitminusI, RhabitplusII, RhabitminusII] = cubic_to_tetr(aC, aT, cT, i, j)


%% Transformation matrices
alpha = aT / aC;  beta = cT / aC;
U = zeros(3,3,3);
U(:,:,1) = [beta  0 0;  0 alpha 0;  0 0 alpha];
U(:,:,2) = [alpha 0 0;  0 beta  0;  0 0 alpha];
U(:,:,3) = [alpha 0 0;  0 alpha 0;  0 0 beta ];

Ui = U(:,:,i);  Uj = U(:,:,j);


%% Calculate twin plane equation components between variants i and j
% First type of twins
kappa = 1;
C = inv(Uj) * (Ui*Ui) * inv(Uj);

[evec,eval]=eig(C);
temp=[transpose(evec(:,1)) eval(1,1);transpose(evec(:,2)) eval(2,2);transpose(evec(:,3)) eval(3,3)];
[temp1,temp2]=sort(temp(:,4));
temp_new=temp(temp2,:);
lam1=temp_new(1,4);  lam2=temp_new(2,4);  lam3=temp_new(3,4);
e1=transpose(temp_new(1,1:3));  e2=transpose(temp_new(2,1:3));  e3=transpose(temp_new(3,1:3));

rho=norm((sqrt(lam3)-sqrt(lam1))/sqrt(lam3-lam1)*(-sqrt(1-lam1)*Uj*e1+kappa*sqrt(lam3-1)*Uj*e3));
aI=rho*(sqrt(lam3*(1-lam1)/(lam3-lam1))*e1+kappa*sqrt(lam1*(lam3-1)/(lam3-lam1))*e3);
nI=unit((sqrt(lam3)-sqrt(lam1))/sqrt(lam3-lam1)*(-sqrt(1-lam1)*Uj*e1+kappa*sqrt(lam3-1)*Uj*e3));
shearI=sqrt(lam3)-sqrt(lam1);

% Second type of twins
kappa = -1;

rho=norm((sqrt(lam3)-sqrt(lam1))/sqrt(lam3-lam1)*(-sqrt(1-lam1)*Uj*e1+kappa*sqrt(lam3-1)*Uj*e3));
aII=rho*(sqrt(lam3*(1-lam1)/(lam3-lam1))*e1+kappa*sqrt(lam1*(lam3-1)/(lam3-lam1))*e3);
nII=unit((sqrt(lam3)-sqrt(lam1))/sqrt(lam3-lam1)*(-sqrt(1-lam1)*Uj*e1+kappa*sqrt(lam3-1)*Uj*e3));
shearII=sqrt(lam3)-sqrt(lam1);


%% Calculate necessary rotations from twinning on variant i
RtwinI = (aI*transpose(nI) + Uj) * inv(Ui);
RtwinII = (aII*transpose(nII) + Uj) * inv(Ui);


%% Calculate habit plane equation components between cubic and tetragonal twin i-j
% First type of twins
delta = dot(aI, Uj*inv(Uj^2-eye(3))*nI);  % requirement #1 for habit plane
lambdaI = 1/2 * (1 - sqrt(1 + 2/delta));  % twin phase fraction
eta = trace(Uj^2) - det(Uj^2) - 2 + norm(aI)^2/(2*delta);  % requirement #2 for habit plane

if delta<=-2 && eta>=0
    C=transpose(Uj+lambdaI*(aI*transpose(nI)))*(Uj+lambdaI*(aI*transpose(nI)));
    [evec,eval]=eig(C);
    temp=[transpose(evec(:,1)) eval(1,1);transpose(evec(:,2)) eval(2,2);transpose(evec(:,3)) eval(3,3)];
    [temp1,temp2]=sort(temp(:,4));
    temp_new=temp(temp2,:);
    lam1=temp_new(1,4);  lam2=temp_new(2,4);  lam3=temp_new(3,4);
    e1=transpose(temp_new(1,1:3));  e2=transpose(temp_new(2,1:3));  e3=transpose(temp_new(3,1:3));
    
    kappa=1;
    m_temp=(sqrt(lam3)-sqrt(lam1))/sqrt(lam3-lam1)*(-sqrt(1-lam1)*e1+kappa*sqrt(lam3-1)*e3);
    rho=norm(m_temp);
    mplusI=m_temp/rho;
    bplusI=rho*(sqrt(lam3*(1-lam1)/(lam3-lam1))*e1+kappa*sqrt(lam1*(lam3-1)/(lam3-lam1))*e3);
    
    kappa=-1;
    m_temp=(sqrt(lam3)-sqrt(lam1))/sqrt(lam3-lam1)*(-sqrt(1-lam1)*e1+kappa*sqrt(lam3-1)*e3);
    rho=norm(m_temp);
    mminusI=m_temp/rho;
    bminusI=rho*(sqrt(lam3*(1-lam1)/(lam3-lam1))*e1+kappa*sqrt(lam1*(lam3-1)/(lam3-lam1))*e3);
    
else  % if habit plane not possible
    bplusI=zeros(3,1);  mplusI=bplusI;  bminusI=bplusI;  mminusI=bplusI;
end


% Second type of twins
delta = dot(aII, Uj*inv(Uj^2-eye(3))*nII);  % requirement #1 for habit plane
lambdaII = 1/2 * (1 - sqrt(1 + 2/delta));  % twin phase fraction
eta = trace(Uj^2) - det(Uj^2) - 2 + norm(aII)^2/(2*delta);  % requirement #2 for habit plane

if delta<=-2 && eta>=0
    C=transpose(Uj+lambdaII*(aI*transpose(nI)))*(Uj+lambdaII*(aI*transpose(nI)));
    [evec,eval]=eig(C);
    temp=[transpose(evec(:,1)) eval(1,1);transpose(evec(:,2)) eval(2,2);transpose(evec(:,3)) eval(3,3)];
    [temp1,temp2]=sort(temp(:,4));
    temp_new=temp(temp2,:);
    lam1=temp_new(1,4);  lam2=temp_new(2,4);  lam3=temp_new(3,4);
    e1=transpose(temp_new(1,1:3));  e2=transpose(temp_new(2,1:3));  e3=transpose(temp_new(3,1:3));
    
    kappa=1;
    m_temp=(sqrt(lam3)-sqrt(lam1))/sqrt(lam3-lam1)*(-sqrt(1-lam1)*e1+kappa*sqrt(lam3-1)*e3);
    rho=norm(m_temp);
    mplusII=m_temp/rho;
    bplusII=rho*(sqrt(lam3*(1-lam1)/(lam3-lam1))*e1+kappa*sqrt(lam1*(lam3-1)/(lam3-lam1))*e3);
    
    kappa=-1;
    m_temp=(sqrt(lam3)-sqrt(lam1))/sqrt(lam3-lam1)*(-sqrt(1-lam1)*e1+kappa*sqrt(lam3-1)*e3);
    rho=norm(m_temp);
    mminusII=m_temp/rho;
    bminusII=rho*(sqrt(lam3*(1-lam1)/(lam3-lam1))*e1+kappa*sqrt(lam1*(lam3-1)/(lam3-lam1))*e3);
    
else  % if habit plane not possible
    bplusII=zeros(3,1); mplusII=bplusII; bminusII=bplusII; mminusII=bplusII;
end


%% Calculate necessary rotations on i and j from habit plane 
RhabitplusI = ((bplusI*transpose(mplusI)) + eye(3)) / (Uj + lambdaI * dyad(aI,nI));
RhabitminusI = ((bminusI*transpose(mminusI)) + eye(3)) / (Uj + lambdaI * dyad(aI,nI));

RhabitplusII = ((bplusII*transpose(mplusII)) + eye(3)) / (Uj + lambdaII * dyad(aII,nII));
RhabitminusII = ((bminusII*transpose(mminusII)) + eye(3)) / (Uj + lambdaII * dyad(aII,nII));