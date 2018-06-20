% Conversion from cubic to tetragonal
%
% INPUTS: cubic and tetragonal lattice parameters, cubic orientation, and
% the twin pair for which you wish to know the tetragonal orientations
%
% Outputs: tetragonal orientations for each correspondence variant (CV) in
% the twin pair of interest
% For each twin pair i-j, there are two types of twins, and for each type
% of twin, there are two types of habit planes (cubic/tetragonal phase
% interfaces). So there are 4 solutions per twin pair. This code outputs
% the orientation for both i and j for each of the 4 solutions, totaling 8
% tetragonal orientations. 
%
% Requires the files tetrOrientations.m, cubic_to_tetr.m, polardecomp.m,
% and unit.m to be in the file path


%% INPUTS
% Cubic lattice parameter
aC = 3.25;

% Tetragonal lattice parameters
aT = 3;  cT = 4;

% Cubic orientation
R_C = [1 0 0;
    0 1 0;
    0 0 1];

TwinPairOfInterest = 1;


%% Transformation matrices
alpha = aT / aC;  beta = cT / aC;
U(:,:,1) = [beta  0 0;  0 alpha 0;  0 0 alpha];
U(:,:,2) = [alpha 0 0;  0 beta  0;  0 0 alpha];
U(:,:,3) = [alpha 0 0;  0 alpha 0;  0 0 beta ];


%% Twinning equation solutions
CVPairs = [1 2;
    2 1;
    2 3;
    3 2;
    1 3;
    3 1];

for ii = TwinPairOfInterest
    % n: twin plane norms
    % m: habit plane norms
    [nI, nII, RtwinI, RtwinII, mplusI, mminusI, mplusII, mminusII, RhabitplusI, RhabitminusI, RhabitplusII, RhabitminusII] ...
        = cubic_to_tetr(aC, aT, cT, CVPairs(ii,1), CVPairs(ii,2));
    
    % Twin type I, Habit solution (+)
    RTetrTmp = tetrOrientations(aC, aT, cT, U(:,:,CVPairs(ii,1)), U(:,:,CVPairs(ii,2)), RtwinI, RhabitplusI);
    for jj = 1 : 2
        % Rotate by cubic orientation
        RTetr(:,:,jj) = R_C * RTetrTmp(:,:,jj);
    end

    % Twin type I, Habit solution (-)
    RTetrTmp = tetrOrientations(aC, aT, cT, U(:,:,CVPairs(ii,1)), U(:,:,CVPairs(ii,2)), RtwinI, RhabitminusI);
    for jj = 1 : 2
        % Rotate by cubic orientation
        RTetr(:,:,jj+2) = R_C * RTetrTmp(:,:,jj);
    end    
    
    % Twin type II, Habit solution (+)
    RTetrTmp = tetrOrientations(aC, aT, cT, U(:,:,CVPairs(ii,1)), U(:,:,CVPairs(ii,2)), RtwinII, RhabitplusII);
    for jj = 1 : 2
        % Rotate by cubic orientation
        RTetr(:,:,jj+4) = R_C * RTetrTmp(:,:,jj);
    end    
    
    % Twin type II, Habit solution (-)
    RTetrTmp = tetrOrientations(aC, aT, cT, U(:,:,CVPairs(ii,1)), U(:,:,CVPairs(ii,2)), RtwinII, RhabitminusII);
    for jj = 1 : 2
        % Rotate by cubic orientation
        RTetr(:,:,jj+6) = R_C * RTetrTmp(:,:,jj);
    end        
end