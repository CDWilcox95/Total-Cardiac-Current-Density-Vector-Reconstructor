function [LL_lung_voxels, LR_lung_voxels, UL_lung_voxels, UR_lung_voxels, ...
          LA_lung_voxels, LP_lung_voxels, UA_lung_voxels, UP_lung_voxels, plot_mesh] ...
          = define_lung_lobes(Jash)

% Jash1 = Jash(1:93,:);
% Jash2 = Jash(94:end,:);

Jash1=Jash(1:303,:);
Jash2=Jash(304:end,:);

Jash1=flip(Jash1);

Jash2=flip(Jash2);

height = 39;
bar_color = NaN(height,size(Jash,2));
plot_mesh = [bar_color ; Jash];

%Josh mesh segmentation Left/right
% horizontal_mid_line = 125; % conformal map
horizontal_mid_line=151;

LL_lung_voxels = unique(Jash1(:,horizontal_mid_line+1:end));  % Lower left
LR_lung_voxels = unique(Jash1(:,1:horizontal_mid_line)); % Lower right
% JM changed it so there are no voxels common to both so volumes are okay
In_LL_and_LR = intersect(LL_lung_voxels,LR_lung_voxels);
Remove_from_LL=In_LL_and_LR(1:2:end);  % Remove odd-indexed entries from LL
Remove_from_LR=In_LL_and_LR(2:2:end); % Remove even-indexed entries from LR
[~,indL] = intersect(LL_lung_voxels,Remove_from_LL);  % ind gives the indices in LL_lung_voxels of what we want to remove
LL_lung_voxels(indL) = [];  % This removes them
[~,indR] = intersect(LR_lung_voxels,Remove_from_LR);  % ind gives the indices in LL_lung_voxels of what we want to remove
LR_lung_voxels(indR) = [];  % This removes them
% Check the intersection now
% in_both = intersect(LL_lung_voxels,LR_lung_voxels);
% size(In_LL_and_LR)
% size(in_both)
% intersect(In_LL_and_LR,in_both)

UL_lung_voxels = unique(Jash2(:,horizontal_mid_line+1:end));
UR_lung_voxels = unique(Jash2(:,1:horizontal_mid_line));
% JM changed it so there are no voxels common to both so volumes are okay
In_UL_and_UR = intersect(UL_lung_voxels,UR_lung_voxels);
Remove_from_UL=In_UL_and_UR(1:2:end);  % Remove odd-indexed entries from LL
Remove_from_UR=In_UL_and_UR(2:2:end); % Remove even-indexed entries from LR
[~,indL] = intersect(UL_lung_voxels,Remove_from_UL);  % ind gives the indices in LL_lung_voxels of what we want to remove
UL_lung_voxels(indL) = [];  % This removes them
[~,indR] = intersect(UR_lung_voxels,Remove_from_UR);  % ind gives the indices in LL_lung_voxels of what we want to remove
UR_lung_voxels(indR) = [];  % This removes them
% Check the intersection now
% in_both = intersect(UL_lung_voxels,UR_lung_voxels);
% size(In_UL_and_UR)
% size(in_both)
% intersect(In_UL_and_UR,in_both)

LL_lung_voxels( LL_lung_voxels == 0) = [];
LR_lung_voxels( LR_lung_voxels == 0) = [];
UL_lung_voxels( UL_lung_voxels == 0) = [];
UR_lung_voxels( UR_lung_voxels == 0) = [];

LL_lung_voxels( LL_lung_voxels < 0) = [];
LR_lung_voxels( LR_lung_voxels < 0) = [];
UL_lung_voxels( UL_lung_voxels < 0) = [];
UR_lung_voxels( UR_lung_voxels < 0) = [];

%Josh mesh segmentation Anterior/Posterior
vertical_mid_line = 45; %conformal map
vertical_mid_line=151;
LA_lung_voxels = unique(Jash1(vertical_mid_line+1:end,:));
LP_lung_voxels = unique(Jash1(1:vertical_mid_line,:));
% JM changed it so there are no voxels common to both so volumes are okay
In_LA_and_LP = intersect(LA_lung_voxels,LP_lung_voxels);
Remove_from_LA=In_LA_and_LP(1:2:end);  % Remove odd-indexed entries from LL
Remove_from_LP=In_LA_and_LP(2:2:end); % Remove even-indexed entries from LR
[~,indL] = intersect(LA_lung_voxels,Remove_from_LA);  % ind gives the indices in LL_lung_voxels of what we want to remove
LA_lung_voxels(indL) = [];  % This removes them
[~,indR] = intersect(LP_lung_voxels,Remove_from_LP);  % ind gives the indices in LL_lung_voxels of what we want to remove
LP_lung_voxels(indR) = [];  % This removes them
% Check the intersection now
% in_both = intersect(LA_lung_voxels,LP_lung_voxels);
% size(In_LA_and_LP)
% size(in_both)
% intersect(In_LA_and_LP,in_both)

UA_lung_voxels = unique(Jash2(vertical_mid_line+3:end,:));
UP_lung_voxels = unique(Jash2(1:vertical_mid_line+2,:));
In_UA_and_UP = intersect(UA_lung_voxels,UP_lung_voxels);
Remove_from_UA=In_UA_and_UP(1:2:end);  % Remove odd-indexed entries from LL
Remove_from_UP=In_UA_and_UP(2:2:end); % Remove even-indexed entries from LR
[~,indL] = intersect(UA_lung_voxels,Remove_from_UA);  % ind gives the indices in LL_lung_voxels of what we want to remove
UA_lung_voxels(indL) = [];  % This removes them
[~,indR] = intersect(UP_lung_voxels,Remove_from_UP);  % ind gives the indices in LL_lung_voxels of what we want to remove
UP_lung_voxels(indR) = [];  % This removes them
% Check the intersection now
% in_both = intersect(UA_lung_voxels,UP_lung_voxels);
% size(In_UA_and_UP)
% size(in_both)
% intersect(In_UA_and_UP,in_both)
% return

LA_lung_voxels( LA_lung_voxels == 0) = [];
LP_lung_voxels( LP_lung_voxels == 0) = [];
UA_lung_voxels( UA_lung_voxels == 0) = [];
UP_lung_voxels( UP_lung_voxels == 0) = [];

LA_lung_voxels( LA_lung_voxels < 0) = [];
LP_lung_voxels( LP_lung_voxels < 0) = [];
UA_lung_voxels( UA_lung_voxels < 0) = [];
UP_lung_voxels( UP_lung_voxels < 0) = [];


end