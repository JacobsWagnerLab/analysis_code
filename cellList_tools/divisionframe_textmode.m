%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% This is not function, this is the command for text mode in MicrobeTracker.
% @auther: Ahmad J Paintdakhi & Setsu Hirano
% @date: July 3 2013
% @copyright 2012-2013 Yale University
%=========================================================================
% ********************** input **********************
%maxframe: Frame to which you want to process mesh detection 
%threshValue: Threshold to find max constriction. 
%             Max constiction is defined as the peak of constriction profile.
%             If the degree of constriction is over than threshValue, this command will
%             pick that frame as division frame.
%fileName1: Set path to location you want to save file and put file name after path.
%           This is for cellList after finding first division.
%fileName2: Same as fileName1. This is for cellList after finding 2nd division.
%fileName3: Same as fileName1. This is for cellList after deleting cell mesh after 2nd division.
%flag: Set flag to 1 if you want to observe profile of each cell, 
%      no analysis will take place if set to 1. Otherwise set to 0.
% ********************** output **********************
%cellList: New mesh will be shown in microbeTracker and saved in the folder you set the path.
%cellIdContainer and frameIdContainer: These two arrays will be saved in the folder you set the path.
%                                      These contain the information of division frames for each cell.
%=========================================================================
%This command allows you to find the frame when each single cell shows max
%degree of constriction for time-lapse datasets. 
% To do so,
% 1. you need to obtain cell mesh without any division (set splitThreshold = 1)
% 2. open textmode window
% 3. put this command in the window
% 4. click run
%Note: You can use this script for sychronous and asynchronous populations.
%      You need to set proper parameters and segmntation before use.
%	   Also splitThreshold should be always set to 1.
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% ----- This part is to find first division frame ------------------------

global cellList p 
p.fitDisplay = 0;
maxframe = 150;
listHistory = [];
threshValue = 0.25;
fileName1 = 'path to location and file name';
flag = 1; 
[cellIdContainer frameIdContainer] = divisionframe1(cellList,threshValue,flag); 
if flag == 0
frameIdContainer = [frameIdContainer' maxframe];
save([fileName1 '_cellAndFrameIdContainerDiv1'],'cellIdContainer','frameIdContainer');
for frame = 1:length(frameIdContainer) - 1
for ii = 1:length(cellIdContainer{frame})
try
lst = forcesplitcellSetsu(frameIdContainer(frame),cellIdContainer{frame}(ii),[]);
catch
end
listHistory = [listHistory lst];
end %for ii = 1:length(cellIdContainer{frame})
if ~isempty(lst) && frame<maxframe
if length(lst) == 1, cellList = CL_removeCell(cell,frameIdContainer(frame),cellList); continue; end
CL_setField(lst(1),frameIdContainer(frame),'timelapse',1,cellList);
CL_setField(lst(2),frameIdContainer(frame),'timelapse',1,cellList);
processTextMode(frameIdContainer(frame),4,listHistory,[1 0 0 0],{0},'',0,0)
processTextMode([frameIdContainer(frame)+1 frameIdContainer(frame+1)],1,listHistory,[1 0 0 0],{0},'',0,0)
end %if ~isempty(lst) && frame<maxframe
end %for frame = 1:length(frameIdContainer) - 1
savemesh(fileName1,[],0,[])
end %if flag ==0

% ----- This part is to find 2nd division frame ------------------------
global cellList
listHistory = [];
threshValue = 0.25;
listSecondValueHistory = [];
frameIdContainer = [];
cellIdContainer =[];
fileName2 = 'path to location and file name';
fileName3 = 'path to location and file name';
flag = 1; 
[cellIdContainer frameIdContainer]= divisionframe2(cellList,threshValue,flag);
if flag ==0
listSecondValueHistory = cell(1,length(frameIdContainer));
frameIdContainer = [frameIdContainer' frameIdContainer(end) + 2];
save([fileName2 '_cellAndFrameIdContainerDiv2'],'cellIdContainer','frameIdContainer');
for frame = 1:length(frameIdContainer) - 1
list2 = [];
for ii = 1:length(cellIdContainer{frame})
try
lst = forcesplitcellSetsu(frameIdContainer(frame),cellIdContainer{frame}(ii),[]);
catch
end
if length(lst)> 1,list2 = [list2 lst(2)];end 
listHistory = [listHistory lst];
end %for ii = 1:length(cellIdContainer{frame})
listSecondValueHistory{frame} = list2;
if ~isempty(lst)
if length(lst) == 1, cellList = CL_removeCell(cell,frameIdContainer(frame),cellList); continue; end
CL_setField(lst(1),frameIdContainer(frame),'timelapse',1,cellList);
CL_setField(lst(2),frameIdContainer(frame),'timelapse',1,cellList);
processTextMode(frameIdContainer(frame),4,lst,[1 0 0 0],{0},'',0,0)
processTextMode([frameIdContainer(frame)+1 frameIdContainer(frame+1)],1,listHistory,[1 0 0 0],{0},'',0,0)
end %if ~isempty(lst)
end %for frame = 1:length(frameIdContainer) - 1
savemesh(fileName2,[],0,[])
save([fileName2 '_cellAndFrameIdContainerDiv2all'],'cellIdContainer','frameIdContainer','listSecondValueHistory');
end

% ----- This part is for deleting cell mesh after [2nd division frame + 2].-----

cellsToDelete = [];
for frame = 1:length(frameIdContainer) - 1
    if length(cellIdContainer{frame}) > 1
          for jjj = 1:length(cellIdContainer{frame})
               cellsToDelete = [cellsToDelete cellIdContainer{frame}(jjj)];
               cellsToDelete = [cellsToDelete listSecondValueHistory{frame}(jjj)];
          end
    else
         cellsToDelete = [cellsToDelete cellIdContainer{frame}];
         cellsToDelete = [cellsToDelete listSecondValueHistory{frame}];
   end
    for ii = frameIdContainer(frame)+3:length(cellList.meshData)
        for iii = 1:length(cellsToDelete)
            if CL_isCell(cellsToDelete(iii), ii, cellList)
               cellList = CL_removeCell(cellsToDelete(iii),ii,cellList);
            end
        end
       
   end
cellsToDelete = [];
end
savemesh(fileName3,[],0,[])