% to do Multiple Particle Tracking using spots in spotList (i.e.detected by
% spotFinderF) and scriptTrackGeneral (from u-Track)

%{
-About-
This scripts allows to do Multiple Particle Tracking using spots in
spotList (i.e, without cellList as detected by spotFinderF) and
scriptTrackGeneral (from u-Track)

-Inputs-
    spotList: spots data (w/o cellList), as detected by spotFinderF
    track_pars: tracking parameters
        track_pars(1)=5 - "quality control"=min length of the trace
        track_pars(2)=2 - single-spot traces with ocasional extra spotes included
        track_pars(3)=0.1 - max fraction of frames with extra spots
    cell: cell ID for a cell from cellList to use

-varargin-
N/A

-Outputs-
movieInfo: utrack-format data for detcted spots
           Orgainzed by frames with 3 fileds for coorinates and
           amplitdue (magnitude) of the spot - xCoord, yCoord and amp:
            movieInfo(frame).xCoord= spot X
            movieInfo(frame).yCoord= spot Y
            movieInfo(frame).amp= spot amplitude

-Example-
  movieInfo=spots_2_movieInfo(cellList,11)
   
-Supplementary-


-Keywords-
spot detection, image analysis, uTrack

-Dependencies-
scriptTrackGeneral
prepare_for_TrackGeneral
Those are key functions in uTrack package

-References-

-Author-
 Ivan Surovtsev, 2016.03.31
 Updated: 2018.05.16, to comply with Code Repo guide
%}

% Set tracking parameters
track_pars=[5,2,0.1]; % track_pars(1)=5 - "quality control"=min length of the trace
% track_pars(2)=2 - single-spot traces with ocasional extra spotes included
% track_pars(3)=0.1 - max fraction of frames with extra spots

% Prepare variables:
traceList_N_all=[];

% prepare parameters for tracking function
prepare_for_TrackGeneral;

% MAIN

% First, convert spots into proper fromat
movieInfo=spotList_2_movieInfo(spotList);
% Now, tracking function call
if ~isempty(movieInfo)
    [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
    % convert into traceList
    [traceList_N,info]=tracksFinal_2_traceList(tracksFinal, track_pars);  
    
    % COMMENT:
    % THIS PART IF PLACING TRACES IN THE CELL IS NOT REQUIRED
    % convert back into cell coordinates
    if ~isempty(traceList_N)
        [traceList_N,trace_in_cell]=traces_2_cell(traceList_N,cellList,cell,track_pars2);
        %   % combine traces from all cells
        %   traceList_N_all=[traceList_N_all,traceList_N];
    end
end

