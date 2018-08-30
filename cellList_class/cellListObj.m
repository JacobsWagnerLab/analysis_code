%{
 
Author: Bradley Parry; Christine Jacobs-Wagner Lab, Yale University
PURPOSE

cellListObj is an interface for Oufti generated cellists. This wrapper
screens the cellList for valid cells and grants the user direct access to
valid cells without the user screening the entire list. By default, a valid
cell is one with a mesh field of size greater than 6. To organize the
cellList, cellListObj internally constructs multiple indices of the
cellList that reference cellList data for straightforward user access.

To maintain backward compatability, cellListObj provides meshData and
cellId which are identical to the fields from the original Oufti cellList.
Objects of this class can be used in any instance requiring an original
Oufti cellList.

INITIALIZATION

cellListWrapper = cellListObj(cellList);
	
cellListObj takes in an Oufti cellList as the sole input argument. At
initialization, the cellListObj traverses the input cellList and
constructs default indexing methods. See addRestriction for more
information on adding conditions that must be met for a cell to be
considered valid and added to the index.

FUNCTION CALLS

addRestriction
USE
    cellListWrapper.addRestriction(@functionHandle)
EXPLANATION
    addRestriction accepts a function handle as an input argument and adds it
    as a condition that must be met for a cell to be valid. The function
    pointed to by functionHandle must return true or false for any input
    argument. functionHandle should take in a single cell from a cellList and
    return logical operators indicating that the restriction has been met and
    the cell is valid (true) or that the cell is not valid (false). Once a new
    restriction condition has been added, the cellListObj will automoatically
    recalculate valid indexing methods.

advanceCell
USE
    cellListWrapper.advanceCell()
    cellListWrapper.advanceCell(N)
EXPLANATION
    advanceCell increments to the next valid cell in the index. If the
    operation was successful, advanceCell returns true. If not (indicating the
    end of valid cells in the cellList has been reached), advanceCell returns
    false. If an input argument is passed, it must be an integer and indicates
    the number of steps forward that should be taken in the index.

calculateForAllCells
USE
    cellListWrapper.calculateForAllCells(@functionHandle)
EXPLANATION
    calculateForAllCells applies the function pointed to by functionHandle to
    all valid cells in the cellList. The output is stored as a cell array. The
    output for each biological cell is returned in its own cell in the cell
    array.

calculateForLineage
USE
    cellListWrapper.calculateForLineage(@functionHandle)
EXPLANATION
    Applies the function pointed to by functionHandle to each cell in a
    lineage. The function pointed to must take in a single Oufti cell; there
    are no restrictions on the output that it may produce.

calculateForAllLineages
USE
    cellListWrapper.calculateForAllLineages(@functionHandle)
EXPLANATION
    Applies the function pointed to by functionHandle to each cell in a
    lineage. The function pointed to must take in a single Oufti cell; there
    are no restrictions on the output that it may produce. Equivalent to
    running calculateForLineage on all cellIds.

cellAtFrameIndex
USE
    cellListWrapper.cellAtFrameIndex(Frame, Index)
EXPLANATION
    Provides redundant access as cellList.meshData{Frame}{Index}. No
    restriction checks are performed to return a cell with this function call
    and Frame and Index is not required to be found in the valid indices array.

cellAtFrameId	
USE
    cellListWrapper.cellAtFrameId(Frame,cellId)
EXPLANATION
    Returns the cell at the indicated frame with the specified cellId. If no
    cell is found, 0 is returned. No restriction checks are performed with this
    function call.

getCurrentCell
USE
    cellListWrapper.getCurrentCell()
EXPLANATION
    Returns the current valid Oufti cell as determined by default restriction
    settings and user supplied restriction rules (if supplied through
    addRestriction). cellListObj internally keeps an index of valid cells and
    allows the user to directly access cells within this list. Use this
    function in conjunction with advanceCell.

getIndex
USE
    cellListWrapper.getIndex
EXPLANATION
    Returns to the user the list of indices (n rows of [Frame, Index])
    determined to be valid, whether by the default method or after the addition
    of some restriction.

getLineage
USE
    cellListWrapper.getLineage(cellId)
EXPLANATION
    Assembles the lineage for the cell corresponding to cellId. A cl_access
    object is returned. Objects of type cl_access are compatible with
    cellListObj and have corresponding functions (see cl_access, below).

getLineageIndex
USE
    cellListWrapper.getLineageIndex
EXPLANATION
    Returns the list of indices for each lineage. The index of the first cellId
    will be in the first cell array. E.g., cellListWrapper.getLineageIndex{1}
    is the index of the lineage corresponding to cellListWrapper.
    uniqueCellIDs(1)

listMethods
USE
    cellListWrapper.listMethods
EXPLANATION
    Prints a list of methods available to the cellList object.

modifyCellList
USE
    cellListWrapper.modifyCellList(@functionHandle)
EXPLANATION
    Modifies each valid cell with the output from the function pointed to by
    function handle. Each cell will be replaced with the output from the
    function; the function pointed to by functionHandle should return all
    aspects of the cell that the user intends to keep.

randomSample
USE
    cellListWrapper.randomSample
EXPLANATION
    Returns a randomly selected cell from the valid cells.

removeRestriction
USE
    cellListWrapper.removeRestriction
EXPLANATION
    Removes any restriction put in place by addRestriction, and resets the
    index of valid cells along with iteration index along the valid cells.

resetIndex
USE
    cellListWrapper.resetIndex
EXPLANATION
    Resets the index of valid cells as well as current location in the index.
    Any sorting or restrictions are eliminated. Any restriction functions are
    not removed and the index will be recalculated with them in place.

sortBy
USE
    cellListWrapper.sortBy(@functionHandle)
    cellListWrapper.sortBy(@functionHandle, direction)
EXPLANATION
    The function pointed to by functionHandle must return a sortable type.
    Cells will be sorted in the direction describted by the optional argument
    direction, which must be either ?ascend?,?descend? If direction is not
    specified sorting will be performed in ascending order.

###########################################################################
Brad Parry
Jacobs-Wagner lab
June 2016
%}

classdef cellListObj < handle
    
    properties
        uniqueCellIDs
        cellId
        meshData
    end
    
    properties(Access = private)
        cellList
        currentIndex
        lineageList
        restrictionHandle
        sortedIndices
        validIndices
    end
    
    methods(Access = private)
        
        %Use restrictionHandle to determine if CL_cell is a valid cell in
        %the cellList
        function out = restrictionFn(obj, CL_cell)
            out = obj.restrictionHandle(CL_cell);
            if ~islogical(out)
                warning('Restriction function does not return logical type')
                out = false;
            end
        end
        
        function getLineageList(obj)
            % constructs indices of cell lineages for all cellIds
            obj.uniqueCellIDs = obj.getUniqueCellIds(obj.cellList);
            obj.lineageList = cell(1,length(obj.uniqueCellIDs));
            
            for C = 1:length(obj.uniqueCellIDs)
                for F = 1:length(obj.cellList.cellId)
                    ix = find(obj.cellList.cellId{F} == obj.uniqueCellIDs(C),1);
                    if ~isempty(ix) && obj.isCell(obj.cellList.meshData{F}{ix(1)})
                        obj.lineageList{C}(end+1,1:2) = [F, ix];
                    end
                end
            end
            
            obj.uniqueCellIDs(cellfun(@isempty, obj.lineageList)) = [];
            obj.lineageList(cellfun(@isempty, obj.lineageList)) = [];
            
        end
        
        %find frame and index of all cells that meet criteria for being
        %a valid cell as well as any criteria provided by the user
        function getValidIndices(obj)
            %basic cocde to begin a calculation over the cellList
            obj.validIndices = [];
            obj.currentIndex = 0;
            %iterate over all frames
            for F = 1:length(obj.cellList.meshData)

                %iterate over all cells
                for C = 1:length(obj.cellList.meshData{F})
                    %check if the cell in the current frame is valid
                    %is the cell empty, does it have a mesh, is the mesh a "real" mesh?
                    %if any user conditions are provided, are they met?
                    if obj.isCell(obj.cellList.meshData{F}{C})
                        obj.validIndices(end+1,1:2) = [F, C];
                    end
                end
            end
        end
        
        %convert an index into a cellId.
        function index = id2index(obj,Frame,CellId)
            %returns NaN if the cellId cannot be found in the frame
            index = NaN;
            
            ix = find(obj.cellList.cellId{Frame} == CellId,1);
            if isempty(ix)
                return
            end
            
            if obj.isCell(Frame,ix)
                index = ix;
            end
        end
        
        %check if a cell in a given frame is a cell
        function answer = isCell(obj, input)
            answer = false;
            if isempty(input)
                return
            end
            if ~isfield(input,'mesh')
                return
            end
            if length(input.mesh) < 6
                return
            end
            if obj.restrictionFn(input) == false
                return
            end
            answer = true;
            return
        end
        
        
        %check if the cell exists with field  NOT IMPLEMENTED
        function answer = isCellWithField(obj,Frame,index,fieldname)
            answer = false;
            if isCell == false
                return
            else 
                if isfield(obj.cellList.meshData{Frame}{index},fieldname)
                    answer = true;
                end
            end
        end
        
        %check that the cell exists with a non-empty field NOT IMPLEMENTED
        function answer = isCellWithNonEmptyField(obj,Frame,index,fieldname)
            answer = isCellWithField(Frame,index,fieldname);
            if answer == false
                return
            else
                if ~isempty(obj.cellList.meshData{Frame}{index}.(fieldname))
                    answer = true;
                end
            end
        end
        
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %       CLASS CONSTRUCTOR
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = cellListObj(cellList,dataType)
            
            if nargin == 0
                error('Class cellListObj requires a single Oufti cellList be provided as an input argument')
            elseif ~isfield(cellList,'meshData') && ~isfield(cellList,'cellId')
                error('Class cellListObj requires an Oufti cellList as input')
            end
            
            if nargin == 1
                dataType = 'timeLapse';
            end
            
            tic
            %initialize and assign the cellList
            obj.cellList = cellList;
            obj.cellId = cellList.cellId;
            obj.meshData = cellList.meshData;
            
            %rather than test for the existence of a restriction handle,
            %the code requires a restriction to always exist. In the
            %construction, initialize a restriction function which always
            %evaluates to true.
            obj.restrictionHandle = @(x) true;
            
            %initialize Frame and Cell indices that are valid
            obj.getValidIndices();
            obj.currentIndex = 0;
            if strcmpi(dataType,'timelapse')
                obj.getLineageList();
            end
            
            t=toc;
            tx = ['wrapper initialized and cellList indexed in ',num2str(t),' secs'];
            disp(tx)
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %       PUBLIC CLASS METHODS
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %add a condition that must evaluate to true for a cell to be valid
        function addRestriction(obj, fn)
            
            if ~isa(fn,'function_handle')
                warning('input argument was not a function handle')
                return
            end
            
            obj.restrictionHandle = fn;
            obj.getValidIndices;
        end
        
        %advance the index of validCells forward
        function out = advanceCell(obj, N)
            if nargin == 1
                N = 1;
            end
            
            out = false;
            if (obj.currentIndex + N) <= size(obj.validIndices,1)
                obj.currentIndex = obj.currentIndex + N;
                out = true;
            else
                obj.resetIndex;
            end
        end
        
        %performs specified calculation on all cells
        function out = calculateForAllCells(obj, fn, inputArgs)
            if ~isa(fn,'function_handle')
                warning('input argument was not a function handle')
                return
            end

            out = cell(1,length(obj.validIndices));
            for k = 1:size(obj.validIndices,1)
                F = obj.validIndices(k,1);
                C = obj.validIndices(k,2);
%                 [F,C]
                if nargin == 3
                    out{1,k} = fn(obj.cellList.meshData{F}{C}, inputArgs);
                else
                    out{1,k} = fn(obj.cellList.meshData{F}{C});
                end
            end
            
        end
        
        %performs specified calculation for all cells in a lineage
        function out = calculateForLineage(obj, cellId, fn, inputArgs)
            out = 0;
            
            if ~isa(fn,'function_handle')
                warning('input argument was not a function handle')
                return
            end
            
            ix = find(obj.uniqueCellIDs == cellId);
            
            if isempty(ix)
                tx = ['No valid lineage found for cell with id ',num2str(cellId)];
                disp(tx)
                return
            end
            
            lineage = cl_access(obj.cellList, obj.lineageList{ix});
            if nargin == 4
                out = lineage.calculateForAllCells(fn, inputArgs);
            else
                out = lineage.calculateForAllCells(fn);
            end
            
        end
        
        %performs specified calculation on all cells in all lineages
        function out = calculateForAllLineages(obj, fn, inputArgs)
            %iterate through uniqueCellIDs
            %   -create cl_access object
            %   -cl_access.calculateForAllCells(fn)
            %   -collect output in cell of cell array, each cellId getting
            %   one cellarray with output organized by frame in sub-cell
            %   arrays
            %   -try and compact the output (convert cells to values if
            %   possible)
            %   
            
            out = [];
            for C = 1:length(obj.uniqueCellIDs)
                if nargin == 3
                    out{end+1} = obj.calculateForLineage(obj.uniqueCellIDs(C), fn, inputArgs);
                else
                    out{end+1} = obj.calculateForLineage(obj.uniqueCellIDs(C), fn);
                end
            end
            
        end
        
        %provides redundant access as cellList.meshData{Frame}{Index}
        function out = cellAtFrameIndex(obj, Frame, Index)
            try
                out = obj.cellList.meshData{Frame}{Index};
            catch
                out = 0;
            end
        end
        
        %return the cell the cell with id specified at Frame
        function out = cellAtFrameId(obj, Frame, cellId)
            out = 0;
            try
                ix = obj.cellList.cellId{Frame} == cellId;
                out = obj.cellList.meshData{Frame}{ix};
            catch

            end
        end
        
        function demo = demograph(obj, signalName, normalization)
            % cl = cellListObj(cellList);
            % signalName = 'signal0';

            if nargin == 2
                nm = @(x) x;
            else
                nm = normalization;
            end

            xhat = 0:0.01:1;
            demo = zeros(length(xhat), size(obj.validIndices,1));
            index = 1;

            for k = 1:size(obj.validIndices,1)
                F = obj.validIndices(k,1);
                C = obj.validIndices(k,2);
                obj.cellList.meshData{F}{C}.mesh

               sl =  cumsum(obj.getStepLength( obj.cellList.meshData{F}{C} ));
               sl = (sl - sl(1)) / (sl(end) - sl(1) );
               segArea = obj.getMeshSegArea( obj.cellList.meshData{F}{C} );
               signal = obj.cellList.meshData{F}{C}.(signalName);

               %get rid of bad values if they exist (Oufti can screw up at the poles)
               kill = signal==0 | segArea==0;
               signal(kill) = [];
               segArea(kill) = [];
               sl(kill) = [];

               signalAreaNorm = signal ./ segArea;
               tmp = interp1(sl, nm(signal), xhat, 'pchip',0);
               demo(:,index) = tmp(:);
               index = index + 1;


        % %                 [F,C]
        %         if nargin == 3
        %             out{1,k} = fn(obj.cellList.meshData{F}{C}, inputArgs);
        %         else
        %             out{1,k} = fn(obj.cellList.meshData{F}{C});
        %         end
            end
        % 
        % 
        % 
        %     xhat = 0:0.01:1;
        %     demo = zeros(length(xhat), size(obj.getIndex,1));
        %     index = 1;
        %     while obj.advanceCell()
        %        sl =  cumsum(getStepLength( obj.getCurrentCell().mesh ));
        %        sl = (sl - sl(1)) / (sl(end) - sl(1) );
        %        segArea = getMeshSegArea( obj.getCurrentCell().mesh );
        %        signal = obj.getCurrentCell().(signalName);
        % 
        %        %get rid of bad values if they exist (Oufti can screw up at the poles)
        %        kill = signal==0 | segArea==0;
        %        signal(kill) = [];
        %        segArea(kill) = [];
        %        sl(kill) = [];
        % 
        %        signalAreaNorm = signal ./ segArea;
        %        tmp = interp1(sl, nm(signal), xhat, 'pchip',0);
        %        demo(:,index) = tmp(:);
        %        index = index + 1;
        %     end
            demo = demo(:,1:(index-1));
        end

        
        
        %Return the current valid cell
        function out = getCurrentCell(obj)
            if obj.currentIndex == 0
                obj.currentIndex = 1;
            end
            F = obj.validIndices(obj.currentIndex,1);
            C = obj.validIndices(obj.currentIndex,2);
            out = obj.cellList.meshData{F}{C};
            out.frame = F;
        end
        
        function out = getCurrentCellIndex(obj)
            out(1) = obj.validIndices(obj.currentIndex,1);
            out(2) = obj.validIndices(obj.currentIndex,2);
            return
        end
        
        % return the valid indices to the user
        function out = getIndex(obj)
            out = obj.validIndices;
        end
        
        %que the lineage associated with the provided cellId
        function lineage = getLineage(obj, cellId)
            lineage = 0;
            ix = find(obj.uniqueCellIDs == cellId);
            
            if isempty(ix)
                tx = ['No valid lineage found for cell with id ',num2str(cellId)];
                disp(tx)
                return
            end
            
            lineage = cl_access(obj.cellList,obj.lineageList{ix});
            
        end
        
        function out = getLineageIndex(obj)
            out = obj.lineageList;
        end

        %directly modify the cellList in place. changes will not be
        %returned, but persist in the cellList structure itself.
        function modifyCellList(obj, modFn, inputArgs)
            for k = 1:size(obj.validIndices,1)
                F = obj.validIndices(k,1);
                C = obj.validIndices(k,2);
                obj.cellList.meshData{F}{C}.frame = F;
                if nargin == 3
                    obj.cellList.meshData{F}{C} = modFn(obj.cellList.meshData{F}{C}, inputArgs);
                else
                    obj.cellList.meshData{F}{C} = modFn(obj.cellList.meshData{F}{C});
                end
            end
        end
        
        function modifyCell(obj, Frame, Cell, fieldname, value )
            obj.cellList.meshData{Frame}{Cell}.(fieldname) = value;
        end
        
        %return a randomly selected valid cell
        function out = randomSample(obj)
            ix = randi( size(obj.validIndices, 1) );
            F = obj.validIndices(ix,1);
            C = obj.validIndices(ix,2);
            
            out = obj.cellList.meshData{F}{C};
        end
        
        %undo any additional restriction that was placed on the valid cell
        %index
        function removeRestriction(obj)
            obj.restrictionHandle = @(x) true;
            obj.resetIndex();
        end

        %reset the currentIndex as well as the list of valid indices.
        function resetIndex(obj)
            obj.getValidIndices()
        end

        %sort valid indices by the function handle targeted by fn.
        %direction should be either 'ascend' or 'descend' and indicates the
        %direction of sorting to occur.
        function sortBy(obj, fn, direction)
            
            if nargin ~= 3
                direction = 'ascend';
            end
            
            values = zeros(1,size(obj.validIndices,1));
            for k = 1:size(obj.validIndices,1)
                F = obj.validIndices(k,1);
                C = obj.validIndices(k,2);
                values(k) = fn(obj.cellList.meshData{F}{C});
            end
            
            [~,ix] = sort(values, direction);
            obj.validIndices = obj.validIndices(ix,:);
        end
        


    end
    
    methods (Static, Access = public)
        
        function listMethods()       
            a = methods('cellListObj');
            b = methods('handle');
            disp( setdiff(a,b) )
        end
        
        function A = getMeshSegArea(oneCell)
            mesh = oneCell.mesh;
            x = cat(2,mesh(1:end-1,[1,3]),...
                mesh(2:end,[3,1]),...
                mesh(1:end-1,1));
            y = cat(2,mesh(1:end-1,[2,4]),...
                mesh(2:end,[4,2]),...
                mesh(1:end-1,2));

            x = double(x);
            y = double(y);

            A = zeros(size(mesh,1)-1,1);
            for k = 1:size(x,1)
                A(k) = polyarea(x(k,:), y(k,:));
            end
        end
        
        function sl = getStepLength(cellMesh)
        cellMesh = cellMesh.mesh;
        x = mean(cellMesh(:,[1,3]),2);
        y = mean(cellMesh(:,[2,4]),2);

        dxy = (diff(x).^2 + diff(y).^2).^(1/2);
        sl = (dxy);

        end
        

        
    end
    
    methods (Static, Access = private)
        
        function cid = getUniqueCellIds(cellList)
            cid = [];
            for F = 1:length(cellList.cellId)
                cid = [cid, cellList.cellId{F}];
            end
            cid = unique(cid);
        end
        
        function answer = isCellArrayUniform(cellArray)
        %     check if any dimension of data in a cell array is uniform. if
        %     so, return the dimension.
        %       -1: not uniform
        %        0: length of each element = 1
        %        1: uniform row count
        %        2: unforrm col count
            answer = -1;

            getUniqueCount = @(x) length(unique(x(x>0)));
            L = cellfun(@length, cellArray);
            R = cellfun(@(x) size(x,1), cellArray);
            C = cellfun(@(x) size(x,2), cellArray);

            if getUniqueCount(L) == 1
                answer = 0;
                return
            end

            Rcount = getUniqueCount(R);
            Ccount = getUniqueCount(C);

            if Rcount == 1 && Ccount ~= 1
                answer = 1;
            elseif Rcount ~= 1 && Ccount == 1
                answer = 2;
            end

        end

%         function cellArray2Matrix()
%             %if cell array inputs are uniform, convert to a matrix
%         end
        
%         function getIndexFromID(obj)
%             ix = find(obj.cellList.cellId{obj.Frame} == obj.cellId,1);
%             if isempty(ix)
%                 obj.cellIndex = NaN;
%             else
%                 obj.cellIndex = ix;
%             end
%         end
   
    end
    
end