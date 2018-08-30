classdef cl_access <  handle
%Author: Bradley Parry; Christine Jacobs-Wagner Lab, Yale University    

    properties(Access = private)
        cellList
        index
        currentIndex
    end
    
    properties
    end
    
    methods(Access = private)
    end
    
    methods
        
        function obj = cl_access(cellList, index)
            %must be initialized with a cellList and an Nby2 matrix
            %specifying frames in the first column and corresponding
            %indices in the second column
            obj.cellList = cellList;
            obj.index = index;
            
            obj.currentIndex = 0;
        end
        
        function out = advanceCell(obj)
            N=1;
            out = false;
            if (obj.currentIndex + N) <= size(obj.index,1)
                obj.currentIndex = obj.currentIndex + N;
                out = true;
            else
                obj.resetIndex;
            end
        end
        
        %perform a calculation on all cells
        function out = calculateForAllCells(obj, fn)
            out = 0;
            if ~isa(fn,'function_handle')
                warning('input argument was not a function handle')
                return
            end
            
            %the output will be arranged as a cell array of N rows
            out = cell(obj.index(end,1) - obj.index(end,1) + 1, 1);
            for k = 1:size(obj.index,1)
                F = obj.index(k,1);
                C = obj.index(k,2);
                out{ obj.index(k,1) - obj.index(1,1) + 1 , 1} = fn(obj.cellList.meshData{F}{C});
            end
            
        end
        
        %allow the user to retrieve a specific frame
        function out = get(Frame)
            out = 0;
            ix = find(obj.index(:,1) == Frame);
            if isempty(ix), return, end
            out = obj.cellList.meshData{Frame}{obj.index(ix,2)};      
        end
        
        function out = getCurrentCell(obj)
            if obj.currentIndex == 0
                obj.currentIndex = 1;
            end

            try
                Frame = obj.index(obj.currentIndex,1);
                Index = obj.index(obj.currentIndex,2);
                out = obj.cellList.meshData{Frame}{Index};
                out.frame = Frame;
            catch
                out = 0;
            end
        end
        
        function resetIndex(obj)
            obj.currentIndex = 0;
        end
    end
    
    methods (Static, Access = public)
        
        function listMethods()
            
            a = methods('cl_access');
            b = methods('handle');
            disp( setdiff(a,b) )

        end
        
    end
    
end