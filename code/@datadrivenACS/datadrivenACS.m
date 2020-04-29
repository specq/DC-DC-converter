classdef datadrivenACS < handle & matlab.mixin.Copyable
    
    properties (Access = public)
        Model = struct('Plant',[],'Frequency',[],'SamplingGrid',[])
        Feedback = struct('objective',[],'constraints',[],'controller',[]...
            ,'parameters',[],'sol',[])
        Feedforward = struct('objective',[],'constraints',[],'controller',...
            [],'parameters',[],'sol',[])
        Logs = datadrivenACS.empty;
    
    end
    
    properties (Access = private)
        M = []; % formated data, internal use
        FB= []; % formated data, internal use
        FF= []; % formated data, internal use
        
        internals = struct('FB',struct('disp',''),...
                           'FF',struct('disp',''));
        
        flag
    end
    
    
    methods (Access = public)
        
        function [this,objs,cons,param_struct] = datadrivenACS(varargin)
            objs = struct('o2W1',[],'o2W2',[],'o2W3',[],...
                'oinfW1',[],'oinfW2',[],'oinfW3',[]);
            cons = struct('cinfW1',[],'cinfW2',[],'cinfW3',[]);% 'c2W1',[],'c2W2',[],'c2W3',[]);
            
            
            this.Feedback.objective = objs;
            this.Feedforward.objective = objs;
            
            this.Feedback.constraints = cons;
            this.Feedforward.constraints = cons;
            
            lpv_struct = struct('enable',false,'theta',[],'additionalPts',[],'fixedDC',1,'K0mat',0);
            
            param_struct = struct('maxIter',20,'tol',1e-3,'exit',0,...
                'c0',0','c1',0,'c2',[],'lpv',lpv_struct);
            
            this.Feedback.controller = struct('order',[],'K0',[],...
                'Fy',tf(1),'Ts',[]);
            
            this.Feedback.parameters = param_struct;
            
            this.Feedforward.controller = struct('order',[],'K0',[],...
                'Fy',tf(1),'Ts',[]);
            this.Feedforward.parameters = param_struct;
            
            this.flag.largeObj =0;
            
        end
        
        function clearLogs(this)
            this.internals.FF.disp = '' ;
            this.internals.FB.disp = '' ;
        end

        function KFB = computeKFB(this)
           KFB = computeK(this,'FB'); 
        end
        function KFB = computeKFF(this)
           KFB = computeK(this,'FF'); 
        end
	function [KFBmat] = solve(this)
		KFBmat = solveFB(this);
		%KFFmat = solveFF(this);
	end
    end

    
end
