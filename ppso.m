function [xopt, fval, exitflag, output] = ppso(func,npars,lb,ub,varargin)

% ppso.m - Find minimum of function using global Particle Swarm
%  Optimization algorithm. Version 1.03.
%
% PROTOTYPE:
%   [xopt, fval, exitflag, output] = ppso(func,npars,lb,ub) or
%   [xopt, fval, exitflag, output] = ppso(func,npars,lb,ub,options) or
%   [xopt, fval, exitflag, output] = ppso(func,npars,lb,ub,options,auxdata)
%
% DESCRIPTION:
% Find minimum of function using a global version of Particle Swarm
% Optimization algorithm, as described in Ref. 1. The values of the
% Cognitive and Social weights are, respectively, cC = 1.49445*r2(0,1) and
% cS = 1.49445*r3(0,1) where r2(0,1) and r3(0,1) are two random numbers
% with uniform distribution between 0 and 1. The Inertial weight is
% cI = 0.5 + 0.5*r1(0,1). A decreasing version of the Inertial weight has
% been tested, but it results in a worst performance of the overall
% algorithm in terms both of speed and final results.
% The oracle penalty method (see Ref. 2) is used for the nonlinear
% constraints handling. l1, l2 or l_inf norm can be chosen for the
% computation of residuals.
% A vectorized version of the code is implemented as well, to speed up the
% simulation if the user-defined function is vectorized too.
% NOTE. No check on the input function, because it is not possible to check
% if a function exists as a nested function in the calling one.
% 
% NOTE. Nonlinear constraints handling does not work properly yet. It has
% been developed though, so the fitness function must return the nonlinear
% constraints, as described in the "INPUT" section.
% 
% To do:
%   check nonlinear constraints handling (not working properly yet!!!)
%   check the code with integer components of the state
%   avoid the check of the whole history when looking for best positions,
%    to improve the speed
%
% INPUT:
%   func           Handle to the function to be used. Note that func must
%                   be a handle, i.e. func = @nameOfTheFunction.
%                   The prototype of the function must be:
%                       [f,c,ceq] = func(x,auxdata);
%                   where
%                       f            : fitness function
%                       c[N,1] < 0   : inequality constraints. It must be a
%                                       column vector if inequality
%                                       constraints are present, an empty
%                                       array otherwise.
%                       ceq[N,1] = 0 : equality constraints. It must be a
%                                       column vector if equality
%                                       constraints are present, an empty
%                                       array otherwise.
%                       x[npars,1]   : input column vector
%                       auxdata : variable containing any auxiliary data
%                   If the vectorized code is used, inputs and outputs of
%                   func must be matrices with the following dimensions:
%                       x[npars,npop]
%                       f[1,npop], c[Nc,npop], ceq[Nceq,npop]
%   npars          Positive integer. Number of unknown parameters.
%   lb             Array of dimension npars containing the lower boundaries
%                   of the unknown parameters. Note, -Inf returns an error.
%   ub             Array of dimension npars containing the upper boundaries
%                   of the unknown parameters. Note, Inf returns an error.
%   Following the optional inputs. If only the i-th optional input has to
%   be specified, put a [] for each of the previous optional inputs not
%   needed.
%   options        Structure containing the options:
%                   npop : Positive integer > 1. Number of particles in the
%                    swarm. Default: 10*npars.
%                   niter : Positive integer. Maximum number of iterations.
%                    Default: 500.
%                   X0 : Initial guess for one of the particle's position.
%                    Default: [].
%                   StallIterLimit : Positive integer. Iterations end when
%                    the relative change in best objective function value
%                    over the last StallIterLimit iterations is less than
%                    options.TolFun. Default: niter.
%                   createNewPop : logical value (0 or 1). Set to 1 if a
%                    new random population shall be created when StallIter
%                    is equal to StallIterLimit/3, set to 0 otherwise.
%                    Default: 1.
%                   TolFun : Nonnegative scalar. Iterations end when the
%                    relative change in best objective function value over
%                    the last StallIterLimit iterations is less than
%                    options.TolFun. Default: 1e-6.
%                   TolCon : Nonnegative scalar. TolCon is used to
%                    determine the feasibility with respect to nonlinear
%                    constraints. Default: 1e-6.
%                   Display : 1 for display the iterations at screen,
%                    0 otherwise. Default: 1.
%                   Plot : 1 for plot of the iterations, 0 otherwise.
%                    Default: 1.
%                   IntCon : Vector of positive integers taking values from
%                    1 to npars. Each value in IntCon represents an x
%                    component that is integer-valued. Default: [].
%                   maxFcount : Positive integer. Iterations end when the
%                    number of function evaluations exceeds
%                    options.maxFcounts. Default: [].
%                   s : User specified settings for random number
%                    generator. Used in order to be able to reproduce
%                    previous results. If no settings are specified,
%                    rng('shuffle') is used. Default: [].
%                   fvalOpt : Optimal known value of the fitness function.
%                    Use it in order to test the algorithm. Default: [].
%                   TerminationErr : Nonnegative scalar. Iterations end
%                    when the error between the computed fitness function
%                    and the known optimal one is less than
%                    options.TerminationErr. options.fvalOpt needed.
%                    Default: [].
%                   Residuals : Value of the norm to be used for the
%                    computation of residuals. One can chose between:
%                    1   : l1 norm
%                           res = Sum(abs(ceq)) + Sum(max(0,c))
%                    2   : l2 norm
%                           res = sqrt(Sum(abs(ceq)^2) + Sum(max(0,c)^2))
%                    Inf : l_inf norm
%                           res = max(abs(ceq),max(0,c))
%                    Default: 1.
%                   Vectorized : Logical value (0 or 1). Set to 1 if the
%                    vectorized version of the code shall be used, set to 0
%                    otherwise. Note that the fitness function shall be
%                    vectorized as well. Default: 0.
%                   Parallel : Logical value (0 or 1). Set to 1 if the
%                    parallelized version of the code shall be used, set to
%                    0 otherwise. If both Vectorized and Parallel are 1,
%                    then the vectorized mode is used. Default: 0.
%                   Save : If 1, save 'iter', 'x' and 'f' at the end of
%                    each iteration in a .mat file. Default: 0.
%                   filename : If options.Save = 1, filename is the name of
%                    the saved file. It must be a character string ending
%                    with '.mat'. Default: 'psoTemp.mat'.
%   auxdata        Variable containing any auxiliary data needed for the
%                   computation of the input function. Note, auxdata can be
%                   an array, a cell or a structure, but it must be a
%                   single variable. Default: [].
%
% OUTPUT:
%   xopt[npars,1]  Best particle located by pso during its iterations.
%   fval           Fitness function evaluated at x, if an optimal solution
%                   has been found. Penalty function evaluated at x,
%                   otherwise.
%   exitflag       Algorithm stopping condition, returned as an integer
%                   identifying the reason the algorithm stopped:
%                    2: Optimal solution found.
%                    1: Feasible solution found. Relative change in the
%                       fitness function less than TolFun.
%                    0: Feasible solution found. Maximum number of
%                       iterations, or function evaluations, reached and
%                       the relative change in the fitness function is
%                       greater than TolFun.
%                   -2: Maximum number of iterations reached, but none of
%                       the particles in the swarm verifies the
%                       constraints.
%   output         Structure containing output:
%                   options: User defined options.
%                   rngstate: Identifier of the random seed used in the
%                    current run. For reproducibilty.
%                   iterations: Number of iterations before the algorithm
%                    stopped.
%                   funccount: Number of functions evaluations.
%                   message: Exit message.
%                   positions: Positions history of each particle of the
%                    swarm. (npars X npop X iterations) matrix.
%                   velocities: Velocity history of each particle of the
%                    swarm. (npars X npop X iterations) matrix.
%                   fitnessFcn: Fitness function history of each particle
%                    of the swarm. (1 X npop X iterations) matrix.
%                   fitnessFcnOpt: Optimal fitness function history along
%                    the swarm. (iterations X 1) matrix.
%                   penalty: penalty function history of each particle of
%                    the swarm. (1 X npop X iterations) matrix.
%                   Oracle: Oracle history of each particle of the swarm.
%                    (1 X npop X iterations) matrix.
%                   residuals: residual history of each particle of the
%                    swarm. (1 X npop X iterations) matrix.
%                   maxconstraint: Maximum constraint excess.
%
% CALLED FUNCTIONS:
%   oracle
%
% REFERENCES:
%   - Pontani, M. and Conway, B., "Particle Swarm Optimization Applied to
%       Space Trajectories", Journal of Guidance, Navigation and Control,
%       Vol. 33, No. 5, 2010, pp. 1429-1441.
%   - Schluter, M. and Gerdts, M., "The oracle penalty method", Journal of
%       Global Optimization, Vol. 47, No. 2, 2010, pp. 293-325.
% 
% AUTHOR:
%	Alessandro Peloni 30/10/2014, MATLAB, ppso.m
% 
% CHANGELOG:
%   04/11/2014, Alessandro Peloni: performances tested with the
%       non-constrained test functions of CEC2005.
%   06/11/2014, Alessandro Peloni: population size and maximum number of
%       iterations added to the optional inputs. Added the possibility to
%       reproduce the results by using the same structure for the random
%       number generator. Second loop vectorized.
%   09/12/2014, Alessandro Peloni: Display at screen and semilog plot fixed
%   13/02/2015, Alessandro Peloni: add the Oracle penalty method for the
%       constraints handling. Only one input function needed, which must
%       return both fitness function and nonlinear constraints.
%   16/02/2015, Alessandro Peloni: Possibility to use the vectorized
%       version of the code added, if the input function is vectorized too.
%   01/05/2015, Alessandro Peloni: Help improved.
%   05/06/2015, Alessandro Peloni: options.Save added.
%   25/06/2015, Alessandro Peloni: Parallelized mode implemented. General
%       improvements.
%   09/04/2016, Alessandro Peloni: General revision.
%   30/06/2016, Alessandro Peloni: Bug related with the "options" entry
%       fixed.
%   29/08/2016, Alessandro Peloni: General revision.
%   02/11/2016, Alessandro Peloni: createNewPop_fLimit removed and fixed an
%       issue related to the creation of a new population.
%
% -------------------------------------------------------------------------
%  Copyright (c) 2016, Alessandro Peloni
%  All rights reserved.
%  Distributed under the BSD 2-Clause License.
% 
%  Please report any bugs/suggestions to the forums at 
%  https://uk.mathworks.com/matlabcentral/fileexchange/58895-ppso
%  
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are met:
%  
%  * Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%  
%  * Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%  
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
%  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
%  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
%  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
%  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% -------------------------------------------------------------------------


% Initialize exitflag
exitflag = [];


%% INITIAL FUNCTION CHECKS

if nargin < 4
    error('ppso:InputCheck','Not enough inputs.');
end

%--------------------------------------------------------------------------
% Check the compulsory inputs
%--------------------------------------------------------------------------
if npars < 1
    error('ppso:InputCheck','There must be at least one parameter.');
end
if ub <= lb
    error('ppso:InputCheck','Lower boundaries must be lower than upper boundaries.');
end
if sum(isinf(-lb)) > 0 || sum(isinf(ub)) > 0
    error('ppso:InputCheck','Lower and upper boundaries cannot be +/-Inf at the moment.');
end
lb = lb(:);
ub = ub(:);


%--------------------------------------------------------------------------
% Check the optional inputs
%--------------------------------------------------------------------------
% Default values
npop = 10*npars;
niter = 500;
X0 = [];
StallIterLimit = niter;
createNewPop = 1;
TolFun = 1e-6;
TolCon = 1e-6;
auxdata = [];
Display = 1;
Plot = 1;
IntCon = [];
maxFcount = [];
s = [];
fvalOpt = [];
TerminationErr = [];
Residuals = 1;
Vectorized = 0;
Parallel = 0;
Save = 0;
filename = 'ppsoTemp.mat';
options = [];

% Input values
switch nargin-4
    case 1
        options = varargin{1};
        if isfield(options,'npop')
            npop = options.npop;
            if npop < 2
                error('ppso:InputCheck','There must be at least 2 particles in the swarm.');
            end
        end
        if isfield(options,'niter')
            niter = options.niter;
            if niter < 1 || mod(niter,1) ~= 0
                error('ppso:InputCheck','The number of iterations must be a positive integer.');
            end
        end
        if isfield(options,'X0')
            X0 = options.X0;
            if sum(size(X0)~=[npars,1]) > 0
                error('ppso:InputCheck','The initial guess X0 must be a column vector with npars rows');
            end
        end
        if isfield(options,'StallIterLimit')
            if options.StallIterLimit < 1
                warning('ppso:InputCheck','StallIterLimit must be a positive integer. Value set to default.')
            else
                StallIterLimit = options.StallIterLimit;
            end
        end
        if isfield(options,'createNewPop')
            switch options.createNewPop
                case 0
                    createNewPop = 0;
                case 1
                    createNewPop = 1;
                otherwise
                    warning('ppso:InputCheck','createNewPop must be a logical value. Value set to default.')
            end
        end
        if isfield(options,'TolFun')
            TolFun = options.TolFun;
        end
        if isfield(options,'TolCon')
            TolCon = options.TolCon;
        end
        if isfield(options,'Display')
            Display = options.Display;
        end
        if isfield(options,'Plot')
            Plot = options.Plot;
        end
        if isfield(options,'IntCon')
            IntCon = options.IntCon;
            if any(IntCon == 0) || any(IntCon > npars)
                warning('ppso:InputCheck','The value of options.IntCon is not valid. Continuing without integer constraints.')
                IntCon = [];
            end
        end
        if isfield(options,'maxFcount')
            maxFcount = options.maxFcount;
        end
        if isfield(options,'s')
            s = options.s;
        end
        if isfield(options,'fvalOpt')
            fvalOpt = options.fvalOpt;
        end
        if isfield(options,'TerminationErr')
            TerminationErr = options.TerminationErr;
        end
        if isfield(options,'Residuals')
            Residuals = options.Residuals;
            if Residuals~=1 && Residuals~=2 && Residuals~=Inf
                warning('ppso:InputCheck','The value of options.Residuals is not valid. Setting options.Residuals to 1.')
                Residuals = 1;
            end
        end
        if isfield(options,'Vectorized')
            Vectorized = options.Vectorized;
        end
        if isfield(options,'Parallel')
            Parallel = options.Parallel;
            if Vectorized
                warning('ppso:InputCheck','Both options.Vectorized and options.Parallel are set to 1. Setting options.Parallel to zero and using the vectorized mode.')
                Parallel = 0;
            end
        end
        if isfield(options,'Save')
            Save = options.Save;
            if isfield(options,'filename')
                filename = options.filename;
            end
        end
    case 2
        if ~isempty(varargin{1})
            options = varargin{1};
            if isfield(options,'npop')
                npop = options.npop;
                if npop < 2
                    error('ppso:InputCheck','There must be at least 2 particles in the swarm.');
                end
            end
            if isfield(options,'niter')
                niter = options.niter;
                if niter < 1 || mod(niter,1) ~= 0
                    error('ppso:InputCheck','The number of iterations must be a positive integer.');
                end
            end
            if isfield(options,'X0')
                X0 = options.X0;
                if sum(size(X0)~=[npars,1]) > 0
                    error('ppso:InputCheck','The initial guess X0 must be a column vector with npars rows');
                end
            end
            if isfield(options,'StallIterLimit')
                if options.StallIterLimit < 1
                    warning('ppso:InputCheck','StallIterLimit must be a positive integer. Value set to default.')
                else
                    StallIterLimit = options.StallIterLimit;
                end
            end
            if isfield(options,'createNewPop')
                switch options.createNewPop
                    case 0
                        createNewPop = 0;
                    case 1
                        createNewPop = 1;
                    otherwise
                        warning('ppso:InputCheck','createNewPop must be a logical value. Value set to default.')
                end
            end
            if isfield(options,'TolFun')
                TolFun = options.TolFun;
            end
            if isfield(options,'TolCon')
                TolCon = options.TolCon;
            end
            if isfield(options,'Display')
                Display = options.Display;
            end
            if isfield(options,'Plot')
                Plot = options.Plot;
            end
            if isfield(options,'IntCon')
                IntCon = options.IntCon;
                if any(IntCon == 0) || any(IntCon > npars)
                    warning('ppso:InputCheck','The value of options.IntCon is not valid. Continuing without integer constraints.')
                    IntCon = [];
                end
            end
            if isfield(options,'maxFcount')
                maxFcount = options.maxFcount;
            end
            if isfield(options,'s')
                s = options.s;
            end
            if isfield(options,'fvalOpt')
                fvalOpt = options.fvalOpt;
            end
            if isfield(options,'TerminationErr')
                TerminationErr = options.TerminationErr;
            end
            if isfield(options,'Residuals')
                Residuals = options.Residuals;
                if Residuals~=1 && Residuals~=2 && Residuals~=Inf
                    warning('ppso:InputCheck','The value of options.Residuals is not valid. Setting options.Residuals to 1.')
                    Residuals = 1;
                end
            end
            if isfield(options,'Vectorized')
                Vectorized = options.Vectorized;
            end
            if isfield(options,'Parallel')
                Parallel = options.Parallel;
                if Vectorized
                    warning('ppso:InputCheck','Both options.Vectorized and options.Parallel are set to 1. Setting options.Parallel to zero and using the vectorized mode.')
                    Parallel = 0;
                end
            end
            if isfield(options,'Save')
                Save = options.Save;
                if isfield(options,'filename')
                    filename = options.filename;
                end
            end
        end
        auxdata = varargin{2};
end


%% INITIALIZATION

% Lower/Upper boundary for the velocity of each particles
d = ub - lb;
d_rep  = repmat(d,1,npop);
lb_rep = repmat(lb,1,npop);
ub_rep = repmat(ub,1,npop);

% Particles of the swarm
x = NaN(npars,npop,niter); % Preallocation for speed
if ~isempty(s)
    s1 = rng(s); % User specified settings for the random number generator, in order to reproduce previous results
else
    s1 = rng('shuffle'); % Random seed for the initial population
end
x(:,:,1) = lb_rep + (ub_rep-lb_rep).*rand(npars,npop); % Initial random particles
if ~isempty(IntCon)
    for i = 1:length(IntCon)
        x(IntCon(i),:,1) = randi([lb(IntCon(i)), ub(IntCon(i))],1,npop); % Initial IntCon parameters have to be integer
    end
end
if ~isempty(X0)
    x(:,randi(npop),1) = X0;
end

% Velocities of the particles of the swarm
v = NaN(npars,npop,niter); % Preallocation for speed
v(:,:,1) = -d_rep + 2*d_rep.*rand(npars,npop); % Initial random velocities
v(IntCon,:,1) = round(v(IntCon,:,1));

% Penalty function
p = NaN(1,npop,niter); % Preallocation for speed

% Fitness function
f = NaN(1,npop,niter); % Preallocation for speed
fopt = NaN(niter,1); % Preallocation for speed
fmean = NaN(niter,1); % Preallocation for speed

% Oracle penalty method
Oracle = NaN(1,npop,niter); % Preallocation for speed
Oracle0 = 1e9; % Initialization of the Oracle parameter with a sufficiently large value
res = NaN(1,npop,niter); % Preallocation for speed
resMin = 1e9; % Set the first resMin to a very high value

% Best position ever visited by each particle
psi = NaN(npars,npop);

% Best and mean penalty function in each generation
popt = NaN(niter,1); % Preallocation for speed
pmean = NaN(niter,1); % Preallocation for speed

% Counter of functions (objective and nonlinear) evaluations
funccount = 0;

% Counter of stall generations
StallIter = 0;


if Plot % Plot fitness function
    figure('Name','Particle Swarm Optimization','NumberTitle','off')
    hold on
end



%% ALGORITHM

if Save
    save(filename,'x','f')
end

for iter = 1:niter
    
    %----------------------------------------------------------------------
    % VECTORIZED MODE
    %----------------------------------------------------------------------
    if Vectorized % Run the input fcn only once (vectorized mode)
        try
            [f(1,:,iter),c,ceq] = func(x(:,:,iter),auxdata);
            if isempty(c)
                c_tmp = false;
            else
                c_tmp = true;
            end
            if isempty(ceq)
                ceq_tmp = false;
            else
                ceq_tmp = true;
            end
        catch
            error('ppso:funcEvaluation','Some error in the user-defined function evaluation occurred.');
        end
        
        
    %----------------------------------------------------------------------
    % PARALLEL MODE
    %----------------------------------------------------------------------
    elseif Parallel % Run the fcn npop times in parallel
        
        if iter == 1 % Only for the first iter, check the structure of the inequality constraints
            try
                [f(1,1,iter),c_tmp,ceq_tmp] = func(x(:,1,1),auxdata);
                
                if isempty(c_tmp)
                    c = [];
                    c_tmp = false;
                else
                    c = NaN(length(c_tmp),npop);
                    c(:,1) = c_tmp;
                    c_tmp = true;
                end
                if isempty(ceq_tmp)
                    ceq = [];
                    ceq_tmp = false;
                else
                    ceq = NaN(length(ceq_tmp),npop);
                    ceq(:,1) = ceq_tmp;
                    ceq_tmp = true;
                end
                
                % Preallocate c and ceq for speed
                if c_tmp && ceq_tmp
                    c = [c_tmp NaN(npars,npop-1)];
                    ceq = [ceq_tmp NaN(npars,npop-1)];
                elseif c_tmp
                    c = [c_tmp NaN(npars,npop-1)];
                elseif ceq_tmp
                    ceq = [ceq_tmp NaN(npars,npop-1)];
                end
                
                parfor particle = 2:npop
                    try
                        if c_tmp && ceq_tmp
                            [f(1,particle,iter),c(:,particle),ceq(:,particle)] = feval(func,x(:,particle,iter),auxdata);
                        elseif c_tmp
                            [f(1,particle,iter),c(:,particle),~] = feval(func,x(:,particle,iter),auxdata);
                        elseif ceq_tmp
                            [f(1,particle,iter),~,ceq(:,particle)] = feval(func,x(:,particle,iter),auxdata);
                        else
                            [f(1,particle,iter),~,~] = feval(func,x(:,particle,iter),auxdata);
                        end
                    catch
                        error('ppso:funcEvaluation','Some error in the user-defined function evaluation occurred.');
                    end
                end
                
            catch
                error('ppso:funcEvaluation','Some error in the user-defined function evaluation occurred.');
            end
            
        else % Other iterations
            
            parfor particle = 1:npop
                try
                    if c_tmp && ceq_tmp
                        [f(1,particle,iter),c(:,particle),ceq(:,particle)] = feval(func,x(:,particle,iter),auxdata);
                    elseif c_tmp
                        [f(1,particle,iter),c(:,particle),~] = feval(func,x(:,particle,iter),auxdata);
                    elseif ceq_tmp
                        [f(1,particle,iter),~,ceq(:,particle)] = feval(func,x(:,particle,iter),auxdata);
                    else
                        [f(1,particle,iter),~,~] = feval(func,x(:,particle,iter),auxdata);
                    end
                catch
                    error('ppso:funcEvaluation','Some error in the user-defined function evaluation occurred.');
                end
            end
            
        end
        
        
    %----------------------------------------------------------------------
    % NORMAL MODE
    %----------------------------------------------------------------------
    else % Non-vectorized, non-parallelized code
        
        
        if iter == 1 % Only for the first iter, check the structure of the inequality constraints
            try
                
                [f(1,1,1),c_tmp,ceq_tmp] = func(x(:,1,1),auxdata);
                
                if isempty(c_tmp)
                    c = [];
                    c_tmp = false;
                else
                    c = NaN(length(c_tmp),npop);
                    c(:,1) = c_tmp;
                    c_tmp = true;
                end
                if isempty(ceq_tmp)
                    ceq = [];
                    ceq_tmp = false;
                else
                    ceq = NaN(length(ceq_tmp),npop);
                    ceq(:,1) = ceq_tmp;
                    ceq_tmp = true;
                end
                
                % Preallocate c and ceq for speed
                if c_tmp && ceq_tmp
                    c = [c_tmp NaN(npars,npop-1)];
                    ceq = [ceq_tmp NaN(npars,npop-1)];
                elseif c_tmp
                    c = [c_tmp NaN(npars,npop-1)];
                elseif ceq_tmp
                    ceq = [ceq_tmp NaN(npars,npop-1)];
                end
                
                for particle = 2:npop
                    try
                        if c_tmp && ceq_tmp
                            [f(1,particle,iter),c(:,particle),ceq(:,particle)] = func(x(:,particle,iter),auxdata);
                        elseif c_tmp
                            [f(1,particle,iter),c(:,particle),~] = func(x(:,particle,iter),auxdata);
                        elseif ceq_tmp
                            [f(1,particle,iter),~,ceq(:,particle)] = func(x(:,particle,iter),auxdata);
                        else
                            [f(1,particle,iter),~,~] = func(x(:,particle,iter),auxdata);
                        end
                    catch
                        error('ppso:funcEvaluation','Some error in the user-defined function evaluation occurred.');
                    end
                    
                end
                
            catch
                error('ppso:funcEvaluation','Some error in the user-defined function evaluation occurred.');
            end
        
        else % Other iterations
            
            for particle = 1:npop
                try
                    if c_tmp && ceq_tmp
                        [f(1,particle,iter),c(:,particle),ceq(:,particle)] = func(x(:,particle,iter),auxdata);
                    elseif c_tmp
                        [f(1,particle,iter),c(:,particle),~] = func(x(:,particle,iter),auxdata);
                    elseif ceq_tmp
                        [f(1,particle,iter),~,ceq(:,particle)] = func(x(:,particle,iter),auxdata);
                    else
                        [f(1,particle,iter),~,~] = func(x(:,particle,iter),auxdata);
                    end
                catch
                    error('ppso:funcEvaluation','Some error in the user-defined function evaluation occurred.');
                end
                
            end
            
        end
        
    end % End of if on Vectorization or Parallelization
    
    funccount = funccount + npop;
    f(1,isnan(f(1,:,iter)),iter) = 1/eps;
    
    if iter == 1
        fopt_loc = f(1,:,1);
        psi = x(:,:,1);
        [fopt(1),I] = min(fopt_loc);
        Y = x(:,I,1);
    end
    
    % Computation of residuals
    if isempty(c) && isempty(ceq) % Non constrained problem
        res(1,:,iter) = zeros(1,npop); % Residuals are manually set to zero
    else
        switch Residuals
            case 1 % l1 norm
                res(1,:,iter) = sum(abs(ceq)) + sum(max(0,c));
            case 2 % l2 norm
                res(1,:,iter) = sqrt(sum(abs(ceq).^2) + sum(max(0,c).^2));
            case Inf % l_inf norm
                if isempty(ceq)
                    res(1,:,iter) = max(0,max(c));
                elseif isempty(c)
                    res(1,:,iter) = max(abs(ceq));
                else
                    res(1,:,iter) = max(max(abs(ceq)),max(0,max(c)));
                end
        end
    end
    
    
    
    if ~isempty(c) || ~isempty(ceq)
        
        % Oracle update
        if iter ~= 1
            tmp = f(1,:,iter-1)<Oracle(1,:,iter-1) & res(1,:,iter-1)<=TolCon;
            Oracle(1,tmp,iter) = f(1,tmp,iter-1);
            Oracle(1,~tmp,iter) = Oracle(1,~tmp,iter-1);
        else
            Oracle(1,:,iter) = Oracle0*ones(1,npop);
        end
        
        % Penalty function
        p(1,:,iter) = oracle(f(1,:,iter),Oracle(1,:,iter),res(1,:,iter),TolCon); % There is no point to compute the penalty function if there are no constraints
        p(1,isnan(p(1,:,iter)),iter) = 1/eps;
        
    end
    
    
    % Determine the minimum residual of the current iteration
    resMin = min(resMin, min(res(1,:,iter)));
    
    
    %----------------------------------------------------------------------
    % Determine the best position ever visited by each particle
    %----------------------------------------------------------------------
    if isempty(c) && isempty(ceq)
        % No nonlinear constraints
        if iter~=1
            l = find(f(1,:,iter) < fopt_loc);
            psi(:,l) = x(:,l,iter);
            fopt_loc(1,l) = f(1,l,iter);
        end
    else
        if resMin > TolCon
            % Nonlinear constraint violation greater than the tolerance
            [~,l] = min(p(1,:,1:iter),[],3);
            tempInd = sub2ind([npop,iter],1:npop,l);
            psi(:,:) = x(:,tempInd); % Historical best position of each particle
        else
            % At least one particle is historically feasible
            feas = find(sum(res(1,:,1:iter)<=TolCon,3)>0); % At least one time is feasible
            infeas = res(1,:,1:iter)>TolCon; % All the infeasibilities
            neverFeas = find(sum(res(1,:,1:iter)>TolCon,3)==iter); % The particles historically never feasible
            f(infeas) = NaN; % All the infeasibilities do not count for the minimum
            [~,lFeas] = min(f(1,:,1:iter),[],3); % Best historical position of each particle
            [~,lNeverFeas] = min(p(1,neverFeas,1:iter),[],3); % Best historical position of each particle which was never feasible
            l(feas) = lFeas(feas);
            l(neverFeas) = lNeverFeas;
            
            tempInd = sub2ind([npop,iter],1:npop,l);
            tempIndFeas = tempInd(feas);
            psi(:,:) = x(:,tempInd); % Historical best position of each particle
            
        end
    end
    
    
    %----------------------------------------------------------------------
    % Determine the global best position ever visited by the entire swarm
    %----------------------------------------------------------------------
    if isempty(c) && isempty(ceq)
        if iter~=1
            [fopt_min,I] = min(fopt_loc);
            if fopt_min < fopt(iter-1)
                fopt(iter) = fopt_min;
                Y = x(:,I,iter);
            else
                fopt(iter) = fopt(iter-1);
            end
        end
    else
        if resMin > TolCon
            [popt(iter),I] = min(p(tempInd));
        else
            [fopt(iter),I] = min(f(tempIndFeas));
            I = tempIndFeas(I);
        end
        Y = x(:,I);
    end
    
    
    %----------------------------------------------------------------------
    % Velocity vector
    %----------------------------------------------------------------------
    % Inertial weight
    r1 = rand(npars,npop);
    cI = (1+r1)/2;
    
    % Cognitive weight
    r2 = rand(npars,npop);
    cC = 1.49445*r2;
    
    % Social weight
    r3 = rand(npars,npop);
    cS = 1.49445*r3;
    
    % Update velocity vector
    v(:,:,iter+1) = cI.*v(:,:,iter) + cC.*(psi-x(:,:,iter)) + cS.*(repmat(Y,1,npop)-x(:,:,iter));
    v(IntCon,:,iter+1) = round(v(IntCon,:,iter+1)); % Velocities of IntCon parameters have to be integer
    
    % Put velocities into bounds
    v_temp = v(:,:,iter+1); % Use a temporary 2D matrix in order to avoid problems with indexing
    % --- Lower boundaries --- %
    temp = v_temp < -d_rep;
    v_temp(temp) = -d_rep(temp);
    v(:,:,iter+1) = v_temp;
    % --- Upper boundaries --- %
    temp = v_temp > d_rep;
    v_temp(temp) = d_rep(temp);
    v(:,:,iter+1) = v_temp;
    
    
    %----------------------------------------------------------------------
    % Update position vector
    %----------------------------------------------------------------------
    % Note. The (iter+1)-th element of the velocity has been used,
    % otherwise the algorithm does not converge.
    x(:,:,iter+1) = x(:,:,iter) + v(:,:,iter+1);
    %         x(:,:,iter+1) = x(:,:,iter) + v(:,:,iter);
    
    % Put positions into bounds
    x_temp = x(:,:,iter+1); % Use a temporary 2D matrix in order to avoid problems with indexing
    v_temp = v(:,:,iter+1); % Use a temporary 2D matrix in order to avoid problems with indexing
    % --- Lower boundaries --- %
    temp = x_temp < lb_rep;
    x_temp(temp) = lb_rep(temp);
    v_temp(temp) = 0;
    x(:,:,iter+1) = x_temp;
    v(:,:,iter+1) = v_temp;
    % --- Upper boundaries --- %
    temp = x_temp > ub_rep;
    x_temp(temp) = ub_rep(temp);
    v_temp(temp) = 0;
    x(:,:,iter+1) = x_temp;
    v(:,:,iter+1) = v_temp;
    
    
    %----------------------------------------------------------------------
    % Check TolFun and stall in generations
    %----------------------------------------------------------------------
    if iter ~= 1
        if resMin <= TolCon
            if abs(fopt(iter)-fopt(iter-1)) <= TolFun*abs(fopt(iter-1))
                StallIter = StallIter + 1;
            else
                StallIter = 0;
            end
        else
            if abs(popt(iter)-popt(iter-1)) <= TolFun*abs(popt(iter-1))
                StallIter = StallIter + 1;
            else
                StallIter = 0;
            end
        end
    else
        StallIter = 0;
    end
    
    
    %----------------------------------------------------------------------
    % Display and plot iterations
    %----------------------------------------------------------------------
    if Display || Plot
        pmean(iter) = mean(p(1,:,iter));
        feas = res(1,:,iter)<=TolCon;
        fmean(iter) = mean(f(1,feas,iter));
    end
    if Display
        if mod(iter,30) == 1 % First line every 20 iterations
            fprintf('\n Iteration     f-count       Best f(x)        Mean f(x)        Residual      Stall iterations\n');
            fprintf('---------------------------------------------------------------------------------------------------\n');
        end
        if resMin <= TolCon
            fprintf('%5.0f         %7.0f    %12.4g    %12.4g    %12.4g    %5.0f\n',iter,funccount,fopt(iter),fmean(iter),resMin,StallIter)
        else
            fprintf('%5.0f         %7.0f    %12.4g    %12.4g    %12.4g    %5.0f\n',iter,funccount,popt(iter),pmean(iter),resMin,StallIter)
        end
    end
    if Plot
        if iter ~= 1 % Do not plot the very first point
            if resMin <= TolCon
                if isempty(fvalOpt) || fvalOpt>=0
                    subplot(2,1,1)
                    semilogy(iter,fopt(iter),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',2)
                    hold on
                    semilogy(iter,fmean(iter),'d','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
                    title('Best and Mean fitness function')
                    
                    subplot(2,1,2)
                    if mod(iter,20)==0
                        hold off
                    end
                    semilogy(iter,fopt(iter),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',2)
                    hold on
                    title('Best fitness function - 20 iterations view')
                else
                    subplot(2,1,1)
                    semilogy(iter,fopt(iter)-fvalOpt,'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',2)
                    hold on
                    semilogy(iter,fmean(iter)-fvalOpt,'d','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
                    title('Best and Mean fitness function (fvalOpt as reference value)')
                    
                    subplot(2,1,2)
                    if mod(iter,20)==0
                        hold off
                    end
                    semilogy(iter,fopt(iter)-fvalOpt,'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',2)
                    hold on
                    title('Best fitness function - 20 iterations view (fvalOpt as reference value)')
                end
                    
            else
                
                if isempty(fvalOpt) || fvalOpt>=0
                    subplot(2,1,1)
                    semilogy(iter,popt(iter),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',2)
                    hold on
                    semilogy(iter,pmean(iter),'d','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
                    title('Best and Mean fitness function')
                    
                    subplot(2,1,2)
                    if mod(iter,20)==0
                        hold off
                    end
                    semilogy(iter,popt(iter),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',2)
                    hold on
                    title('Best fitness function - 20 iterations view')
                else
                    subplot(2,1,1)
                    semilogy(iter,popt(iter)-fvalOpt,'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',2)
                    hold on
                    semilogy(iter,pmean(iter)-fvalOpt,'d','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
                    title('Best and Mean fitness function (fvalOpt as reference value)')
                    
                    subplot(2,1,2)
                    if mod(iter,20)==0
                        hold off
                    end
                    semilogy(iter,popt(iter)-fvalOpt,'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',2)
                    hold on
                    title('Best fitness function - 20 iterations view (fvalOpt as reference value)')
                end
            end
            drawnow % to update online the plots
            
        end
    end
    
    
    %----------------------------------------------------------------------
    % Break loop criteria
    %----------------------------------------------------------------------
    
    % Break the loop if StallIterLimit is reached
    if resMin<=TolCon && StallIter == StallIterLimit
        message = ['Optimization terminated: average change in the fitness function ' ...
            'less than options.TolFun over options.StallIterLimit number of consecutive iterations.'];
        exitflag = 1;
        break
    end
    
    % Break the loop if maximum number of function evaluations is reached
    if ~isempty(maxFcount) && funccount > maxFcount
        if iter ~= 1 && resMin <= TolCon && ...
                abs(fopt(iter)-fopt(iter-1)) < TolFun*abs(fopt(iter-1))
            message = ['Optimization terminated: maximum number of function evaluation reached, '...
                'feasible solution found and relative change in the fitness function less than TolFun.'];
            exitflag = 1;
            break
        elseif iter ~= 1 && resMin > TolCon
            message = ['Optimization terminated: maximum number of function evaluation reached '...
                'but nonlinear constraints are not satisfied.'];
            exitflag = -2;
            break
        else
            message = ['Optimization terminated: maximum number of function evaluations reached ' ...
                'and feasible solution found.'];
            exitflag = 0;
            break
        end
    end
    
    % Break the loop if the error between the objective function and the
    % known optimum is less than TerminationErr
    if ~isempty(fvalOpt) && ~isempty(TerminationErr) && fopt(iter)-fvalOpt < TerminationErr
        message = 'Optimization terminated: optimal solution found.';
        exitflag = 2;
        break
    end
    
    
    
    %----------------------------------------------------------------------
    % New population if createNewPop==1 & StallIter>StallIterLimit/3
    %----------------------------------------------------------------------
    if createNewPop
        if StallIter>StallIterLimit/3
            
            StallIter = 0;
            
            % Position
            x(:,1,iter+1) = Y;
            x(:,2:end,iter+1) = lb_rep(:,2:end) + (ub_rep(:,2:end)-lb_rep(:,2:end)).*rand(npars,npop-1); % Random particles
            if ~isempty(IntCon)
                for i = 1:length(IntCon)
                    x(IntCon(i),2:end,iter+1) = randi([lb(IntCon(i)), ub(IntCon(i))],1,npop-1); % IntCon parameters have to be integer
                end
            end
            
            % Velocity
            v(:,:,iter+1) = -d_rep + 2*d_rep.*rand(npars,npop); % Random velocities
            v(IntCon,:,iter+1) = round(v(IntCon,:,iter+1));
            
        end
    end % end of createNewPop loop
    
    
    %----------------------------------------------------------------------
    % Save iter, x and f into a .mat file if required
    %----------------------------------------------------------------------
    if Save
        % Saving with '-v6' version creates a larger file, because no
        % compression is used. However, the saving process is faster.
%         save(filename,'iter','x','f','-v6')
        if iter > 1
            if resMin <= TolCon
                if fopt(iter) < fopt(iter-1) % Save only when a better value for f is found
                    save(filename,'iter','x','f','-append')
                end
            else
                if popt(iter) < popt(iter-1) % Save only when a better value for p is found
                    save(filename,'iter','x','p','-append')
                end
            end
        end
    end
    
end % end of iterations loop


%--------------------------------------------------------------------------
% Delete the extra data stored in preallocation
%--------------------------------------------------------------------------
if iter < niter
    x(:,:,iter+1:end) = [];
    v(:,:,iter+1:end) = [];
    f(:,:,iter+1:end) = [];
    p(:,:,iter+1:end) = [];
    fopt(iter+1:end,:) = [];
    Oracle(:,:,iter+1:end) = [];
    res(:,:,iter+1:end) = [];
end


% if Plot
%     legend(ax1,'Best fitness function', 'Mean fitness function')
%     legend(ax2,'Best fitness function only')
%     figure('Name','Particle Swarm Optimization - semilog view','NumberTitle','off')
%     if Jopt(iter) > 0
%         semilogy(1:iter,Jopt,'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',2)
%         legend('Best fitness function - semilog view')
%     else
%         semilogy(1:iter,Jopt-Jopt(iter),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',2)
%         legend('Best traslated fitness function - semilog view')
%     end
% end



%% WORK ON OUTPUTS

%--------------------------------------------------------------------------
% Optimal solution
%--------------------------------------------------------------------------
if resMin <= TolCon
    fval = fopt(iter);
else
    fval = popt(iter);
end
xopt = Y;

if c_tmp || ceq_tmp
    [~,c_f,ceq_f] = func(xopt,auxdata);
    funccount = funccount + 1;
end


%--------------------------------------------------------------------------
% Analyze the solution found
%--------------------------------------------------------------------------
if iter == niter && isempty(exitflag)
    if resMin > TolCon % No feasible point found
        message = 'Optimization terminated: no feasible point has been found.';
        exitflag = -2;
    else
        if abs(fopt(iter)-fopt(iter-1)) < TolFun*abs(fopt(iter-1))
            message = ['Optimization terminated: maximum number of iterations reached, ' ...
                'feasible solution found and relative change in the fitness function less than TolFun.'];
            exitflag = 1;
        else
            message = ['Optimization terminated: maximum number of iterations reached ' ...
                'and feasible solution found.'];
            exitflag = 0;
        end
    end
end

disp(' ')
disp(message)


%--------------------------------------------------------------------------
% Output structure
%--------------------------------------------------------------------------
output.options = options;
output.rngstate = s1;
output.iterations = iter;
output.funccount = funccount;
output.message = message;
output.positions = x;
output.velocities = v;
output.fitnessFcn = f;
output.fitnessFcnOpt = fopt;
output.penalty = p;
output.Oracle = Oracle;
output.residuals = res;
if c_tmp && ceq_tmp
    cMax = max(0,max(c_f));
    output.maxconstraint = max(cMax,max(abs(ceq_f)));
elseif c_tmp
    cMax = max(0,max(c_f));
    output.maxconstraint = cMax;
elseif ceq_tmp
    output.maxconstraint = max(abs(ceq_f));
else
    output.maxconstraint = 0;
end


return

