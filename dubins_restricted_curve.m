%DUBINS_CURVE   Find the restricted Dubins path (shortest curve) between 
%   two points when restricted by orientation at all time.
%   PATH = DUBINS_CURVE(P1, P2, r, psi, delta, stepsize) finds the shortest
%   curve that connects two points in the Euclidean plane with a constraint
%   of the orientation of the path. The start and finish orientations P1 
%   and P2 are defined as [x, y, theta]. The psi and delta corresponds
%   to the restricted angles. Psi being the resitricted orientation and 
%   delta the range, with the range of restricted angle being
%   [ psi - delta/2, psi + delta/2]
%   The turning radius (r), stepsize and delta will be defined 
%   automatically if their value is <= 0. 
%
%   The output PATH is an [mx3] array consisting of m rows of [x, y, theta] values. 
%
%   PATH = DUBINS_CURVE(P1, P2, r, psi, delta, stepsize, quiet) performs the same as above,
%   however if quiet == true, then no plots of command window output will be
%   generated. Ommitting this input will result in quiet = false/0.
%
%   This function handles the interface to dubins_restricted_core.m to give a more
%   intuitive tool for finding the Dubins path. 

% Reference:
%       https://github.com/AndrewWalker/Dubins-Curves#shkel01
%       Shkel, A. M. and Lumelsky, V. (2001). "Classification of the Dubins
%                  set". Robotics and Autonomous Systems 34 (2001) 179�V202
%       https://github.com/UlysseVautier/MATLAB/Dubins-Curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Source: Andrew Walker, Ulysse Vautier
% Date: 2019.01.01
% contact: ulysse.vautier [at] gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2019, Ulysse Vautier                                             % 
%                                                                                %
% Permission is hereby granted, free of charge, to any person obtaining a copy   %
% of this software and associated documentation files (the "Software"), to deal  %
% in the Software without restriction, including without limitation the rights   %
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      %
% copies of the Software, and to permit persons to whom the Software is          %  
% furnished to do so, subject to the following conditions:                       %
%                                                                                %
% The above copyright notice and this permission notice shall be included in     %
% all copies or substantial portions of the Software.                            %
%                                                                                %
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     %
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       %
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    %
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         %
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  %
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN      %
% THE SOFTWARE.                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ path ] = dubins_restricted_curve( p1, p2, r, psi, delta, stepsize, quiet )
%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % there are 8 types of dubin's curve, only one will have minimum cost
    % LSLSL = 1;
	% LSRSL = 2;
	% RSLSR = 3;
	% RSRSR = 4;
	% RSLSL = 5;
    % RSRSL = 6;
	% LSLSR = 7;
    % LSRSR = 8;
    
    % The three segment types a path can be made up of
    % L_SEG = 1;
    % S_SEG = 2;
    % R_SEG = 3;

    % The segment types for each of the Path types
    %{
    DIRDATA = [ L_SEG, S_SEG, L_SEG, S_SEG, L_SEG ;...
                L_SEG, S_SEG, R_SEG, S_SEG, L_SEG ;...
                R_SEG, S_SEG, L_SEG, S_SEG, R_SEG ;...
                R_SEG, S_SEG, R_SEG, S_SEG, R_SEG ;...
                R_SEG, S_SEG, L_SEG, S_SEG, L_SEG ;...
                R_SEG, S_SEG, R_SEG, S_SEG, L_SEG ;...
                L_SEG, S_SEG, L_SEG, S_SEG, R_SEG ;...
                L_SEG, S_SEG, R_SEG, S_SEG, R_SEG ]; 
    %}
      
    % the return parameter from dubins_core
    % param.p_init = p1;              % the initial configuration
    % param.SEG_param = [0, 0, 0];    % the lengths of the three segments
    % param.angle = [0, 0];          % angles of departure and arrival from C1 and C3
    % param.r = r;                    % model forward velocity / model angular velocity turning radius
    % param.type = -1;                % path type. one of LSL, LSR, ... 
    %param.achievable = 0;           % which paths are achievable
    %%%%%%%%%%%%%%%%%%%%%%%%% END DEFINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Handle inputs.
    if nargin < 4
        error('Function requires at least four inputs.');
    end
    if nargin < 5
        delta = 0.1;
    end
    if nargin < 6
        stepsize = 0;
    end
    if nargin < 7 
        quiet = 0;  %Default/undefined is not quiet
    end
    
    if ~quiet
        close(findobj('type','figure','name','Dubins curve'));
        tic;
    end
    
    % main function
    param = dubins_restricted_core(p1, p2, r, psi, delta);
    if stepsize <= 0
        stepsize = dubins_length(param)/1000;
    end
    path = dubins_path_sample_many(param, stepsize);
    
    % plot if not quiet
    if ~quiet
        disp(param)
        disp('dubins calculation time'); toc;
        % plotting
        tic;    % most of the time is spent on plotting
        figure('name','Dubins curve');
        DrawPacman(p1, p2, r, psi, delta);
        plot(path(:,1), path(:,2), 'LineWidth' , 2); axis equal; hold on
        scatter(p1(1), p1(2), 45, '*','r','LineWidth',1); hold on;
        scatter(p2(1), p2(2), 45, 'square','b','LineWidth',1); hold on;
        text(p1(1), p1(2),'start','HorizontalAlignment','center');
        text(p2(1), p2(2),'end','VerticalAlignment','top');
        disp('plot drawing time'); toc;
    end
end

function path = dubins_path_sample_many( param, stepsize)
    if param.STATUS < 0
        path = 0;
        return
    end
    length = dubins_length(param);
    path = -1 * ones(floor(length/stepsize), 3);
    x = 0;
    i = 1;
    while x <= length
        path(i, :) = dubins_path_sample( param, x );
        x = x + stepsize;
        i = i + 1;
    end
    return
end

function l = dubins_length(param)
    l = param.SEG_param(1);
    for i=2:length(param.SEG_param)
        l = l + param.SEG_param(i);
    end
    l = l * param.r;
end


%{
 * Calculate the configuration along the path, using the parameter t
 *
 * @param path - an initialised path
 * @param t    - a length measure, where 0 <= t < dubins_path_length(path)
 * @param q    - the configuration result
 * @returns    - -1 if 't' is not in the correct range
%}
function end_pt = dubins_path_sample(param, t)
    if( t < 0 || t >= dubins_length(param) || param.STATUS < 0)
        end_pt = -1;
        return;
    end

    % tprime is the normalised variant of the parameter t
    tprime = t / param.r;

    % In order to take rho != 1 into account this function needs to be more complex
    % than it would be otherwise. The transformation is done in five stages.
    %
    % 1. translate the components of the initial configuration to the origin
    % 2. generate the target configuration
    % 3. transform the target configuration
    %      scale the target configuration
    %      translate the target configration back to the original starting point
    %      normalise the target configurations angular component

    % The translated initial configuration
    p_init = [0, 0, param.p_init(3) ];
    
    %%%%%%%%%%%%%%%%%%%%%%%%% DEFINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The three segment types a path can be made up of
    L_SEG = 1;
    S_SEG = 2;
    R_SEG = 3;

    % The segment types for each of the Path types
    DIRDATA = [ L_SEG, S_SEG, L_SEG, S_SEG, L_SEG ;...
                L_SEG, S_SEG, R_SEG, S_SEG, L_SEG ;...
                R_SEG, S_SEG, L_SEG, S_SEG, R_SEG ;...
                R_SEG, S_SEG, R_SEG, S_SEG, R_SEG ;...
                R_SEG, S_SEG, L_SEG, S_SEG, L_SEG ;...
                R_SEG, S_SEG, R_SEG, S_SEG, L_SEG ;...
                L_SEG, S_SEG, L_SEG, S_SEG, R_SEG ;...
                L_SEG, S_SEG, R_SEG, S_SEG, R_SEG ]; 
            
    DIRDATASTANDARD = [ L_SEG, S_SEG, L_SEG ;...
                R_SEG, S_SEG, R_SEG ;...
                R_SEG, S_SEG, L_SEG ;...
                L_SEG, S_SEG, R_SEG ]; 
    %%%%%%%%%%%%%%%%%%%%%%%%% END DEFINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generate the target configuration
    
    types = 0;
    if(length(param.SEG_param) < 5)
        types = DIRDATASTANDARD(round(param.type/2), :);
    else
        types = DIRDATA(param.type, :);
    end
    
    param1 = 0;
    mid_pt1 = p_init;
    for i=1:length(param.SEG_param)-1
        param1 = [param1 param.SEG_param(i)];
        mid_pt1 = [mid_pt1 ;dubins_segment( param1(i+1), mid_pt1(i,:), types(i) )];
    end
    
    end_pt = dubins_segment( tprime-sum(param1(1:end)), mid_pt1(end,:),  types(end) );
    for i=1:length(param.SEG_param)-1
        if( tprime < sum(param1(2:i+1)) ) 
            end_pt = dubins_segment( tprime-sum(param1(1:i)), mid_pt1(i,:),  types(i) );
            break;
        end
    end
    
    %{
    types = DIRDATA(param.type, :);
    param1 = param.SEG_param(1);
    param2 = param.SEG_param(2);
    param3 = param.SEG_param(3);
    param4 = param.SEG_param(4);
    mid_pt1 = dubins_segment( param1, p_init, types(1) );
    mid_pt2 = dubins_segment( param2, mid_pt1,  types(2) );
    mid_pt3 = dubins_segment( param3, mid_pt2,  types(3) );
    mid_pt4 = dubins_segment( param4, mid_pt3,  types(4) );
    
    % Actual calculation of the position of tprime within the curve
    if( tprime < sum(param1(2:2)) ) 
        end_pt = dubins_segment( tprime, mid_pt1(1,:),  types(1) );
    elseif( tprime < sum(param1(2:3)) ) 
        end_pt = dubins_segment( tprime-sum(param1(1:2)), mid_pt1(2,:),  types(2) );
    elseif( tprime < sum(param1(2:4)) ) 
        end_pt = dubins_segment( tprime-sum(param1(1:3)), mid_pt1(3,:),  types(3) );
    elseif( tprime < sum(param1(2:5)) ) 
        end_pt = dubins_segment( tprime-sum(param1(1:4)), mid_pt1(4,:),  types(4) );
    else 
        end_pt = dubins_segment( tprime-sum(param1(1:5)), mid_pt1(5,:),  types(5) );
    end
    %}
    
    % scale the target configuration, translate back to the original starting point
    end_pt(1) = end_pt(1) * param.r + param.p_init(1);
    end_pt(2) = end_pt(2) * param.r + param.p_init(2);
    end_pt(3) = mod(end_pt(3), 2*pi);
    return;
end

%{
 returns the parameter of certain location according to an inititalpoint,
 segment type, and its corresponding parameter
%}
function seg_end = dubins_segment(seg_param, seg_init, seg_type)
    L_SEG = 1;
    S_SEG = 2;
    R_SEG = 3;
    if( seg_type == L_SEG ) 
        seg_end(1) = seg_init(1) + sin(seg_init(3)+seg_param) - sin(seg_init(3));
        seg_end(2) = seg_init(2) - cos(seg_init(3)+seg_param) + cos(seg_init(3));
        seg_end(3) = seg_init(3) + seg_param;
    elseif( seg_type == R_SEG )
        seg_end(1) = seg_init(1) - sin(seg_init(3)-seg_param) + sin(seg_init(3));
        seg_end(2) = seg_init(2) + cos(seg_init(3)-seg_param) - cos(seg_init(3));
        seg_end(3) = seg_init(3) - seg_param;
    elseif( seg_type == S_SEG ) 
        seg_end(1) = seg_init(1) + cos(seg_init(3)) * seg_param;
        seg_end(2) = seg_init(2) + sin(seg_init(3)) * seg_param;
        seg_end(3) = seg_init(3);
    end
end

