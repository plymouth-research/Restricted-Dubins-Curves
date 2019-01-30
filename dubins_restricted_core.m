% This function will find a dubins curve that connect two points when
% restricted to a certain angle
% Input: p1 and p2 are two row vector, e.g. [x, y, theta], that defines a 
% 2-D point and p1ing/ending direction.
% r is the turning rate of the vehicle.
% psi is the restricted angle
% delta is the range of restricted angle : [psi-delta/2,psi+delta/2]
% Output: the function returns 4 field:
%       type: defines one of the 8 shape of the CSCSC curve
%       r: turning radius, also the scaling factor for the curve paramater
%       SEG_param: defines the parameter of the three segments in row vector
% Reference:
%       https://github.com/UlysseVautier/MATLAB/Dubins-Curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Source: Ulysse Vautier
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


function [ param ] = dubins_restricted_core( p1, p2, r, psi, delta )
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
    
    % the return parameter
    param.p_init = p1;              % the initial configuration
    param.SEG_param = [0, 0, 0, 0, 0];    % the lengths of the three segments
    param.angle = [0, 0];          % angles of departure and arrival from C1 and C3
    param.r = r;                    % model forward velocity / model angular velocity turning radius
    param.type = -1;                % path type. one of LSL, LSR, ... 
    param.STATUS = 0;
    param.achievable = 0;           % which paths are achievable

    %%%%%%%%%%%%%%%%%%%%%%%%% p1 %%%%%%%%%%%%%%%%%%%%%%%%%
    % First, basic properties and normalization of the problem
    dx = p2(1) - p1(1);
    dy = p2(2) - p1(2);
    D = sqrt( dx^2 + dy^2 );
    d = D / r;                  % distance is shrunk by r, this make lengh calculation very easy
    if( r <= 0 )
        param.STATUS = -1;
        return;
    end
    
    theta = mod(atan2( dy, dx ), 2*pi);
    alpha = mod((p1(3) - theta), 2*pi);
    beta  = mod((p2(3) - theta), 2*pi);
    
    %% Is Dubins Achievable
    % Helpers
    u.left = 1;
    u.right = -1;
    restriction = wrapTo2Pi([psi-delta,psi+delta]);
    
    % Circles Origins
    Ol = [p1(1:2);p2(1:2)] + [r*cos(p1(3)+u.left*pi/2), r*sin(p1(3)+u.left*pi/2); r*cos(p2(3)+u.left*pi/2), r*sin(p2(3)+u.left*pi/2)];
    Or = [p1(1:2);p2(1:2)] + [r*cos(p1(3)+u.right*pi/2), r*sin(p1(3)+u.right*pi/2); r*cos(p2(3)+u.right*pi/2), r*sin(p2(3)+u.right*pi/2)];
    
    % Bitangent points
    %Outerleft = zeros(1,2);
    %Outerright = zeros(1,2);
    Innerleft = 2*pi*ones(1,2);
    Innerright = 2*pi*ones(1,2);
    
    Outerleft = wrapTo2Pi(atan2(Ol(2,2)-Ol(1,2), Ol(2,1)-Ol(1,1)));
    Outerleft(2) = Outerleft(1);

    if(norm(Or(2,:)-Ol(1,:)) > 2*r)
        Innerleft = wrapTo2Pi(atan2(Or(2,2)-Ol(1,2), Or(2,1)-Ol(1,1))-u.left*acos(2*r/norm(Or(2,:)-Ol(1,:))));
        Innerleft(2) = wrapTo2Pi(Innerleft(1)-pi);
    end

    Outerright = wrapTo2Pi(atan2(Or(2,2)-Or(1,2), Or(2,1)-Or(1,1)));
    Outerright(2) = Outerright(1);

    if(norm(Ol(2,:)-Or(1,:)) > 2*r)
        Innerright = wrapTo2Pi(atan2(Ol(2,2)-Or(1,2), Ol(2,1)-Or(1,1))-u.right*acos(2*r/norm(Ol(2,:)-Or(1,:))));
        Innerright(2) = wrapTo2Pi(Innerright(1)-pi);
    end
    
    param.achievable = is_dubins_achievable(p1, p2, Ol, Or, Outerleft, Outerright, Innerleft, Innerright, restriction, r);

    disp(param.achievable)
    %% Second, we find all possible curves
    if(Innerright(1) ~= 2*pi)
        Innerright = wrapTo2Pi(Innerright + [pi/2,-pi/2]);
    end
    if(Innerleft(1) ~= 2*pi)
        Innerleft = wrapTo2Pi(Innerleft + [pi/2,-pi/2]);
    end
    
    best_word = -1;
    best_cost = -1;
    if(~param.achievable.LSL)
        test_param(1,:) = dubins_CSCSC(p1, p2, r, Outerleft, restriction, u.left, u.left, u.left);
        test_param(2,:) = dubins_CSCSC(p1, p2, r, Outerleft, restriction, u.left, u.right, u.left);
    else
        test_param(1,:) = dubins_LSL(alpha,beta,d);
        test_param(2,:).SEG = -1;
    end
    
    if(~param.achievable.RSR)
        test_param(3,:) = dubins_CSCSC(p1, p2, r, Outerright, restriction, u.right, u.left, u.right);
        test_param(4,:) = dubins_CSCSC(p1, p2, r, Outerright, restriction, u.right, u.right, u.right);
    else
        test_param(3,:) = dubins_RSR(alpha,beta,d);
        test_param(4,:).SEG = -1;
    end
    if(Innerright(1) ~= 2*pi)
        if(~param.achievable.RSL)
            test_param(5,:) = dubins_CSCSC(p1, p2, r, Innerright, restriction, u.right, u.left, u.left);
            test_param(6,:) = dubins_CSCSC(p1, p2, r, Innerright, restriction, u.right, u.right, u.left);
        else
            test_param(5,:) = dubins_RSL(alpha,beta,d);
            test_param(6,:).SEG = -1;
        end
    else
        test_param(5,:).SEG = -1;
        test_param(6,:).SEG = -1;
    end
        
    
    if(Innerleft(1) ~= 2*pi)
        if(~param.achievable.LSR)
            test_param(7,:) = dubins_CSCSC(p1, p2, r, Innerleft, restriction, u.left, u.left, u.right);
            test_param(8,:) = dubins_CSCSC(p1, p2, r, Innerleft, restriction, u.left, u.right, u.right);
        else
            test_param(7,:) = dubins_LSR(alpha,beta,d);
            test_param(8,:).SEG = -1;
        end
    else
        test_param(7,:).SEG = -1;
        test_param(8,:).SEG = -1;
    end
    
    for i = 1:1:8
        if(test_param(i).SEG ~= -1) 
            cost = sum(test_param(i,:).SEG);
            if(cost < best_cost) || (best_cost == -1)
                best_word = i;
                best_cost = cost;
                param.SEG_param = test_param(i,:).SEG;
                param.angle = test_param(i,:).angle;
                param.type = i;
            end
        end
    end

    if(best_word == -1) 
        param.STATUS = -2;             % NO PATH
        return;
    else
        return;
    end
end

function param = dubins_CSCSC(p1, p2, r, bitangent, restriction, C1,C2,C3)
    dx = p2(1) + r*cos(p2(3)+C3*pi/2) - (p1(1)  + r*cos(p1(3)+C1*pi/2));
    dy = p2(2) + r*sin(p2(3)+C3*pi/2) - (p1(2)  + r*sin(p1(3)+C1*pi/2));
    
    aangle = bitangent(1);
    eangle = bitangent(2);
    
    [rd,ric1] = min(wrapTo2Pi(C1*(restriction - p1(3))));
    [rd,ric2] = min(wrapTo2Pi(C3*(p2(3)-restriction)));
    
    [a,ai] = min([wrapTo2Pi(C1*(bitangent(1) - p1(3))), wrapTo2Pi(C1*(restriction(ric1)-p1(3)))]);
    if(ai == 2 || wrapTo2Pi(C1*(bitangent(1) - p1(3))) == wrapTo2Pi(C1*(restriction(ric1)-p1(3))))
        aangle = restriction(ric1);
    end
    [e,ei] = min([wrapTo2Pi(C3*(p2(3) - bitangent(2))),wrapTo2Pi(C3*(p2(3)-restriction(ric2)))]);
    if(ei == 2 || wrapTo2Pi(C3*(p2(3) - bitangent(2))) == wrapTo2Pi(C3*(p2(3)-restriction(ric2))))
        eangle = restriction(ric2);
    end
    
    c = mod(C1*C2*C3*(eangle-aangle), 2*pi);
    if(c > mod(C1*C2*C3*(restriction((-C2+3)/2) - aangle),2*pi) || c == 0)
        param.SEG(1) = -1;
        param.angle = [0 0];
        return;
    end
    %if(c < pi)
    %    param.SEG(1) = -1;
    %    param.angle = [0 0];
    %    return;
    %end
    
    eangleb = eangle-C3*pi/2;
    aangleb = aangle+C1*pi/2;
    b = (-sin(eangleb).*dx+cos(eangleb).*dy)./sin(aangleb-eangleb);
    d = -1./cos(eangleb).*(cos(aangleb).*b-dx);
    
    i = 0;
    if(b <0 && b > -1000)
        while(b <0 && i < 1000)
            aangle = wrapTo2Pi(aangle - C1*0.001);
            eangleb = eangle-C3*pi/2;
            aangleb = aangle+C1*pi/2;
            b = (-sin(eangleb).*dx+cos(eangleb).*dy)./sin(aangleb-eangleb);
            d = -1./cos(eangleb).*(cos(aangleb).*b-dx);
            i = i+1;
        end
    end
    
    i = 0;
    if(d <0 && b > -1000)
        while(d <0 && i < 1000)
            eangle = wrapTo2Pi(eangle + C3*0.001);
            eangleb = eangle-C3*pi/2;
            aangleb = aangle+C1*pi/2;
            b = (-sin(eangleb).*dx+cos(eangleb).*dy)./sin(aangleb-eangleb);
            d = -1./cos(eangleb).*(cos(aangleb).*b-dx);
            i = i+1;
        end
    end
    
    if(b < 0 || d < 0)
        param.SEG(1) = -1;
        param.angle = [0 0];
        return;
    end        
    
    L = 2*abs(tan((aangle-eangle)/2));
    param.SEG(1) = wrapTo2Pi(C1*(aangle - p1(3))); 
    param.SEG(2) = d/r+L; 
    param.SEG(3) = c;
    param.SEG(4) = b/r+L; 
    param.SEG(5) = wrapTo2Pi(C3*(p2(3) - eangle));
    param.angle(1) = aangle;
    param.angle(2) = eangle;
end

function param = is_dubins_achievable(p1,p2, Ol, Or, Outerleft, Outerright, Innerleft, Innerright, restriction, r)
    param.LSL = -1;
    param.RSL = -1;
    param.RSR = -1;
    param.LSR = -1;
    
    u.left = 1;
    u.right = -1;

    Lstart = wrapTo2Pi(p1(3)-u.left*pi/2);
    Lfinal = wrapTo2Pi(p2(3)-u.left*pi/2);
    Rstart = wrapTo2Pi(p1(3)-u.right*pi/2);
    Rfinal = wrapTo2Pi(p2(3)-u.right*pi/2);

    Lwindminstart = wrapTo2Pi(restriction(1)-u.left*pi/2);
    Lwindminfinal = wrapTo2Pi(restriction(2)-u.left*pi/2);
    Rwindminstart = wrapTo2Pi(restriction(2)-u.right*pi/2);
    Rwindminfinal = wrapTo2Pi(restriction(1)-u.right*pi/2);
    
    Outerleft = Outerleft - [u.left*pi/2, u.left*pi/2];
    Outerright = Outerright - [u.right*pi/2, u.right*pi/2];

    %LSL
    Lts = mod(-u.left*(Lstart - Outerleft(1)),2*pi) < mod(-u.left*(Lstart - Lwindminstart),2*pi);
    Ltf = mod(-u.left*(Outerleft(2) - Lfinal),2*pi) < mod(-u.left*(Lwindminfinal - Lfinal),2*pi);
    param.LSL = Lts && Ltf;
    clear Lts Ltf

    %LSR
    if(norm(Or(2,:)-Ol(1,:)) > 2*r)
        Lts = mod(-u.left*(Lstart - Innerleft(1)),2*pi) < mod(-u.left*(Lstart - Lwindminstart),2*pi);
        Rtf = mod(-u.right*(Innerleft(2) - Rfinal),2*pi) < mod(-u.right*(Rwindminfinal - Rfinal),2*pi);
        param.LSR= Lts && Rtf;
        clear Lts Rtf
    else
        param.LSR = 0;
    end

    %RSR
    Rts = mod(-u.right*(Rstart - Outerright(1)),2*pi) < mod(-u.right*(Rstart - Rwindminstart),2*pi);
    Rtf = mod(-u.right*(Outerright(2) - Rfinal),2*pi) < mod(-u.right*(Rwindminfinal - Rfinal),2*pi);
    param.RSR = Rts && Rtf;
    clear Rts Rtf

    %RSL
    if(norm(Ol(2,:)-Or(1,:)) > 2*r)
        Rts = mod(-u.right*(Rstart - Innerright(1)),2*pi) < mod(-u.right*(Rstart - Rwindminstart),2*pi);
        Ltf = mod(-u.left*(Innerright(2) - Lfinal),2*pi) < mod(-u.left*(Lwindminfinal - Lfinal),2*pi);
        param.RSL = Rts && Ltf;
        clear Rts Ltf
    else
        param.RSL = 0;
    end
end

function param = dubins_LSL(alpha, beta, d)
    tmp0 = d + sin(alpha) - sin(beta);
    p_squared = 2 + (d*d) -(2*cos(alpha - beta)) + (2*d*(sin(alpha) - sin(beta)));
    if( p_squared < 0 )
        param = -1;
        return;
    else
        tmp1 = atan2( (cos(beta)-cos(alpha)), tmp0 );
        t = mod((-alpha + tmp1 ), 2*pi);
        p = sqrt( p_squared );
        q = mod((beta - tmp1 ), 2*pi);
        param.SEG(1) = t; 
        param.SEG(2) = p; 
        param.SEG(3) = q;
        param.angle = [0,0];
        return ;
    end
end
function param = dubins_LSR(alpha, beta, d)
    p_squared = -2 + (d*d) + (2*cos(alpha - beta)) + (2*d*(sin(alpha)+sin(beta)));
    if( p_squared < 0 )
        param = -1;
        return;
    else
        p    = sqrt( p_squared );
        tmp2 = atan2( (-cos(alpha)-cos(beta)), (d+sin(alpha)+sin(beta)) ) - atan2(-2.0, p);
        t    = mod((-alpha + tmp2), 2*pi);
        q    = mod(( -mod((beta), 2*pi) + tmp2 ), 2*pi);
        param.SEG(1) = t; 
        param.SEG(2) = p; 
        param.SEG(3) = q;
        param.angle = [0,0];
        return ;
    end
end
function param = dubins_RSL(alpha, beta, d)
    p_squared = (d*d) -2 + (2*cos(alpha - beta)) - (2*d*(sin(alpha)+sin(beta)));
    if( p_squared< 0 ) 
        param = -1; 
        return;
    else
        p    = sqrt( p_squared );
        tmp2 = atan2( (cos(alpha)+cos(beta)), (d-sin(alpha)-sin(beta)) ) - atan2(2.0, p);
        t    = mod((alpha - tmp2), 2*pi);
        q    = mod((beta - tmp2), 2*pi);
        param.SEG(1) = t;
        param.SEG(2) = p; 
        param.SEG(3) = q;
        param.angle = [0,0];
        return ;
    end
end
function param = dubins_RSR(alpha, beta, d)
    tmp0 = d-sin(alpha)+sin(beta);
    p_squared = 2 + (d*d) -(2*cos(alpha - beta)) + (2*d*(sin(beta)-sin(alpha)));
    if( p_squared < 0 )
        param = -1; 
        return;
    else
        tmp1 = atan2( (cos(alpha)-cos(beta)), tmp0 );
        t = mod(( alpha - tmp1 ), 2*pi);
        p = sqrt( p_squared );
        q = mod(( -beta + tmp1 ), 2*pi);
        param.SEG(1) = t; 
        param.SEG(2) = p; 
        param.SEG(3) = q;
        param.angle = [0,0];
        return;
    end
end
function param = dubins_RLR(alpha, beta, d)
    tmp_rlr = (6. - d*d + 2*cos(alpha - beta) + 2*d*(sin(alpha)-sin(beta))) / 8.;
    if( abs(tmp_rlr) > 1)
        param = -1; 
        return;
    else
        p = mod(( 2*pi - acos( tmp_rlr ) ), 2*pi);
        t = mod((alpha - atan2( cos(alpha)-cos(beta), d-sin(alpha)+sin(beta) ) + mod(p/2, 2*pi)), 2*pi);
        q = mod((alpha - beta - t + mod(p, 2*pi)), 2*pi);
        param.SEG(1) = t;
        param.SEG(2) = p;
        param.SEG(3) = q;
        param.angle = [0,0];
        
        return;
    end
end
function param = dubins_LRL(alpha, beta, d)
    tmp_lrl = (6. - d*d + 2*cos(alpha - beta) + 2*d*(- sin(alpha) + sin(beta))) / 8.;
    if( abs(tmp_lrl) > 1)
        param = -1; return;
    else
        p = mod(( 2*pi - acos( tmp_lrl ) ), 2*pi);
        t = mod((-alpha - atan2( cos(alpha)-cos(beta), d+sin(alpha)-sin(beta) ) + p/2), 2*pi);
        q = mod((mod(beta, 2*pi) - alpha -t + mod(p, 2*pi)), 2*pi);
        param.SEG(1) = t;
        param.SEG(2) = p;
        param.SEG(3) = q;
        param.angle = [0,0];
        return;
    end
end

