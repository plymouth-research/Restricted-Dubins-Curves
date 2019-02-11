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


function [ param, all ] = dubins_restricted_core( p1, p2, r, psi, delta )
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
    
    % Bitangent angles
    %Outerleft = zeros(1,2);
    %Outerright = zeros(1,2);
    Innerleft = 2*pi*ones(1,2);
    Innerright = 2*pi*ones(1,2);
    
    Outerleft = wrapTo2Pi(atan2(Ol(2,2)-Ol(1,2), Ol(2,1)-Ol(1,1)))-u.left*pi/2;
    Outerleft(2) = Outerleft(1);
    
    TangentsL = wrapTo2Pi(atan2(Ol(2,2)-Ol(1,2), Ol(2,1)-Ol(1,1)))+acos(norm(Ol(2,:)-Ol(1,:))/(4*r));
    TangentsL(2) = wrapTo2Pi(atan2(Ol(2,2)-Ol(1,2), Ol(2,1)-Ol(1,1)))-pi-acos(norm(Ol(2,:)-Ol(1,:))/(4*r));
    
    TangentsR = wrapTo2Pi(atan2(Or(2,2)-Or(1,2), Or(2,1)-Or(1,1)))-acos(norm(Or(2,:)-Or(1,:))/(4*r));
    TangentsR(2) = wrapTo2Pi(atan2(Or(2,2)-Or(1,2), Or(2,1)-Or(1,1)))+pi+acos(norm(Or(2,:)-Or(1,:))/(4*r));

    if(norm(Or(2,:)-Ol(1,:)) > 2*r)
        Innerleft = wrapTo2Pi(atan2(Or(2,2)-Ol(1,2), Or(2,1)-Ol(1,1))-u.left*acos(2*r/norm(Or(2,:)-Ol(1,:))));
        Innerleft(2) = wrapTo2Pi(Innerleft(1)+pi);
    end

    Outerright = wrapTo2Pi(atan2(Or(2,2)-Or(1,2), Or(2,1)-Or(1,1)))-u.right*pi/2;
    Outerright(2) = Outerright(1);

    if(norm(Ol(2,:)-Or(1,:)) > 2*r)
        Innerright = wrapTo2Pi(atan2(Ol(2,2)-Or(1,2), Ol(2,1)-Or(1,1))-u.right*acos(2*r/norm(Ol(2,:)-Or(1,:))));
        Innerright(2) = wrapTo2Pi(Innerright(1)+pi);
    end
    
    param.achievable = is_dubins_achievable(p1, p2, Ol, Or, Outerleft, Outerright, Innerleft, Innerright, TangentsL, TangentsR, restriction, r);

    disp(param.achievable)
    %% Calculate all possible curves if not achievable
    
    test_param(9,:).SEG = -1;
    test_param(10,:).SEG = -1;
    test_param(11,:).SEG = -1;
    test_param(12,:).SEG = -1;
    test_param(9,:).angle = -1;
    test_param(10,:).angle = -1;
    test_param(11,:).angle = -1;
    test_param(12,:).angle = -1;
    % Compute path for each configuration
    if(~param.achievable.LSL && ~param.achievable.LRL)
        test_param(1,:) = dubins_CSCSC2(p1, p2, r, Outerleft, psi, delta, u.left, u.left, u.left);
        test_param(2,:) = dubins_CSCSC2(p1, p2, r, Outerleft, psi, delta, u.left, u.right, u.left);
    else
        if(norm(Or(2,:)-Ol(1,:)) > 4*r)
            test_param(1,:) = dubins_LSL(alpha,beta,d);
            test_param(2,:).SEG = -1;
        else
            test_param(1,:).SEG = -1;
            test_param(2,:).SEG = -1;
            test_param(10,:) = dubins_LRL(alpha,beta,d);
        end
    end
    
    if(~param.achievable.RSR && ~param.achievable.RLR)
        test_param(3,:) = dubins_CSCSC2(p1, p2, r, Outerright, psi, delta, u.right, u.left, u.right);
        test_param(4,:) = dubins_CSCSC2(p1, p2, r, Outerright, psi, delta, u.right, u.right, u.right);
            test_param(11,:).SEG = -1;
            test_param(12,:).SEG = -1;
    else
        if(norm(Or(2,:)-Ol(1,:)) > 4*r)
            test_param(3,:) = dubins_RSR(alpha,beta,d);
            test_param(4,:).SEG = -1;
        else
            test_param(3,:).SEG = -1;
            test_param(4,:).SEG = -1;
            test_param(12,:) = dubins_RLR(alpha,beta,d);
        end
    end
    
    if(norm(Ol(2,:)-Or(1,:)) > 2*r)
        if(~param.achievable.RSL)
            test_param(5,:) = dubins_CSCSC2(p1, p2, r, Innerright, psi, delta, u.right, u.left, u.left);
            test_param(6,:) = dubins_CSCSC2(p1, p2, r, Innerright, psi, delta, u.right, u.right, u.left);
        else
            test_param(5,:) = dubins_RSL(alpha,beta,d);
            test_param(6,:).SEG = -1;
        end
    else
        test_param(5,:).SEG = -1;
        test_param(6,:).SEG = -1;
    end
        
    
    if(norm(Or(2,:)-Ol(1,:)) > 2*r)
        if(~param.achievable.LSR)
            test_param(7,:) = dubins_CSCSC2(p1, p2, r, Innerleft, psi, delta, u.left, u.left, u.right);
            test_param(8,:) = dubins_CSCSC2(p1, p2, r, Innerleft, psi, delta, u.left, u.right, u.right);
        else
            test_param(7,:) = dubins_LSR(alpha,beta,d);
            test_param(8,:).SEG = -1;
        end
    else
        test_param(7,:).SEG = -1;
        test_param(8,:).SEG = -1;
    end
    
    % Find the shortest one
    best_word = -1;
    best_cost = -1;
    for i = 1:1:12
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
    
    for i = 1:1:12
        all(i).SEG_param = test_param(i,:).SEG;
        all(i).angle = test_param(i,:).angle;
        all(i).type = i;
        all(i).p_init = p1;
        all(i).r = r;
        all(i).STATUS = 0;
        all(i).achievable = 0;
    end

    if(best_word == -1) 
        param.STATUS = -2;             % NO PATH
        return;
    else
        return;
    end
end

function param = dubins_CSCSC2(p1, p2, r, bitangent, psi, delta, C1,C2,C3)
    dx = p2(1) + r*cos(p2(3)+C3*pi/2) - (p1(1)  + r*cos(p1(3)+C1*pi/2));
    dy = p2(2) + r*sin(p2(3)+C3*pi/2) - (p1(2)  + r*sin(p1(3)+C1*pi/2));

    % Tranforming to circle frame
    theta1b = wrapTo2Pi(p1(3)-C1*pi/2);
    theta2b = wrapTo2Pi(p2(3)-C3*pi/2);

    %domain1 = [theta1b, restriction((-C1+3)/2)-C1*pi/2];
    %domain2 = [theta2b, restriction((C3+3)/2)-C3*pi/2];
    domain1 = [wrapToPi(theta1b), wrapToPi(wrapTo2Pi(psi)-C1*wrapTo2Pi(delta+pi/2))];
    domain2 = [wrapToPi(wrapTo2Pi(psi)+C3*wrapTo2Pi(delta-pi/2))-C3*2*pi, wrapToPi(theta2b)];

    %if(C1*domain1(1) > C1*domain1(2))
    %    domain1(1) = domain1(1)-2*pi;
    %end
    %if(domain2(1) > domain2(2))
    %    domain2(1) = domain2(1)-2*pi;
    %end

    domain1d1 = wrapTo2Pi(domain1(1):C1*0.01:domain1(2));
    domain1d2 = wrapTo2Pi(domain2(1):C3*0.01:domain2(2));

    d1 = repelem(domain1d1,size(domain1d2,2));
    d2 = repmat(domain1d2,[1,size(domain1d1,2)]);

    d1C2 = wrapTo2Pi(d1+(-C1*C2+1)/2*pi);
    d2C2 = wrapTo2Pi(d2+(-C2*C3+1)/2*pi);
    c = wrapTo2Pi(C2*(d2C2-d1C2));

    d1c = d1(c <= wrapTo2Pi(C2*(wrapTo2Pi(wrapTo2Pi(psi)-C2*wrapTo2Pi(delta+pi/2)) - d1C2)));
    d2c = d2(c <= wrapTo2Pi(C2*(wrapTo2Pi(wrapTo2Pi(psi)-C2*wrapTo2Pi(delta+pi/2)) - d1C2)));

    d1C2 = wrapTo2Pi(d1c+(-C1*C2+1)/2*pi);
    d2C2 = wrapTo2Pi(d2c+(-C2*C3+1)/2*pi);
    c = wrapTo2Pi(C2*(d2C2-d1C2));

    [b,d]=distTT2(d1c',d2c',dx,dy,r,c<pi,C1,C3);
    bp = b(b>=0&d>=0);
    dp = d(b>=0&d>=0);

    d1p = d1c(b>=0&d>=0);
    d2p = d2c(b>=0&d>=0);

    d1s = abs(wrapToPi(bitangent(1))-wrapToPi(d1p));
    d2s = abs(wrapToPi(d2p) - wrapToPi(bitangent(2)));
    [a,ai]=min(d1s);
    [e,ei]=min(d2s);
    [c,ci]=min(bp+dp);
    %%
    if(size(ai,2) == 0 || size(ei,2) == 0)
        param.SEG(1) = -1;
        param.angle = [0 0];
        return;
    end
    a = [d1p(ai(1)), d2p(ai(1)); d1p(ei(1)), d2p(ei(1))];
    tta = [bp(ai(1)),dp(ai(1));bp(ei(1)),dp(ei(1))];
    toc

    [anglen,anglei] = min([sum(tta(1,:)),sum(tta(2,:))]);
    aangle = a(anglei,1);
    eangle = a(anglei,2);
    bdist = tta(anglei,1);
    ddist = tta(anglei,2);

    %aangle = d1p(ci);
    %eangle = d2p(ci);
    %bdist = bp(ci);
    %ddist = dp(ci);
    
    aangleC2 = wrapTo2Pi(aangle+(-C1*C2+1)/2*pi);
    eangleC2 = wrapTo2Pi(eangle+(-C2*C3+1)/2*pi);
    cangle = wrapTo2Pi(C2*(eangleC2-aangleC2));
    
    param.SEG(1) = wrapTo2Pi(C1*(aangle - theta1b)); 
    param.SEG(2) = bdist/r;
    param.SEG(3) = cangle;
    param.SEG(4) = ddist/r;
    param.SEG(5) = wrapTo2Pi(C3*(theta2b - eangle));
    param.angle(1) = aangle;
    param.angle(2) = eangle;
end

function param = dubins_CSCSC(p1, p2, r, bitangent, restriction, C1,C2,C3)
    % Measure distances between circles
    dx = p2(1) + r*cos(p2(3)+C3*pi/2) - (p1(1)  + r*cos(p1(3)+C1*pi/2));
    dy = p2(2) + r*sin(p2(3)+C3*pi/2) - (p1(2)  + r*sin(p1(3)+C1*pi/2));
    
    % Tranforming to circle frame
    theta1b = p1(3)-C1*pi/2;
    theta2b = p2(3)-C3*pi/2;
    
    % Compute angle on C1 and C3, closest to the bitangent
    aangle = bitangent(1);
    eangle = bitangent(2);
    %[rd,ric1] = min(wrapTo2Pi(C1*(restriction - p1(3))));
    %[rd,ric2] = min(wrapTo2Pi(C3*(p2(3)-restriction)));
    
    [a,ai] = min([wrapTo2Pi(C1*(bitangent(1)-theta1b)), wrapTo2Pi(C1*(restriction((-C1+3)/2)-C1*pi/2-theta1b))]);%, wrapTo2Pi(-C1*(bitangent(1)-theta1b))]);
    %if(ai == 2 || wrapTo2Pi(C1*(bitangent(1) - theta1b)) == wrapTo2Pi(C1*(restriction((-C1+3)/2)-C1*pi/2-theta1b)))
    %    aangle = restriction(ric1);
    %end
    if(ai == 1)
        aangle = bitangent(1);
    elseif(ai == 2)
        aangle = wrapTo2Pi(restriction((-C1+3)/2)-C1*pi/2);
    else
        aangle = theta1b;
    end
    [e,ei] = min([wrapTo2Pi(C3*(theta2b - bitangent(2))),wrapTo2Pi(C3*(theta2b-restriction((C3+3)/2)+C3*pi/2))]);%,wrapTo2Pi(-C3*(theta2b - bitangent(2)))]);
    %if(ei == 2 || wrapTo2Pi(C3*(theta2b - bitangent(2))) == wrapTo2Pi(C3*(theta2b-restriction((-C3+3)/2)+C3*pi/2)))
    %    eangle = restriction(ric2);
    %end
    if(ei == 1)
        eangle = bitangent(2);
    elseif(ei == 2)
        eangle = wrapTo2Pi(restriction((C3+3)/2)-C3*pi/2);
    %else
    %    eangle = theta2b;
    end
    
    % Check if C2 is in the no-go zone
    c = wrapTo2Pi(C2*(eangle-aangle));
    bitangentC2 = aangle + (-C1*C2+1)/2*pi;
    if(c > wrapTo2Pi(C2*(restriction((-C2+3)/2)-C2*pi/2 - bitangentC2)))
        param.SEG(1) = -1;
        param.angle = [0 0];
        return;
    end
    
    % Calculate Segments algebraic distances (based on the intersection)
    [b,d] = distTT(aangle+C1*pi/2,eangle-C3*pi/2,dx,dy,r);
    
    % Fix domains of intersection
    i = 0;
    if(b <0 && b > -10000)
        while(b <0 && i < 1000)
            aangle = wrapTo2Pi(aangle - C1*0.001);
            %eangle = wrapTo2Pi(eangle + C3*0.001);
            [b,d] = distTT(aangle-C1*pi/2,eangle-C3*pi/2,dx,dy,r);
            i = i+1;
        end
    end
    
    i = 0;
    if(d <0 && d > -1000)
        while(d <0 && i < 1000)
            eangle = wrapTo2Pi(eangle + C3*0.001);
            %aangle = wrapTo2Pi(aangle - C1*0.001);
            [b,d] = distTT(aangle-C1*pi/2,eangle-C3*pi/2,dx,dy,r);
            i = i+1;
        end
    end
    
    % Check if intersecting
    if(b < 0 || d < 0)
        param.SEG(1) = -1;
        param.angle = [0 0];
        return;
    end        
    
    % Added length to path
    L = 2*abs(tan((aangle-eangle)/2));
    %L = 0;
    
    param.SEG(1) = wrapTo2Pi(C1*(aangle - theta1b)); 
    param.SEG(2) = d/r+L; 
    param.SEG(3) = c;
    param.SEG(4) = b/r+L; 
    param.SEG(5) = wrapTo2Pi(C3*(theta2b - eangle));
    param.angle(1) = aangle;
    param.angle(2) = eangle;
end

function param = is_dubins_achievable(p1,p2, Ol, Or, Outerleft, Outerright, Innerleft, Innerright, TangentsL, TangentsR, restriction, r)
    param.LSL = -1;
    param.RSL = -1;
    param.RSR = -1;
    param.LSR = -1;
    param.LRL = -1;
    param.RLR = -1;
    
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

    %LSL
    Lts = wrapTo2Pi(-u.left*(Lstart - Outerleft(1))) < wrapTo2Pi(-u.left*(Lstart - Lwindminstart));
    Ltf = wrapTo2Pi(-u.left*(Outerleft(2) - Lfinal)) < wrapTo2Pi(-u.left*(Lwindminfinal - Lfinal));
    param.LSL = Lts && Ltf;
    clear Lts Ltf

    %LSR
    if(norm(Or(2,:)-Ol(1,:)) > 2*r)
        Lts = wrapTo2Pi(-u.left*(Lstart - Innerleft(1))) < wrapTo2Pi(-u.left*(Lstart - Lwindminstart));
        Rtf = wrapTo2Pi(-u.right*(Innerleft(2) - Rfinal)) < wrapTo2Pi(-u.right*(Rwindminfinal - Rfinal));
        param.LSR= Lts && Rtf;
        clear Lts Rtf
    else
        param.LSR = 0;
    end

    %RSR
    Rts = wrapTo2Pi(-u.right*(Rstart - Outerright(1))) < wrapTo2Pi(-u.right*(Rstart - Rwindminstart));
    Rtf = wrapTo2Pi(-u.right*(Outerright(2) - Rfinal)) < wrapTo2Pi(-u.right*(Rwindminfinal - Rfinal));
    param.RSR = Rts && Rtf;
    clear Rts Rtf

    %RSL
    if(norm(Ol(2,:)-Or(1,:)) > 2*r)
        Rts = wrapTo2Pi(-u.right*(Rstart - Innerright(1))) < wrapTo2Pi(-u.right*(Rstart - Rwindminstart));
        Ltf = wrapTo2Pi(-u.left*(Innerright(2) - Lfinal)) < wrapTo2Pi(-u.left*(Lwindminfinal - Lfinal));
        param.RSL = Rts && Ltf;
        clear Rts Ltf
    else
        param.RSL = 0;
    end

    %LRL
    if(norm(Ol(2,:)-Ol(1,:)) <= 4*r)
        Lts = wrapTo2Pi(-u.left*(Lstart - TangentsL(1))) < wrapTo2Pi(-u.left*(Lstart - Lwindminstart));
        Ltf = wrapTo2Pi(-u.left*(TangentsL(2) - Lfinal)) < wrapTo2Pi(-u.left*(Lwindminfinal - Lfinal));
        param.LRL = Lts && Ltf;
        clear Lts Ltf
    else
        param.LRL = 0;
    end

    %RLR
    if(norm(Or(2,:)-Or(1,:)) <= 4*r)
        Rts = wrapTo2Pi(-u.right*(Rstart - TangentsR(1))) < wrapTo2Pi(-u.right*(Rstart - Rwindminstart));
        Rtf = wrapTo2Pi(-u.right*(TangentsR(2) - Rfinal)) < wrapTo2Pi(-u.right*(Rwindminfinal - Rfinal));
        param.RLR = Rts && Rtf;
        clear Rts Rtf
    else
        param.RLR = 0;
    end
end

function [b,d] = distTT2(theta1,theta2, dx, dy,r,s,C1,C3)
    theta1b = theta1+C1*pi/2;
    theta2b = theta2+C3*pi/2;
    dx = dx + r*(cos(theta2) - cos(theta1));
    dy = dy + r*(sin(theta2) - sin(theta1));
    L = r*abs(tan((theta1b-theta2b)/2));
    
    flag = (cos(theta1b) ~= 0);
    theta1b1 = theta1b(flag);
    theta2b1 = theta2b(flag);
    dx1 = dx(flag);
    dy1 = dy(flag);
    b(flag) = (sin(theta1b1).*dx1-cos(theta1b1).*dy1)./sin(theta2b1-theta1b1);
    d(flag) = -1./cos(theta1b1).*(cos(theta2b1).*b(flag)'+dx1);
    
    flag = (cos(theta2b) ~= 0);
    theta1b1 = theta1b(flag);
    theta2b1 = theta2b(flag);
    dx1 = dx(flag);
    dy1 = dy(flag);
    b(flag) = (-sin(theta2b1).*dx1+cos(theta2b1).*dy1)./sin(theta1b1-theta2b1);
    d(flag) = -1./cos(theta2b1).*(cos(theta1b1).*b(flag)'-dx1);
    
    flag = (sin(theta1b) ~= 0);
    theta1b1 = theta1b(flag);
    theta2b1 = theta2b(flag);
    dx1 = dx(flag);
    dy1 = dy(flag);
    b(flag) = (cos(theta1b1).*dy1-sin(theta1b1).*dx1)./sin(theta1b1-theta2b1);
    d(flag) = -1./sin(theta1b1).*(sin(theta2b1).*b(flag)'+dy1);
    
    flag = (sin(theta2b) ~= 0);
    theta1b1 = theta1b(flag);
    theta2b1 = theta2b(flag);
    dx1 = dx(flag);
    dy1 = dy(flag);
    b(flag) = (-cos(theta2b1).*dy1+sin(theta2b1).*dx1)./sin(theta2b1-theta1b1);
    d(flag) = -1./sin(theta2b1).*(sin(theta1b1).*b(flag)'-dy1);
    
    
    b = b + (-2*s+1).*L';
    d = d + (-2*s+1).*L';
end

function [b,d] = distTT(theta1b,theta2b, dx, dy,r)
    %b = (-sin(eangle).*dx+cos(eangle).*dy)./sin(aangle-eangle);
    %d = -1./cos(eangle).*(cos(aangle).*b-dx);
    if(abs(sin(theta2b-theta1b) ) < 10e-10)
        b = -1;
        d = -1;
        return;
    end
    if(abs(cos(theta2b)) >  10e-10)
        b = (sin(theta1b).*dx-cos(theta1b).*dy)./sin(theta2b-theta1b);
        d = -1./cos(theta1b).*(cos(theta2b).*b+dx);
    elseif(abs(sin(theta2b)) >  10e-10)
        b = -(cos(theta1b).*dy-sin(theta1b).*dx)./sin(theta1b-theta2b);
        d = 1./sin(theta1b).*(sin(theta2b).*b+dy);
    end
    %b = -(cos(theta1b).*dy-sin(theta1b).*dx)./sin(theta1b-theta2b)
    %d = 1./sin(theta1b).*(sin(theta2b).*b+dy)
    %b = -(-sin(theta2b).*dx+cos(theta2b).*dy)./sin(theta1b-theta2b)
    %d = -1./cos(theta2b).*(cos(theta1b).*b-dx)
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

