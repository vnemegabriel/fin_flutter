function [N_res, M_res, Q_res] = CalcularEsfuerzosCompletos(mesh, geometry, material, D_global, evalPoints)

    nEle = size(mesh.connect, 1);
    nEval = size(evalPoints, 2);
    t = geometry.t; 
    E = material.E; 
    nu = material.nu;
    G = E / (2*(1+nu));
    
    Dm = (E*t/(1-nu^2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
    Db = (E*t^3/(12*(1-nu^2))) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
    Ds = (5/6) * G * t * eye(2);
    
    N_res = zeros(3, nEval, nEle); 
    M_res = zeros(3, nEval, nEle); 
    Q_res = zeros(2, nEval, nEle);
    
    isQuad = (size(mesh.connect, 2) == 4);
    
    for iEle = 1:nEle
        numNodes = size(mesh.connect, 2);
        nodesElem = mesh.nodes(mesh.connect(iEle, :), :);
        [T_full, nodesLoc] = getTransformationMatrix(nodesElem);
        
        eleDofs = zeros(1, numNodes*6);
        for iNode = 1:numNodes
            eleDofs((iNode-1)*6+1 : iNode*6) = (mesh.connect(iEle, iNode)-1)*6 + (1:6);
        end
        D_local = T_full * D_global(eleDofs);
        
        x = nodesLoc(:,1); y = nodesLoc(:,2);
        
        if isQuad
            %QLQL
            Bg1_1 = zeros(1, 24); Bg1_3 = zeros(1, 24);
            Bg2_2 = zeros(1, 24); Bg2_4 = zeros(1, 24);
            Bg1_1(3)=-0.5; Bg1_1(9)=0.5; Bg1_1(4)=-0.25*(y(2)-y(1)); Bg1_1(10)=-0.25*(y(2)-y(1)); Bg1_1(5)=0.25*(x(2)-x(1)); Bg1_1(11)=0.25*(x(2)-x(1));
            Bg1_3(21)=-0.5; Bg1_3(15)=0.5; Bg1_3(22)=-0.25*(y(3)-y(4)); Bg1_3(16)=-0.25*(y(3)-y(4)); Bg1_3(23)=0.25*(x(3)-x(4)); Bg1_3(17)=0.25*(x(3)-x(4));
            Bg2_2(9)=-0.5; Bg2_2(15)=0.5; Bg2_2(10)=-0.25*(y(3)-y(2)); Bg2_2(16)=-0.25*(y(3)-y(2)); Bg2_2(11)=0.25*(x(3)-x(2)); Bg2_2(17)=0.25*(x(3)-x(2));
            Bg2_4(3)=-0.5; Bg2_4(21)=0.5; Bg2_4(4)=-0.25*(y(4)-y(1)); Bg2_4(22)=-0.25*(y(4)-y(1)); Bg2_4(5)=0.25*(x(4)-x(1)); Bg2_4(23)=0.25*(x(4)-x(1));
        end
        for p = 1:nEval
            xi = evalPoints(1, p);
            eta = evalPoints(2, p);
            
            if isQuad
                dN_nat = 0.25 * [-(1-eta),  (1-eta),  (1+eta), -(1+eta);
                                 -(1-xi),  -(1+xi),   (1+xi),   (1-xi)];
                J = dN_nat * nodesLoc;
                dN_xy = J \ dN_nat;
                
                Bm = zeros(3, 24); Bb = zeros(3, 24);
                for k = 1:4
                    col = 6*(k-1);
                    Bm(1, col+1) = dN_xy(1,k);
                    Bm(2, col+2) = dN_xy(2,k);
                    Bm(3, col+1) = dN_xy(2,k); Bm(3, col+2) = dN_xy(1,k);
                    
                    Bb(1, col+5) =  dN_xy(1,k); Bb(2, col+4) = -dN_xy(2,k);
                    Bb(3, col+5) =  dN_xy(2,k); Bb(3, col+4) = -dN_xy(1,k);
                end
                
                Bg1 = 0.5*(1 - eta)*Bg1_1 + 0.5*(1 + eta)*Bg1_3;
                Bg2 = 0.5*(1 - xi)*Bg2_4  + 0.5*(1 + xi)*Bg2_2;
                Bs = J \ [Bg1; Bg2];
            end
            
            N_res(:, p, iEle) = Dm * Bm * D_local;
            M_res(:, p, iEle) = Db * Bb * D_local;
            Q_res(:, p, iEle) = Ds * Bs * D_local;
        end
    end
end