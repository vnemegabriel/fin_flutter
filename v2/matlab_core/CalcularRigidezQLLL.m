function [Ke_global] = CalcularRigidezQLLL(nodes3D, geometry, material, integrationType)
    % nodes3D: [4x3] Coordenadas globales de los 4 nodos del elemento
    % geometry.t: espesor
    % material.E, material.nu: Propiedades elásticas
    % integrationType: 'selective' (recomendado) o 'full'

    % 1. Parámetros y constantes
    E = material.E; 
    nu = material.nu; 
    t = geometry.t;
    Ke_local = zeros(24, 24);
    
    % 2. Sistema Local y Transformación
    [T, nodesLocal] = getTransformationMatrix(nodes3D);
    

    % 3. Matrices Constitutivas 

    % Membrana (Plane Stress)
    Dm = (E*t/(1-nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    
    % Flexión (Mindlin)
    Db = (E*t^3/(12*(1-nu^2))) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    
    % Corte (k=5/6)
    G = E / (2*(1+nu));
    Ds = (5/6) * G * t * eye(2);

    % INTEGRACIÓN DE MEMBRANA Y FLEXIÓN (2x2 full) 
    [gp, gw] = getGaussPoints(2);
    Area_elem = 0; % Para calcular la rigidez de drilling luego
    
    for i = 1:2    %Integracion Full en flexion
        for j = 1:2
            [Bm, Bb, detJ] = getMatricesMB(nodesLocal, gp(i), gp(j));
            WJ = detJ * gw(i) * gw(j);
            Area_elem = Area_elem + WJ;

            % Mapeo para Membrana (u, v) -> DOFs 1, 2 de cada nodo
            idxM = [1, 2, 7, 8, 13, 14, 19, 20];
            Ke_local(idxM, idxM) = Ke_local(idxM, idxM) + Bm' * Dm * Bm * WJ;
            
            % Mapeo para Flexión (th_x, th_y) -> DOFs 4, 5 de cada nodo
            idxB = [4, 5, 10, 11, 16, 17, 22, 23];
            Ke_local(idxB, idxB) = Ke_local(idxB, idxB) + Bb' * Db * Bb * WJ;
        end
    end

    % INTEGRACIÓN DE CORTE (Selectiva o Full)
    if strcmpi(integrationType, 'full')
        nGaussS = 2;
    else
        nGaussS = 1;
    end
    
    [gpS, gwS] = getGaussPoints(nGaussS);
    for i = 1:nGaussS
        for j = 1:nGaussS
            [Bs, detJ] = getMatrixS(nodesLocal, gpS(i), gpS(j));
            WJ = detJ * gwS(i) * gwS(j);
            % Mapeo para Corte (w, th_x, th_y) -> DOFs 3, 4, 5 de cada nodo
            idxS = [3, 4, 5, 9, 10, 11, 15, 16, 17, 21, 22, 23];
            Ke_local(idxS, idxS) = Ke_local(idxS, idxS) + Bs' * Ds * Bs * WJ;
        end
    end

    % RIGIDEZ DE DRILLING 
    % Evita singularidad en la rotación normal al plano local
    k_drill = 1e-3 * E * t * Area_elem; % Ver justificacion en Onate
    idxD = [6, 12, 18, 24]; % DOFs th_z local
    for d = 1:4
        Ke_local(idxD(d), idxD(d)) = Ke_local(idxD(d), idxD(d)) + k_drill;
    end

    % 4. Transformación al Sistema Global
    Ke_global = T' * Ke_local * T;
end


%% Funciones Extras

function [gp, gw] = getGaussPoints(n)
    % n: número de puntos por dirección (1 o 2)
    if n == 1
        gp = 0; gw = 2;
    elseif n == 2
        gp = [-1/sqrt(3), 1/sqrt(3)];
        gw = [1, 1];
    end
end

function [N, dN] = shapeQ4(xi, eta)
    % N: Funciones de forma
    % dN: Derivadas respecto a xi (fila 1) y eta (fila 2)
    N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    
    dN = 0.25 * [-(1-eta),  (1-eta), (1+eta), -(1+eta);  % dN/dxi
                 -(1-xi),  -(1+xi),  (1+xi),  (1-xi)];  % dN/deta
end

function [Bm, Bb, detJ] = getMatricesMB(nodesLoc, xi, eta)
    [~, dN_nat] = shapeQ4(xi, eta);
    J = dN_nat * nodesLoc;
    detJ = det(J);
    dN = inv(J) * dN_nat;
    
    % Membrana (Igual que antes)
    Bm = zeros(3, 8);
    for i = 1:4
        col = (i-1)*2 + 1;
        Bm(:, col:col+1) = [dN(1,i), 0; 0, dN(2,i); dN(2,i), dN(1,i)];
    end
    
    % Flexión: Relaciona rot_x y rot_y con curvaturas
    % kappa_x = d(rot_y)/dx
    % kappa_y = -d(rot_x)/dy
    Bb = zeros(3, 8);
    for i = 1:4
        col = (i-1)*2 + 1;
        Bb(:, col:col+1) = [0,       dN(1,i);
                           -dN(2,i), 0;
                           -dN(1,i), dN(2,i)];
    end
end

function [Bs, detJ] = getMatrixS(nodesLoc, xi, eta)
    % Relaciona {w1, th_x1, th_y1, ...} con {gamma_xz, gamma_yz}
    [N, dN_nat] = shapeQ4(xi, eta);
    J = dN_nat * nodesLoc;
    detJ = det(J);
    dN = inv(J) * dN_nat;
    
    Bs = zeros(2, 12); % 4 nodos * 3 DOFs (w, th_x, th_y)
    for i = 1:4
        col = (i-1)*3 + 1;
        Bs(:, col:col+2) = [dN(1,i),    0,    N(i);
                            dN(2,i), -N(i),    0 ];
    end
end

