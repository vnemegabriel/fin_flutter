function [Cb, Cs] = getElementConstitutive(material, t)
    % Extracción de propiedades del material
    E = material.E;
    G = material.E/(2*(1+material.nu));
    nu = material.nu;
    
    % Término de rigidez de deformación plana
    Ep = E / (1 - nu^2); 
    
    % Módulo de corte modificado con factor de corrección de Mindlin (5/6)
    Gp = (5/6) * G; 
    
    % Matriz constitutiva de flexión generalizada (Cb)
    Cb = (t^3 / 12) * [Ep,      nu*Ep,   0;
                       nu*Ep,   Ep,      0;
                       0,       0,       G];
                   
    % Matriz constitutiva de corte transversal generalizada (Cs)
    Cs = t * [Gp, 0;
              0,  Gp];
end