function R_finale = ElemRotationMatrix(order, varargin)

    
    % Definizione delle funzioni anonime per le rotazioni attorno agli assi
    % order: specifica l'ordine di applicazione delle rotazioni 'zxz', 'zy', 'x' ecc.. 
    % varagin: Se si vuole una valutazione numerica inserisci gli angoli
    % necessari alla rotazione. Ad esempio: ('xyz', pi/2, pi, pi/3) o
    % ('xy', pi/2, pi/3) ecc..
    check = true;

    if length(varargin) ~= length(order)
        error('Numero di angoli (%d) non corrisponde alla lunghezza dell''ordine (%d)', ...
            length(varargin), length(order));
    end

    for i = 1:length(order)
        if isa(varargin{i}, 'sym')
            check=false;
        end
    end
         

    RotX = @(alpha) [1, 0, 0; 
                     0, cos(alpha), -sin(alpha); 
                     0, sin(alpha), cos(alpha)];
                 
    RotY = @(beta) [cos(beta), 0, sin(beta); 
                     0, 1, 0; 
                     -sin(beta), 0, cos(beta)];
                 
    RotZ = @(gamma) [cos(gamma), -sin(gamma), 0; 
                     sin(gamma), cos(gamma), 0; 
                     0, 0, 1];
    
    
    % Inizializza R_finale con la matrice identità
    R_finale = eye(3);
    
    % Applica le rotazioni in base all'ordine specificato
    for i = 1:length(order)
        angle = varargin{i};
        switch upper(order(i))
            case 'X'
                R_finale = R_finale * RotX(angle);
            case 'Y'
                R_finale = R_finale * RotY(angle);
            case 'Z'
                R_finale = R_finale * RotZ(angle);
            otherwise
                error('Asse non valido: %s. Usa solo X, Y o Z.', order(i));
        end
    end
    if check
        check_rotation_matrix(R_finale)
    end

end
